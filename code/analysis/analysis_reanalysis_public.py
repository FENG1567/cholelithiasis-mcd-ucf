#!/usr/bin/env python3
"""Frozen, local, streaming reanalysis for MCD attribution-gap revision.

Reads NCHS public-use files without modifying them.  It stores only aggregate,
de-identified statistics in the assigned revision folder and uses at most eight
top-level worker processes (no nested parallelism).
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, multiprocessing as mp, os, platform, subprocess
import statistics, sys, time, zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

OUT = Path("results")
NCHS = Path("data/raw")
XW_PATH = Path("config/schema_crosswalk.csv")
ONTOLOGY_PATH = Path("data/derived/ucd_destination_ontology.csv")
YEARS = list(range(1999, 2025))
REC_LEN = {y: (440 if y < 2003 else 488 if y < 2013 else 490 if y < 2020 else 817) for y in YEARS}
SEVEN_ZIP = Path("7z")
DISEASES = {
    "I21": {"I21"}, "C18": {"C18"}, "E10_E14": {"E10", "E11", "E12", "E13", "E14"},
    "N18": {"N18"}, "A40_A41": {"A40", "A41"},
}
# The candidate-owned, mutually exclusive 113-cause-recode destination ontology.
# K81–K83 exact UCD prefixes supersede this recode map when classifying records.
EXPECTED_DESTINATIONS = {"CIRCULATORY","DIABETES_NUTRITIONAL","DIGESTIVE_OTHER","EXTERNAL","INFECTIOUS_PARASITIC","NEOPLASMS","OTHER","OTHER_GALLBLADDER_BILIARY","RESPIRATORY"}
# Pre-specified sparse-cell rebin for fixed-weight standardization: 0-14 and
# 15-24 are combined before any outcome calculation, never reweighted by year.
AGE_ORDER = ["0-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75-84", "85+"]
SEX_ORDER = ["Male", "Female"]
RACE40_GROUPS = ("White", "Black", "AIAN", "Asian", "NHOPI", "Multiple")
RACE40_GROUP_MAP = {
    **{"01": "White", "02": "Black", "03": "AIAN"},
    **{f"{i:02d}": "Asian" for i in range(4, 11)},
    **{f"{i:02d}": "NHOPI" for i in range(11, 15)},
    **{f"{i:02d}": "Multiple" for i in range(15, 41)},
}

def norm(v: str) -> str:
    return (v or "").strip().upper().replace(".", "")

def age_years(s: str):
    s = (s or "").strip()
    if len(s) < 3 or not s[0].isdigit() or s[0] == "9": return None
    unit, value = int(s[0]), s[1:4] if len(s) >= 4 else s[1:3]
    if not value.isdigit(): return None
    n = int(value)
    if unit in (0, 1): return float(n)
    if unit == 2: return n / 12.0
    if unit in (3, 4, 5, 6): return 0.0
    return None

def age_group(s: str) -> str:
    a = age_years(s)
    if a is None: return "unknown"
    for lo, hi, label in [(0,14,"0-14"),(15,24,"15-24"),(25,34,"25-34"),(35,44,"35-44"),
                          (45,54,"45-54"),(55,64,"55-64"),(65,74,"65-74"),(75,84,"75-84")]:
        if lo <= a <= hi: return label
    return "85+" if a >= 85 else "unknown"

def sex_group(s: str) -> str:
    return "Male" if s.strip() in ("M", "1") else "Female" if s.strip() in ("F", "2") else "unknown"

def pos(row: dict, name: str):
    s = (row.get(name) or "").strip()
    if not s or s == "MISSING" or "-" not in s: return None
    parts=s.split("-")
    if len(parts)!=2 or not all(x.isdigit() for x in parts): return None
    a, b = parts; return int(a)-1, int(b)

def slc(rec: bytes, p):
    return rec[p[0]:p[1]].decode("ascii", "replace") if p else ""

def int_or_zero(s: str) -> int:
    try: return int(s.strip())
    except (ValueError, AttributeError): return 0

def axis_codes(rec: bytes, p, n: int, width: int):
    if not p or n <= 0: return []
    out=[]
    for i in range(min(n, 20)):
        off = p[0] + i*width
        out.append(norm(rec[off:off+4].decode("ascii", "replace") if width == 5 else rec[off+2:off+6].decode("ascii", "replace")))
    return out

def entity_k80_part(rec: bytes, p, n: int):
    if not p: return "entity_absent"
    for i in range(min(n, 20)):
        off = p[0] + i*7
        code = norm(rec[off+2:off+6].decode("ascii", "replace"))
        if code.startswith("K80"):
            q = rec[off:off+1].decode("ascii", "replace").strip()
            return "Part II" if q == "6" else "Part I" if q in {"1","2","3","4","5"} else "entity_unknown"
    return "entity_no_K80"

def k80_subtype(codes):
    """Pre-specified mutually exclusive priority, never file-order dependent."""
    found={f"K80.{c[3]}" for c in codes if len(c)>=4 and c.startswith("K80") and c[3] in "0123458"}
    for label in ("K80.3","K80.0","K80.4","K80.1","K80.5","K80.2","K80.8"):
        if label in found: return label
    return "K80_other_or_unspecified"

def k80_labels(codes):
    labels=sorted({f"K80.{c[3]}" for c in codes if len(c)>=4 and c.startswith("K80") and c[3] in "0123458"})
    return labels or ["K80_other_or_unspecified"]

def severity(subtype: str):
    if subtype == "K80.0": return "acute_cholecystitis"
    if subtype == "K80.3": return "cholangitis"
    if subtype in {"K80.1","K80.4"}: return "other_cholecystitis"
    if subtype in {"K80.2","K80.5"}: return "without_cholecystitis"
    return "other_or_unspecified"

def complexity(n: int):
    return "1" if n == 1 else "2" if n == 2 else "3-4" if 3 <= n <= 4 else "5-9" if 5 <= n <= 9 else "10+" if n >= 10 else "unknown"

def k80_code_count(codes):
    n=len(k80_labels(codes)); return "1" if n==1 else "2" if n==2 else "3+"

def bridged(row, rec):
    r5 = slc(rec, row["p"]["race_recode5"]).strip()
    if r5 in {"1","2","3","4"}: return {"1":"White","2":"Black","3":"AIAN","4":"API"}[r5]
    rd = slc(rec, row["p"]["race_detail"]).strip()
    if rd in {"01","1"}: return "White"
    if rd in {"02","2"}: return "Black"
    if rd in {"03","3"}: return "AIAN"
    return "API" if rd else "missing"

def race40_group(row, rec):
    """Collapse official Race Recode 40 (489–490) into frozen reporting groups."""
    code = slc(rec, row["p"]["race_recode40"]).strip()
    return RACE40_GROUP_MAP.get(code, "missing")

def single_race(row, rec):
    return race40_group(row, rec)

def load_crosswalk():
    ans={}
    with XW_PATH.open(encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f):
            r["year"] = int(r["year"]); r["p"] = {k: pos(r,k) for k in r if k not in {"year","p"}}
            ans[r["year"]]=r
    return ans

def load_dest_map():
    out={}
    with ONTOLOGY_PATH.open(encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f):
            key=int(str(r.get("recode_113", "")).strip())
            destination=(r.get("destination") or "").strip()
            if key in out: raise ValueError(f"duplicate 113-cause recode ontology key: {key}")
            out[key]=destination
    expected_keys=set(range(1,136))
    if set(out)!=expected_keys:
        missing=sorted(expected_keys-set(out)); extra=sorted(set(out)-expected_keys)
        raise ValueError(f"ontology must contain each 113-cause recode 1–135 exactly once; missing={missing}; extra={extra}")
    if set(out.values()) != EXPECTED_DESTINATIONS:
        raise ValueError(f"ontology must use exactly the frozen nine destination categories; got={sorted(set(out.values()))}")
    return out

def source_for(year: int):
    d=NCHS/"raw"/str(year); extracted=d/"extracted"
    files=[p for p in extracted.glob("*") if p.is_file()] if extracted.exists() else []
    if files: return ("file", max(files, key=lambda x:x.stat().st_size), "")
    zips=list(d.glob("*.zip"))
    if not zips: raise FileNotFoundError(f"no source for {year}")
    z=max(zips,key=lambda x:x.stat().st_size)
    with zipfile.ZipFile(z) as zz: member=max((i for i in zz.infolist() if not i.is_dir()),key=lambda i:i.file_size).filename
    return ("zip",z,member)

def make_year_result(year: int):
    xw=load_crosswalk(); destmap=load_dest_map(); row=xw[year]; mode, src, member=source_for(year)
    blank=lambda: {"A":0,"B":0,"B_official":0}
    agg={"year":year,"source_mode":mode,"source_path":str(src),"zip_member":member,"record_length_expected":REC_LEN[year],
         "total_records":0,"bad_length":0,"residents":0,"nonresidents":0,"main":blank(),"a_prime":0,
         "cells":defaultdict(blank),"bridged":defaultdict(blank),"single_race":defaultdict(blank),"bridged_cells":defaultdict(lambda: defaultdict(blank)),"single_race_cells":defaultdict(lambda: defaultdict(blank)),
         "groups":{g:defaultdict(blank) for g in ("subtype","severity","entity_part","complexity","k80_code_count","multilabel_subtype")},
         "dest_raw":Counter(),"dest_cells":defaultdict(Counter),"other_ucd":Counter(),"missing":Counter(),
         "field_available":{},"cross":{d:{"M":0,"U":0,"U_official":0,"cells":defaultdict(lambda:{"M":0,"U":0})} for d in DISEASES}}
    for field in ("age_detail","sex","race_recode5","race_recode6","race_recode40","race_detail","ucd_code","record_n","record_block","entity_n","entity_block"):
        agg["field_available"][field]=int(row["p"].get(field) is not None)
    def bump(d, a, u):
        if a: d["A"] += 1
        if u and a: d["B"] += 1
        if u: d["B_official"] += 1
    def process(fh):
        for lineno, line in enumerate(fh,1):
            rec=line.rstrip(b"\r\n"); agg["total_records"]+=1
            if len(rec)!=REC_LEN[year]: agg["bad_length"]+=1; continue
            resident=slc(rec,row["p"]["resident_status"]).strip() in {"1","2","3"}
            if not resident: agg["nonresidents"]+=1; continue
            agg["residents"]+=1
            rn=int_or_zero(slc(rec,row["p"]["record_n"])); en=int_or_zero(slc(rec,row["p"]["entity_n"]))
            codes=axis_codes(rec,row["p"]["record_block"],rn,5); roots={c[:3] for c in codes if c}
            ucd=norm(slc(rec,row["p"]["ucd_code"])); uroot=ucd[:3]
            ag=age_group(slc(rec,row["p"]["age_detail"])); sx=sex_group(slc(rec,row["p"]["sex"])); cell=f"{ag}|{sx}"
            a="K80" in roots; u=uroot=="K80"
            if a or u:
                bump(agg["main"],a,u); agg["a_prime"] += int(a or u)
                if ag == "unknown": agg["missing"]["age_group_unknown"] += 1
                if sx == "unknown": agg["missing"]["sex_unknown"] += 1
                if not ucd: agg["missing"]["ucd_blank"] += 1
                if a:
                    br=bridged(row,rec); sr=single_race(row,rec)
                    if year >= 2018 and sr == "missing":
                        agg["missing"]["race_recode40"] += 1
                    bump(agg["cells"][cell],True,u); bump(agg["bridged"][br],True,u); bump(agg["single_race"][sr],True,u)
                    if ag in AGE_ORDER and sx in SEX_ORDER:
                        bump(agg["bridged_cells"][br][cell],True,u); bump(agg["single_race_cells"][sr][cell],True,u)
                    st=k80_subtype(codes); part=entity_k80_part(rec,row["p"]["entity_block"],en)
                    for group,label in (("subtype",st),("severity",severity(st)),("entity_part",part),("complexity",complexity(rn)),("k80_code_count",k80_code_count(codes))): bump(agg["groups"][group][label],True,u)
                    for label in k80_labels(codes): bump(agg["groups"]["multilabel_subtype"][label],True,u)
                    if not u:
                        rec113=int_or_zero(slc(rec,row["p"]["ucd_113"]))
                        # Exact K81–K83 UCD prefixes take precedence over the
                        # 113-cause recode (including recode 96). Missing/invalid
                        # recodes remain in OTHER so the published partition is
                        # always closed over the frozen nine categories.
                        exact = "OTHER_GALLBLADDER_BILIARY" if uroot in {"K81","K82","K83"} else destmap.get(rec113,"OTHER")
                        agg["dest_raw"][exact]+=1; agg["dest_cells"][cell][exact]+=1
                        if exact=="OTHER": agg["other_ucd"][ucd or "<blank>"]+=1
            for disease, roots_needed in DISEASES.items():
                m=bool(roots & roots_needed); uu=uroot in roots_needed
                if m or uu:
                    d=agg["cross"][disease]
                    d["M"]+=int(m); d["U"]+=int(m and uu); d["U_official"]+=int(uu)
                    if ag in AGE_ORDER and sx in SEX_ORDER and m:
                        d["cells"][cell]["M"]+=1; d["cells"][cell]["U"]+=int(uu)
    if mode=="file":
        with src.open("rb") as f: process(f)
    else:
        try:
            with zipfile.ZipFile(src) as z:
                with z.open(member) as f: process(f)
        except NotImplementedError:
            if not SEVEN_ZIP.exists(): raise RuntimeError(f"ZIP compression unsupported by Python and 7-Zip unavailable: {src}")
            # 7-Zip writes the member bytes to stdout; no archive member is extracted to disk.
            proc=subprocess.Popen([str(SEVEN_ZIP),"e","-so",str(src),member],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            assert proc.stdout is not None
            process(proc.stdout)
            stderr=proc.stderr.read().decode("utf-8","replace") if proc.stderr else ""
            if proc.wait()!=0: raise RuntimeError(f"7-Zip stream failure for {src}::{member}: {stderr[-1000:]}")
    # convert defaultdicts to ordinary dicts for stable IPC/JSON
    agg["cells"]=dict(agg["cells"]); agg["bridged"]=dict(agg["bridged"]); agg["single_race"]=dict(agg["single_race"]); agg["bridged_cells"]={k:dict(v) for k,v in agg["bridged_cells"].items()}; agg["single_race_cells"]={k:dict(v) for k,v in agg["single_race_cells"].items()}
    agg["groups"]={k:dict(v) for k,v in agg["groups"].items()}; agg["dest_raw"]=dict(agg["dest_raw"]); agg["dest_cells"]={k:dict(v) for k,v in agg["dest_cells"].items()}; agg["other_ucd"]=dict(agg["other_ucd"])
    for d in agg["cross"].values(): d["cells"]=dict(d["cells"])
    return agg

def wilson(k,n,z=1.959963984540054):
    if n<=0:return (float("nan"),float("nan"),float("nan"))
    p=k/n; den=1+z*z/n; c=p+z*z/(2*n); w=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n)); return p,(c-w)/den,(c+w)/den

def direct(year_agg, weights, numerator="B", denominator="A"):
    usable=[]
    for cell,w in weights.items():
        q=year_agg.get(cell,{})
        n=q.get(denominator,0); k=q.get(numerator,0)
        if n>0 and w>0: usable.append((w,k,n,cell))
    support=sum(x[0] for x in usable)
    if support<=0:return {"p":float("nan"),"lo":float("nan"),"hi":float("nan"),"var":float("nan"),"support":0.0,"cells":0}
    p=sum((w/support)*(k/n) for w,k,n,_ in usable); var=sum((w/support)**2*(k/n)*(1-k/n)/n for w,k,n,_ in usable)
    se=math.sqrt(max(var,0)); return {"p":p,"lo":max(0,p-1.96*se),"hi":min(1,p+1.96*se),"var":var,"support":support,"cells":len(usable)}

def weights_from(annual, year_filter, item_key="cells", denominator="A", support_years=YEARS):
    c=Counter()
    for y in year_filter:
        for cell,q in annual[y][item_key].items():
            if cell.split("|")[0] in AGE_ORDER and cell.split("|")[1] in SEX_ORDER: c[cell]+=q.get(denominator,0)
    # Fixed reference distribution is restricted before estimation to cells with
    # positive denominator in every target year; no annual renormalization.
    common={cell for cell in c if all(annual[y][item_key].get(cell,{}).get(denominator,0)>0 for y in support_years)}
    c=Counter({cell:n for cell,n in c.items() if cell in common})
    total=sum(c.values()); return {k:v/total for k,v in c.items()} if total else {}

def safe_float(x): return "" if not math.isfinite(x) else x
def write_csv(path:Path, rows, fields=None):
    path.parent.mkdir(parents=True,exist_ok=True)
    fields=fields or (list(rows[0].keys()) if rows else [])
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore");w.writeheader();w.writerows(rows)

def annual_table(annual):
    rows=[]
    for y in YEARS:
        q=annual[y]["main"]; p,lo,hi=wilson(q["B"],q["A"])
        rows.append({"year":y,"A_record_axis_K80":q["A"],"B_main_A_and_UCD_K80":q["B"],"B_official_all_UCD_K80":q["B_official"],"A_prime_A_union_B_official":annual[y]["a_prime"],"gap_A_minus_B":q["A"]-q["B"],"UCF":p,"UCF_CI_lo":lo,"UCF_CI_hi":hi,"orphan_UCD_K80":q["B_official"]-q["B"]})
    return rows

def apply_prespecified_age_rebin(annual):
    """Merge 0-14/15-24 into 0-24 in all age-sex aggregation modules."""
    def rekey(cell):
        age,sex=cell.split("|"); return ("0-24" if age in {"0-14","15-24"} else age)+"|"+sex
    def merge_stats(target, source):
        for cell,q in source.items():
            key=rekey(cell); dest=target.setdefault(key,{})
            for field,val in q.items(): dest[field]=dest.get(field,0)+val
    for q in annual.values():
        merged={};merge_stats(merged,q["cells"]);q["cells"]=merged
        for key in ("bridged_cells","single_race_cells"):
            outer={}
            for race,cells in q[key].items():
                outer[race]={};merge_stats(outer[race],cells)
            q[key]=outer
        for disease,d in q["cross"].items():
            merged={};merge_stats(merged,d["cells"]);d["cells"]=merged
        dc={}
        for cell,counts in q["dest_cells"].items():
            key=rekey(cell); dest=dc.setdefault(key,{})
            for cat,n in counts.items(): dest[cat]=dest.get(cat,0)+n
        q["dest_cells"]=dc

def standardized_tables(annual):
    schemes={"primary_2018_2024":weights_from(annual,range(2018,2025)),"sensitivity_1999":weights_from(annual,[1999]),"sensitivity_full_period":weights_from(annual,YEARS)}
    weightrows=[]; rows=[]; details={}
    for scheme,w in schemes.items():
        for cell,value in w.items(): weightrows.append({"scheme":scheme,"age_sex_cell":cell,"weight":value})
        out={}
        for y in YEARS:
            z=direct(annual[y]["cells"],w);out[y]=z
            rows.append({"scheme":scheme,"year":y,"std_UCF":safe_float(z["p"]),"CI_lo":safe_float(z["lo"]),"CI_hi":safe_float(z["hi"]),"variance":safe_float(z["var"]),"effective_weight_sum":z["support"],"supported_cells":z["cells"]})
        details[scheme]=out
    contrasts=[]
    for scheme,out in details.items():
        base=out[1999]; end=out[2024]
        rd=end["p"]-base["p"]; se=math.sqrt(end["var"]+base["var"])
        rr=end["p"]/base["p"] if base["p"]>0 else float("nan")
        lrr=math.sqrt(end["var"]/end["p"]**2+base["var"]/base["p"]**2) if base["p"]>0 and end["p"]>0 else float("nan")
        for y in (1999,2015,2024):
            z=out[y]; contrasts.append({"scheme":scheme,"contrast_or_year":str(y),"std_UCF":safe_float(z["p"]),"CI_lo":safe_float(z["lo"]),"CI_hi":safe_float(z["hi"]),"RD": "","RD_CI_lo":"","RD_CI_hi":"","RR":"","RR_CI_lo":"","RR_CI_hi":""})
        contrasts.append({"scheme":scheme,"contrast_or_year":"2024_minus_1999","std_UCF":"","CI_lo":"","CI_hi":"","RD":safe_float(rd),"RD_CI_lo":safe_float(rd-1.96*se),"RD_CI_hi":safe_float(rd+1.96*se),"RR":safe_float(rr),"RR_CI_lo":safe_float(rr*math.exp(-1.96*lrr)) if math.isfinite(lrr) else "","RR_CI_hi":safe_float(rr*math.exp(1.96*lrr)) if math.isfinite(lrr) else ""})
    return weightrows,rows,contrasts,details

def interaction_diagnostics(annual):
    """Aggregated binomial LR/deviance diagnostics for the requested interactions.

    The saturated three-way cell estimator is algebraically the empirical
    cell-proportion marginal standardization reported in the companion table.
    """
    try: import numpy as np
    except Exception as exc: return [{"model":"unavailable","status":"dependency_unavailable","detail":repr(exc),"n_cells":"","df_model":"","deviance":"","lr_vs_main":"","df_diff":"","p_value":""}]
    data=[]
    for y in YEARS:
        for cell,q in annual[y]["cells"].items():
            ag,sx=cell.split("|")
            if ag in AGE_ORDER and sx in SEX_ORDER and q["A"]>0:data.append((str(y),ag,sx,q["A"],q["B"]))
    levels=[sorted({r[i] for r in data}) for i in range(3)]
    def one(v, lev): return [int(v==x) for x in lev[1:]]
    def matrix(kind):
        rows=[]
        for y,a,s,_,_ in data:
            yy,aa,ss=one(y,levels[0]),one(a,levels[1]),one(s,levels[2]); x=[1]+yy+aa+ss
            if kind in {"year_x_age","all_pairwise","saturated_year_x_age_x_sex"}: x += [u*v for u in yy for v in aa]
            if kind in {"year_x_sex","all_pairwise","saturated_year_x_age_x_sex"}: x += [u*v for u in yy for v in ss]
            if kind in {"age_x_sex","all_pairwise","saturated_year_x_age_x_sex"}: x += [u*v for u in aa for v in ss]
            if kind=="saturated_year_x_age_x_sex": x += [u*v*w for u in yy for v in aa for w in ss]
            rows.append(x)
        return np.asarray(rows,dtype=float)
    n=np.asarray([r[3] for r in data],float);k=np.asarray([r[4] for r in data],float); prop=k/n
    def gammaincc(a,x):
        """Numerical Recipes regularized upper incomplete gamma, no SciPy dependency."""
        if x<=0:return 1.0
        gln=math.lgamma(a)
        if x<a+1:
            ap=a; summ=1/a; delta=summ
            for _ in range(10000):
                ap+=1;delta*=x/ap;summ+=delta
                if abs(delta)<abs(summ)*1e-14:break
            return max(0.0,min(1.0,1-summ*math.exp(-x+a*math.log(x)-gln)))
        b=x+1-a; c=1e300; d=1/b; h=d
        for i in range(1,10000):
            an=-i*(i-a); b+=2; d=an*d+b
            if abs(d)<1e-300:d=1e-300
            c=b+an/c
            if abs(c)<1e-300:c=1e-300
            d=1/d; delta=d*c; h*=delta
            if abs(delta-1)<1e-14:break
        return max(0.0,min(1.0,math.exp(-x+a*math.log(x)-gln)*h))
    def fit(X):
        eta=np.full(len(n),math.log((k.sum()+0.5)/(n.sum()-k.sum()+0.5)))
        for _ in range(100):
            mu=1/(1+np.exp(-np.clip(eta,-30,30))); w=np.maximum(n*mu*(1-mu),1e-9); z=eta+(prop-mu)/np.maximum(mu*(1-mu),1e-9)
            beta=np.linalg.lstsq(X*np.sqrt(w)[:,None],z*np.sqrt(w),rcond=1e-10)[0]; new=X@beta
            if np.max(np.abs(new-eta))<1e-8: eta=new;break
            eta=new
        mu=np.clip(1/(1+np.exp(-np.clip(eta,-30,30))),1e-12,1-1e-12); term=np.zeros(len(n)); mk=k>0; mn=(n-k)>0
        term[mk]+=k[mk]*np.log(k[mk]/(n[mk]*mu[mk])); term[mn]+=(n[mn]-k[mn])*np.log((n[mn]-k[mn])/(n[mn]*(1-mu[mn])))
        return float(2*term.sum()),int(np.linalg.matrix_rank(X))
    kinds=["main_effects","year_x_age","year_x_sex","age_x_sex","all_pairwise","saturated_year_x_age_x_sex"]
    fitted={};out=[]
    for kind in kinds:
        try:
            dev,rank=fit(matrix(kind));fitted[kind]=(dev,rank);out.append({"model":kind,"status":"fit_numpy_irls","detail":"aggregated binomial IRLS; LR p uses regularized upper incomplete gamma","n_cells":len(data),"df_model":rank-1,"df_residual":len(data)-rank,"deviance":dev,"deviance_per_df":dev/(len(data)-rank) if len(data)>rank else "","lr_vs_main":"","df_diff":"","p_value":""})
        except Exception as exc:out.append({"model":kind,"status":"fit_failed","detail":repr(exc),"n_cells":len(data),"df_model":"","df_residual":"","deviance":"","deviance_per_df":"","lr_vs_main":"","df_diff":"","p_value":""})
    if "main_effects" in fitted:
        d0,r0=fitted["main_effects"]
        for r in out:
            if r["model"] in fitted and r["model"]!="main_effects":
                d,rank=fitted[r["model"]]; lr=max(0,d0-d);dfd=rank-r0;r["lr_vs_main"]=lr;r["df_diff"]=dfd;r["p_value"]=gammaincc(dfd/2,lr/2) if dfd>0 else 1.0
    out.append({"model":"empirical_marginal_equivalence","status":"by_construction","detail":"Fixed-weight sum of observed year×age×sex binomial proportions equals saturated three-way cell marginal estimate; annual values in interaction_empirical_marginal_sensitivity.csv.","n_cells":len(data),"df_model":"","df_residual":"","deviance":"","deviance_per_df":"","lr_vs_main":"","df_diff":"","p_value":""})
    return out

def decomposition(annual, group):
    q0=annual[1999]["groups"][group];q1=annual[2024]["groups"][group]; keys=sorted(set(q0)|set(q1)); a0=sum(v["A"] for v in q0.values());a1=sum(v["A"] for v in q1.values()); rows=[]
    composition=selection=0.0
    for k in keys:
        x0=q0.get(k,{"A":0,"B":0});x1=q1.get(k,{"A":0,"B":0});w0=x0["A"]/a0 if a0 else 0;w1=x1["A"]/a1 if a1 else 0;p0=x0["B"]/x0["A"] if x0["A"] else 0;p1=x1["B"]/x1["A"] if x1["A"] else 0
        comp=.5*(p0+p1)*(w1-w0); sel=.5*(w0+w1)*(p1-p0);composition+=comp;selection+=sel
        for year,x,w,p in ((1999,x0,w0,p0),(2024,x1,w1,p1)):
            pp,lo,hi=wilson(x["B"],x["A"]);rows.append({"group":group,"stratum":k,"year":year,"A":x["A"],"B":x["B"],"UCF":safe_float(pp),"UCF_CI_lo":safe_float(lo),"UCF_CI_hi":safe_float(hi),"composition_weight":w,"kitagawa_composition_component":comp,"kitagawa_selection_component":sel})
    total=(annual[2024]["main"]["B"]/annual[2024]["main"]["A"])-(annual[1999]["main"]["B"]/annual[1999]["main"]["A"])
    summary={"group":group,"UCF_change_2024_minus_1999":total,"composition_component":composition,"selection_component":selection,"decomposition_residual":total-composition-selection}
    return rows,summary

def destination_tables(annual):
    cats=sorted(set().union(*(set(annual[y]["dest_raw"]) for y in YEARS)))
    weights=weights_from(annual,range(2018,2025),item_key="dest_cells",denominator=None) if False else None
    wc=Counter()
    for y in range(2018,2025):
        for cell,d in annual[y]["dest_cells"].items():
            if cell.split("|")[0] in AGE_ORDER and cell.split("|")[1] in SEX_ORDER: wc[cell]+=sum(d.values())
    common={cell for cell in wc if all(sum(annual[y]["dest_cells"].get(cell,{}).values())>0 for y in YEARS)}
    wc=Counter({cell:n for cell,n in wc.items() if cell in common})
    total=sum(wc.values());weights={k:v/total for k,v in wc.items()} if total else {}
    rows=[]; other=[]
    for y in YEARS:
        den=sum(annual[y]["dest_raw"].values()); cellagg={c:{"A":sum(d.values()),"B":d.get(c,0)} for c in cats for _,d in annual[y]["dest_cells"].items()}
        # direct needs per-cell denominator; build for each category
        for c in cats:
            per={cell:{"A":sum(d.values()),"B":d.get(c,0)} for cell,d in annual[y]["dest_cells"].items()}
            z=direct(per,weights); n=annual[y]["dest_raw"].get(c,0);p,lo,hi=wilson(n,den)
            rows.append({"year":y,"destination":c,"gap_denominator_A_minus_B":den,"raw_n":n,"raw_probability":safe_float(p),"raw_CI_lo":safe_float(lo),"raw_CI_hi":safe_float(hi),"std_probability":safe_float(z["p"]),"std_CI_lo":safe_float(z["lo"]),"std_CI_hi":safe_float(z["hi"]),"std_variance":safe_float(z["var"]),"effective_weight_sum":z["support"]})
        for code,n in sorted(annual[y]["other_ucd"].items()): other.append({"year":y,"destination":"OTHER","ucd":code,"raw_n":n,"gap_denominator":den,"raw_probability":n/den if den else ""})
    return rows,other,[{"scheme":"destination_primary_2018_2024_gap_age_sex","age_sex_cell":k,"weight":v} for k,v in weights.items()]

def destination_contrasts(destrows):
    out=[]
    for cat in sorted({r["destination"] for r in destrows}):
        r0=next(r for r in destrows if r["destination"]==cat and r["year"]==1999);r1=next(r for r in destrows if r["destination"]==cat and r["year"]==2024)
        p0=float(r0["std_probability"]);p1=float(r1["std_probability"]);v0=float(r0["std_variance"]);v1=float(r1["std_variance"]);rd=p1-p0;se=math.sqrt(v0+v1);rr=p1/p0 if p0>0 else float("nan"); sl=math.sqrt(v0/p0**2+v1/p1**2) if p0>0 and p1>0 else float("nan")
        out.append({"destination":cat,"std_probability_1999":p0,"CI_lo_1999":r0["std_CI_lo"],"CI_hi_1999":r0["std_CI_hi"],"std_probability_2024":p1,"CI_lo_2024":r1["std_CI_lo"],"CI_hi_2024":r1["std_CI_hi"],"RD_2024_minus_1999":rd,"RD_CI_lo":rd-1.96*se,"RD_CI_hi":rd+1.96*se,"RR_2024_over_1999":safe_float(rr),"RR_CI_lo":safe_float(rr*math.exp(-1.96*sl)) if math.isfinite(sl) else "","RR_CI_hi":safe_float(rr*math.exp(1.96*sl)) if math.isfinite(sl) else "","inference_note":"Fixed common age-sex weights; independent-year stratified-binomial approximation. Destination categories are compositional and correlated across categories; each single-category contrast is valid as reported, but no between-category independence is assumed."})
    return out

def race_tables(annual):
    rows=[]
    for regime, window, key in (("bridged",range(1999,2021),"bridged"),("single_race",range(2018,2025),"single_race")):
        categories=list(RACE40_GROUPS) if regime == "single_race" else sorted(set().union(*(set(annual[y][key]) for y in window)))
        # One common fixed age-sex standard population per regime, across all races.
        wc=Counter()
        for y in window:
            for cat in categories:
                for cell,q in annual[y][key+"_cells"].get(cat,{}).items(): wc[cell]+=q["A"]
        denom=sum(wc.values()); w={c:n/denom for c,n in wc.items()} if denom else {}
        for cat in categories:
            for y in window:
                q=annual[y][key].get(cat,{"A":0,"B":0});p,lo,hi=wilson(q["B"],q["A"])
                z=direct(annual[y][key+"_cells"].get(cat,{}),w)
                estimable=abs(z["support"]-1.0)<1e-12
                rows.append({"regime":regime,"year":y,"race":cat,"A":q["A"],"B":q["B"],"UCF":safe_float(p),"CI_lo":safe_float(lo),"CI_hi":safe_float(hi),"std_UCF":safe_float(z["p"]) if estimable else "","std_CI_lo":safe_float(z["lo"]) if estimable else "","std_CI_hi":safe_float(z["hi"]) if estimable else "","effective_weight_sum":z["support"],"standardization_status":"common pooled age-sex distribution among record-axis coded deaths within regime" if estimable else "not_estimable: common regime weight has zero-denominator cell; crude result retained"})
    return rows

def cross_disease(annual):
    rows=[]; standardized=[]; contrasts=[]; weights=[]
    for disease in DISEASES:
        wc=Counter()
        for y in range(2018,2025):
            for cell,q in annual[y]["cross"][disease]["cells"].items(): wc[cell]+=q["M"]
        common={cell for cell in wc if all(annual[y]["cross"][disease]["cells"].get(cell,{}).get("M",0)>0 for y in YEARS)}
        wc=Counter({cell:n for cell,n in wc.items() if cell in common})
        den=sum(wc.values()); w={c:n/den for c,n in wc.items()} if den else {}
        for cell,value in w.items(): weights.append({"disease":disease,"scheme":"own_disease_2018_2024_age_sex_mention_distribution","age_sex_cell":cell,"weight":value})
        zby={}
        for y in YEARS:
            q=annual[y]["cross"][disease];p,lo,hi=wilson(q["U"],q["M"])
            rows.append({"disease":disease,"year":y,"M_record_axis_mention":q["M"],"U_main_M_and_UCD":q["U"],"U_official_all_UCD":q["U_official"],"UCF":safe_float(p),"CI_lo":safe_float(lo),"CI_hi":safe_float(hi),"orphan_UCD":q["U_official"]-q["U"]})
            z=direct(q["cells"],w,numerator="U",denominator="M"); zby[y]=z
            standardized.append({"disease":disease,"scheme":"own_disease_2018_2024_age_sex_mention_distribution","year":y,"std_UCF":safe_float(z["p"]),"CI_lo":safe_float(z["lo"]),"CI_hi":safe_float(z["hi"]),"variance":safe_float(z["var"]),"effective_weight_sum":z["support"],"supported_cells":z["cells"]})
        for y in (1999,2015,2024):
            z=zby[y];contrasts.append({"disease":disease,"contrast_or_year":str(y),"std_UCF":safe_float(z["p"]),"CI_lo":safe_float(z["lo"]),"CI_hi":safe_float(z["hi"]),"RD":"","RD_CI_lo":"","RD_CI_hi":"","RR":"","RR_CI_lo":"","RR_CI_hi":""})
        a,b=zby[1999],zby[2024];rd=b["p"]-a["p"];se=math.sqrt(a["var"]+b["var"]);rr=b["p"]/a["p"] if a["p"]>0 else float("nan");lrr=math.sqrt(a["var"]/a["p"]**2+b["var"]/b["p"]**2) if a["p"]>0 and b["p"]>0 else float("nan")
        contrasts.append({"disease":disease,"contrast_or_year":"2024_minus_1999","std_UCF":"","CI_lo":"","CI_hi":"","RD":safe_float(rd),"RD_CI_lo":safe_float(rd-1.96*se),"RD_CI_hi":safe_float(rd+1.96*se),"RR":safe_float(rr),"RR_CI_lo":safe_float(rr*math.exp(-1.96*lrr)) if math.isfinite(lrr) else "","RR_CI_hi":safe_float(rr*math.exp(1.96*lrr)) if math.isfinite(lrr) else ""})
    return rows,standardized,contrasts,weights

def issues(annual, std_details):
    main1999=annual[1999]["main"];main2024=annual[2024]["main"]
    return [
      {"issue_id":"R1_standardization","requirement":"K80 age-sex direct standardization and sensitivity weights","status":"resolved","evidence":"tables/k80_standardized_annual.csv; tables/k80_standardized_contrasts.csv"},
      {"issue_id":"R2_interactions","requirement":"year×age/year×sex/age×sex sensitivity","status":"resolved","evidence":"tables/interaction_model_diagnostics.csv; tables/interaction_empirical_marginal_sensitivity.csv"},
      {"issue_id":"R3_composition","requirement":"subtype/severity/Part I-II/complexity composition and decomposition","status":"resolved","evidence":"tables/k80_composition_ucf.csv; tables/k80_kitagawa_decomposition.csv"},
      {"issue_id":"R4_destination","requirement":"gap denominator, raw and standardized destination probabilities, OTHER detail","status":"resolved","evidence":"tables/ucd_destination_annual.csv; tables/ucd_destination_other_ucd.csv"},
      {"issue_id":"R5_covid_medcoder","requirement":"segment only if supported","status":"descriptive_only","evidence":"tables/era_break_evidence.csv (COVID period has two annual observations; no stable segmented trend)"},
      {"issue_id":"R6_cross_disease","requirement":"preselected cross-disease benchmark","status":"resolved","evidence":"tables/cross_disease_ucf_annual.csv; tables/cross_disease_standardized_annual.csv; tables/cross_disease_standardized_contrasts.csv"},
      {"issue_id":"R7_axis_inconsistency","requirement":"record-axis/UCD reconciliation count and A-prime sensitivity","status":"resolved","evidence":"tables/k80_annual_main.csv; no record-level examples are released"},
      {"issue_id":"R8_race_missingness","requirement":"Race Recode 40, annual flow, 2021 audit, missingness","status":"resolved","evidence":"tables/race_ucf_annual.csv; tables/missingness_and_field_audit.csv; evidence/NCHS_coding_and_field_evidence.md"},
    ]

def tests(annual, stdweights, stdrows, destrows, crossrows, crossstd, racestd):
    results=[]
    def add(name, passed, detail):results.append({"test":name,"passed":int(bool(passed)),"detail":detail})
    add("main_B_subset_A",all(v["main"]["B"]<=v["main"]["A"] for v in annual.values()),"main B is defined as A∩B_official")
    add("gap_identity",all(v["main"]["A"]-v["main"]["B"]>=0 for v in annual.values()),"annual gap nonnegative")
    for scheme in {r["scheme"] for r in stdweights}:
        s=sum(r["weight"] for r in stdweights if r["scheme"]==scheme);add(f"weight_sum_{scheme}",abs(s-1)<1e-12,f"sum={s}")
    direct_support=[r for r in stdrows+destrows+crossstd+racestd if "effective_weight_sum" in r and r["effective_weight_sum"]!=""]
    for r in direct_support:
        key="_".join(str(r.get(k,"")) for k in ("scheme","regime","race","disease","destination","year") if r.get(k,"")!="")
        flagged=str(r.get("standardization_status","")).startswith("not_estimable")
        add(f"direct_effective_weight_sum_{key}",abs(float(r["effective_weight_sum"])-1.0)<1e-12 or flagged,f"effective_weight_sum={r['effective_weight_sum']}; {'explicitly not estimable' if flagged else 'must equal 1 for an estimable fixed-weight result'}")
    for y in YEARS:
        s=sum(float(r["std_probability"]) for r in destrows if r["year"]==y and r["std_probability"]!=""); add(f"destination_probability_sum_{y}",abs(s-1)<1e-8,f"sum={s}")
        raw=[r for r in destrows if r["year"]==y]; add(f"destination_mutually_exclusive_{y}",sum(int(r["raw_n"]) for r in raw)==int(raw[0]["gap_denominator_A_minus_B"]) if raw else False,f"raw category total={sum(int(r['raw_n']) for r in raw) if raw else 0}")
        add(f"destination_ontology_{y}",set(r["destination"] for r in raw)==EXPECTED_DESTINATIONS,"exactly the nine frozen observable destination categories are represented (zero counts permitted)")
    finite=True
    for r in stdrows+destrows+crossrows+crossstd:
        for k in ("CI_lo","CI_hi","std_CI_lo","std_CI_hi"):
            if k in r and r[k] != "": finite &= math.isfinite(float(r[k])) and float(r[k])<=float(r[k.replace("lo","hi")]) if "lo" in k and k.replace("lo","hi") in r else math.isfinite(float(r[k]))
    add("ci_finite_ordered",finite,"all reported probability CIs finite and ordered")
    add("cross_B_subset_A",all(r["U_main_M_and_UCD"]<=r["M_record_axis_mention"] for r in crossrows),"each disease U≤M")
    expected_race40={"01":"White","02":"Black","03":"AIAN",**{f"{i:02d}":"Asian" for i in range(4,11)},**{f"{i:02d}":"NHOPI" for i in range(11,15)},**{f"{i:02d}":"Multiple" for i in range(15,41)}}
    add("race40_field_present_2018_2024",all(annual[y]["field_available"].get("race_recode40")==1 for y in range(2018,2025)),"official Race Recode 40 at positions 489-490 is present in every single-race year")
    add("race40_mapping_stable",RACE40_GROUP_MAP==expected_race40,"01 White; 02 Black; 03 AIAN; 04-10 Asian; 11-14 NHOPI; 15-40 Multiple")
    race_rows=[r for r in racestd if r["regime"]=="single_race" and 2018 <= int(r["year"]) <= 2024]
    race_close=True
    for y in range(2018,2025):
        rr=[r for r in race_rows if int(r["year"])==y]
        race_close &= set(r["race"] for r in rr)==set(RACE40_GROUPS)
        race_close &= sum(int(r["A"]) for r in rr)==int(annual[y]["main"]["A"])
        race_close &= sum(int(r["B"]) for r in rr)==int(annual[y]["main"]["B"])
    add("race40_yearly_category_closure",race_close,"six Race Recode 40 reporting groups close to annual K80 A and B in 2018-2024, including 2021")
    add("race40_no_missing_among_A",all(annual[y]["missing"].get("race_recode40",0)==0 for y in range(2018,2025)),"Race Recode 40 is nonblank for every A record in the single-race regime")
    return results

def manifest():
    files=[]
    for folder in (OUT/"analysis",OUT/"tables",OUT/"logs",OUT/"reports",OUT/"scripts"):
        for p in sorted(folder.rglob("*")):
            if p.is_file() and p != OUT/"manifests"/"analysis_manifest.json":
                h=hashlib.sha256(p.read_bytes()).hexdigest();files.append({"path":str(p.relative_to(OUT)).replace("\\","/"),"bytes":p.stat().st_size,"sha256":h})
    obj={"generated_utc":datetime.now(timezone.utc).isoformat(),"manifest_scope":"assigned analysis outputs; manifest file itself excluded to avoid self-reference","files":files}
    (OUT/"manifests"/"analysis_manifest.json").write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding="utf-8")

def main():
    global OUT, NCHS, XW_PATH, ONTOLOGY_PATH, SEVEN_ZIP
    ap=argparse.ArgumentParser();ap.add_argument("--workers",type=int,default=8);ap.add_argument("--postprocess",action="store_true",help="rebuild outputs from a locally generated annual aggregate JSON without rescanning inputs");ap.add_argument("--output-root",default="results",help="directory for generated aggregate outputs");ap.add_argument("--nchs-data-root",default="data/raw",help="directory containing officially obtained NCHS annual public-use files");ap.add_argument("--schema-crosswalk",default="config/schema_crosswalk.csv",help="official year-specific layout crosswalk supplied by the reproducing researcher");ap.add_argument("--ontology",default="data/derived/ucd_destination_ontology.csv",help="public nine-category destination ontology");ap.add_argument("--seven-zip",default="7z",help="optional 7-Zip executable for unsupported archive compression");args=ap.parse_args();OUT=Path(args.output_root);NCHS=Path(args.nchs_data_root);XW_PATH=Path(args.schema_crosswalk);ONTOLOGY_PATH=Path(args.ontology);SEVEN_ZIP=Path(args.seven_zip);workers=max(1,min(8,args.workers,len(YEARS)))
    t=time.time()
    if args.postprocess:
        annual={int(k):v for k,v in json.loads((OUT/"analysis"/"annual_aggregates.json").read_text(encoding="utf-8")).items()}
    else:
        ctx=mp.get_context("spawn")
        with ctx.Pool(workers) as pool: results=pool.map(make_year_result,YEARS)
        annual={r["year"]:r for r in results}
        (OUT/"analysis"/"annual_aggregates.json").write_text(json.dumps(annual,ensure_ascii=False,indent=2),encoding="utf-8")
    apply_prespecified_age_rebin(annual)
    write_csv(OUT/"tables"/"k80_annual_main.csv",annual_table(annual))
    weights,stdrows,contrasts,details=standardized_tables(annual);write_csv(OUT/"tables"/"k80_standardization_weights.csv",weights);write_csv(OUT/"tables"/"k80_standardized_annual.csv",stdrows);write_csv(OUT/"tables"/"k80_standardized_contrasts.csv",contrasts)
    # Empirical marginal values are the saturated year×age×sex cell estimand; fitted interaction diagnostics are reported separately.
    introws=[dict(r,method="saturated_empirical_year_age_sex_marginal") for r in stdrows];write_csv(OUT/"tables"/"interaction_empirical_marginal_sensitivity.csv",introws)
    write_csv(OUT/"tables"/"interaction_model_diagnostics.csv",interaction_diagnostics(annual))
    comp=[];decomp=[]
    for g in ("subtype","severity","entity_part","complexity"):
        a,b=decomposition(annual,g);comp+=a;decomp.append(b)
    write_csv(OUT/"tables"/"k80_composition_ucf.csv",comp);write_csv(OUT/"tables"/"k80_kitagawa_decomposition.csv",decomp)
    multi=[]
    for y in YEARS:
        for g in ("k80_code_count","multilabel_subtype"):
            for label,q in annual[y]["groups"][g].items():
                p,lo,hi=wilson(q["B"],q["A"]);multi.append({"analysis":"multi_label_sensitivity" if g=="multilabel_subtype" else "k80_code_count","year":y,"stratum":label,"A_records":q["A"],"B_records":q["B"],"UCF":safe_float(p),"CI_lo":safe_float(lo),"CI_hi":safe_float(hi),"note":"multilabel subtype rows are non-mutually-exclusive; do not sum" if g=="multilabel_subtype" else "mutually exclusive count of distinct valid K80 four-character codes"})
    write_csv(OUT/"tables"/"k80_multilabel_sensitivity.csv",multi)
    dest,other,dw=destination_tables(annual);write_csv(OUT/"tables"/"ucd_destination_annual.csv",dest);write_csv(OUT/"tables"/"ucd_destination_contrasts.csv",destination_contrasts(dest));write_csv(OUT/"tables"/"ucd_destination_other_ucd.csv",other);write_csv(OUT/"tables"/"ucd_destination_standardization_weights.csv",dw)
    race=race_tables(annual);write_csv(OUT/"tables"/"race_ucf_annual.csv",race);cross,crossstd,crosscon,crossw=cross_disease(annual);write_csv(OUT/"tables"/"cross_disease_ucf_annual.csv",cross);write_csv(OUT/"tables"/"cross_disease_standardized_annual.csv",crossstd);write_csv(OUT/"tables"/"cross_disease_standardized_contrasts.csv",crosscon);write_csv(OUT/"tables"/"cross_disease_standardization_weights.csv",crossw)
    flow=[];missing=[]
    for y in YEARS:
        q=annual[y];flow.append({"year":y,"source_mode":q["source_mode"],"total_records":q["total_records"],"bad_length":q["bad_length"],"resident_records":q["residents"],"A":q["main"]["A"],"B":q["main"]["B"],"gap":q["main"]["A"]-q["main"]["B"]})
        for field,avail in q["field_available"].items():
            note = ""
            if field == "race_recode40" and 2018 <= y <= 2024:
                note = "official Race Recode 40 at positions 489-490; single-race descriptive regime"
            elif field == "race_recode5" and y <= 2020:
                note = "bridged Race Recode 5 at position 450; bridged descriptive regime"
            elif field == "race_recode6" and y >= 2022:
                note = "Race Recode 6 at position 450; not used for the 2018-2024 Race Recode 40 module"
            missing.append({"year":y,"field":field,"available_in_layout":avail,"missing_among_A":q["missing"].get(field,0),"note":note})
    write_csv(OUT/"tables"/"cohort_annual_flow.csv",flow);write_csv(OUT/"tables"/"missingness_and_field_audit.csv",missing)
    eras=[]
    for name, yrs, reason in (("COVID_2020_2021",[2020,2021],"two annual points; segmented trend not statistically stable"),("MedCoder_2022_2024",[2022,2023,2024],"three annual points only; descriptive slope shown, no causal interruption claim")):
        vals=[annual[y]["main"]["B"]/annual[y]["main"]["A"] for y in yrs]; eras.append({"period":name,"years":"-".join(map(str,yrs)),"n_annual_points":len(yrs),"mean_UCF":statistics.mean(vals),"min_UCF":min(vals),"max_UCF":max(vals),"status":"descriptive_only","reason":reason})
    write_csv(OUT/"tables"/"era_break_evidence.csv",eras)
    primary=[r for r in stdrows if r["scheme"]=="primary_2018_2024"]; peak=max(primary,key=lambda r:float(r["std_UCF"]))
    write_csv(OUT/"tables"/"k80_2015_context.csv",[{"year":r["year"],"std_UCF":r["std_UCF"],"CI_lo":r["CI_lo"],"CI_hi":r["CI_hi"],"is_2015_pre_requested_reference":int(r["year"]==2015),"is_descriptive_maximum_over_annual_series":int(r["year"]==peak["year"]),"interpretation":"2015 was requested in the frozen contrast set; annual maximum is descriptive only and no post-hoc peak hypothesis test is performed."} for r in primary if r["year"] in (1999,2015,2024,peak["year"])])
    write_csv(OUT/"reports"/"analysis_issue_resolution.csv",issues(annual,details))
    testrows=tests(annual,weights,stdrows,dest,cross,crossstd,race);write_csv(OUT/"logs"/"verification_results.csv",testrows);passed=all(r["passed"] for r in testrows)
    report=["# Statistical reanalysis report","",f"Core definition totals: A={sum(annual[y]['main']['A'] for y in YEARS)}, B={sum(annual[y]['main']['B'] for y in YEARS)}, gap={sum(annual[y]['main']['A']-annual[y]['main']['B'] for y in YEARS)}.",f"The official-UCD reconciliation count was {sum(annual[y]['main']['B_official']-annual[y]['main']['B'] for y in YEARS)}; individual record examples are intentionally not written or released.","COVID and MedCoder windows are descriptive only.","",f"Runtime: Python {sys.version.split()[0]}; platform {platform.platform()}; elapsed {time.time()-t:.1f} seconds."]
    (OUT/"reports"/"statistical_analysis_report.md").write_text("\n".join(report)+"\n",encoding="utf-8")
    (OUT/"logs"/"analysis_run.json").write_text(json.dumps({"finished_utc":datetime.now(timezone.utc).isoformat(),"workers":workers,"elapsed_seconds":time.time()-t,"tests_passed":passed,"python":sys.version,"platform":platform.platform()},ensure_ascii=False,indent=2),encoding="utf-8")
    manifest()
    print(json.dumps({"A":sum(annual[y]["main"]["A"] for y in YEARS),"B":sum(annual[y]["main"]["B"] for y in YEARS),"tests_passed":passed,"workers":workers},ensure_ascii=False));return 0 if passed else 2

if __name__=="__main__": raise SystemExit(main())
