#!/usr/bin/env python3
"""Build a current-candidate-bound robustness matrix for the K80 UCF study.

This script is intentionally independent of every prior robustness-matrix
artifact.  It reads the final candidate aggregate/tables and, for definitions
not retained in that aggregate, scans the same 26 NCHS public-use source files.
Only de-identified age-by-sex counts leave the record scanner.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import multiprocessing as mp
import os
import platform
import subprocess
import sys
import time
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

ROOT = Path("results")
SCRIPT = Path(__file__).resolve()
CORE_SCRIPT = Path(__file__).with_name("analysis_reanalysis_public.py")
AGGREGATE = ROOT / "analysis" / "annual_aggregates.json"
INPUT_LEDGER = ROOT / "analysis" / "input_availability.csv"
WEIGHTS_TABLE = ROOT / "tables" / "k80_standardization_weights.csv"
ANNUAL_TABLE = ROOT / "tables" / "k80_standardized_annual.csv"
CONTRAST_TABLE = ROOT / "tables" / "k80_standardized_contrasts.csv"
OUTPUT_CSV = ROOT / "tables" / "robustness_matrix_current.csv"
OUTPUT_LOG = ROOT / "logs" / "robustness_current_run.json"
OUTPUT_EVIDENCE = ROOT / "evidence" / "robustness_current_evidence.md"

YEARS = tuple(range(1999, 2025))
REC_LEN = {y: (440 if y < 2003 else 488 if y < 2013 else 490 if y < 2020 else 817) for y in YEARS}
AGE_ORDER = ("0-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75-84", "85+")
SEX_ORDER = ("Male", "Female")
SEVEN_ZIP = Path("7z")
Z95 = 1.96
LAYOUT_VERSION = (
    "NCHS-MCD-core-v2: 1999-2002 resident=20,sex=59,age=64-66,ucd=142-145,"
    "entity_n=160-161,entity=162-301/7,record_n=338-339,record=341-440/5; "
    "2003-2024 resident=20,sex=69,age=70-73,ucd=146-149,entity_n=163-164,"
    "entity=165-304/7,record_n=341-342,record=344-443/5"
)

# Import only the current candidate's pure record-parsing helpers.  The module's
# output-building functions are never called; in particular no legacy staging
# data or robustness artifact is read.
sys.path.insert(0, str(CORE_SCRIPT.parent))
import analysis_reanalysis_public as core  # noqa: E402


def sha256_file(path: Path, block: int = 8 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(block)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def canonical_bundle_hash(items: Iterable[tuple[str, str]]) -> str:
    payload = "\n".join(f"{name}:{digest.lower()}" for name, digest in sorted(items))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def load_nchs_inputs() -> dict[int, dict[str, Any]]:
    rows = [r for r in read_csv(INPUT_LEDGER) if r.get("component") == "NCHS_MCD"]
    out: dict[int, dict[str, Any]] = {}
    for row in rows:
        year = int(row["year"])
        source = Path(row["path"])
        if year in out:
            raise RuntimeError(f"duplicate NCHS input row for {year}")
        if not source.is_file():
            raise FileNotFoundError(f"missing NCHS input for {year}: {source}")
        out[year] = {
            "year": year,
            "access_mode": row["access_mode"],
            "path": source,
            "zip_member": row.get("zip_member", ""),
            "expected_bytes": int(row["bytes"]),
            "expected_sha256": row["sha256"].lower(),
        }
    if set(out) != set(YEARS):
        raise RuntimeError(f"NCHS years incomplete: {sorted(set(YEARS) - set(out))}")
    return out


def _blank_endpoint() -> dict[str, Any]:
    return {"A": 0, "B": 0, "cells": defaultdict(lambda: {"A": 0, "B": 0})}


def _bump(endpoint: dict[str, Any], cell: str, denominator: bool, numerator: bool) -> None:
    if denominator:
        endpoint["A"] += 1
        endpoint["cells"][cell]["A"] += 1
        if numerator:
            endpoint["B"] += 1
            endpoint["cells"][cell]["B"] += 1


def _layout(year: int) -> dict[str, tuple[int, int]]:
    """Return zero-based slices from the official year-specific NCHS layouts."""
    if year <= 2002:
        return {
            "resident": (19, 20), "sex": (58, 59), "age": (63, 66), "ucd": (141, 145),
            "entity_n": (159, 161), "entity": (161, 301),
            "record_n": (337, 339), "record": (340, 440),
        }
    return {
        "resident": (19, 20), "sex": (68, 69), "age": (69, 73), "ucd": (145, 149),
        "entity_n": (162, 164), "entity": (164, 304),
        "record_n": (340, 342), "record": (343, 443),
    }


def _age_sex_cell(rec: bytes, layout: dict[str, tuple[int, int]]) -> str:
    age = core.age_group(rec[slice(*layout["age"])].decode("ascii", "replace"))
    if age in {"0-14", "15-24"}:
        age = "0-24"
    sex = core.sex_group(rec[slice(*layout["sex"])].decode("ascii", "replace"))
    return f"{age}|{sex}"


def _scan_stream(year: int, fh: Any, source_hasher: Any | None) -> dict[str, Any]:
    layout = _layout(year)
    record = _blank_endpoint()
    entity = _blank_endpoint()
    expanded = _blank_endpoint()
    total = bad_length = residents = nonresidents = 0
    for line in fh:
        if source_hasher is not None:
            source_hasher.update(line)
        total += 1
        rec = line.rstrip(b"\r\n")
        if len(rec) != REC_LEN[year]:
            bad_length += 1
            continue
        if rec[slice(*layout["resident"])].decode("ascii", "replace").strip() not in {"1", "2", "3"}:
            nonresidents += 1
            continue
        residents += 1
        cell = _age_sex_cell(rec, layout)
        rn = core.int_or_zero(rec[slice(*layout["record_n"])].decode("ascii", "replace"))
        en = core.int_or_zero(rec[slice(*layout["entity_n"])].decode("ascii", "replace"))
        record_codes = core.axis_codes(rec, layout["record"], rn, 5)
        entity_codes = core.axis_codes(rec, layout["entity"], en, 7)
        record_roots = {code[:3] for code in record_codes if code}
        entity_roots = {code[:3] for code in entity_codes if code}
        ucd_root = core.norm(rec[slice(*layout["ucd"])].decode("ascii", "replace"))[:3]
        record_a = "K80" in record_roots
        entity_a = "K80" in entity_roots
        expanded_a = bool(record_roots.intersection({"K80", "K81", "K82", "K83"}))
        _bump(record, cell, record_a, ucd_root == "K80")
        _bump(entity, cell, entity_a, ucd_root == "K80")
        _bump(expanded, cell, expanded_a, ucd_root in {"K80", "K81", "K82", "K83"})
    for endpoint in (record, entity, expanded):
        endpoint["cells"] = dict(endpoint["cells"])
    return {
        "year": year,
        "total_records": total,
        "bad_length": bad_length,
        "residents": residents,
        "nonresidents": nonresidents,
        "record_k80": record,
        "entity_k80": entity,
        "expanded_k80_k83": expanded,
    }


def scan_year(spec: dict[str, Any]) -> dict[str, Any]:
    year = int(spec["year"])
    source = Path(spec["path"])
    mode = spec["access_mode"]
    started = time.time()
    actual_size = source.stat().st_size
    if mode == "extracted":
        h = hashlib.sha256()
        with source.open("rb") as fh:
            result = _scan_stream(year, fh, h)
        actual_hash = h.hexdigest()
    elif mode == "zip_stream":
        actual_hash = sha256_file(source)
        member = spec["zip_member"]
        try:
            with zipfile.ZipFile(source) as archive:
                with archive.open(member) as fh:
                    result = _scan_stream(year, fh, None)
        except NotImplementedError:
            if not SEVEN_ZIP.is_file():
                raise RuntimeError(f"unsupported ZIP compression and 7-Zip missing for {year}")
            proc = subprocess.Popen(
                [str(SEVEN_ZIP), "e", "-so", str(source), member],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            assert proc.stdout is not None
            result = _scan_stream(year, proc.stdout, None)
            stderr = proc.stderr.read().decode("utf-8", "replace") if proc.stderr else ""
            if proc.wait() != 0:
                raise RuntimeError(f"7-Zip stream failure for {year}: {stderr[-1000:]}")
    else:
        raise RuntimeError(f"unsupported NCHS access mode for {year}: {mode}")
    result.update(
        {
            "source_basename": source.name,
            "access_mode": mode,
            "zip_member": spec.get("zip_member", ""),
            "source_bytes": actual_size,
            "source_sha256": actual_hash,
            "expected_bytes": spec["expected_bytes"],
            "expected_sha256": spec["expected_sha256"],
            "source_binding_pass": bool(
                actual_size == spec["expected_bytes"] and actual_hash == spec["expected_sha256"]
            ),
            "elapsed_seconds": time.time() - started,
        }
    )
    return result


def load_aggregate() -> dict[int, dict[str, Any]]:
    raw = json.loads(AGGREGATE.read_text(encoding="utf-8"))
    annual = {int(year): value for year, value in raw.items()}
    if set(annual) != set(YEARS):
        raise RuntimeError("final aggregate does not contain exactly 1999-2024")
    return annual


def rebin_aggregate_cells(cells: dict[str, dict[str, int]]) -> dict[str, dict[str, int]]:
    out: dict[str, dict[str, int]] = {}
    for cell, counts in cells.items():
        age, sex = cell.split("|")
        key = f"{'0-24' if age in {'0-14', '15-24'} else age}|{sex}"
        target = out.setdefault(key, {"A": 0, "B": 0, "B_official": 0})
        for field in ("A", "B", "B_official"):
            target[field] += int(counts.get(field, 0))
    return out


def load_weights() -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = defaultdict(dict)
    for row in read_csv(WEIGHTS_TABLE):
        out[row["scheme"]][row["age_sex_cell"]] = float(row["weight"])
    return dict(out)


def pooled_weights(
    annual: dict[int, dict[str, Any]], endpoint: str, pool_years: Iterable[int],
    allowed_ages: set[str] | None = None,
) -> dict[str, float]:
    counts: Counter[str] = Counter()
    for year in pool_years:
        for cell, q in annual[year][endpoint]["cells"].items():
            age, sex = cell.split("|")
            if age in AGE_ORDER and sex in SEX_ORDER and (allowed_ages is None or age in allowed_ages):
                counts[cell] += int(q["A"])
    common = {
        cell
        for cell in counts
        if all(int(annual[year][endpoint]["cells"].get(cell, {}).get("A", 0)) > 0 for year in YEARS)
    }
    counts = Counter({cell: n for cell, n in counts.items() if cell in common})
    total = sum(counts.values())
    if total <= 0:
        raise RuntimeError(f"no common support for {endpoint}")
    return {cell: count / total for cell, count in counts.items()}


def restricted_weights(weights: dict[str, float], allowed_ages: set[str]) -> dict[str, float]:
    kept = {cell: w for cell, w in weights.items() if cell.split("|")[0] in allowed_ages}
    total = sum(kept.values())
    if total <= 0:
        raise RuntimeError("restricted standard has zero weight")
    return {cell: w / total for cell, w in kept.items()}


def direct(cells: dict[str, dict[str, int]], weights: dict[str, float]) -> dict[str, float]:
    usable = []
    for cell, weight in weights.items():
        q = cells.get(cell, {})
        n = int(q.get("A", 0))
        k = int(q.get("B", 0))
        if n > 0 and weight > 0:
            usable.append((weight, k, n))
    support = sum(item[0] for item in usable)
    if not math.isclose(support, 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise RuntimeError(f"incomplete fixed-standard support: {support}")
    p = sum(weight * (k / n) for weight, k, n in usable)
    variance = sum((weight**2) * (k / n) * (1.0 - k / n) / n for weight, k, n in usable)
    se = math.sqrt(max(variance, 0.0))
    return {
        "p": p,
        "lo": max(0.0, p - Z95 * se),
        "hi": min(1.0, p + Z95 * se),
        "variance": variance,
        "weight_sum": support,
    }


def wilson(k: int, n: int) -> dict[str, float]:
    if n <= 0:
        raise RuntimeError("Wilson interval requires n > 0")
    p = k / n
    den = 1.0 + Z95 * Z95 / n
    center = p + Z95 * Z95 / (2.0 * n)
    width = Z95 * math.sqrt(p * (1.0 - p) / n + Z95 * Z95 / (4.0 * n * n))
    return {"p": p, "lo": (center - width) / den, "hi": (center + width) / den,
            "variance": p * (1.0 - p) / n, "weight_sum": float("nan")}


def contrast(first: dict[str, float], last: dict[str, float]) -> dict[str, float]:
    rd = last["p"] - first["p"]
    se = math.sqrt(first["variance"] + last["variance"])
    return {"rd": rd, "lo": rd - Z95 * se, "hi": rd + Z95 * se}


def table_reference() -> tuple[dict[tuple[str, int], dict[str, float]], dict[str, dict[str, float]]]:
    annual: dict[tuple[str, int], dict[str, float]] = {}
    for row in read_csv(ANNUAL_TABLE):
        annual[(row["scheme"], int(row["year"]))] = {
            "p": float(row["std_UCF"]), "lo": float(row["CI_lo"]),
            "hi": float(row["CI_hi"]), "variance": float(row["variance"]),
            "weight_sum": float(row["effective_weight_sum"]),
        }
    contrasts: dict[str, dict[str, float]] = {}
    for row in read_csv(CONTRAST_TABLE):
        if row["contrast_or_year"] == "2024_minus_1999":
            contrasts[row["scheme"]] = {
                "rd": float(row["RD"]), "lo": float(row["RD_CI_lo"]), "hi": float(row["RD_CI_hi"])
            }
    return annual, contrasts


def row_from_estimates(
    analysis_id: str, label: str, axis: str, code_set: str, denominator: str,
    numerator: str, standardization: str, weight_definition: str,
    first: dict[str, float], last: dict[str, float], rd: dict[str, float],
    counts: dict[int, tuple[int, int]], input_hash: str, source: str,
    ci_method: str,
) -> dict[str, Any]:
    return {
        "analysis_id": analysis_id,
        "analysis_label": label,
        "mention_axis": axis,
        "code_set": code_set,
        "denominator_definition": denominator,
        "numerator_definition": numerator,
        "standardization": standardization,
        "standard_weight_definition": weight_definition,
        "ci_method": ci_method,
        "estimate_1999": first["p"],
        "ci_lower_1999": first["lo"],
        "ci_upper_1999": first["hi"],
        "estimate_2024": last["p"],
        "ci_lower_2024": last["lo"],
        "ci_upper_2024": last["hi"],
        "rd_2024_minus_1999": rd["rd"],
        "rd_ci_lower": rd["lo"],
        "rd_ci_upper": rd["hi"],
        "A_1999": counts[1999][0],
        "B_1999": counts[1999][1],
        "A_2024": counts[2024][0],
        "B_2024": counts[2024][1],
        "weight_sum": "" if math.isnan(first.get("weight_sum", float("nan"))) else first["weight_sum"],
        "input_sha256": input_hash,
        "source_binding": source,
        "execution_status": "PASS",
    }


def checkpoint() -> int:
    _, contrasts = table_reference()
    current = contrasts["primary_2018_2024"]["rd"]
    print(json.dumps({
        "checkpoint": "PASS",
        "current_primary_rd": current,
        "legacy_0.1080_used": False,
        "matches_current_candidate": abs(current - 0.1046499068926785) < 1e-12,
    }, ensure_ascii=False))
    return 0 if abs(current - 0.1046499068926785) < 1e-12 and abs(current - 0.1080) > 1e-4 else 2


def main() -> int:
    global ROOT, AGGREGATE, INPUT_LEDGER, WEIGHTS_TABLE, ANNUAL_TABLE, CONTRAST_TABLE, OUTPUT_CSV, OUTPUT_LOG, OUTPUT_EVIDENCE, SEVEN_ZIP
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--checkpoint", action="store_true")
    parser.add_argument("--analysis-root", default="results", help="output directory generated by analysis_reanalysis_public.py")
    parser.add_argument("--seven-zip", default="7z")
    args = parser.parse_args()
    ROOT = Path(args.analysis_root); AGGREGATE = ROOT / "analysis" / "annual_aggregates.json"; INPUT_LEDGER = ROOT / "analysis" / "input_availability.csv"; WEIGHTS_TABLE = ROOT / "tables" / "k80_standardization_weights.csv"; ANNUAL_TABLE = ROOT / "tables" / "k80_standardized_annual.csv"; CONTRAST_TABLE = ROOT / "tables" / "k80_standardized_contrasts.csv"; OUTPUT_CSV = ROOT / "tables" / "robustness_matrix_current.csv"; OUTPUT_LOG = ROOT / "logs" / "robustness_current_run.json"; OUTPUT_EVIDENCE = ROOT / "evidence" / "robustness_current_evidence.md"; SEVEN_ZIP = Path(args.seven_zip)
    if args.checkpoint:
        return checkpoint()

    started = time.time()
    start_utc = datetime.now(timezone.utc).isoformat()
    workers = max(1, min(8, int(args.workers), len(YEARS)))
    input_specs = load_nchs_inputs()
    aggregate = load_aggregate()
    aggregate_cells = {year: rebin_aggregate_cells(aggregate[year]["cells"]) for year in YEARS}
    final_weights = load_weights()
    reference_annual, reference_contrasts = table_reference()

    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=workers) as pool:
        scan_results = pool.map(scan_year, [input_specs[year] for year in YEARS])
    scanned = {int(result["year"]): result for result in scan_results}

    raw_hash_items = [(f"NCHS_MCD_{year}", scanned[year]["source_sha256"]) for year in YEARS]
    raw_bundle_hash = canonical_bundle_hash(raw_hash_items)
    aggregate_hash = sha256_file(AGGREGATE)
    table_hash_items = [
        ("annual_aggregates.json", aggregate_hash),
        ("k80_standardization_weights.csv", sha256_file(WEIGHTS_TABLE)),
        ("k80_standardized_annual.csv", sha256_file(ANNUAL_TABLE)),
        ("k80_standardized_contrasts.csv", sha256_file(CONTRAST_TABLE)),
    ]
    current_bundle_hash = canonical_bundle_hash(table_hash_items + raw_hash_items)

    verifications: list[dict[str, Any]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        verifications.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})

    check("all_26_sources_bound", all(scanned[y]["source_binding_pass"] for y in YEARS),
          f"bound={sum(int(scanned[y]['source_binding_pass']) for y in YEARS)}/26")
    check("all_record_lengths_valid", all(scanned[y]["bad_length"] == 0 for y in YEARS),
          f"bad_length_total={sum(scanned[y]['bad_length'] for y in YEARS)}")
    check("scan_total_matches_final_aggregate",
          all(scanned[y]["total_records"] == int(aggregate[y]["total_records"]) for y in YEARS),
          f"scanned={sum(scanned[y]['total_records'] for y in YEARS)}; aggregate={sum(int(aggregate[y]['total_records']) for y in YEARS)}")
    check("resident_total_matches_final_aggregate",
          all(scanned[y]["residents"] == int(aggregate[y]["residents"]) for y in YEARS),
          f"scanned={sum(scanned[y]['residents'] for y in YEARS)}; aggregate={sum(int(aggregate[y]['residents']) for y in YEARS)}")
    check("record_k80_counts_match_final_aggregate",
          all(scanned[y]["record_k80"]["A"] == int(aggregate[y]["main"]["A"]) and
              scanned[y]["record_k80"]["B"] == int(aggregate[y]["main"]["B"]) for y in YEARS),
          "annual record-axis A/B exact comparison")

    # Verify that final table weights themselves are reproducible from the final aggregate.
    derived_weight_schemes = {
        "primary_2018_2024": pooled_weights(
            {y: {"record_k80": {"cells": aggregate_cells[y]}} for y in YEARS},
            "record_k80", range(2018, 2025)),
        "sensitivity_1999": pooled_weights(
            {y: {"record_k80": {"cells": aggregate_cells[y]}} for y in YEARS},
            "record_k80", [1999]),
        "sensitivity_full_period": pooled_weights(
            {y: {"record_k80": {"cells": aggregate_cells[y]}} for y in YEARS},
            "record_k80", YEARS),
    }
    for scheme, weights in final_weights.items():
        check(f"weight_sum_{scheme}", math.isclose(sum(weights.values()), 1.0, abs_tol=1e-12),
              f"sum={sum(weights.values()):.17g}")
        derived = derived_weight_schemes[scheme]
        same = set(weights) == set(derived) and all(abs(weights[c] - derived[c]) < 1e-15 for c in weights)
        check(f"weights_reproduced_{scheme}", same, f"cells={len(weights)}")

    endpoint_counts = {
        endpoint: {year: (scanned[year][endpoint]["A"], scanned[year][endpoint]["B"]) for year in YEARS}
        for endpoint in ("record_k80", "entity_k80", "expanded_k80_k83")
    }
    for endpoint, counts in endpoint_counts.items():
        check(f"B_subset_A_{endpoint}", all(b <= a for a, b in counts.values()), "all 26 annual B<=A")

    rows: list[dict[str, Any]] = []
    primary_weights = final_weights["primary_2018_2024"]

    # Current primary, and two alternative fixed standards, are emitted exactly
    # as frozen in the current candidate after independent recomputation checks.
    scheme_specs = [
        ("ROBUST-01", "Primary current definition", "primary_2018_2024",
         "2018-2024 pooled record-axis K80 mention age×sex distribution; cells require positive A in every 1999-2024 year"),
        ("ROBUST-07", "1999 standard population", "sensitivity_1999",
         "1999 record-axis K80 mention age×sex distribution; cells require positive A in every 1999-2024 year"),
        ("ROBUST-08", "Full-period standard population", "sensitivity_full_period",
         "1999-2024 pooled record-axis K80 mention age×sex distribution; cells require positive A in every 1999-2024 year"),
    ]
    for analysis_id, label, scheme, weight_def in scheme_specs:
        recomputed_first = direct(aggregate_cells[1999], final_weights[scheme])
        recomputed_last = direct(aggregate_cells[2024], final_weights[scheme])
        recomputed_rd = contrast(recomputed_first, recomputed_last)
        ref_first = reference_annual[(scheme, 1999)]
        ref_last = reference_annual[(scheme, 2024)]
        ref_rd = reference_contrasts[scheme]
        values_match = all(
            abs(a - b) < 1e-12
            for a, b in [
                (recomputed_first["p"], ref_first["p"]), (recomputed_first["lo"], ref_first["lo"]),
                (recomputed_first["hi"], ref_first["hi"]), (recomputed_last["p"], ref_last["p"]),
                (recomputed_last["lo"], ref_last["lo"]), (recomputed_last["hi"], ref_last["hi"]),
                (recomputed_rd["rd"], ref_rd["rd"]), (recomputed_rd["lo"], ref_rd["lo"]),
                (recomputed_rd["hi"], ref_rd["hi"]),
            ]
        )
        check(f"current_table_recomputed_{scheme}", values_match, "point estimates and 95% CIs")
        rows.append(row_from_estimates(
            analysis_id, label, "record axis", "K80", "resident deaths with ≥1 record-axis K80 code",
            "denominator deaths whose UCD is K80", "fixed direct age×sex standardization", weight_def,
            ref_first, ref_last, ref_rd,
            {y: (int(aggregate[y]["main"]["A"]), int(aggregate[y]["main"]["B"])) for y in YEARS},
            current_bundle_hash, "current final aggregate + current final standardization tables + same-version raw NCHS bundle",
            "within-cell binomial variance; independent-year normal 95% CI; limits for proportions clipped to [0,1]",
        ))

    # A-prime is an intentionally crude denominator sensitivity.  B is the
    # official UCD-K80 count and is therefore a subset of A union UCD-K80.
    aprime_counts = {
        year: (int(aggregate[year]["a_prime"]), int(aggregate[year]["main"]["B_official"]))
        for year in YEARS
    }
    aprime_first = wilson(aprime_counts[1999][1], aprime_counts[1999][0])
    aprime_last = wilson(aprime_counts[2024][1], aprime_counts[2024][0])
    aprime_rd = contrast(aprime_first, aprime_last)
    rows.append(row_from_estimates(
        "ROBUST-02", "A-prime crude denominator sensitivity", "record axis plus UCD union", "K80",
        "resident deaths with record-axis K80 OR UCD K80 (A′)", "denominator deaths whose official UCD is K80",
        "crude (not standardized)", "not applicable; annual A′ denominator",
        aprime_first, aprime_last, aprime_rd, aprime_counts, aggregate_hash,
        "current final aggregate", "Wilson 95% CI for annual proportions; independent-year binomial normal 95% CI for RD",
    ))
    check("B_subset_A_A_prime", all(b <= a for a, b in aprime_counts.values()), "all 26 annual B<=A′")

    scanned_annual = {year: scanned[year] for year in YEARS}
    entity_weights = pooled_weights(scanned_annual, "entity_k80", range(2018, 2025))
    expanded_weights = pooled_weights(scanned_annual, "expanded_k80_k83", range(2018, 2025))
    check("weight_sum_entity_axis", math.isclose(sum(entity_weights.values()), 1.0, abs_tol=1e-12),
          f"sum={sum(entity_weights.values()):.17g}; cells={len(entity_weights)}")
    check("weight_sum_expanded_endpoint", math.isclose(sum(expanded_weights.values()), 1.0, abs_tol=1e-12),
          f"sum={sum(expanded_weights.values()):.17g}; cells={len(expanded_weights)}")

    for analysis_id, label, endpoint, axis, code_set, denom, numer, weights, weight_def in [
        ("ROBUST-03", "Entity-axis K80 sensitivity", "entity_k80", "entity axis", "K80",
         "resident deaths with ≥1 entity-axis K80 code", "denominator deaths whose UCD is K80",
         entity_weights, "2018-2024 pooled entity-axis K80 mention age×sex distribution; cells require positive A in every 1999-2024 year"),
        ("ROBUST-04", "Expanded biliary-code sensitivity", "expanded_k80_k83", "record axis", "K80-K83",
         "resident deaths with ≥1 record-axis K80-K83 code", "denominator deaths whose UCD is K80-K83",
         expanded_weights, "2018-2024 pooled record-axis K80-K83 mention age×sex distribution; cells require positive A in every 1999-2024 year"),
    ]:
        first = direct(scanned[1999][endpoint]["cells"], weights)
        last = direct(scanned[2024][endpoint]["cells"], weights)
        rd = contrast(first, last)
        rows.append(row_from_estimates(
            analysis_id, label, axis, code_set, denom, numer, "fixed direct age×sex standardization",
            weight_def, first, last, rd, endpoint_counts[endpoint], raw_bundle_hash,
            "same-version raw NCHS bundle rescanned with current candidate parser",
            "within-cell binomial variance; independent-year normal 95% CI; limits for proportions clipped to [0,1]",
        ))

    age_specs = [
        ("ROBUST-05", "Age ≥25 years", {"25-34", "35-44", "45-54", "55-64", "65-74", "75-84", "85+"}),
        ("ROBUST-06", "Age ≥65 years", {"65-74", "75-84", "85+"}),
    ]
    for analysis_id, label, ages in age_specs:
        weights = restricted_weights(primary_weights, ages)
        check(f"weight_sum_{analysis_id}", math.isclose(sum(weights.values()), 1.0, abs_tol=1e-12),
              f"sum={sum(weights.values()):.17g}; cells={len(weights)}")
        first = direct(scanned[1999]["record_k80"]["cells"], weights)
        last = direct(scanned[2024]["record_k80"]["cells"], weights)
        rd = contrast(first, last)
        counts: dict[int, tuple[int, int]] = {}
        for year in YEARS:
            cells = scanned[year]["record_k80"]["cells"]
            counts[year] = (
                sum(int(q["A"]) for cell, q in cells.items() if cell.split("|")[0] in ages),
                sum(int(q["B"]) for cell, q in cells.items() if cell.split("|")[0] in ages),
            )
        rows.append(row_from_estimates(
            analysis_id, label, "record axis", "K80",
            f"resident deaths aged {'≥25' if analysis_id == 'ROBUST-05' else '≥65'} with ≥1 record-axis K80 code",
            "denominator deaths whose UCD is K80", "fixed direct age×sex standardization within age restriction",
            "current primary 2018-2024 pooled record-axis K80 age×sex weights restricted to eligible ages and renormalized once",
            first, last, rd, counts, current_bundle_hash,
            "same-version raw NCHS bundle + current final primary weights",
            "within-cell binomial variance; independent-year normal 95% CI; limits for proportions clipped to [0,1]",
        ))
        check(f"B_subset_A_{analysis_id}", all(b <= a for a, b in counts.values()), "all 26 annual B<=A")

    # Sort into the intended display order after building from logically grouped sources.
    rows.sort(key=lambda row: row["analysis_id"])
    for row in rows:
        numeric = [
            float(row["estimate_1999"]), float(row["ci_lower_1999"]), float(row["ci_upper_1999"]),
            float(row["estimate_2024"]), float(row["ci_lower_2024"]), float(row["ci_upper_2024"]),
            float(row["rd_2024_minus_1999"]), float(row["rd_ci_lower"]), float(row["rd_ci_upper"]),
        ]
        ordered = (
            numeric[1] <= numeric[0] <= numeric[2]
            and numeric[4] <= numeric[3] <= numeric[5]
            and numeric[7] <= numeric[6] <= numeric[8]
        )
        complete = all(math.isfinite(value) for value in numeric)
        subset = int(row["B_1999"]) <= int(row["A_1999"]) and int(row["B_2024"]) <= int(row["A_2024"])
        check(f"row_integrity_{row['analysis_id']}", ordered and complete and subset,
              f"CI_ordered={ordered}; finite={complete}; B<=A={subset}")

    primary_row = next(row for row in rows if row["analysis_id"] == "ROBUST-01")
    primary_match = all(
        abs(float(primary_row[field]) - reference) < 1e-12
        for field, reference in [
            ("estimate_1999", reference_annual[("primary_2018_2024", 1999)]["p"]),
            ("ci_lower_1999", reference_annual[("primary_2018_2024", 1999)]["lo"]),
            ("ci_upper_1999", reference_annual[("primary_2018_2024", 1999)]["hi"]),
            ("estimate_2024", reference_annual[("primary_2018_2024", 2024)]["p"]),
            ("ci_lower_2024", reference_annual[("primary_2018_2024", 2024)]["lo"]),
            ("ci_upper_2024", reference_annual[("primary_2018_2024", 2024)]["hi"]),
            ("rd_2024_minus_1999", reference_contrasts["primary_2018_2024"]["rd"]),
            ("rd_ci_lower", reference_contrasts["primary_2018_2024"]["lo"]),
            ("rd_ci_upper", reference_contrasts["primary_2018_2024"]["hi"]),
        ]
    )
    check("primary_point_and_ci_exact_current_candidate", primary_match,
          f"RD={primary_row['rd_2024_minus_1999']:.15f}")
    check("legacy_base_rd_not_used", abs(float(primary_row["rd_2024_minus_1999"]) - 0.1080) > 1e-4,
          "legacy 0.1080 is neither an input nor the current primary estimate")
    check("endpoint_support_1999_2024", len(rows) >= 8 and all(int(r["A_1999"]) > 0 and int(r["A_2024"]) > 0 for r in rows),
          f"rows={len(rows)}")

    fields = [
        "analysis_id", "analysis_label", "mention_axis", "code_set", "denominator_definition",
        "numerator_definition", "standardization", "standard_weight_definition", "ci_method",
        "estimate_1999", "ci_lower_1999", "ci_upper_1999", "estimate_2024", "ci_lower_2024",
        "ci_upper_2024", "rd_2024_minus_1999", "rd_ci_lower", "rd_ci_upper", "A_1999", "B_1999",
        "A_2024", "B_2024", "weight_sum", "input_sha256", "source_binding", "execution_status",
    ]
    all_pass = all(item["status"] == "PASS" for item in verifications)
    if not all_pass:
        for row in rows:
            row["execution_status"] = "FAIL"
    write_csv(OUTPUT_CSV, rows, fields)
    table_hash = sha256_file(OUTPUT_CSV)

    evidence_lines = [
        "# Current-candidate robustness evidence",
        "",
        f"**Execution status:** {'PASS' if all_pass else 'FAIL'}  ",
        f"**Candidate aggregate SHA256:** `{aggregate_hash}`  ",
        f"**Raw NCHS bundle SHA256:** `{raw_bundle_hash}`  ",
        f"**Current-bound input bundle SHA256:** `{current_bundle_hash}`  ",
        f"**Output matrix SHA256:** `{table_hash}`",
        "",
        "## Why this replaces the legacy matrix",
        "",
        f"The current primary fixed-standard RD is `{float(primary_row['rd_2024_minus_1999']):.15f}` "
        "(2024 minus 1999), reproduced from the final candidate aggregate and final candidate weights. "
        "The former approximate value 0.1080 was not read or used. No prior robustness-matrix file is an input.",
        "",
        "## Frozen definitions",
        "",
        "- Primary: record-axis K80 mention denominator; K80 UCD numerator; final 2018-2024 fixed age×sex standard.",
        "- A-prime: crude union of record-axis K80 or official K80 UCD; numerator is official K80 UCD.",
        "- Entity-axis: at least one entity-axis K80 code; numerator additionally has K80 UCD.",
        "- Expanded boundary: at least one record-axis K80-K83 code; numerator additionally has K80-K83 UCD.",
        "- Age restrictions: record-axis K80 endpoint, with the final primary standard restricted to eligible age cells and renormalized once.",
        "- Alternative standards: the final 1999 and full-period fixed standards, without annual renormalization.",
        "- COVID years and the MedCoder window are not endpoint sensitivities here; they remain descriptive institutional windows only.",
        "",
        "## Inference",
        "",
        "Fixed-standard estimates use within-cell binomial variances and independent-year normal 95% intervals; annual proportion limits are clipped to [0,1]. "
        "A-prime annual intervals are Wilson intervals, with an independent-year binomial normal interval for the RD. These are design-based descriptive "
        "contrasts of recorded UCD selection, not causal effects or clinical adjudication.",
        "",
        "## Verification",
        "",
        "| Check | Status | Detail |",
        "|---|---:|---|",
    ]
    evidence_lines.extend(f"| {v['check']} | {v['status']} | {v['detail']} |" for v in verifications)
    evidence_lines.extend([
        "",
        f"Passed: `{sum(v['status'] == 'PASS' for v in verifications)}/{len(verifications)}`.",
        "",
        "Only de-identified aggregate counts were retained; no line-level identifiers or examples were written.",
    ])
    OUTPUT_EVIDENCE.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_EVIDENCE.write_text("\n".join(evidence_lines) + "\n", encoding="utf-8")

    log = {
        "task": "current-robustness-reanalysis",
        "status": "PASS" if all_pass else "FAIL",
        "started_utc": start_utc,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.time() - started,
        "workers": workers,
        "python": sys.version,
        "platform": platform.platform(),
        "layout_version": LAYOUT_VERSION,
        "code": {
            "analysis_robustness_current.py": sha256_file(SCRIPT),
            "analysis_reanalysis.py": sha256_file(CORE_SCRIPT),
        },
        "candidate_inputs": {
            "annual_aggregates.json": aggregate_hash,
            "input_availability.csv": sha256_file(INPUT_LEDGER),
            "k80_standardization_weights.csv": sha256_file(WEIGHTS_TABLE),
            "k80_standardized_annual.csv": sha256_file(ANNUAL_TABLE),
            "k80_standardized_contrasts.csv": sha256_file(CONTRAST_TABLE),
            "raw_nchs_bundle_sha256": raw_bundle_hash,
            "current_bound_bundle_sha256": current_bundle_hash,
        },
        "raw_sources": [
            {
                "logical_id": f"NCHS_MCD_{year}",
                "year": year,
                "access_mode": scanned[year]["access_mode"],
                "source_basename": scanned[year]["source_basename"],
                "zip_member": scanned[year]["zip_member"],
                "bytes": scanned[year]["source_bytes"],
                "sha256": scanned[year]["source_sha256"],
                "expected_sha256_match": scanned[year]["source_binding_pass"],
                "records_scanned": scanned[year]["total_records"],
                "resident_records": scanned[year]["residents"],
                "bad_length": scanned[year]["bad_length"],
                "elapsed_seconds": scanned[year]["elapsed_seconds"],
            }
            for year in YEARS
        ],
        "records_scanned": sum(scanned[y]["total_records"] for y in YEARS),
        "resident_records": sum(scanned[y]["residents"] for y in YEARS),
        "result_rows": len(rows),
        "primary_rd": float(primary_row["rd_2024_minus_1999"]),
        "legacy_0_1080_used": False,
        "verification": {
            "passed": sum(v["status"] == "PASS" for v in verifications),
            "failed": sum(v["status"] == "FAIL" for v in verifications),
            "total": len(verifications),
            "checks": verifications,
        },
        "outputs": {
            "robustness_matrix_current.csv": table_hash,
            "robustness_current_evidence.md": sha256_file(OUTPUT_EVIDENCE),
        },
        "legacy_dependencies": [],
        "institutional_windows": "COVID and MedCoder windows retained as descriptive only; not labeled endpoint sensitivities.",
    }
    OUTPUT_LOG.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_LOG.write_text(json.dumps(log, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    summary = {
        "status": log["status"],
        "elapsed_seconds": log["elapsed_seconds"],
        "workers": workers,
        "records_scanned": log["records_scanned"],
        "resident_records": log["resident_records"],
        "result_rows": len(rows),
        "primary_rd": log["primary_rd"],
        "verification": f"{log['verification']['passed']}/{log['verification']['total']} PASS",
        "outputs": [str(OUTPUT_CSV), str(OUTPUT_LOG), str(OUTPUT_EVIDENCE)],
    }
    print(json.dumps(summary, ensure_ascii=False))
    return 0 if all_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
