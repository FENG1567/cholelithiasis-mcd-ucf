#!/usr/bin/env python3
"""Independent full-scan sensitivity analysis for Figure 3 definitions.

Reads the NCHS public-use fixed-width Multiple Cause of Death files using the
frozen schema crosswalk and the frozen resident/record-axis K80 cohort
definition. Only aggregate, de-identified outputs are written.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import multiprocessing as mp
import platform
import subprocess
import sys
import time
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

OUT = Path("results")
NCHS = Path("data/raw")
SCHEMA_PATH = Path("config/schema_crosswalk.csv")
SEVEN_ZIP = Path("7z")
YEARS = list(range(1999, 2025))
REC_LEN = {y: (440 if y < 2003 else 488 if y < 2013 else 490 if y < 2020 else 817) for y in YEARS}
ENTITY_STATUSES = ("Part_I_only", "Part_II_only", "Both", "entity_no_K80", "entity_unknown")
PART_STRATEGIES = ("Part_I_priority", "Part_II_priority", "Both_separate", "Both_excluded")
COMPLEXITY_COUNTERS = ("record_n", "unique_record_axis_icd")
COMPLEXITY_LEVELS = ("1", "2", "3-4", "5-9", "10+", "unknown")


def norm(value: str) -> str:
    return (value or "").strip().upper().replace(".", "")


def pos(row: dict, name: str):
    value = (row.get(name) or "").strip()
    if not value or value == "MISSING" or "-" not in value:
        return None
    left, right = value.split("-", 1)
    if not left.isdigit() or not right.isdigit():
        return None
    return int(left) - 1, int(right)


def slc(record: bytes, position):
    return record[position[0]:position[1]].decode("ascii", "replace") if position else ""


def int_or_zero(value: str) -> int:
    try:
        return int(value.strip())
    except (AttributeError, TypeError, ValueError):
        return 0


def axis_codes(record: bytes, position, n: int):
    if not position or n <= 0:
        return []
    values = []
    for i in range(min(n, 20)):
        offset = position[0] + i * 5
        code = norm(record[offset:offset + 4].decode("ascii", "replace"))
        if code:
            values.append(code)
    return values


def entity_k80_status(record: bytes, position, n: int):
    """Classify K80 across every entity position, not only the first one."""
    if not position or n <= 0:
        return "entity_no_K80"
    parts = set()
    for i in range(min(n, 20)):
        offset = position[0] + i * 7
        code = norm(record[offset + 2:offset + 6].decode("ascii", "replace"))
        if not code.startswith("K80"):
            continue
        indicator = record[offset:offset + 1].decode("ascii", "replace").strip()
        if indicator in {"1", "2", "3", "4", "5"}:
            parts.add("I")
        elif indicator == "6":
            parts.add("II")
        else:
            parts.add("UNKNOWN")
    if parts == {"I"}:
        return "Part_I_only"
    if parts == {"II"}:
        return "Part_II_only"
    if parts == {"I", "II"}:
        return "Both"
    if parts:
        return "entity_unknown"
    return "entity_no_K80"


def strategy_group(status: str, strategy: str):
    if strategy == "Part_I_priority":
        if status in {"Part_I_only", "Both"}:
            return "Part_I"
        if status == "Part_II_only":
            return "Part_II"
        return status
    if strategy == "Part_II_priority":
        if status in {"Part_II_only", "Both"}:
            return "Part_II"
        if status == "Part_I_only":
            return "Part_I"
        return status
    if strategy == "Both_separate":
        return status
    if strategy == "Both_excluded":
        return "excluded_Both" if status == "Both" else status
    raise ValueError(strategy)


def complexity_label(count: int) -> str:
    if count == 1:
        return "1"
    if count == 2:
        return "2"
    if 3 <= count <= 4:
        return "3-4"
    if 5 <= count <= 9:
        return "5-9"
    if count >= 10:
        return "10+"
    return "unknown"


def load_crosswalk():
    result = {}
    with SCHEMA_PATH.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            year = int(row["year"])
            row["year"] = year
            row["p"] = {key: pos(row, key) for key in row if key not in {"year", "p"}}
            result[year] = row
    return result


def source_for(year: int):
    folder = NCHS / "raw" / str(year)
    extracted = folder / "extracted"
    files = [path for path in extracted.glob("*") if path.is_file()] if extracted.exists() else []
    if files:
        return "file", max(files, key=lambda path: path.stat().st_size), ""
    archives = list(folder.glob("*.zip"))
    if not archives:
        raise FileNotFoundError(f"no NCHS source for {year}: {folder}")
    archive = max(archives, key=lambda path: path.stat().st_size)
    with zipfile.ZipFile(archive) as handle:
        member = max((item for item in handle.infolist() if not item.is_dir()), key=lambda item: item.file_size).filename
    return "zip_stream", archive, member


def open_source(mode, source: Path, member: str):
    if mode == "file":
        return source.open("rb"), None
    try:
        archive = zipfile.ZipFile(source)
        return archive.open(member), archive
    except NotImplementedError:
        if not SEVEN_ZIP.exists():
            raise RuntimeError(f"ZIP compression unsupported and 7-Zip unavailable: {source}")
        process = subprocess.Popen([str(SEVEN_ZIP), "e", "-so", str(source), member],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if process.stdout is None:
            raise RuntimeError(f"7-Zip did not expose stdout: {source}")
        return process.stdout, process


def close_source(handle, owner):
    try:
        handle.close()
    finally:
        if owner is None:
            return
        if isinstance(owner, zipfile.ZipFile):
            owner.close()
            return
        stderr = owner.stderr.read().decode("utf-8", "replace") if owner.stderr else ""
        code = owner.wait()
        if code != 0:
            raise RuntimeError(f"7-Zip stream failed exit={code}: {stderr[-1000:]}")


def blank():
    return {"A": 0, "B": 0}


def new_year_result(year: int):
    return {
        "year": year,
        "source_mode": "",
        "source_path": "",
        "zip_member": "",
        "record_length_expected": REC_LEN[year],
        "total_records": 0,
        "bad_length": 0,
        "residents": 0,
        "nonresidents": 0,
        "A": 0,
        "B": 0,
        "entity_status": {name: blank() for name in ENTITY_STATUSES},
        "strategies": {name: defaultdict(blank) for name in PART_STRATEGIES},
        "complexity": {name: defaultdict(blank) for name in COMPLEXITY_COUNTERS},
        "complexity_count_difference": 0,
        "complexity_difference_by_delta": Counter(),
    }


def scan_year(year: int):
    crosswalk = load_crosswalk()
    row = crosswalk[year]
    mode, source, member = source_for(year)
    result = new_year_result(year)
    result["source_mode"] = mode
    result["source_path"] = str(source)
    result["zip_member"] = member

    def process(handle):
        for line in handle:
            record = line.rstrip(b"\r\n")
            result["total_records"] += 1
            if len(record) != REC_LEN[year]:
                result["bad_length"] += 1
                continue
            resident = slc(record, row["p"]["resident_status"]).strip() in {"1", "2", "3"}
            if not resident:
                result["nonresidents"] += 1
                continue
            result["residents"] += 1
            record_n = int_or_zero(slc(record, row["p"]["record_n"]))
            entity_n = int_or_zero(slc(record, row["p"]["entity_n"]))
            codes = axis_codes(record, row["p"]["record_block"], record_n)
            has_record_k80 = "K80" in {code[:3] for code in codes if code}
            ucd = norm(slc(record, row["p"]["ucd_code"]))
            has_ucd_k80 = ucd[:3] == "K80"
            if not has_record_k80:
                continue

            result["A"] += 1
            result["B"] += int(has_ucd_k80)
            status = entity_k80_status(record, row["p"]["entity_block"], entity_n)
            if status not in result["entity_status"]:
                status = "entity_unknown"
            result["entity_status"][status]["A"] += 1
            result["entity_status"][status]["B"] += int(has_ucd_k80)
            for strategy in PART_STRATEGIES:
                group = strategy_group(status, strategy)
                result["strategies"][strategy][group]["A"] += 1
                result["strategies"][strategy][group]["B"] += int(has_ucd_k80)

            raw_count = record_n
            unique_count = len(set(codes))
            raw_label = complexity_label(raw_count)
            unique_label = complexity_label(unique_count)
            result["complexity"]["record_n"][raw_label]["A"] += 1
            result["complexity"]["record_n"][raw_label]["B"] += int(has_ucd_k80)
            result["complexity"]["unique_record_axis_icd"][unique_label]["A"] += 1
            result["complexity"]["unique_record_axis_icd"][unique_label]["B"] += int(has_ucd_k80)
            if raw_count != unique_count:
                result["complexity_count_difference"] += 1
                result["complexity_difference_by_delta"][str(unique_count - raw_count)] += 1

    handle, owner = open_source(mode, source, member)
    try:
        process(handle)
    finally:
        close_source(handle, owner)
    result["entity_status"] = {key: dict(value) for key, value in result["entity_status"].items()}
    result["strategies"] = {key: {group: dict(value) for group, value in groups.items()} for key, groups in result["strategies"].items()}
    result["complexity"] = {key: {group: dict(value) for group, value in groups.items()} for key, groups in result["complexity"].items()}
    result["complexity_difference_by_delta"] = dict(result["complexity_difference_by_delta"])
    return result


def wilson(k, n, z=1.959963984540054):
    if n <= 0:
        return math.nan, math.nan, math.nan
    p = k / n
    den = 1 + z * z / n
    centre = p + z * z / (2 * n)
    width = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return p, max(0.0, (centre - width) / den), min(1.0, (centre + width) / den)


def kitagawa(groups_by_year: dict, start: int, end: int):
    g0, g1 = groups_by_year[start], groups_by_year[end]
    den0 = sum(value["A"] for value in g0.values())
    den1 = sum(value["A"] for value in g1.values())
    p0_total = sum(value["B"] for value in g0.values()) / den0 if den0 else math.nan
    p1_total = sum(value["B"] for value in g1.values()) / den1 if den1 else math.nan
    composition = selection = 0.0
    entries = []
    for group in sorted(set(g0) | set(g1)):
        a0, b0 = g0.get(group, blank())["A"], g0.get(group, blank())["B"]
        a1, b1 = g1.get(group, blank())["A"], g1.get(group, blank())["B"]
        w0, w1 = (a0 / den0 if den0 else 0.0), (a1 / den1 if den1 else 0.0)
        p0, p1 = (b0 / a0 if a0 else 0.0), (b1 / a1 if a1 else 0.0)
        comp = 0.5 * (p0 + p1) * (w1 - w0)
        sel = 0.5 * (w0 + w1) * (p1 - p0)
        composition += comp
        selection += sel
        entries.append({"group": group, "A_start": a0, "B_start": b0, "A_end": a1, "B_end": b1,
                        "start_weight": w0, "end_weight": w1, "start_UCF": p0, "end_UCF": p1,
                        "composition_component": comp, "selection_component": sel})
    total = p1_total - p0_total
    return {"start_den": den0, "end_den": den1, "start_UCF": p0_total, "end_UCF": p1_total,
            "total_change": total, "composition_total": composition, "selection_total": selection,
            "residual": total - composition - selection, "entries": entries}


def finite_or_blank(value):
    return "" if not math.isfinite(value) else value


def make_rows(annual):
    fields = ["analysis", "scenario", "comparison", "start_year", "end_year", "year", "group",
              "A_records", "B_records", "denominator_records", "UCF", "UCF_CI_lo", "UCF_CI_hi",
              "start_weight", "end_weight", "start_UCF", "end_UCF", "composition_component",
              "selection_component", "total_change", "composition_total", "selection_total",
              "decomposition_residual", "excluded_A_records", "excluded_B_records", "metric",
              "metric_value", "note"]
    rows = []
    def append(**kwargs):
        row = {field: "" for field in fields}
        row.update(kwargs)
        rows.append(row)

    for year in YEARS:
        q = annual[year]
        p, lo, hi = wilson(q["B"], q["A"])
        append(analysis="overall", scenario="primary_record_axis_K80", comparison="annual", year=year,
               group="all_A_records", A_records=q["A"], B_records=q["B"], denominator_records=q["A"],
               UCF=finite_or_blank(p), UCF_CI_lo=finite_or_blank(lo), UCF_CI_hi=finite_or_blank(hi),
               note="A=resident record-axis K80; B=A and UCD root K80; crude UCF.")

    for year in YEARS:
        q = annual[year]
        for group in ENTITY_STATUSES:
            c = q["entity_status"].get(group, blank())
            p, lo, hi = wilson(c["B"], c["A"])
            append(analysis="entity_status", scenario="raw_entity_axis_status", comparison="annual", year=year,
                   group=group, A_records=c["A"], B_records=c["B"], denominator_records=c["A"],
                   UCF=finite_or_blank(p), UCF_CI_lo=finite_or_blank(lo), UCF_CI_hi=finite_or_blank(hi),
                   note="Exhaustive entity-axis K80 status among A.")

    strategy_groups = {}
    for strategy in PART_STRATEGIES:
        by_year = {}
        for year in YEARS:
            groups = annual[year]["strategies"][strategy]
            excluded_a = groups.get("excluded_Both", blank())["A"]
            excluded_b = groups.get("excluded_Both", blank())["B"]
            included = {group: c for group, c in groups.items() if group != "excluded_Both"}
            by_year[year] = included
            den = sum(c["A"] for c in included.values())
            for group, c in sorted(included.items()):
                p, lo, hi = wilson(c["B"], c["A"])
                note = {
                    "Part_I_priority": "Both assigned to Part I; all A retained.",
                    "Part_II_priority": "Both assigned to Part II; all A retained.",
                    "Both_separate": "Part I-only, Part II-only, Both and entity-no-K80 are separate; all A retained.",
                    "Both_excluded": "Both records excluded from Part-stratified denominator; other A retained.",
                }[strategy]
                append(analysis="part_rule", scenario=strategy, comparison="annual", year=year, group=group,
                       A_records=c["A"], B_records=c["B"], denominator_records=den, UCF=finite_or_blank(p),
                       UCF_CI_lo=finite_or_blank(lo), UCF_CI_hi=finite_or_blank(hi),
                       excluded_A_records=excluded_a, excluded_B_records=excluded_b, note=note)
            b = sum(c["B"] for c in included.values())
            p, lo, hi = wilson(b, den)
            append(analysis="part_rule_overall", scenario=strategy, comparison="annual", year=year,
                   group="included_all", A_records=den, B_records=b, denominator_records=den,
                   UCF=finite_or_blank(p), UCF_CI_lo=finite_or_blank(lo), UCF_CI_hi=finite_or_blank(hi),
                   excluded_A_records=excluded_a, excluded_B_records=excluded_b,
                   note="Overall UCF over strategy-defined included A records.")
        strategy_groups[strategy] = by_year

    for strategy, by_year in strategy_groups.items():
        for start, end, label in ((1999, 2024, "1999_vs_2024"), (2018, 2024, "2018_vs_2024")):
            result = kitagawa(by_year, start, end)
            for e in result["entries"]:
                append(analysis="part_rule_kitagawa", scenario=strategy, comparison=label,
                       start_year=start, end_year=end, group=e["group"], A_records=e["A_start"],
                       B_records=e["B_start"], denominator_records=result["start_den"],
                       start_weight=e["start_weight"], end_weight=e["end_weight"], start_UCF=e["start_UCF"],
                       end_UCF=e["end_UCF"], composition_component=e["composition_component"],
                       selection_component=e["selection_component"], total_change=result["total_change"],
                       composition_total=result["composition_total"], selection_total=result["selection_total"],
                       decomposition_residual=result["residual"],
                       note="Kitagawa two-way descriptive accounting; components are not causal mediation.")
            append(analysis="part_rule_kitagawa_total", scenario=strategy, comparison=label, start_year=start,
                   end_year=end, group="TOTAL", A_records=result["start_den"],
                   B_records=sum(e["B_start"] for e in result["entries"]), denominator_records=result["start_den"],
                   start_UCF=result["start_UCF"], end_UCF=result["end_UCF"], total_change=result["total_change"],
                   composition_total=result["composition_total"], selection_total=result["selection_total"],
                   decomposition_residual=result["residual"], note="TOTAL row is endpoint UCF change and identity check.")

    complexity_groups = {}
    for counter in COMPLEXITY_COUNTERS:
        by_year = {}
        for year in YEARS:
            groups = annual[year]["complexity"][counter]
            by_year[year] = groups
            for group in COMPLEXITY_LEVELS:
                c = groups.get(group, blank())
                p, lo, hi = wilson(c["B"], c["A"])
                note = ("record_n: raw record-axis ICD position count." if counter == "record_n"
                        else "unique_record_axis_icd: distinct nonblank normalized four-character record-axis ICD codes.")
                append(analysis="complexity", scenario=counter, comparison="annual", year=year, group=group,
                       A_records=c["A"], B_records=c["B"], denominator_records=c["A"], UCF=finite_or_blank(p),
                       UCF_CI_lo=finite_or_blank(lo), UCF_CI_hi=finite_or_blank(hi), note=note)
        complexity_groups[counter] = by_year

    for counter, by_year in complexity_groups.items():
        for start, end, label in ((1999, 2024, "1999_vs_2024"), (2018, 2024, "2018_vs_2024")):
            result = kitagawa(by_year, start, end)
            for e in result["entries"]:
                append(analysis="complexity_kitagawa", scenario=counter, comparison=label,
                       start_year=start, end_year=end, group=e["group"], A_records=e["A_start"],
                       B_records=e["B_start"], denominator_records=result["start_den"],
                       start_weight=e["start_weight"], end_weight=e["end_weight"], start_UCF=e["start_UCF"],
                       end_UCF=e["end_UCF"], composition_component=e["composition_component"],
                       selection_component=e["selection_component"], total_change=result["total_change"],
                       composition_total=result["composition_total"], selection_total=result["selection_total"],
                       decomposition_residual=result["residual"],
                       note="Kitagawa two-way descriptive accounting; components are not causal mediation.")
            append(analysis="complexity_kitagawa_total", scenario=counter, comparison=label, start_year=start,
                   end_year=end, group="TOTAL", A_records=result["start_den"],
                   B_records=sum(e["B_start"] for e in result["entries"]), denominator_records=result["start_den"],
                   start_UCF=result["start_UCF"], end_UCF=result["end_UCF"], total_change=result["total_change"],
                   composition_total=result["composition_total"], selection_total=result["selection_total"],
                   decomposition_residual=result["residual"], note="TOTAL row is endpoint UCF change and identity check.")

    for year in YEARS:
        q = annual[year]
        append(analysis="complexity_counter_audit", scenario="record_n_vs_unique_record_axis_icd",
               comparison="annual", year=year, group="all_A_records", A_records=q["A"],
               metric="count_records_with_record_n_neq_unique_count",
               metric_value=q["complexity_count_difference"],
               note="The count objects differ when their integer counts differ; target total is 46 across 1999-2024.")
    return fields, rows


def definitions_rows():
    fields = ["definition_id", "domain", "rule_or_group", "level", "field_or_codes", "operation",
              "threshold_or_priority", "denominator_role", "notes", "source_layout"]
    rows = []
    def add(*values):
        rows.append(dict(zip(fields, values)))
    add("COHORT-A", "cohort", "A", "record", "record-axis roots include K80", "resident record with >=1 record-axis ICD code whose first 3 characters are K80", "record_n positions 1..20; K80 root", "A denominator", "Record-axis coding cohort, not literal certificate mentions.", "record_n + record_block") 
    add("COHORT-B", "cohort", "B", "record", "A and UCD root K80", "A record whose single UCD code begins K80", "UCD field exact root K80", "B numerator", "B is nested in A by construction.", "ucd_code")
    add("COHORT-RES", "cohort", "resident", "record", "resident_status in {1,2,3}", "retain resident deaths only", "codes 1,2,3", "analysis universe", "Malformed lengths and nonresidents are excluded before A/B classification.", "resident_status")
    add("SEV-1", "severity", "acute_cholecystitis", "K80 subtype", "K80.0", "mutually exclusive subtype priority", "priority 2 (after K80.3)", "within-A stratum", "Exact four-character K80.0; the established priority is K80.3, K80.0, K80.4, K80.1, K80.5, K80.2, K80.8.", "record-axis K80 codes")
    add("SEV-2", "severity", "cholangitis", "K80 subtype", "K80.3", "mutually exclusive subtype priority", "priority 1", "within-A stratum", "Exact four-character K80.3; highest priority in the established subtype order.", "record-axis K80 codes")
    add("SEV-3", "severity", "other_cholecystitis", "K80 subtype", "K80.1; K80.4", "mutually exclusive subtype priority", "K80.4 priority 3; K80.1 priority 4", "within-A stratum", "K80.1 and K80.4 are combined after applying the established order.", "record-axis K80 codes")
    add("SEV-4", "severity", "without_cholecystitis", "K80 subtype", "K80.2; K80.5", "mutually exclusive subtype priority", "K80.5 priority 5; K80.2 priority 6", "within-A stratum", "K80.2 and K80.5 are combined after applying the established order.", "record-axis K80 codes")
    add("SEV-5", "severity", "other_or_unspecified", "K80 subtype", "K80.8; other/blank K80 suffix", "fallback after specified subtype codes", "priority 7/fallback", "within-A stratum", "Sparse or unspecified K80 records retained; exact priority order is K80.3, K80.0, K80.4, K80.1, K80.5, K80.2, K80.8.", "record-axis K80 codes")
    add("ENT-1", "entity_part", "Part_I_only", "entity status", "K80 in Part I and not Part II", "evaluate all entity positions", "status set={I}", "within-A stratum", "Part assignment is from entity axis, not record-axis order.", "entity_n + entity_block; 7-character positions")
    add("ENT-2", "entity_part", "Part_II_only", "entity status", "K80 in Part II and not Part I", "evaluate all entity positions", "status set={II}", "within-A stratum", "Part II indicator is 6.", "entity_n + entity_block; 7-character positions")
    add("ENT-3", "entity_part", "Both", "entity status", "K80 in Part I and Part II", "evaluate all entity positions", "status set={I,II}", "within-A stratum or exclusion", "Key dual-location sensitivity target; expected total 337.", "entity_n + entity_block; 7-character positions")
    add("ENT-4", "entity_part", "entity_no_K80", "entity status", "no entity-axis K80 among A", "evaluate all entity positions", "empty K80 set", "within-A stratum", "Record-axis K80 can occur without entity-axis K80; expected total 2.", "entity_n + entity_block; 7-character positions")
    add("ENT-5", "entity_part", "entity_unknown", "entity status", "K80 with unknown Part indicator", "retain explicit unknown status", "unknown indicator", "within-A stratum", "Should be zero under the scanned public-use layouts.", "entity_n + entity_block; 7-character positions")
    add("PART-1", "part_rule", "Part_I_priority", "strategy", "Both -> Part_I", "assign Both to Part I; retain all A", "Part I precedence", "all A retained", "Part II-only remains Part II; entity-no-K80 is retained.", "ENT-1..ENT-5")
    add("PART-2", "part_rule", "Part_II_priority", "strategy", "Both -> Part_II", "assign Both to Part II; retain all A", "Part II precedence", "all A retained", "Part I-only remains Part I; entity-no-K80 is retained.", "ENT-1..ENT-5")
    add("PART-3", "part_rule", "Both_separate", "strategy", "Part I-only; Part II-only; Both; entity_no_K80", "retain four mutually exclusive statuses", "no priority", "all A retained", "Both displayed as a separate stratum.", "ENT-1..ENT-5")
    add("PART-4", "part_rule", "Both_excluded", "strategy", "Both excluded", "remove Both from numerator and denominator; retain other A", "exclude status=Both", "A/B denominator excludes Both", "Sensitivity only; primary cohort unchanged.", "ENT-1..ENT-5")
    add("CMP-1", "complexity", "record_n", "certificate complexity", "record_n", "count raw record-axis ICD positions", "1; 2; 3-4; 5-9; 10+; unknown", "within-A stratum", "Original manuscript counter; bounded to 20 NCHS record positions.", "record_n + record_block")
    add("CMP-2", "complexity", "unique_record_axis_icd", "certificate complexity", "distinct normalized nonblank record-axis ICD codes", "deduplicate first record_n codes, then count", "1; 2; 3-4; 5-9; 10+; unknown", "within-A stratum", "Sensitivity for repeated ICD codes.", "record_n + record_block")
    add("CMP-3", "complexity", "unknown", "certificate complexity", "nonpositive or invalid count", "assign unknown; retain denominator", "count <=0 or unavailable", "within-A stratum", "No record silently dropped for complexity failure.", "record_n + record_block")
    add("KIT-1", "decomposition", "Kitagawa", "accounting", "w_g and p_g", "0.5*(p0+p1)*(w1-w0) + 0.5*(w0+w1)*(p1-p0)", "1999/2024 and 2018/2024", "descriptive identity", "Components are not causal mediation; residual must be <1e-10.", "aggregate A/B by group")
    return fields, rows


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def validate(annual, rows):
    checks = []
    def check(name, passed, detail):
        checks.append({"name": name, "passed": bool(passed), "detail": detail})
    total_a = sum(q["A"] for q in annual.values())
    total_b = sum(q["B"] for q in annual.values())
    check("all_years_present", set(annual) == set(YEARS), sorted(annual))
    check("bad_length_zero", sum(q["bad_length"] for q in annual.values()) == 0, sum(q["bad_length"] for q in annual.values()))
    check("resident_A_B_totals", total_a == 51084 and total_b == 27514, {"A": total_a, "B": total_b})
    status_total = sum(sum(c["A"] for c in q["entity_status"].values()) for q in annual.values())
    check("entity_status_closure", status_total == total_a, {"entity_A": status_total, "A": total_a})
    status_counts = {key: sum(q["entity_status"][key]["A"] for q in annual.values()) for key in ENTITY_STATUSES}
    check("known_entity_status_counts", status_counts == {"Part_I_only": 30489, "Part_II_only": 20256, "Both": 337, "entity_no_K80": 2, "entity_unknown": 0}, status_counts)
    diff_total = sum(q["complexity_count_difference"] for q in annual.values())
    by_delta = {delta: sum(q["complexity_difference_by_delta"].get(delta, 0) for q in annual.values()) for delta in sorted(set(delta for q in annual.values() for delta in q["complexity_difference_by_delta"]))}
    check("known_complexity_counter_difference", diff_total == 46, {"difference_records": diff_total, "by_delta": by_delta})
    check("complexity_record_n_closure", all(sum(c["A"] for c in q["complexity"]["record_n"].values()) == q["A"] for q in annual.values()), "record_n categories close to A in every year")
    check("complexity_unique_closure", all(sum(c["A"] for c in q["complexity"]["unique_record_axis_icd"].values()) == q["A"] for q in annual.values()), "unique categories close to A in every year")
    for strategy in PART_STRATEGIES:
        for year in YEARS:
            groups = annual[year]["strategies"][strategy]
            included = sum(c["A"] for g, c in groups.items() if g != "excluded_Both")
            expected = annual[year]["A"] - groups.get("excluded_Both", blank())["A"]
            check(f"{strategy}_{year}_closure", included == expected, {"included_A": included, "expected_A": expected})
    for row in rows:
        if row.get("analysis", "").endswith("_kitagawa_total"):
            residual = float(row["decomposition_residual"])
            check(f"{row['analysis']}_{row['scenario']}_{row['comparison']}_identity", abs(residual) < 1e-10, residual)
    return checks


def input_manifest():
    records = []
    for year in YEARS:
        mode, source, member = source_for(year)
        records.append({"year": year, "mode": mode, "path": str(source), "member": member,
                        "bytes": source.stat().st_size, "sha256": sha256(source)})
    records.append({"year": "schema", "mode": "read_only", "path": str(SCHEMA_PATH), "member": "",
                    "bytes": SCHEMA_PATH.stat().st_size, "sha256": sha256(SCHEMA_PATH)})
    return records


def build_evidence(annual, checks, meta):
    status_a = {key: sum(q["entity_status"][key]["A"] for q in annual.values()) for key in ENTITY_STATUSES}
    status_b = {key: sum(q["entity_status"][key]["B"] for q in annual.values()) for key in ENTITY_STATUSES}
    diff_total = sum(q["complexity_count_difference"] for q in annual.values())
    strategy_groups = {}
    for strategy in PART_STRATEGIES:
        strategy_groups[strategy] = {}
        for year in YEARS:
            groups = annual[year]["strategies"][strategy]
            strategy_groups[strategy][year] = {group: counts for group, counts in groups.items() if group != "excluded_Both"}
    complexity_groups = {
        counter: {year: annual[year]["complexity"][counter] for year in YEARS}
        for counter in COMPLEXITY_COUNTERS
    }
    endpoint_summaries = []
    for label, groups_by_scenario in (("Part rule", strategy_groups), ("Complexity counter", complexity_groups)):
        for scenario, by_year in groups_by_scenario.items():
            for start, end, comparison in ((1999, 2024, "1999-2024"), (2018, 2024, "2018-2024")):
                endpoint_summaries.append((label, scenario, comparison, kitagawa(by_year, start, end)))
    lines = [
        "# Figure 3 composition-sensitivity evidence", "",
        "## Scope and cohort", "",
        "Independent full scan of NCHS public-use MCD fixed-width files for 1999-2024. The frozen schema crosswalk is used; resident statuses 1/2/3 are retained; A is a resident record with a record-axis K80 root; B is A with UCD root K80. No individual-level records are written.", "",
        f"Scanned records: **{sum(q['total_records'] for q in annual.values()):,}**; residents: **{sum(q['residents'] for q in annual.values()):,}**; A: **{sum(q['A'] for q in annual.values()):,}**; B: **{sum(q['B'] for q in annual.values()):,}**; malformed lengths: **{sum(q['bad_length'] for q in annual.values())}**.", "",
        "## Entity-axis Part audit", "", "| Entity status | A records | B records | UCF |", "|---|---:|---:|---:|"
    ]
    for key in ENTITY_STATUSES:
        p = status_b[key] / status_a[key] if status_a[key] else math.nan
        lines.append(f"| {key} | {status_a[key]:,} | {status_b[key]:,} | {p:.10f} |")
    lines += [
        "", f"The target counts are reproduced: Part I-only 30,489; Part II-only 20,256; Both 337; entity-axis no K80 2; entity-unknown {status_a['entity_unknown']}. These exhaustive statuses close to A=51,084.", "",
        "## Part-rule sensitivity", "",
        "The sensitivity table contains annual UCFs for Part I priority, Part II priority, Both separate, and Both excluded, as well as 1999-2024 and 2018-2024 Kitagawa rows. Both-excluded removes Both from numerator and denominator only for that sensitivity; it does not alter the primary cohort.", "",
        "Each endpoint identity is UCF change = composition component + within-stratum component + residual. The residual tolerance is 1e-10. These are descriptive accounting terms, not causal mediation.", "",
        "## Endpoint Kitagawa totals", "", "| Sensitivity | Endpoint | UCF change | Composition | Within-stratum | Residual |", "|---|---|---:|---:|---:|---:|"
    ]
    for label, scenario, comparison, result in endpoint_summaries:
        lines.append(f"| {label}: {scenario} | {comparison} | {result['total_change']:.10f} | {result['composition_total']:.10f} | {result['selection_total']:.10f} | {result['residual']:.2e} |")
    part_1999 = [result for label, scenario, comparison, result in endpoint_summaries if label == "Part rule" and comparison == "1999-2024"]
    comp_1999 = [result for label, scenario, comparison, result in endpoint_summaries if label == "Complexity counter" and comparison == "1999-2024"]
    part_change_spread = max(result["total_change"] for result in part_1999) - min(result["total_change"] for result in part_1999)
    part_composition_spread = max(result["composition_total"] for result in part_1999) - min(result["composition_total"] for result in part_1999)
    complexity_composition_spread = max(result["composition_total"] for result in comp_1999) - min(result["composition_total"] for result in comp_1999)
    lines += [
        "", f"Across Part rules at 1999-2024, the maximum spread in total UCF change is {part_change_spread:.10f} ({part_change_spread * 100:.4f} percentage points); the maximum spread in the composition component is {part_composition_spread:.10f} ({part_composition_spread * 100:.4f} percentage points). The two complexity counters have identical total UCF change, while their composition components differ by {complexity_composition_spread:.10f} ({complexity_composition_spread * 100:.4f} percentage points).", "",
        "## Complexity-counter sensitivity", "",
        f"record_n counts raw record-axis ICD positions. unique_record_axis_icd deduplicates nonblank normalized four-character record-axis ICD codes before applying the same bins: 1, 2, 3-4, 5-9, 10+, and unknown. The counters differ for **{diff_total:,}** A records; annual counts and deltas are in the run JSON. Both category systems close to A in every year.", "",
        "## Verification", "", "| Check | Result | Detail |", "|---|---|---|"
    ]
    for item in checks:
        lines.append(f"| {item['name']} | {'PASS' if item['passed'] else 'FAIL'} | {json.dumps(item['detail'], ensure_ascii=False)} |")
    lines += [
        "", "## Reproduction and hashes", "",
        f"Command: python code/analysis/analysis_composition_sensitivity_public.py --workers {meta['workers']}", 
        f"Runtime UTC: {meta['finished_utc']}; Python: {meta['python'].splitlines()[0]}; platform: {meta['platform']}.", "",
        "Final SHA-256 values for the script, both tables, and this evidence are stored in logs/composition_sensitivity_run.json. The analysis is limited to recorded coding/attribution patterns; it is not clinical adjudication and A-B is not an error count.", ""
    ]
    return "\n".join(lines)


def main():
    global OUT, NCHS, SCHEMA_PATH, SEVEN_ZIP
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output-root", default="results")
    parser.add_argument("--nchs-data-root", default="data/raw")
    parser.add_argument("--schema-crosswalk", default="config/schema_crosswalk.csv")
    parser.add_argument("--seven-zip", default="7z")
    args = parser.parse_args()
    OUT = Path(args.output_root); NCHS = Path(args.nchs_data_root); SCHEMA_PATH = Path(args.schema_crosswalk); SEVEN_ZIP = Path(args.seven_zip)
    workers = max(1, min(8, args.workers, len(YEARS)))
    started = time.time()
    if workers == 1:
        results = [scan_year(year) for year in YEARS]
    else:
        context = mp.get_context("spawn")
        with context.Pool(workers) as pool:
            results = pool.map(scan_year, YEARS)
    annual = {result["year"]: result for result in results}
    fields, rows = make_rows(annual)
    checks = validate(annual, rows)
    if not all(item["passed"] for item in checks):
        raise RuntimeError("validation failed: " + repr([item for item in checks if not item["passed"]]))
    definition_fields, definitions = definitions_rows()
    sensitivity_path = OUT / "tables" / "k80_composition_sensitivity.csv"
    definitions_path = OUT / "tables" / "k80_composition_definitions.csv"
    log_path = OUT / "logs" / "composition_sensitivity_run.json"
    evidence_path = OUT / "evidence" / "composition_sensitivity_evidence.md"
    write_csv(sensitivity_path, fields, rows)
    write_csv(definitions_path, definition_fields, definitions)
    finished = datetime.now(timezone.utc).isoformat()
    meta = {
        "analysis": "independent_k80_composition_sensitivity",
        "finished_utc": finished,
        "workers": workers,
        "years": YEARS,
        "record_lengths": REC_LEN,
        "python": sys.version,
        "platform": platform.platform(),
        "elapsed_seconds_scan_and_write": time.time() - started,
        "cohort_definition": "resident_status in {1,2,3}; A=record-axis K80; B=A and UCD root K80",
        "entity_status_counts_A": {key: sum(q["entity_status"][key]["A"] for q in annual.values()) for key in ENTITY_STATUSES},
        "entity_status_counts_B": {key: sum(q["entity_status"][key]["B"] for q in annual.values()) for key in ENTITY_STATUSES},
        "complexity_counter_difference_records": sum(q["complexity_count_difference"] for q in annual.values()),
        "complexity_counter_difference_by_delta": {delta: sum(q["complexity_difference_by_delta"].get(delta, 0) for q in annual.values()) for delta in sorted(set(delta for q in annual.values() for delta in q["complexity_difference_by_delta"]))},
        "annual_scan_summary": {str(year): {"source_mode": q["source_mode"], "source_path": q["source_path"], "zip_member": q["zip_member"],
            "total_records": q["total_records"], "bad_length": q["bad_length"], "residents": q["residents"], "nonresidents": q["nonresidents"],
            "A": q["A"], "B": q["B"], "entity_status": q["entity_status"], "complexity_count_difference": q["complexity_count_difference"],
            "complexity_difference_by_delta": q["complexity_difference_by_delta"]} for year, q in sorted(annual.items())},
        "input_files": input_manifest(),
        "checks": checks,
        "outputs": {"sensitivity_table": str(sensitivity_path.relative_to(OUT)).replace("\\\\", "/"),
                    "definitions_table": str(definitions_path.relative_to(OUT)).replace("\\\\", "/"),
                    "evidence": str(evidence_path.relative_to(OUT)).replace("\\\\", "/")},
    }
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    evidence_path.parent.mkdir(parents=True, exist_ok=True)
    evidence_path.write_text(build_evidence(annual, checks, meta), encoding="utf-8")
    meta["script_sha256"] = sha256(Path(__file__).resolve())
    meta["final_output_hashes"] = {"sensitivity_table": sha256(sensitivity_path), "definitions_table": sha256(definitions_path), "evidence": sha256(evidence_path)}
    meta["log_sha256_before_final_log_update"] = sha256(log_path)
    meta["elapsed_seconds_total"] = time.time() - started
    log_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"A": sum(q["A"] for q in annual.values()), "B": sum(q["B"] for q in annual.values()),
        "entity_status_counts_A": meta["entity_status_counts_A"], "complexity_counter_difference_records": meta["complexity_counter_difference_records"],
        "checks_passed": sum(item["passed"] for item in checks), "checks_total": len(checks), "script_sha256": meta["script_sha256"],
        "output_hashes": meta["final_output_hashes"], "log_sha256": sha256(log_path)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
