#!/usr/bin/env python3
"""Public release verifier for the aggregate reproducibility package.

Run from the repository root: python verify_release.py
Use --write-manifest only after intentional release-content changes.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "docs" / "FILE_MANIFEST.csv"
TEXT_SUFFIXES = {".py", ".r", ".md", ".txt", ".csv", ".json", ".xml", ".cff", ".yml", ".yaml"}
EXCLUDED_FROM_MANIFEST = {"docs/FILE_MANIFEST.csv"}
EXPECTED_DESTINATIONS = {
    "CIRCULATORY", "DIABETES_NUTRITIONAL", "DIGESTIVE_OTHER", "EXTERNAL",
    "INFECTIOUS_PARASITIC", "NEOPLASMS", "OTHER", "OTHER_GALLBLADDER_BILIARY", "RESPIRATORY",
}
REQUIRED = {
    "README.md", "CITATION.cff", "LICENSE", ".gitignore", "verify_release.py",
    "docs/REPRODUCIBILITY.md", "docs/DATA_AVAILABILITY.md", "docs/DATA_DICTIONARY.md",
    "code/analysis/ANALYSIS_SPECIFICATION.md", "code/figures/revision_figures.R",
    "environment/requirements.txt", "environment/R_environment.md", "config/paths.example.json",
    "data/derived/Table1_Cohort_and_estimands.csv", "data/derived/Table2_Main_estimates.csv",
}
FORBIDDEN_CONTENT = [
    (re.compile("master" "2333", re.I), "password marker"),
    (re.compile("C:" + chr(92) * 2 + "Users" + chr(92) * 2 + "|C:/" + "Users/|" + chr(92) * 4 + "Users" + chr(92) * 2, re.I), "private absolute path"),
    (re.compile("internal_" "private|record_axis_" "codes|orphan_" "examples", re.I), "private audit material"),
    (re.compile("__py" "cache__", re.I), "cache marker"),
    (re.compile(r"\b(?:password|passwd)\s*[:=]", re.I), "credential assignment"),
    (re.compile(r"\b(?:[0-9]{1,3}\.){3}[0-9]{1,3}\b"), "IP address"),
]


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def release_files() -> list[Path]:
    return sorted(
        p
        for p in ROOT.rglob("*")
        if p.is_file()
        and not (p.relative_to(ROOT).parts and p.relative_to(ROOT).parts[0] == ".git")
        and rel(p) not in EXCLUDED_FROM_MANIFEST
    )


def write_manifest() -> None:
    rows = [{"path": rel(p), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in release_files()]
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(rows)


def check(label: str, condition: bool, detail: str, checks: list[dict]) -> None:
    checks.append({"check": label, "passed": bool(condition), "detail": detail})


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-manifest", action="store_true", help="write docs/FILE_MANIFEST.csv, then verify it")
    args = parser.parse_args()
    if args.write_manifest:
        write_manifest()

    checks: list[dict] = []
    missing = sorted(p for p in REQUIRED if not (ROOT / p).is_file())
    check("required_release_files", not missing, f"missing={missing}", checks)

    main_figs = [ROOT / "data" / "figure_source_data" / f"Figure{i}_source_data.csv" for i in range(1, 7)]
    supp_figs = [ROOT / "data" / "figure_source_data" / f"Supplementary_Figure_S{i}_source_data.csv" for i in range(1, 4)]
    missing_figs = [rel(p) for p in main_figs + supp_figs if not p.is_file()]
    check("figure_source_data_coverage", not missing_figs, f"missing={missing_figs}; main=6; supplementary=3", checks)

    xmls = list((ROOT / "queries" / "cdc_wonder_xml").glob("*.xml")) if (ROOT / "queries" / "cdc_wonder_xml").exists() else []
    check("wonder_query_coverage", len(xmls) >= 12 and (ROOT / "queries" / "cdc_wonder_xml" / "query_manifest.csv").is_file(), f"xml_count={len(xmls)}", checks)

    annual_path = ROOT / "data" / "derived" / "k80_annual_main.csv"
    try:
        annual = list(csv.DictReader(annual_path.open(encoding="utf-8-sig", newline="")))
        a_total = sum(int(r["A_record_axis_K80"]) for r in annual)
        b_total = sum(int(r["B_main_A_and_UCD_K80"]) for r in annual)
        identity = all(int(r["gap_A_minus_B"]) == int(r["A_record_axis_K80"]) - int(r["B_main_A_and_UCD_K80"]) for r in annual)
        check("core_totals", a_total == 51084 and b_total == 27514, f"A={a_total}; B={b_total}", checks)
        check("annual_gap_identity", identity and len(annual) == 26, f"years={len(annual)}", checks)
    except Exception as exc:
        check("annual_tables_readable", False, repr(exc), checks)

    destination_path = ROOT / "data" / "derived" / "ucd_destination_annual.csv"
    try:
        destination = list(csv.DictReader(destination_path.open(encoding="utf-8-sig", newline="")))
        closure = True
        years = sorted({int(r["year"]) for r in destination})
        for year in years:
            rows = [r for r in destination if int(r["year"]) == year]
            closure &= {r["destination"] for r in rows} == EXPECTED_DESTINATIONS
            closure &= abs(sum(float(r["std_probability"]) for r in rows) - 1.0) < 1e-8
            closure &= sum(int(r["raw_n"]) for r in rows) == int(rows[0]["gap_denominator_A_minus_B"])
        check("destination_partition_and_probability", closure and len(years) == 26, f"years={len(years)}; categories=9", checks)
    except Exception as exc:
        check("destination_table_readable", False, repr(exc), checks)

    text_hits: list[str] = []
    for path in release_files():
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        try:
            content = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for pattern, name in FORBIDDEN_CONTENT:
            if pattern.search(content):
                text_hits.append(f"{rel(path)}: {name}")
    check("public_boundary_scan", not text_hits, f"hits={text_hits}", checks)

    if not MANIFEST.is_file():
        check("manifest_present", False, "docs/FILE_MANIFEST.csv is absent", checks)
    else:
        rows = list(csv.DictReader(MANIFEST.open(encoding="utf-8-sig", newline="")))
        expected = {row["path"]: row for row in rows}
        actual = {rel(p): p for p in release_files()}
        mismatch = sorted(set(expected).symmetric_difference(actual))
        mismatch += sorted(path for path, p in actual.items() if path in expected and (expected[path]["sha256"] != sha256(p) or int(expected[path]["bytes"]) != p.stat().st_size))
        check("sha256_manifest", not mismatch, f"entries={len(rows)}; mismatches={mismatch}", checks)

    report = {"passed": all(item["passed"] for item in checks), "checks": checks}
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
