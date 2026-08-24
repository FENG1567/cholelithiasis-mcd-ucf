#!/usr/bin/env python3
"""Joint-bootstrap inference extension for the frozen K80 endpoint analyses.

This script intentionally reads only the frozen aggregate file and existing
point-estimate tables.  It neither rescans raw mortality files nor alters the
primary A/B/UCF estimands.  The unit sampled for Kitagawa analyses is the
joint (stratum, B) contingency table within each endpoint year; the unit for
destination analyses is the full nine-category multinomial vector within each
age-sex cell.  Thus composition and conditional UCF, or all destination
categories, are re-estimated jointly on every bootstrap replicate.
"""

from __future__ import annotations

import csv
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path("results")
ANNUAL_PATH = ROOT / "analysis" / "annual_aggregates.json"
DECOMP_POINT_PATH = ROOT / "tables" / "k80_kitagawa_decomposition.csv"
DEST_POINT_PATH = ROOT / "tables" / "ucd_destination_contrasts.csv"
DEST_WEIGHT_PATH = ROOT / "tables" / "ucd_destination_standardization_weights.csv"
OUT_DECOMP = ROOT / "tables" / "k80_kitagawa_decomposition_uncertainty.csv"
OUT_DEST = ROOT / "tables" / "ucd_destination_contrasts_simultaneous.csv"
OUT_REPORT = ROOT / "reports" / "inference_extension_report.md"
OUT_LOG = ROOT / "logs" / "inference_extension_run.json"
OUT_MAIN_TABLE = ROOT / "tables" / "Table2_Main_estimates.csv"

SEED = 20260823
REPS = 100_000
ALPHA = 0.05
GROUPS = ("subtype", "severity", "entity_part", "complexity")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def finite_float(value: float) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("non-finite output")
    return value


def signed_points(value: float) -> str:
    """Format a probability difference as signed percentage points."""
    points = 100.0 * float(value)
    return ("+" if points >= 0 else "−") + f"{abs(points):.2f} percentage points"


def update_main_estimates(dest_rows: list[dict[str, object]]) -> None:
    """Bind the publication summary table to the current simultaneous output."""
    rows = read_csv(OUT_MAIN_TABLE)
    by_destination = {str(row["destination"]): row for row in dest_rows}
    targets = {
        "CIRCULATORY destination RD": "CIRCULATORY",
        "OTHER destination RD": "OTHER",
        "DIGESTIVE_OTHER destination RD": "DIGESTIVE_OTHER",
    }
    updated = set()
    for row in rows:
        destination = targets.get(row.get("Analysis", ""))
        if destination is None:
            continue
        source = by_destination[destination]
        row["Estimate"] = signed_points(float(source["RD_2024_minus_1999_point"]))
        row["95% CI or component"] = "FWER simultaneous 95% CI {} to {}".format(
            signed_points(float(source["simultaneous_95_ci_lo"])).replace(" percentage points", ""),
            signed_points(float(source["simultaneous_95_ci_hi"])).replace(" percentage points", ""),
        )
        updated.add(destination)
    if updated != set(targets.values()):
        raise ValueError(f"main estimates destination rows missing: {sorted(set(targets.values()) - updated)}")
    write_csv(OUT_MAIN_TABLE, rows, ["Analysis", "Estimate", "95% CI or component", "Interpretation"])


def percentile(values: np.ndarray, q: float) -> float:
    return finite_float(np.quantile(values, q))


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def kitagawa_from_counts(a0: np.ndarray, b0: np.ndarray, a1: np.ndarray, b1: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized symmetric Kitagawa components for one or more replicates."""
    total_a0 = a0.sum(axis=1)
    total_a1 = a1.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        w0 = a0 / total_a0[:, None]
        w1 = a1 / total_a1[:, None]
        p0 = np.divide(b0, a0, out=np.zeros_like(b0, dtype=float), where=a0 > 0)
        p1 = np.divide(b1, a1, out=np.zeros_like(b1, dtype=float), where=a1 > 0)
    composition = 0.5 * ((p0 + p1) * (w1 - w0)).sum(axis=1)
    within = 0.5 * ((w0 + w1) * (p1 - p0)).sum(axis=1)
    crude = b1.sum(axis=1) / total_a1 - b0.sum(axis=1) / total_a0
    return crude, composition, within


def group_joint_probabilities(annual: dict, year: str, group: str, strata: list[str]) -> tuple[int, np.ndarray]:
    data = annual[year]["groups"][group]
    n_total = int(annual[year]["main"]["A"])
    counts: list[int] = []
    for stratum in strata:
        values = data.get(stratum, {"A": 0, "B": 0})
        a = int(values["A"])
        b = int(values["B"])
        if b < 0 or a < b:
            raise ValueError("invalid group counts: {} {} {}".format(year, group, stratum))
        counts.extend((b, a - b))
    if sum(counts) != n_total:
        raise ValueError("group counts do not sum to A: {} {}".format(year, group))
    probabilities = np.asarray(counts, dtype=float) / n_total
    return n_total, probabilities


def bootstrap_kitagawa(annual: dict, point_rows: dict[str, dict[str, str]], rng: np.random.Generator) -> tuple[list[dict[str, object]], dict[str, object]]:
    output: list[dict[str, object]] = []
    diagnostics: dict[str, object] = {}
    for group in GROUPS:
        strata = sorted(set(annual["1999"]["groups"][group]) | set(annual["2024"]["groups"][group]))
        n0, prob0 = group_joint_probabilities(annual, "1999", group, strata)
        n1, prob1 = group_joint_probabilities(annual, "2024", group, strata)
        draw0 = rng.multinomial(n0, prob0, size=REPS)
        draw1 = rng.multinomial(n1, prob1, size=REPS)
        b0, a0_not_b = draw0[:, 0::2], draw0[:, 1::2]
        b1, a1_not_b = draw1[:, 0::2], draw1[:, 1::2]
        a0, a1 = b0 + a0_not_b, b1 + a1_not_b
        crude_boot, composition_boot, within_boot = kitagawa_from_counts(a0, b0, a1, b1)
        closure = crude_boot - composition_boot - within_boot

        point_a0 = np.asarray([[int(annual["1999"]["groups"][group].get(s, {"A": 0})["A"]) for s in strata]], dtype=float)
        point_b0 = np.asarray([[int(annual["1999"]["groups"][group].get(s, {"B": 0})["B"]) for s in strata]], dtype=float)
        point_a1 = np.asarray([[int(annual["2024"]["groups"][group].get(s, {"A": 0})["A"]) for s in strata]], dtype=float)
        point_b1 = np.asarray([[int(annual["2024"]["groups"][group].get(s, {"B": 0})["B"]) for s in strata]], dtype=float)
        crude_point, composition_point, within_point = kitagawa_from_counts(point_a0, point_b0, point_a1, point_b1)
        baseline = point_rows[group]
        existing_composition = float(baseline["composition_component"])
        existing_within = float(baseline["selection_component"])
        existing_crude = float(baseline["UCF_change_2024_minus_1999"])
        point_residual = finite_float(crude_point[0] - composition_point[0] - within_point[0])
        if max(abs(composition_point[0] - existing_composition), abs(within_point[0] - existing_within), abs(crude_point[0] - existing_crude)) > 1e-12:
            raise ValueError("Kitagawa point estimate does not reproduce frozen table for {}".format(group))

        covariance = np.cov(np.vstack((composition_boot, within_boot)), ddof=1)
        row = {
            "group": group,
            "bootstrap_reps": REPS,
            "seed": SEED,
            "n_1999": n0,
            "n_2024": n1,
            "crude_change_point": finite_float(crude_point[0]),
            "composition_point": finite_float(composition_point[0]),
            "composition_bootstrap_se": finite_float(np.std(composition_boot, ddof=1)),
            "composition_percentile_ci_lo": percentile(composition_boot, ALPHA / 2),
            "composition_percentile_ci_hi": percentile(composition_boot, 1 - ALPHA / 2),
            "within_stratum_point": finite_float(within_point[0]),
            "within_stratum_bootstrap_se": finite_float(np.std(within_boot, ddof=1)),
            "within_stratum_percentile_ci_lo": percentile(within_boot, ALPHA / 2),
            "within_stratum_percentile_ci_hi": percentile(within_boot, 1 - ALPHA / 2),
            "component_covariance": finite_float(covariance[0, 1]),
            "component_correlation": finite_float(covariance[0, 1] / math.sqrt(covariance[0, 0] * covariance[1, 1])),
            "component_sum_point": finite_float(composition_point[0] + within_point[0]),
            "component_sum_bootstrap_se": finite_float(np.std(composition_boot + within_boot, ddof=1)),
            "component_sum_percentile_ci_lo": percentile(composition_boot + within_boot, ALPHA / 2),
            "component_sum_percentile_ci_hi": percentile(composition_boot + within_boot, 1 - ALPHA / 2),
            "point_closure_residual": point_residual,
            "bootstrap_closure_mean_residual": finite_float(np.mean(closure)),
            "bootstrap_closure_max_abs_residual": finite_float(np.max(np.abs(closure))),
            "frozen_point_max_abs_difference": finite_float(max(abs(composition_point[0] - existing_composition), abs(within_point[0] - existing_within), abs(crude_point[0] - existing_crude))),
            "method": "Endpoint-year joint multinomial bootstrap over (stratum, UCD-K80 status); symmetric Kitagawa re-estimated jointly each replicate; percentile 95% CI.",
        }
        output.append(row)
        diagnostics[group] = {
            "strata": strata,
            "point_closure_residual": row["point_closure_residual"],
            "bootstrap_closure_max_abs_residual": row["bootstrap_closure_max_abs_residual"],
            "component_covariance": row["component_covariance"],
        }
    return output, diagnostics


def destination_weights() -> dict[str, float]:
    rows = read_csv(DEST_WEIGHT_PATH)
    weights = {row["age_sex_cell"]: float(row["weight"]) for row in rows}
    if abs(sum(weights.values()) - 1.0) > 1e-12:
        raise ValueError("destination weights do not sum to one")
    return weights


def endpoint_destination_probability(annual: dict, year: str, categories: list[str], weights: dict[str, float]) -> np.ndarray:
    value = np.zeros(len(categories), dtype=float)
    for cell, weight in weights.items():
        counts = annual[year]["dest_cells"].get(cell)
        if counts is None:
            raise ValueError("missing destination cell {} in {}".format(cell, year))
        total = sum(int(counts.get(category, 0)) for category in categories)
        if total <= 0:
            raise ValueError("zero destination denominator in weighted cell {} {}".format(cell, year))
        value += weight * np.asarray([int(counts.get(category, 0)) / total for category in categories])
    if abs(value.sum() - 1.0) > 1e-12:
        raise ValueError("point destination probabilities do not sum to one")
    return value


def bootstrap_destination_probabilities(annual: dict, year: str, categories: list[str], weights: dict[str, float], rng: np.random.Generator) -> np.ndarray:
    result = np.zeros((REPS, len(categories)), dtype=float)
    for cell, weight in weights.items():
        counts = annual[year]["dest_cells"].get(cell)
        vector = np.asarray([int(counts.get(category, 0)) for category in categories], dtype=int)
        total = int(vector.sum())
        if total <= 0:
            raise ValueError("zero destination denominator in weighted cell {} {}".format(cell, year))
        draw = rng.multinomial(total, vector.astype(float) / total, size=REPS)
        result += weight * (draw.astype(float) / total)
    return result


def normal_two_sided_p(z: float) -> float:
    return finite_float(math.erfc(abs(z) / math.sqrt(2.0)))


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    m = len(p_values)
    order = np.argsort(p_values)
    sorted_p = p_values[order]
    adjusted = sorted_p * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    q = np.empty(m, dtype=float)
    q[order] = np.minimum(adjusted, 1.0)
    return q


def bootstrap_destinations(annual: dict, frozen_rows: dict[str, dict[str, str]], rng: np.random.Generator) -> tuple[list[dict[str, object]], dict[str, object]]:
    categories = sorted(frozen_rows)
    weights = destination_weights()
    point_1999 = endpoint_destination_probability(annual, "1999", categories, weights)
    point_2024 = endpoint_destination_probability(annual, "2024", categories, weights)
    point_rd = point_2024 - point_1999
    for index, category in enumerate(categories):
        if abs(point_rd[index] - float(frozen_rows[category]["RD_2024_minus_1999"])) > 1e-12:
            raise ValueError("destination point estimate does not reproduce frozen table: {}".format(category))

    boot_1999 = bootstrap_destination_probabilities(annual, "1999", categories, weights, rng)
    boot_2024 = bootstrap_destination_probabilities(annual, "2024", categories, weights, rng)
    boot_rd = boot_2024 - boot_1999
    probability_sum_error = max(float(np.max(np.abs(boot_1999.sum(axis=1) - 1.0))), float(np.max(np.abs(boot_2024.sum(axis=1) - 1.0))))
    rd_sum_error = float(np.max(np.abs(boot_rd.sum(axis=1))))
    covariance = np.cov(boot_rd, rowvar=False, ddof=1)
    se = np.sqrt(np.diag(covariance))
    if np.any(se <= 0) or not np.all(np.isfinite(se)):
        raise ValueError("invalid destination bootstrap standard error")
    centered_t = (boot_rd - point_rd[None, :]) / se[None, :]
    max_t_critical = percentile(np.max(np.abs(centered_t), axis=1), 1 - ALPHA)
    p_values = np.asarray([normal_two_sided_p(point_rd[i] / se[i]) for i in range(len(categories))])
    q_values = benjamini_hochberg(p_values)
    cov_row_sums = covariance.sum(axis=1)
    cov_col_sums = covariance.sum(axis=0)
    output: list[dict[str, object]] = []
    covariance_columns = ["covariance_{}".format(category) for category in categories]
    for i, category in enumerate(categories):
        frozen = frozen_rows[category]
        row: dict[str, object] = {
            "destination": category,
            "bootstrap_reps": REPS,
            "seed": SEED,
            "std_probability_1999_point": finite_float(point_1999[i]),
            "std_probability_2024_point": finite_float(point_2024[i]),
            "RD_2024_minus_1999_point": finite_float(point_rd[i]),
            "nominal_ci_lo_frozen": finite_float(float(frozen["RD_CI_lo"])),
            "nominal_ci_hi_frozen": finite_float(float(frozen["RD_CI_hi"])),
            "bootstrap_percentile_ci_lo": percentile(boot_rd[:, i], ALPHA / 2),
            "bootstrap_percentile_ci_hi": percentile(boot_rd[:, i], 1 - ALPHA / 2),
            "bootstrap_se": finite_float(se[i]),
            "normal_reference_two_sided_p": finite_float(p_values[i]),
            "bh_fdr_q": finite_float(q_values[i]),
            "simultaneous_95_ci_lo": finite_float(point_rd[i] - max_t_critical * se[i]),
            "simultaneous_95_ci_hi": finite_float(point_rd[i] + max_t_critical * se[i]),
            "simultaneous_excludes_zero": bool((point_rd[i] - max_t_critical * se[i] > 0) or (point_rd[i] + max_t_critical * se[i] < 0)),
            "max_t_critical_95": finite_float(max_t_critical),
            "covariance_row_sum": finite_float(cov_row_sums[i]),
            "covariance_column_sum": finite_float(cov_col_sums[i]),
            "point_estimate_abs_difference_frozen": finite_float(abs(point_rd[i] - float(frozen["RD_2024_minus_1999"]))),
            "inference_note": "Joint stratified multinomial bootstrap preserves nine-category compositional dependence; 95% FWER max-t interval across all nine destination RDs. Nominal frozen CI is retained separately and is not used for category selection.",
        }
        for column, value in zip(covariance_columns, covariance[i]):
            row[column] = finite_float(value)
        output.append(row)
    diagnostics = {
        "categories": categories,
        "max_t_critical_95": finite_float(max_t_critical),
        "max_probability_sum_abs_error": finite_float(probability_sum_error),
        "max_rd_sum_abs_error": finite_float(rd_sum_error),
        "max_abs_covariance_row_sum": finite_float(np.max(np.abs(cov_row_sums))),
        "max_abs_covariance_column_sum": finite_float(np.max(np.abs(cov_col_sums))),
    }
    return output, diagnostics


def make_report(decomp_rows: list[dict[str, object]], dest_rows: list[dict[str, object]], diagnostics: dict[str, object]) -> str:
    retained = [row["destination"] for row in dest_rows if row["simultaneous_excludes_zero"]]
    lines = [
        "# Joint inference extension for decomposition and destination contrasts",
        "",
        "## Scope",
        "",
        "This extension leaves the frozen A/B/UCF definitions and all primary point estimates unchanged. It adds endpoint-year bootstrap uncertainty for each of four alternative Kitagawa partitions and joint compositional inference for the nine mutually exclusive destination RDs. The four Kitagawa partitions are separate descriptive accounting decompositions and must not be added together or interpreted as causal mediation effects.",
        "",
        "## Bootstrap design",
        "",
        "- Random seed: `{}`; replicates: `{:,}` for each analysis.".format(SEED, REPS),
        "- Kitagawa: in each endpoint year, the observed joint `(stratum, UCD-K80 status)` table was resampled multinomially conditional on the observed number of A records. Composition weights and stratum-specific UCFs were recomputed together before applying the symmetric identity. The reported component covariance comes from the same replicate distribution; the component sum was checked against the resampled crude endpoint difference on every replicate.",
        "- Destinations: within every fixed-weight age-sex cell and endpoint year, the full nine-category destination vector was resampled multinomially. Standardized probabilities were then recomputed with the frozen gap weights. This preserves the probability-sum constraint and the negative cross-category covariance. Simultaneous intervals use the 95th percentile of the joint max absolute studentized bootstrap deviation (max-t), producing family-wise 95% coverage for the nine RD intervals.",
        "- Nominal intervals retained in the destination output are the pre-existing fixed-weight independent-year intervals for reference only. They are not used to select categories. `normal_reference_two_sided_p` and BH q values are supplementary, derived from the joint-bootstrap standard errors; the FWER max-t intervals are the primary multiplicity-controlled inference.",
        "",
        "## Kitagawa results",
        "",
        "| Partition | Composition (95% percentile CI) | Within-stratum (95% percentile CI) | Covariance | Max bootstrap closure residual |",
        "|---|---|---|---:|---:|",
    ]
    for row in decomp_rows:
        lines.append("| {group} | {composition_point:.6f} ({composition_percentile_ci_lo:.6f}, {composition_percentile_ci_hi:.6f}) | {within_stratum_point:.6f} ({within_stratum_percentile_ci_lo:.6f}, {within_stratum_percentile_ci_hi:.6f}) | {component_covariance:.8g} | {bootstrap_closure_max_abs_residual:.3g} |".format(**row))
    lines += [
        "",
        "## Destination results",
        "",
        "The joint max-t critical value was `{:.6f}`. Destinations whose simultaneous 95% RD interval excludes zero: **{}**.".format(diagnostics["destinations"]["max_t_critical_95"], ", ".join(retained) if retained else "none"),
        "",
        "| Destination | RD | FWER simultaneous 95% CI | BH q | Excludes 0 simultaneously? |",
        "|---|---:|---|---:|---|",
    ]
    for row in dest_rows:
        lines.append("| {destination} | {RD_2024_minus_1999_point:.6f} | ({simultaneous_95_ci_lo:.6f}, {simultaneous_95_ci_hi:.6f}) | {bh_fdr_q:.6g} | {simultaneous_excludes_zero} |".format(**row))
    lines += [
        "",
        "## Constraint checks and limitations",
        "",
        "- Maximum absolute replicate probability-sum error: `{:.3g}`; maximum absolute replicate RD-sum error: `{:.3g}`.".format(diagnostics["destinations"]["max_probability_sum_abs_error"], diagnostics["destinations"]["max_rd_sum_abs_error"]),
        "- Maximum absolute destination RD covariance row/column sums: `{:.3g}` / `{:.3g}`; zero is expected from the compositional constraint up to floating-point precision.".format(diagnostics["destinations"]["max_abs_covariance_row_sum"], diagnostics["destinations"]["max_abs_covariance_column_sum"]),
        "- These are conditional resampling inferences for the observed endpoint record distributions and frozen standard weights. They quantify sampling stability of descriptive certificate/coding patterns; they do not create causal, clinical-validity, or biological-substitution interpretations.",
        "",
        "Machine-readable outputs: `tables/k80_kitagawa_decomposition_uncertainty.csv`, `tables/ucd_destination_contrasts_simultaneous.csv`, and `logs/inference_extension_run.json`.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    global ROOT, ANNUAL_PATH, DECOMP_POINT_PATH, DEST_POINT_PATH, DEST_WEIGHT_PATH, OUT_DECOMP, OUT_DEST, OUT_REPORT, OUT_LOG, OUT_MAIN_TABLE
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-root", default="results", help="output directory generated by analysis_reanalysis_public.py")
    args = parser.parse_args()
    ROOT = Path(args.analysis_root); ANNUAL_PATH = ROOT / "analysis" / "annual_aggregates.json"; DECOMP_POINT_PATH = ROOT / "tables" / "k80_kitagawa_decomposition.csv"; DEST_POINT_PATH = ROOT / "tables" / "ucd_destination_contrasts.csv"; DEST_WEIGHT_PATH = ROOT / "tables" / "ucd_destination_standardization_weights.csv"; OUT_DECOMP = ROOT / "tables" / "k80_kitagawa_decomposition_uncertainty.csv"; OUT_DEST = ROOT / "tables" / "ucd_destination_contrasts_simultaneous.csv"; OUT_REPORT = ROOT / "reports" / "inference_extension_report.md"; OUT_LOG = ROOT / "logs" / "inference_extension_run.json"; OUT_MAIN_TABLE = ROOT / "tables" / "Table2_Main_estimates.csv"
    annual = json.loads(ANNUAL_PATH.read_text(encoding="utf-8"))
    point_decomp = {row["group"]: row for row in read_csv(DECOMP_POINT_PATH)}
    point_dest = {row["destination"]: row for row in read_csv(DEST_POINT_PATH)}
    if set(point_decomp) != set(GROUPS):
        raise ValueError("unexpected Kitagawa groups")
    if len(point_dest) != 9:
        raise ValueError("expected nine destination categories")
    rng = np.random.default_rng(SEED)
    decomp_rows, decomp_diag = bootstrap_kitagawa(annual, point_decomp, rng)
    dest_rows, dest_diag = bootstrap_destinations(annual, point_dest, rng)
    decomp_fields = list(decomp_rows[0])
    destination_fields = list(dest_rows[0])
    write_csv(OUT_DECOMP, decomp_rows, decomp_fields)
    write_csv(OUT_DEST, dest_rows, destination_fields)
    update_main_estimates(dest_rows)
    diagnostics = {
        "seed": SEED,
        "bootstrap_reps": REPS,
        "input_sha256": {
            "analysis/annual_aggregates.json": file_sha256(ANNUAL_PATH),
            "tables/k80_kitagawa_decomposition.csv": file_sha256(DECOMP_POINT_PATH),
            "tables/ucd_destination_contrasts.csv": file_sha256(DEST_POINT_PATH),
            "tables/ucd_destination_standardization_weights.csv": file_sha256(DEST_WEIGHT_PATH),
        },
        "kitagawa": decomp_diag,
        "destinations": dest_diag,
        "output_sha256": {
            "tables/k80_kitagawa_decomposition_uncertainty.csv": file_sha256(OUT_DECOMP),
            "tables/ucd_destination_contrasts_simultaneous.csv": file_sha256(OUT_DEST),
            "tables/Table2_Main_estimates.csv": file_sha256(OUT_MAIN_TABLE),
        },
        "verification": {
            "all_finite": True,
            "kitagawa_points_reproduce_frozen": True,
            "destination_points_reproduce_frozen": True,
            "destination_joint_compositional_bootstrap": True,
            "fwer_simultaneous_max_t": True,
        },
    }
    OUT_REPORT.write_text(make_report(decomp_rows, dest_rows, diagnostics), encoding="utf-8")
    diagnostics["output_sha256"]["reports/inference_extension_report.md"] = file_sha256(OUT_REPORT)
    OUT_LOG.write_text(json.dumps(diagnostics, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
