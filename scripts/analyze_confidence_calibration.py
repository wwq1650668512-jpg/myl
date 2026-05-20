from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.web.service import GutPredictionService


BENCHMARK_RESULTS_CSV = ROOT / "predictions" / "evaluation" / "fusion_comparison" / "case_based_results.csv"
DEFAULT_OUTPUT_DIR = ROOT / "predictions" / "evaluation" / "confidence"

FIXED_CASES: list[dict[str, str]] = [
    {
        "drug_name": "Rifaximin",
        "smiles": "CC1=C(C(=O)NC(=C1C)C2=CC3=C(C=C2O)C(=O)C4=C(C3=O)C(=C(C=C4OC)O)OC)C",
    },
    {
        "drug_name": "Vancomycin",
        "smiles": "C[C@H]1[C@H]([C@@](C[C@@H](O1)O[C@@H]2[C@H]([C@@H]([C@H](O[C@H]2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)[C@H]([C@H](C(=O)N[C@H](C(=O)N[C@H]5C(=O)N[C@@H]7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)[C@H](NC(=O)[C@H]([C@@H](C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)[C@@H](CC(C)C)NC)O)Cl)CO)O)O)(C)N)O",
    },
    {
        "drug_name": "Lubiprostone",
        "smiles": "CC(C)CCCC(C)C1CCC2C1(CCC3C2CCC4=CC(=O)CCC34C)C",
    },
    {
        "drug_name": "Metronidazole",
        "smiles": "CN1C=NC(=N1)COC2=NC(=O)N(C2=O)C",
    },
]

STATUS_NUMERIC = {"FAIL": 0.0, "PARTIAL": 0.5, "PASS": 1.0}
SEVERE_WARNINGS = {"over-suppression", "core-butyrate-suppression", "ecology-risk", "drug-profile-conflict"}


def _safe_float(value: object) -> float | None:
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(numeric):
        return None
    return float(numeric)


def _predict_fixed_case_confidence() -> pd.DataFrame:
    service = GutPredictionService()
    rows: list[dict[str, object]] = []
    for case in FIXED_CASES:
        payload = service.predict_custom_drug(drug_name=case["drug_name"], smiles=case["smiles"])
        warnings = payload.get("warning_flags", [])
        rows.append(
            {
                "Drug": case["drug_name"],
                "confidence_score": _safe_float(payload.get("confidence_score")),
                "confidence_tier": str(payload.get("confidence_tier") or ""),
                "warning_flags": json.dumps(warnings, ensure_ascii=False),
                "warning_count": len(warnings) if isinstance(warnings, list) else 0,
                "severe_warning_count": len([item for item in warnings if item in SEVERE_WARNINGS])
                if isinstance(warnings, list)
                else 0,
                "confidence_explanation": str(payload.get("confidence_explanation") or ""),
            }
        )
    return pd.DataFrame(rows)


def _load_benchmark_results(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Benchmark results not found: {path}")
    frame = pd.read_csv(path)
    required = {"assertion_id", "Drug", "Status"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"Benchmark results missing columns: {missing}")
    return frame.loc[:, ["assertion_id", "Drug", "Status"]].copy()


def _build_calibration_summary(assertion_level: pd.DataFrame) -> pd.DataFrame:
    tier_order = ["high", "medium", "low"]
    distribution = (
        assertion_level.groupby(["confidence_tier", "Status"], dropna=False)
        .size()
        .unstack(fill_value=0)
        .reindex(tier_order, fill_value=0)
        .reindex(columns=["PASS", "PARTIAL", "FAIL"], fill_value=0)
        .rename_axis("confidence_tier")
        .reset_index()
    )
    distribution["total"] = distribution[["PASS", "PARTIAL", "FAIL"]].sum(axis=1)
    for status in ["PASS", "PARTIAL", "FAIL"]:
        distribution[f"{status.lower()}_rate"] = distribution[status] / distribution["total"].replace(0, pd.NA)

    distribution_rows = distribution.assign(
        row_type="tier_distribution",
        metric="",
        value=pd.NA,
        n=distribution["total"],
    )[
        [
            "row_type",
            "confidence_tier",
            "PASS",
            "PARTIAL",
            "FAIL",
            "total",
            "pass_rate",
            "partial_rate",
            "fail_rate",
            "metric",
            "value",
            "n",
        ]
    ]

    pearson = assertion_level["confidence_score"].corr(assertion_level["correctness_score"], method="pearson")
    spearman = assertion_level["confidence_score"].corr(assertion_level["correctness_score"], method="spearman")

    correlation_rows = pd.DataFrame(
        [
            {
                "row_type": "correlation",
                "confidence_tier": "",
                "PASS": pd.NA,
                "PARTIAL": pd.NA,
                "FAIL": pd.NA,
                "total": pd.NA,
                "pass_rate": pd.NA,
                "partial_rate": pd.NA,
                "fail_rate": pd.NA,
                "metric": "pearson_confidence_vs_correctness",
                "value": pearson,
                "n": len(assertion_level),
            },
            {
                "row_type": "correlation",
                "confidence_tier": "",
                "PASS": pd.NA,
                "PARTIAL": pd.NA,
                "FAIL": pd.NA,
                "total": pd.NA,
                "pass_rate": pd.NA,
                "partial_rate": pd.NA,
                "fail_rate": pd.NA,
                "metric": "spearman_confidence_vs_correctness",
                "value": spearman,
                "n": len(assertion_level),
            },
        ]
    )
    return pd.concat([distribution_rows, correlation_rows], ignore_index=True)


def _build_alignment_summary(assertion_level: pd.DataFrame) -> pd.DataFrame:
    work = assertion_level.copy()
    work["low_or_warning"] = work["confidence_tier"].eq("low") | (work["warning_count"] > 0)
    work["no_severe_warning"] = work["severe_warning_count"].eq(0)

    rows: list[dict[str, object]] = []
    for status in ["FAIL", "PARTIAL", "PASS"]:
        subset = work[work["Status"] == status]
        if subset.empty:
            continue
        rows.append(
            {
                "row_type": "status_profile",
                "check_id": f"{status}_profile",
                "benchmark_status": status,
                "n_cases": int(len(subset)),
                "match_count": pd.NA,
                "match_rate": pd.NA,
                "low_confidence_rate": float(subset["confidence_tier"].eq("low").mean()),
                "any_warning_rate": float((subset["warning_count"] > 0).mean()),
                "severe_warning_rate": float((subset["severe_warning_count"] > 0).mean()),
                "definition": "diagnostic rates per benchmark status",
            }
        )

    fail_subset = work[work["Status"] == "FAIL"]
    fail_matched = fail_subset["low_or_warning"] if not fail_subset.empty else pd.Series(dtype=bool)
    rows.append(
        {
            "row_type": "alignment_check",
            "check_id": "fail_should_have_low_or_warning",
            "benchmark_status": "FAIL",
            "n_cases": int(len(fail_subset)),
            "match_count": int(fail_matched.sum()) if not fail_subset.empty else 0,
            "match_rate": float(fail_matched.mean()) if not fail_subset.empty else pd.NA,
            "low_confidence_rate": pd.NA,
            "any_warning_rate": pd.NA,
            "severe_warning_rate": pd.NA,
            "definition": "confidence_tier == low OR warning_count > 0",
        }
    )

    pass_subset = work[work["Status"] == "PASS"]
    pass_matched = pass_subset["no_severe_warning"] if not pass_subset.empty else pd.Series(dtype=bool)
    rows.append(
        {
            "row_type": "alignment_check",
            "check_id": "pass_should_not_have_severe_warning",
            "benchmark_status": "PASS",
            "n_cases": int(len(pass_subset)),
            "match_count": int(pass_matched.sum()) if not pass_subset.empty else 0,
            "match_rate": float(pass_matched.mean()) if not pass_subset.empty else pd.NA,
            "low_confidence_rate": pd.NA,
            "any_warning_rate": pd.NA,
            "severe_warning_rate": pd.NA,
            "definition": "severe_warning_count == 0",
        }
    )

    return pd.DataFrame(rows)


def _fmt(value: object, digits: int = 3) -> str:
    if value is None:
        return "N/A"
    if isinstance(value, str):
        return value
    numeric = _safe_float(value)
    if numeric is None or math.isnan(numeric):
        return "N/A"
    return f"{numeric:.{digits}f}"


def _write_analysis_markdown(
    *,
    output_path: Path,
    assertion_level: pd.DataFrame,
    calibration_summary: pd.DataFrame,
    alignment_summary: pd.DataFrame,
) -> None:
    tier_table = calibration_summary[calibration_summary["row_type"] == "tier_distribution"].copy()
    corr_table = calibration_summary[calibration_summary["row_type"] == "correlation"].copy()
    fail_check = alignment_summary[alignment_summary["check_id"] == "fail_should_have_low_or_warning"]
    pass_check = alignment_summary[alignment_summary["check_id"] == "pass_should_not_have_severe_warning"]

    lines: list[str] = []
    lines.append("# Confidence Calibration Analysis")
    lines.append("")
    lines.append("## Scope")
    lines.append("")
    lines.append(f"- Benchmark assertions: {len(assertion_level)}")
    lines.append(f"- Fixed test drugs: {assertion_level['Drug'].nunique()}")
    lines.append("- Correctness mapping: PASS=1.0, PARTIAL=0.5, FAIL=0.0")
    lines.append("- Severe warning set: over-suppression, core-butyrate-suppression, ecology-risk, drug-profile-conflict")
    lines.append("")
    lines.append("## Tier Distribution (PASS/PARTIAL/FAIL)")
    lines.append("")
    if tier_table.empty:
        lines.append("- No tier distribution rows.")
    else:
        for _, row in tier_table.iterrows():
            lines.append(
                f"- {row['confidence_tier']}: PASS={int(row['PASS'])}, PARTIAL={int(row['PARTIAL'])}, FAIL={int(row['FAIL'])}, "
                f"pass_rate={_fmt(row['pass_rate'])}, fail_rate={_fmt(row['fail_rate'])}"
            )
    lines.append("")
    lines.append("## Correlation With Benchmark Correctness")
    lines.append("")
    if corr_table.empty:
        lines.append("- Correlation rows unavailable.")
    else:
        for _, row in corr_table.iterrows():
            lines.append(f"- {row['metric']}: {_fmt(row['value'])} (n={int(row['n'])})")
    lines.append("")
    lines.append("## Benchmark Alignment Checks")
    lines.append("")
    if fail_check.empty:
        lines.append("- FAIL check unavailable.")
    else:
        row = fail_check.iloc[0]
        lines.append(
            f"- FAIL 对齐（低置信或有 warning）: {int(row['match_count'])}/{int(row['n_cases'])} = {_fmt(row['match_rate'])}"
        )
    if pass_check.empty:
        lines.append("- PASS check unavailable.")
    else:
        row = pass_check.iloc[0]
        lines.append(
            f"- PASS 对齐（无严重 warning）: {int(row['match_count'])}/{int(row['n_cases'])} = {_fmt(row['match_rate'])}"
        )
    lines.append("")
    lines.append("## Quick Read")
    lines.append("")
    lines.append("- 如果 low tier 的 FAIL 占比高、high tier 的 PASS 占比高，说明分层有区分度。")
    lines.append("- 如果相关系数为正，说明 confidence_score 趋势与 benchmark correctness 同向。")
    lines.append("- 对齐检查可直接监控“FAIL 是否被提醒、PASS 是否不过度报警”。")
    lines.append("")
    output_path.write_text("\n".join(lines), encoding="utf-8")


def run(benchmark_results_csv: Path, output_dir: Path) -> dict[str, str]:
    benchmark = _load_benchmark_results(benchmark_results_csv)
    confidence_by_drug = _predict_fixed_case_confidence()
    assertion_level = benchmark.merge(confidence_by_drug, on="Drug", how="left", validate="many_to_one")
    assertion_level = assertion_level[assertion_level["Status"].isin({"PASS", "PARTIAL", "FAIL"})].copy()
    assertion_level["correctness_score"] = assertion_level["Status"].map(STATUS_NUMERIC)

    calibration_summary = _build_calibration_summary(assertion_level)
    alignment_summary = _build_alignment_summary(assertion_level)

    output_dir.mkdir(parents=True, exist_ok=True)
    calibration_csv = output_dir / "confidence_calibration_summary.csv"
    alignment_csv = output_dir / "confidence_benchmark_alignment.csv"
    analysis_md = output_dir / "confidence_calibration_analysis.md"

    calibration_summary.to_csv(calibration_csv, index=False)
    alignment_summary.to_csv(alignment_csv, index=False)
    _write_analysis_markdown(
        output_path=analysis_md,
        assertion_level=assertion_level,
        calibration_summary=calibration_summary,
        alignment_summary=alignment_summary,
    )

    return {
        "confidence_calibration_summary.csv": str(calibration_csv),
        "confidence_benchmark_alignment.csv": str(alignment_csv),
        "confidence_calibration_analysis.md": str(analysis_md),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run confidence calibration/alignment analysis using fixed benchmark cases.")
    parser.add_argument("--benchmark-results", type=Path, default=BENCHMARK_RESULTS_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    result = run(args.benchmark_results, args.output_dir)
    print(json.dumps({"status": "ok", **result}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
