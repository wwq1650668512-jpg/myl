from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.web.service import GutPredictionService


DEFAULT_DRUG_QUERIES = ["Rifaximin", "Vancomycin hydrochloride"]


def _safe_float(value: object) -> float:
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    return float(numeric) if pd.notna(numeric) else float("nan")


def _extract_disease_table(profile: dict[str, object]) -> pd.DataFrame:
    rows = profile.get("candidate_diseases", [])
    if not rows:
        return pd.DataFrame(
            columns=[
                "disease_name",
                "disease_score_raw_only",
                "disease_score_mechanism",
                "mechanism_delta",
            ]
        )
    frame = pd.DataFrame(rows)
    frame["disease_score_raw_only"] = pd.to_numeric(frame["disease_score_raw_only"], errors="coerce")
    frame["disease_score_mechanism"] = pd.to_numeric(frame["disease_score_mechanism"], errors="coerce")
    frame["mechanism_delta"] = pd.to_numeric(frame["mechanism_delta"], errors="coerce")
    frame["raw_rank"] = frame["disease_score_raw_only"].rank(method="min", ascending=False).astype(int)
    frame["mechanism_rank"] = frame["disease_score_mechanism"].rank(method="min", ascending=False).astype(int)
    frame["rank_shift"] = frame["raw_rank"] - frame["mechanism_rank"]

    def _mechanism_value(row: pd.Series, key: str) -> float:
        payload = row.get("mechanism_scores")
        if isinstance(payload, dict):
            return _safe_float(payload.get(key))
        return float("nan")

    for key in [
        "anti_inflammatory_score",
        "pro_inflammatory_score",
        "butyrate_support_score",
        "barrier_protection_score",
        "toxin_risk_score",
        "mucus_degradation_score",
        "pathobiont_load",
        "competition_vs_crossfeeding_proxy",
        "mechanism_benefit_score",
        "mechanism_risk_score",
        "mechanism_balance_score",
    ]:
        frame[key] = frame.apply(lambda row: _mechanism_value(row, key), axis=1)

    return frame.sort_values("disease_score_mechanism", ascending=False)


def _write_drug_markdown(
    output_path: Path,
    drug_query: str,
    profile: dict[str, object],
    table: pd.DataFrame,
) -> None:
    lines: list[str] = []
    lines.append(f"# Case Study: {drug_query}")
    lines.append("")
    drug_meta = profile.get("drug", {}) if isinstance(profile.get("drug"), dict) else {}
    lines.append(f"- Prestwick: {drug_meta.get('prestwick_id', 'N/A')}")
    lines.append(f"- Chemical name: {drug_meta.get('chemical_name', 'N/A')}")
    lines.append(f"- Candidate diseases: {len(table)}")
    lines.append("")
    if table.empty:
        lines.append("No candidate disease rows were produced.")
        output_path.write_text("\n".join(lines), encoding="utf-8")
        return

    view_columns = [
        "disease_name",
        "disease_score_raw_only",
        "disease_score_mechanism",
        "mechanism_delta",
        "raw_rank",
        "mechanism_rank",
        "rank_shift",
        "anti_inflammatory_score",
        "pro_inflammatory_score",
        "butyrate_support_score",
        "barrier_protection_score",
        "toxin_risk_score",
        "mucus_degradation_score",
        "pathobiont_load",
        "competition_vs_crossfeeding_proxy",
        "mechanism_balance_score",
    ]
    top = table.loc[:, [column for column in view_columns if column in table.columns]].head(8)
    lines.append("## Top Diseases (Mechanism Layer)")
    lines.append("")
    lines.append(top.to_markdown(index=False))
    lines.append("")

    top_row = table.iloc[0]
    lines.append(f"## Top Disease: {top_row['disease_name']}")
    lines.append("")
    contributor_payload = top_row.get("mechanism_top_contributors", {})
    if isinstance(contributor_payload, dict):
        for mechanism_key in [
            "anti_inflammatory_score",
            "pro_inflammatory_score",
            "toxin_risk_score",
            "pathobiont_load",
        ]:
            rows = contributor_payload.get(mechanism_key, [])
            if not rows:
                continue
            lines.append(f"### {mechanism_key}")
            lines.append("")
            contributor_table = pd.DataFrame(rows).head(5)
            lines.append(contributor_table.to_markdown(index=False))
            lines.append("")

    output_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate disease-ranking ablation for mechanism-layer scoring.")
    parser.add_argument(
        "--integrated-predictions",
        default=ROOT / "predictions/step2/baseline_scaffold_v1_83/predictions.csv",
        type=Path,
        help="Integrated Step1+Step2 prediction table used by the web service.",
    )
    parser.add_argument(
        "--drug-query",
        action="append",
        default=[],
        help="Drug query used for case study. Repeatable. Defaults to Rifaximin and Vancomycin hydrochloride.",
    )
    parser.add_argument(
        "--output-dir",
        default=ROOT / "predictions/evaluation/mechanism_layer_case_study",
        type=Path,
        help="Directory for case-study outputs.",
    )
    args = parser.parse_args()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    queries = [item.strip() for item in args.drug_query if item and item.strip()] or list(DEFAULT_DRUG_QUERIES)

    service = GutPredictionService(integrated_predictions_path=args.integrated_predictions)
    summary_rows: list[dict[str, object]] = []
    for query in queries:
        profile = service.get_drug_profile(query)
        table = _extract_disease_table(profile)
        safe_name = "".join(ch.lower() if ch.isalnum() else "_" for ch in query).strip("_") or "drug"

        table_path = output_dir / f"{safe_name}_disease_ablation.csv"
        table.to_csv(table_path, index=False)
        json_path = output_dir / f"{safe_name}_profile.json"
        json_path.write_text(json.dumps(profile, ensure_ascii=False, indent=2), encoding="utf-8")
        md_path = output_dir / f"{safe_name}_case_study.md"
        _write_drug_markdown(md_path, query, profile, table)

        if not table.empty:
            top = table.iloc[0]
            summary_rows.append(
                {
                    "drug_query": query,
                    "prestwick_id": profile.get("drug", {}).get("prestwick_id") if isinstance(profile.get("drug"), dict) else None,
                    "top_disease": top.get("disease_name"),
                    "top_disease_score_raw_only": _safe_float(top.get("disease_score_raw_only")),
                    "top_disease_score_mechanism": _safe_float(top.get("disease_score_mechanism")),
                    "top_disease_mechanism_delta": _safe_float(top.get("mechanism_delta")),
                    "top_mechanism_balance_score": _safe_float(top.get("mechanism_balance_score")),
                }
            )

    summary = pd.DataFrame(summary_rows)
    summary_path = output_dir / "case_study_summary.csv"
    summary.to_csv(summary_path, index=False)
    payload = {
        "output_dir": str(output_dir),
        "n_case_studies": int(len(summary)),
        "summary_path": str(summary_path),
    }
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
