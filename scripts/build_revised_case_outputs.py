from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
FUSION_DIR = ROOT / "predictions" / "evaluation" / "fusion_comparison"
BEFORE_RESULTS_PATH = FUSION_DIR / "revised_case_based_results_before_realism_fix.csv"
AFTER_RESULTS_PATH = FUSION_DIR / "case_based_results.csv"
PROFILE_EVIDENCE_PATH = FUSION_DIR / "profile_evidence.json"

RANKING_SUMMARY_PATH = FUSION_DIR / "ranking_benchmark_summary.csv"
ECOLOGY_SUMMARY_PATH = FUSION_DIR / "ecology_benchmark_summary.csv"
REVISED_RESULTS_PATH = FUSION_DIR / "revised_case_based_results.csv"
REVISED_SUMMARY_MD_PATH = FUSION_DIR / "revised_case_based_summary.md"

RANKING_ASSERTIONS = {"RFX_01", "VAN_02", "MTZ_01"}


def _status_rank(status: str) -> int:
    mapping = {"FAIL": 0, "PARTIAL": 1, "PASS": 2}
    return mapping.get(str(status).strip().upper(), -1)


def _change_label(before: str, after: str) -> str:
    b = _status_rank(before)
    a = _status_rank(after)
    if a > b:
        return "improved"
    if a < b:
        return "regressed"
    return "unchanged"


def _load_weighted_table(drug_file: str) -> pd.DataFrame:
    table = pd.read_csv(FUSION_DIR / drug_file)
    sub = table[table["fusion_mode"] == "weighted_0.3_0.7"].copy()
    return sub.sort_values("rank").reset_index(drop=True)


def _rank_of(table: pd.DataFrame, patterns: list[str]) -> int | None:
    mask = table["disease_name"].astype(str).str.contains("|".join(patterns), case=False, regex=True, na=False)
    matched = table.loc[mask]
    if matched.empty:
        return None
    return int(matched["rank"].min())


def _ecology_check_rows(before: pd.DataFrame, profile: dict[str, object]) -> pd.DataFrame:
    before_status = {row["assertion_id"]: row["Status"] for _, row in before.iterrows()}
    before_notes = {row["assertion_id"]: row.get("Key conflicting output", "") for _, row in before.iterrows()}

    rows: list[dict[str, object]] = []

    # Rifaximin: core butyrate inhibition.
    rifa = dict(profile.get("Rifaximin", {}))
    rifa_rows = rifa.get("butyrate_rows", []) or []
    strong = 0
    for row in rifa_rows:
        label = str(row.get("predicted_effect_label", "")).strip().lower()
        prob = pd.to_numeric(pd.Series([row.get("predicted_inhibit_probability")]), errors="coerce").fillna(0.0).iloc[0]
        score = pd.to_numeric(pd.Series([row.get("predicted_effect_score")]), errors="coerce").fillna(0.0).iloc[0]
        if label == "inhibit" and (prob >= 0.5 or score < 0.0):
            strong += 1
    if strong <= 2:
        rifa_status = "PASS"
    elif strong <= 3:
        rifa_status = "PARTIAL"
    else:
        rifa_status = "FAIL"
    rows.append(
        {
            "check_id": "ECO_RFX_BUTYRATE",
            "drug_name": "Rifaximin",
            "check_name": "core_butyrate_not_strongly_inhibited",
            "status_before": before_status.get("RFX_02", "NA"),
            "status_after": rifa_status,
            "change": _change_label(before_status.get("RFX_02", "FAIL"), rifa_status),
            "details_before": before_notes.get("RFX_02", ""),
            "details_after": f"strong_core_butyrate_inhibit={strong}/{max(len(rifa_rows),1)}",
        }
    )

    # Lubiprostone: broad antimicrobial pattern.
    lubi = dict(profile.get("Lubiprostone", {}))
    panel_size = int(lubi.get("panel_size", 0) or 0)
    inhibit_count = int(dict(lubi.get("step1_counts", {})).get("inhibit", 0))
    inhibit_fraction = float(inhibit_count / max(panel_size, 1))
    if inhibit_fraction <= 0.45:
        lubi_status = "PASS"
    elif inhibit_fraction <= 0.60:
        lubi_status = "PARTIAL"
    else:
        lubi_status = "FAIL"
    rows.append(
        {
            "check_id": "ECO_LUB_BROAD_ANTIMICROBIAL",
            "drug_name": "Lubiprostone",
            "check_name": "not_broad_antimicrobial",
            "status_before": before_status.get("LUB_02", "NA"),
            "status_after": lubi_status,
            "change": _change_label(before_status.get("LUB_02", "FAIL"), lubi_status),
            "details_before": before_notes.get("LUB_02", ""),
            "details_after": f"inhibit_fraction={inhibit_fraction:.3f}, inhibit_count={inhibit_count}, panel_size={panel_size}",
        }
    )

    # Vancomycin: high ecological risk should be visible.
    vanco = dict(profile.get("Vancomycin", {}))
    vanco_risk = float(pd.to_numeric(pd.Series([vanco.get("ecological_risk_score")]), errors="coerce").fillna(0.0).iloc[0])
    if vanco_risk >= 0.15:
        vanco_status = "PASS"
    elif vanco_risk >= 0.12:
        vanco_status = "PARTIAL"
    else:
        vanco_status = "FAIL"
    rows.append(
        {
            "check_id": "ECO_VANCO_HIGH_RISK",
            "drug_name": "Vancomycin",
            "check_name": "high_ecological_risk_visible",
            "status_before": before_status.get("VAN_01", "NA"),
            "status_after": vanco_status,
            "change": _change_label(before_status.get("VAN_01", "FAIL"), vanco_status),
            "details_before": before_notes.get("VAN_01", ""),
            "details_after": f"ecological_risk_score={vanco_risk:.3f}",
        }
    )

    # Metronidazole: benefit + risk coexist.
    metro = dict(profile.get("Metronidazole", {}))
    metro_table = _load_weighted_table("metronidazole.csv")
    infect_rank = _rank_of(metro_table, [r"abscess", r"脓肿", r"\bibd\b", r"炎症性肠病", r"\bcd\b", r"克罗恩", r"fistula", r"肛瘘"])
    metro_panel = int(metro.get("panel_size", 0) or 0)
    metro_inhibit = int(dict(metro.get("step1_counts", {})).get("inhibit", 0))
    metro_inhibit_frac = float(metro_inhibit / max(metro_panel, 1))
    metro_risk = float(pd.to_numeric(pd.Series([metro.get("ecological_risk_score")]), errors="coerce").fillna(0.0).iloc[0])
    benefit_ok = infect_rank is not None and infect_rank <= 4
    risk_ok = metro_risk >= 0.15 or metro_inhibit_frac >= 0.30
    if benefit_ok and risk_ok:
        metro_status = "PASS"
    elif benefit_ok or risk_ok:
        metro_status = "PARTIAL"
    else:
        metro_status = "FAIL"
    # Before status synthesized from MTZ_01 + MTZ_02.
    b1 = str(before_status.get("MTZ_01", "FAIL"))
    b2 = str(before_status.get("MTZ_02", "FAIL"))
    if b1 == "PASS" and b2 == "PASS":
        metro_before = "PASS"
    elif b1 == "FAIL" and b2 == "FAIL":
        metro_before = "FAIL"
    else:
        metro_before = "PARTIAL"
    rows.append(
        {
            "check_id": "ECO_MTZ_BENEFIT_PLUS_RISK",
            "drug_name": "Metronidazole",
            "check_name": "benefit_and_risk_coexist",
            "status_before": metro_before,
            "status_after": metro_status,
            "change": _change_label(metro_before, metro_status),
            "details_before": f"MTZ_01={b1}, MTZ_02={b2}",
            "details_after": f"infection_rank={infect_rank}, ecological_risk={metro_risk:.3f}, inhibit_fraction={metro_inhibit_frac:.3f}",
        }
    )

    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build revised ranking/ecology benchmark summaries.")
    parser.add_argument("--before", type=Path, default=BEFORE_RESULTS_PATH)
    parser.add_argument("--after", type=Path, default=AFTER_RESULTS_PATH)
    parser.add_argument("--profile", type=Path, default=PROFILE_EVIDENCE_PATH)
    parser.add_argument("--ranking-out", type=Path, default=RANKING_SUMMARY_PATH)
    parser.add_argument("--ecology-out", type=Path, default=ECOLOGY_SUMMARY_PATH)
    parser.add_argument("--revised-results-out", type=Path, default=REVISED_RESULTS_PATH)
    parser.add_argument("--summary-md-out", type=Path, default=REVISED_SUMMARY_MD_PATH)
    args = parser.parse_args()

    before = pd.read_csv(args.before)
    after = pd.read_csv(args.after)
    profile_payload = json.loads(args.profile.read_text(encoding="utf-8"))

    merged = before.loc[:, ["assertion_id", "Status"]].merge(
        after,
        on="assertion_id",
        how="right",
        suffixes=("_before", ""),
    )
    merged = merged.rename(columns={"Status_before": "Status_before", "Status": "Status_after"})
    merged["Status_before"] = merged["Status_before"].fillna("NA")
    merged["Status_after"] = merged["Status_after"].fillna("NA")
    merged["Status_change"] = [
        _change_label(before_status, after_status)
        for before_status, after_status in zip(merged["Status_before"], merged["Status_after"])
    ]
    merged.to_csv(args.revised_results_out, index=False)

    ranking = merged[merged["assertion_id"].isin(RANKING_ASSERTIONS)].copy()
    ranking["layer"] = "ranking"
    ranking = ranking.loc[
        :,
        ["layer", "assertion_id", "Drug", "Status_before", "Status_after", "Status_change", "Key conflicting output"],
    ].rename(columns={"Key conflicting output": "details_after"})
    args.ranking_out.parent.mkdir(parents=True, exist_ok=True)
    ranking.to_csv(args.ranking_out, index=False)

    ecology = _ecology_check_rows(before=before, profile=profile_payload)
    ecology.to_csv(args.ecology_out, index=False)

    def _count(df: pd.DataFrame, col: str) -> str:
        vc = df[col].value_counts().to_dict()
        return f"PASS={vc.get('PASS',0)}, PARTIAL={vc.get('PARTIAL',0)}, FAIL={vc.get('FAIL',0)}"

    lines: list[str] = []
    lines.append("# Revised Case-Based Summary")
    lines.append("")
    lines.append("## Ranking Layer")
    lines.append(f"- Before: {_count(ranking, 'Status_before')}")
    lines.append(f"- After: {_count(ranking, 'Status_after')}")
    lines.append("")
    improved_ranking = ranking[ranking["Status_change"] == "improved"]
    if improved_ranking.empty:
        lines.append("- FAIL -> PASS: none")
    else:
        lines.append("- Improved assertions:")
        for _, row in improved_ranking.iterrows():
            lines.append(f"  - {row['assertion_id']} ({row['Drug']}): {row['Status_before']} -> {row['Status_after']}")
    lines.append("")

    lines.append("## Ecology Layer")
    lines.append(f"- Before: {_count(ecology, 'status_before')}")
    lines.append(f"- After: {_count(ecology, 'status_after')}")
    lines.append("")
    improved_ecology = ecology[ecology["change"] == "improved"]
    if improved_ecology.empty:
        lines.append("- FAIL -> PASS: none")
    else:
        lines.append("- Improved checks:")
        for _, row in improved_ecology.iterrows():
            lines.append(f"  - {row['check_id']} ({row['drug_name']}): {row['status_before']} -> {row['status_after']}")
    lines.append("")

    lines.append("## Remaining Failures")
    remaining_ranking = merged[merged["Status_after"] == "FAIL"]
    remaining_ecology = ecology[ecology["status_after"] == "FAIL"]
    if remaining_ranking.empty and remaining_ecology.empty:
        lines.append("- none")
    else:
        for _, row in remaining_ranking.iterrows():
            lines.append(f"- Ranking {row['assertion_id']} ({row['Drug']}): {row.get('Key conflicting output', '')}")
        for _, row in remaining_ecology.iterrows():
            lines.append(f"- Ecology {row['check_id']} ({row['drug_name']}): {row.get('details_after', '')}")
    lines.append("")

    args.summary_md_out.write_text("\n".join(lines), encoding="utf-8")

    print(
        {
            "status": "ok",
            "ranking_summary": str(args.ranking_out),
            "ecology_summary": str(args.ecology_out),
            "revised_results": str(args.revised_results_out),
            "revised_summary_md": str(args.summary_md_out),
        }
    )


if __name__ == "__main__":
    main()
