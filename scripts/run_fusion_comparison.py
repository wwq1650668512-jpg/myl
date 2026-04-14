from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.mechanism_layer import fuse_disease_scores
from gut_drug_microbiome.mechanism_layer import fusion_weights
from gut_drug_microbiome.web.service import GutPredictionService


FUSION_MODES = ["raw_only", "mechanism_only", "weighted_0.3_0.7"]
DRUG_PROFILE_PRIORS = {
    "rifaximin": "eubiotic_modulator",
    "vancomycin": "disruptive_antibiotic",
    "lubiprostone": "host_secretagogue",
    "metronidazole": "contextual_antimicrobial",
}

# Step 0 fixed inputs (must not be modified).
FIXED_CASES = [
    {
        "drug_name": "Rifaximin",
        "file_stem": "rifaximin",
        "smiles": "CC1=C(C(=O)NC(=C1C)C2=CC3=C(C=C2O)C(=O)C4=C(C3=O)C(=C(C=C4OC)O)OC)C",
    },
    {
        "drug_name": "Vancomycin",
        "file_stem": "vancomycin",
        "smiles": "C[C@H]1[C@H]([C@@](C[C@@H](O1)O[C@@H]2[C@H]([C@@H]([C@H](O[C@H]2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)[C@H]([C@H](C(=O)N[C@H](C(=O)N[C@H]5C(=O)N[C@@H]7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)[C@H](NC(=O)[C@H]([C@@H](C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)[C@@H](CC(C)C)NC)O)Cl)CO)O)O)(C)N)O",
    },
    {
        "drug_name": "Lubiprostone",
        "file_stem": "lubiprostone",
        "smiles": "CC(C)CCCC(C)C1CCC2C1(CCC3C2CCC4=CC(=O)CCC34C)C",
    },
    {
        "drug_name": "Metronidazole",
        "file_stem": "metronidazole",
        "smiles": "CN1C=NC(=N1)COC2=NC(=O)N(C2=O)C",
    },
]


def _coerce_float(value: object, default: float = 0.0) -> float:
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(numeric):
        return float(default)
    return float(numeric)


def _candidate_df(profile: dict[str, object]) -> pd.DataFrame:
    rows = profile.get("candidate_diseases", [])
    if not rows:
        raise RuntimeError("candidate_diseases is empty; cannot run fusion comparison.")

    frame = pd.DataFrame(rows)
    frame["raw_score"] = pd.to_numeric(frame.get("disease_score_raw_only"), errors="coerce").fillna(0.0)
    frame["mechanism_score"] = frame.get("mechanism_scores", pd.Series([{}] * len(frame))).map(
        lambda payload: _coerce_float(payload.get("mechanism_balance_score"), default=0.0)
        if isinstance(payload, dict)
        else 0.0
    )
    for key in [
        "anti_inflammatory_score",
        "pro_inflammatory_score",
        "butyrate_support_score",
        "barrier_protection_score",
        "toxin_risk_score",
        "mucus_degradation_score",
        "pathobiont_load",
        "competition_vs_crossfeeding_proxy",
    ]:
        frame[key] = frame.get("mechanism_scores", pd.Series([{}] * len(frame))).map(
            lambda payload, k=key: _coerce_float(payload.get(k), default=0.0) if isinstance(payload, dict) else 0.0
        )
    return frame


def _mode_scored_table(candidate_df: pd.DataFrame, fusion_mode: str) -> pd.DataFrame:
    return _mode_scored_table_with_gating(
        candidate_df=candidate_df,
        fusion_mode=fusion_mode,
        drug_name="",
        drug_profile="unknown",
        microbe_effect_strength=0.0,
        drug_butyrate_support_score=0.0,
    )


def _drug_profile(drug_name: str) -> str:
    return DRUG_PROFILE_PRIORS.get(str(drug_name).strip().lower(), "unknown")


def _is_infection_context(disease_name: str) -> bool:
    return _contains_any(
        disease_name,
        [r"abscess", r"脓肿", r"感染", r"fistula", r"瘘", r"\bibd\b", r"炎症性肠病", r"\bcd\b", r"克罗恩", r"colitis"],
    )


def _is_inflammatory_context(disease_name: str) -> bool:
    return _contains_any(
        disease_name,
        [
            r"\bcd\b",
            r"克罗恩",
            r"\bibd\b",
            r"炎症性肠病",
            r"\bcrc\b",
            r"结直肠癌",
            r"\bra\b",
            r"类风湿",
            r"\bsle\b",
            r"红斑狼疮",
            r"\bms\b",
            r"自身免疫",
            r"肝炎",
        ],
    )


def _is_host_symptom_context(disease_name: str) -> bool:
    return _contains_any(disease_name, [r"constipation", r"便秘", r"\bibs\b", r"肠易激", r"diarrhea", r"腹泻"])


def _compute_butyrate_support_from_profile(profile_payload: dict[str, object]) -> float:
    rows = profile_payload.get("butyrate_rows", [])
    if not isinstance(rows, list) or not rows:
        return 0.0
    supports: list[float] = []
    for row in rows:
        if not isinstance(row, dict):
            continue
        label = str(row.get("predicted_effect_label", "")).strip().lower()
        inhibit_probability = _coerce_float(row.get("predicted_inhibit_probability"), default=0.5)
        effect_score = _coerce_float(row.get("predicted_effect_score"), default=0.0)
        if label == "inhibit":
            support = max(0.0, 1.0 - inhibit_probability)
        elif label == "promote":
            support = min(1.0, 0.5 + max(effect_score, 0.0))
        else:
            support = max(0.0, 0.7 - max(0.0, -effect_score))
        supports.append(float(support))
    return float(sum(supports) / max(len(supports), 1))


def _compute_row_ecological_risk(frame: pd.DataFrame, microbe_effect_strength: float) -> pd.Series:
    pro = pd.to_numeric(frame.get("pro_inflammatory_score"), errors="coerce").fillna(0.0)
    toxin = pd.to_numeric(frame.get("toxin_risk_score"), errors="coerce").fillna(0.0)
    pathobiont = pd.to_numeric(frame.get("pathobiont_load"), errors="coerce").fillna(0.0)
    mucus = pd.to_numeric(frame.get("mucus_degradation_score"), errors="coerce").fillna(0.0)
    anti = pd.to_numeric(frame.get("anti_inflammatory_score"), errors="coerce").fillna(0.0)
    butyrate = pd.to_numeric(frame.get("butyrate_support_score"), errors="coerce").fillna(0.0)
    barrier = pd.to_numeric(frame.get("barrier_protection_score"), errors="coerce").fillna(0.0)

    risk_component = 0.40 * pro + 0.25 * toxin + 0.20 * pathobiont + 0.15 * mucus
    protection_component = 0.20 * anti + 0.20 * butyrate + 0.15 * barrier
    ecology = (risk_component - protection_component).clip(lower=0.0)
    ecology = (ecology + 0.35 * float(microbe_effect_strength)).clip(lower=0.0, upper=1.0)
    return ecology


def _mode_scored_table_with_gating(
    candidate_df: pd.DataFrame,
    fusion_mode: str,
    drug_name: str,
    drug_profile: str,
    microbe_effect_strength: float,
    drug_butyrate_support_score: float,
) -> pd.DataFrame:
    raw_weight, mechanism_weight = fusion_weights(fusion_mode)
    scored = candidate_df.copy()
    scored["final_score_base"] = [
        fuse_disease_scores(raw_score=raw, mechanism_score=mechanism, fusion_mode=fusion_mode)
        for raw, mechanism in zip(scored["raw_score"], scored["mechanism_score"])
    ]
    scored["drug_profile"] = str(drug_profile)
    scored["microbe_effect_strength"] = float(microbe_effect_strength)
    scored["drug_butyrate_support_score"] = float(drug_butyrate_support_score)
    scored["ecological_risk_score"] = _compute_row_ecological_risk(scored, microbe_effect_strength=microbe_effect_strength)

    gating_multipliers: list[float] = []
    gating_reasons: list[str] = []
    for _, row in scored.iterrows():
        disease_name = str(row.get("disease_name", ""))
        row_butyrate = _coerce_float(row.get("butyrate_support_score"), default=0.0)
        row_ecology_risk = _coerce_float(row.get("ecological_risk_score"), default=0.0)
        multiplier = 1.0
        reasons: list[str] = []

        if drug_profile == "eubiotic_modulator":
            butyrate_drop_threshold = max(0.06, 0.7 * float(drug_butyrate_support_score))
            if row_butyrate < butyrate_drop_threshold and _is_inflammatory_context(disease_name):
                multiplier *= 0.5
                reasons.append(f"rifaximin_butyrate_drop<{butyrate_drop_threshold:.3f}")

        elif drug_profile == "disruptive_antibiotic":
            dynamic_threshold = max(0.12, 0.85 * _coerce_float(scored["ecological_risk_score"].median(), default=0.12))
            if row_ecology_risk >= dynamic_threshold:
                multiplier *= 0.6
                reasons.append(f"vancomycin_high_ecological_risk>={dynamic_threshold:.3f}")

        elif drug_profile == "host_secretagogue":
            if float(microbe_effect_strength) >= 0.45 and not _is_host_symptom_context(disease_name):
                multiplier *= 0.3
                reasons.append("lubiprostone_high_microbe_effect_strength")

        elif drug_profile == "contextual_antimicrobial":
            if _is_infection_context(disease_name):
                multiplier *= 1.15
                reasons.append("metronidazole_infection_context_bonus")
            multiplier *= max(0.35, 1.0 - 0.8 * row_ecology_risk)
            reasons.append("metronidazole_ecological_risk_penalty")

        gating_multipliers.append(float(multiplier))
        gating_reasons.append("|".join(reasons) if reasons else "none")

    scored["gating_multiplier"] = gating_multipliers
    scored["gating_reason"] = gating_reasons
    scored["final_score"] = scored["final_score_base"] * scored["gating_multiplier"]

    scored = scored.sort_values(["final_score", "raw_score", "mechanism_score"], ascending=[False, False, False]).reset_index(drop=True)
    scored["rank"] = scored.index + 1
    scored["fusion_mode"] = fusion_mode
    scored["raw_contribution_share"] = float(raw_weight)
    scored["mechanism_contribution_share"] = float(mechanism_weight)
    scored["raw_weighted_component"] = float(raw_weight) * scored["raw_score"]
    scored["mechanism_weighted_component"] = float(mechanism_weight) * scored["mechanism_score"]
    return scored


def _contains_any(text: str, patterns: list[str]) -> bool:
    value = str(text or "").lower()
    return any(re.search(pattern, value, flags=re.IGNORECASE) for pattern in patterns)


def _best_rank(df: pd.DataFrame, patterns: list[str]) -> int | None:
    mask = df["disease_name"].map(lambda value: _contains_any(str(value), patterns))
    matched = df.loc[mask]
    if matched.empty:
        return None
    return int(matched["rank"].min())


def _top3_mean(df: pd.DataFrame, column: str) -> float:
    return float(pd.to_numeric(df.head(3)[column], errors="coerce").fillna(0.0).mean())


def _build_rule_row(
    fusion_mode: str,
    drug_name: str,
    rule_id: str,
    description: str,
    passed: bool,
    details: str,
) -> dict[str, object]:
    return {
        "record_type": "rule",
        "fusion_mode": fusion_mode,
        "drug_name": drug_name,
        "rule_id": rule_id,
        "rule_description": description,
        "status": "PASS" if passed else "FAIL",
        "passed": bool(passed),
        "details": details,
    }


def _sanity_benchmark(
    per_drug_modes: dict[str, dict[str, pd.DataFrame]],
    inhibit_counts: dict[str, int],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    antibiotics = ["Rifaximin", "Vancomycin", "Metronidazole"]
    antibiotic_mean_inhibit = float(sum(inhibit_counts.get(name, 0) for name in antibiotics) / max(len(antibiotics), 1))

    for mode in FUSION_MODES:
        # Rifaximin rules
        rifa_raw = per_drug_modes["Rifaximin"]["raw_only"]
        rifa = per_drug_modes["Rifaximin"][mode]
        target_patterns = [r"\bibs\b", r"肠易激", r"便秘", r"constipation"]
        rank_mode = _best_rank(rifa, target_patterns)
        rank_raw = _best_rank(rifa_raw, target_patterns)
        rifa_rule_1 = rank_mode is not None and (
            rank_mode <= 3 or (rank_raw is not None and rank_mode < rank_raw)
        )
        rows.append(
            _build_rule_row(
                mode,
                "Rifaximin",
                "RIFA_1",
                "IBS/便秘排名应上升或进入前3",
                rifa_rule_1,
                f"rank_mode={rank_mode}, rank_raw={rank_raw}",
            )
        )

        baseline_butyrate = _top3_mean(rifa_raw, "butyrate_support_score")
        mode_butyrate = _top3_mean(rifa, "butyrate_support_score")
        rifa_rule_2 = mode_butyrate >= 0.8 * baseline_butyrate
        rows.append(
            _build_rule_row(
                mode,
                "Rifaximin",
                "RIFA_2",
                "butyrate_support_score 不应显著下降",
                rifa_rule_2,
                f"top3_mode={mode_butyrate:.4f}, top3_raw={baseline_butyrate:.4f}",
            )
        )

        # Vancomycin rules
        vanco = per_drug_modes["Vancomycin"][mode]
        top3 = vanco.head(3)["final_score"].tolist()
        top1 = float(top3[0]) if top3 else 0.0
        high_count = int((vanco["final_score"] >= 0.9 * top1).sum()) if top1 > 0 else 0
        concentration = float(top1 / max(sum(top3), 1e-8)) if top3 else 0.0
        vanco_rule_1 = high_count <= 3 and concentration >= 0.34
        rows.append(
            _build_rule_row(
                mode,
                "Vancomycin",
                "VANCO_1",
                "不应多个疾病都高分（top3应相对集中）",
                vanco_rule_1,
                f"high_count={high_count}, concentration={concentration:.4f}",
            )
        )
        cd_rank = _best_rank(vanco, [r"克罗恩", r"\bcd\b", r"crohn"])
        vanco_rule_2 = cd_rank is None or cd_rank > 3
        rows.append(
            _build_rule_row(
                mode,
                "Vancomycin",
                "VANCO_2",
                "CD 排名不应进入前3",
                vanco_rule_2,
                f"cd_rank={cd_rank}",
            )
        )

        # Lubiprostone rules
        lubi = per_drug_modes["Lubiprostone"][mode]
        lubi_inhibit = int(inhibit_counts.get("Lubiprostone", 0))
        lubi_rule_1 = lubi_inhibit <= 0.7 * antibiotic_mean_inhibit
        rows.append(
            _build_rule_row(
                mode,
                "Lubiprostone",
                "LUBI_1",
                "抑制数量应明显低于抗生素",
                lubi_rule_1,
                f"lubi_inhibit={lubi_inhibit}, antibiotic_mean_inhibit={antibiotic_mean_inhibit:.2f}",
            )
        )
        barrier_mean = _top3_mean(lubi, "barrier_protection_score")
        antimicrobial_mean = 0.5 * (_top3_mean(lubi, "pro_inflammatory_score") + _top3_mean(lubi, "toxin_risk_score"))
        lubi_rule_2 = barrier_mean >= antimicrobial_mean
        rows.append(
            _build_rule_row(
                mode,
                "Lubiprostone",
                "LUBI_2",
                "机制层应偏 barrier/secretion 而非 antimicrobial",
                lubi_rule_2,
                f"barrier_top3={barrier_mean:.4f}, antimicrobial_top3={antimicrobial_mean:.4f}",
            )
        )

        # Metronidazole rules
        metro = per_drug_modes["Metronidazole"][mode]
        infect_rank = _best_rank(metro, [r"abscess", r"脓肿", r"\bibd\b", r"炎症性肠病", r"colitis"])
        metro_rule_1 = infect_rank is not None and infect_rank <= 4
        rows.append(
            _build_rule_row(
                mode,
                "Metronidazole",
                "METRO_1",
                "感染/炎症相关疾病应有一定信号",
                metro_rule_1,
                f"infection_or_ibd_rank={infect_rank}",
            )
        )
        pro_mean = _top3_mean(metro, "pro_inflammatory_score")
        toxin_mean = _top3_mean(metro, "toxin_risk_score")
        metro_rule_2 = (pro_mean >= 0.08) or (toxin_mean >= 0.05)
        rows.append(
            _build_rule_row(
                mode,
                "Metronidazole",
                "METRO_2",
                "应存在生态风险信号（pro_inflammatory 或 toxin）",
                metro_rule_2,
                f"pro_top3={pro_mean:.4f}, toxin_top3={toxin_mean:.4f}",
            )
        )

    summary = pd.DataFrame(rows)
    mode_summary_rows = []
    for mode, group in summary.groupby("fusion_mode", dropna=False):
        pass_count = int(group["passed"].sum())
        total = int(len(group))
        mode_summary_rows.append(
            {
                "record_type": "mode_summary",
                "fusion_mode": mode,
                "drug_name": "__ALL__",
                "rule_id": "PASS_COUNT",
                "rule_description": "mode-level pass count summary",
                "status": f"{pass_count}/{total}",
                "passed": pass_count == total,
                "details": f"pass_count={pass_count}, total={total}, pass_rate={pass_count/max(total,1):.3f}",
            }
        )
    return pd.concat([summary, pd.DataFrame(mode_summary_rows)], ignore_index=True)


def _build_summary_markdown(
    output_path: Path,
    combined_table: pd.DataFrame,
    sanity: pd.DataFrame,
    per_drug_modes: dict[str, dict[str, pd.DataFrame]],
) -> None:
    mode_rule = sanity[sanity["record_type"] == "rule"].copy()
    mode_agg = (
        mode_rule.groupby("fusion_mode", dropna=False)["passed"]
        .agg(pass_count="sum", total="count")
        .reset_index()
    )
    mode_agg["pass_rate"] = mode_agg["pass_count"] / mode_agg["total"]
    mode_agg = mode_agg.sort_values(["pass_count", "pass_rate"], ascending=[False, False]).reset_index(drop=True)
    best_mode = str(mode_agg.iloc[0]["fusion_mode"]) if not mode_agg.empty else "N/A"

    lines: list[str] = []
    lines.append("# Fusion Comparison Summary")
    lines.append("")
    lines.append("## 1. 三种融合方式的整体差异")
    lines.append("")
    diff = (
        combined_table.groupby("fusion_mode", dropna=False)["final_score"]
        .agg(mean_final_score="mean", std_final_score="std", max_final_score="max")
        .reset_index()
        .sort_values("mean_final_score", ascending=False)
    )
    lines.append(diff.to_markdown(index=False))
    lines.append("")

    lines.append("## 2. Sanity benchmark 最优模式")
    lines.append("")
    lines.append(mode_agg.to_markdown(index=False))
    lines.append("")
    lines.append(f"- 最佳模式（PASS 数最多）：`{best_mode}`")
    lines.append("")

    lines.append("## 3. 关键案例分析")
    lines.append("")
    for drug_name in ["Rifaximin", "Vancomycin"]:
        rows = []
        for mode in FUSION_MODES:
            table = per_drug_modes[drug_name][mode]
            top = table.iloc[0]
            rows.append(
                {
                    "fusion_mode": mode,
                    "top_disease": top["disease_name"],
                    "top_final_score": float(top["final_score"]),
                    "top_raw_score": float(top["raw_score"]),
                    "top_mechanism_score": float(top["mechanism_score"]),
                }
            )
        lines.append(f"### {drug_name}")
        lines.append("")
        lines.append(pd.DataFrame(rows).to_markdown(index=False))
        lines.append("")

    lines.append("## 4. 结论")
    lines.append("")
    lines.append(f"- 当前最优融合方式：`{best_mode}`（按 sanity PASS 数量）。")
    lines.append("- 建议下一步尝试：按药物类型做 gating（antibiotic vs modulator）或使用非线性融合（例如 piecewise / sigmoid gating）。")
    lines.append("")

    output_path.write_text("\n".join(lines), encoding="utf-8")


def run() -> dict[str, object]:
    output_dir = ROOT / "predictions" / "evaluation" / "fusion_comparison"
    output_dir.mkdir(parents=True, exist_ok=True)

    service = GutPredictionService()
    per_drug_modes: dict[str, dict[str, pd.DataFrame]] = {}
    inhibit_counts: dict[str, int] = {}
    profile_evidence: dict[str, dict[str, object]] = {}
    combined_rows: list[pd.DataFrame] = []

    for case in FIXED_CASES:
        drug_name = case["drug_name"]
        smiles = case["smiles"]
        file_stem = case["file_stem"]

        prediction = service.predict_custom_drug(drug_name=drug_name, smiles=smiles)
        profile = prediction.get("profile", {})
        if not isinstance(profile, dict):
            raise RuntimeError(f"{drug_name}: invalid profile payload.")

        aggregated = profile.get("aggregated", {})
        step1_counts = aggregated.get("step1_counts", {}) if isinstance(aggregated, dict) else {}
        inhibit_counts[drug_name] = int(step1_counts.get("inhibit", 0))

        panel_effect = pd.DataFrame(profile.get("panel_effect_microbes", []))
        if panel_effect.empty:
            butyrate_rows: list[dict[str, object]] = []
            inhibit_prob_ge_07 = 0
            inhibit_prob_ge_05 = 0
            effect_score_le_02 = 0
            effect_score_ge_02 = 0
        else:
            panel_effect["_inh"] = pd.to_numeric(panel_effect.get("predicted_inhibit_probability"), errors="coerce").fillna(0.0)
            panel_effect["_eff"] = pd.to_numeric(panel_effect.get("predicted_effect_score"), errors="coerce").fillna(0.0)
            butyrate_mask = panel_effect.get("microbe_label", pd.Series("", index=panel_effect.index)).astype(str).str.contains(
                r"Faecalibacterium prausnitzii|Roseburia|Eubacterium rectale",
                case=False,
                regex=True,
                na=False,
            )
            butyrate_rows = (
                panel_effect.loc[
                    butyrate_mask,
                    ["microbe_label", "predicted_effect_label", "predicted_inhibit_probability", "predicted_effect_score"],
                ]
                .to_dict(orient="records")
            )
            inhibit_prob_ge_07 = int((panel_effect["_inh"] >= 0.7).sum())
            inhibit_prob_ge_05 = int((panel_effect["_inh"] >= 0.5).sum())
            effect_score_le_02 = int((panel_effect["_eff"] <= -0.2).sum())
            effect_score_ge_02 = int((panel_effect["_eff"] >= 0.2).sum())

        profile_evidence[drug_name] = {
            "session_id": str(prediction.get("session_id", "")),
            "step1_counts": step1_counts,
            "mean_predicted_inhibit_probability": aggregated.get("mean_predicted_inhibit_probability")
            if isinstance(aggregated, dict)
            else None,
            "mean_predicted_effect_score": aggregated.get("mean_predicted_effect_score") if isinstance(aggregated, dict) else None,
            "butyrate_rows": butyrate_rows,
            "panel_size": int(len(panel_effect)),
            "inhibit_prob_ge_0_7": inhibit_prob_ge_07,
            "inhibit_prob_ge_0_5": inhibit_prob_ge_05,
            "effect_score_le_-0_2": effect_score_le_02,
            "effect_score_ge_0_2": effect_score_ge_02,
        }

        panel_size = int(len(panel_effect))
        inhibit_count = int(step1_counts.get("inhibit", 0))
        microbe_effect_strength = float(inhibit_count / max(panel_size, 1))
        drug_profile = _drug_profile(drug_name)
        drug_butyrate_support_score = _compute_butyrate_support_from_profile(profile_evidence[drug_name])

        base = _candidate_df(profile)
        mode_tables: dict[str, pd.DataFrame] = {}
        for mode in FUSION_MODES:
            mode_tables[mode] = _mode_scored_table_with_gating(
                candidate_df=base,
                fusion_mode=mode,
                drug_name=drug_name,
                drug_profile=drug_profile,
                microbe_effect_strength=microbe_effect_strength,
                drug_butyrate_support_score=drug_butyrate_support_score,
            )

        drug_ecological_risk_score = float(mode_tables["weighted_0.3_0.7"]["ecological_risk_score"].mean())
        profile_evidence[drug_name]["drug_profile"] = drug_profile
        profile_evidence[drug_name]["microbe_effect_strength"] = microbe_effect_strength
        profile_evidence[drug_name]["drug_butyrate_support_score"] = drug_butyrate_support_score
        profile_evidence[drug_name]["ecological_risk_score"] = drug_ecological_risk_score

        rank_compare = pd.DataFrame({"disease_name": base["disease_name"]})
        for mode in FUSION_MODES:
            rank_compare[f"rank_{mode}"] = mode_tables[mode].set_index("disease_name")["rank"].reindex(rank_compare["disease_name"]).values

        long_rows = []
        for mode in FUSION_MODES:
            table = mode_tables[mode].copy()
            merged = table.merge(rank_compare, on="disease_name", how="left")
            long_rows.append(merged)
            combined_rows.append(merged.assign(drug_name=drug_name))
        drug_output = pd.concat(long_rows, ignore_index=True)
        drug_output = drug_output.sort_values(["fusion_mode", "rank", "disease_name"]).reset_index(drop=True)

        keep_columns = [
            "fusion_mode",
            "disease_name",
            "raw_score",
            "mechanism_score",
            "final_score_base",
            "final_score",
            "rank",
            "rank_raw_only",
            "rank_mechanism_only",
            "rank_weighted_0.3_0.7",
            "raw_contribution_share",
            "mechanism_contribution_share",
            "raw_weighted_component",
            "mechanism_weighted_component",
            "anti_inflammatory_score",
            "pro_inflammatory_score",
            "butyrate_support_score",
            "barrier_protection_score",
            "toxin_risk_score",
            "mucus_degradation_score",
            "pathobiont_load",
            "competition_vs_crossfeeding_proxy",
            "ecological_risk_score",
            "drug_profile",
            "microbe_effect_strength",
            "drug_butyrate_support_score",
            "gating_multiplier",
            "gating_reason",
        ]
        existing_columns = [column for column in keep_columns if column in drug_output.columns]
        drug_path = output_dir / f"{file_stem}.csv"
        drug_output.loc[:, existing_columns].to_csv(drug_path, index=False)

        per_drug_modes[drug_name] = mode_tables

    combined = pd.concat(combined_rows, ignore_index=True)
    sanity = _sanity_benchmark(per_drug_modes=per_drug_modes, inhibit_counts=inhibit_counts)
    sanity_path = output_dir / "sanity_summary.csv"
    sanity.to_csv(sanity_path, index=False)

    summary_path = output_dir / "fusion_comparison_summary.md"
    _build_summary_markdown(
        output_path=summary_path,
        combined_table=combined,
        sanity=sanity,
        per_drug_modes=per_drug_modes,
    )

    profile_evidence_path = output_dir / "profile_evidence.json"
    profile_evidence_path.write_text(json.dumps(profile_evidence, ensure_ascii=False, indent=2), encoding="utf-8")

    return {
        "output_dir": str(output_dir),
        "drug_files": [f"{item['file_stem']}.csv" for item in FIXED_CASES],
        "sanity_summary": str(sanity_path),
        "summary_markdown": str(summary_path),
        "profile_evidence": str(profile_evidence_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run fixed-case fusion mode comparison and sanity benchmark."
    )
    parser.parse_args()
    try:
        result = run()
    except Exception as exc:  # noqa: BLE001
        print(
            json.dumps(
                {
                    "status": "failed",
                    "error": str(exc),
                },
                ensure_ascii=False,
                indent=2,
            )
        )
        raise SystemExit(1) from exc
    print(json.dumps({"status": "ok", **result}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
