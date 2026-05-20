from __future__ import annotations

import pandas as pd

from gut_drug_microbiome.mechanism_layer import attach_action_signals
from gut_drug_microbiome.mechanism_layer import combine_disease_scores
from gut_drug_microbiome.mechanism_layer import compute_mechanism_layer
from gut_drug_microbiome.mechanism_layer import fuse_disease_scores
from gut_drug_microbiome.mechanism_layer import fusion_weights
from gut_drug_microbiome.mechanism_layer import infer_microbe_trait_priors


def _build_base_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "nt_code": "NT_BEN",
                "microbe_label": "Faecalibacterium prausnitzii",
                "species_label": "Faecalibacterium prausnitzii",
                "genus": "Faecalibacterium",
                "display_step1_predicted_effect_score": 0.8,
                "display_step1_predicted_inhibit_probability": 0.1,
                "predicted_promote_probability_refined": 0.9,
                "predicted_cross_feeding_reference_flag": True,
                "predicted_enzyme_step1_promote_support_score": 0.7,
                "predicted_enzyme_step1_inhibit_risk_score": 0.0,
                "predicted_mechanism_support_score": 0.6,
                "predicted_parent_depletion_fraction": -0.2,
                "applicability_flag": True,
                "relation_weight": 1.0,
            },
            {
                "nt_code": "NT_RISK",
                "microbe_label": "Escherichia coli",
                "species_label": "Escherichia coli",
                "genus": "Escherichia",
                "display_step1_predicted_effect_score": -0.7,
                "display_step1_predicted_inhibit_probability": 0.9,
                "predicted_promote_probability_refined": 0.1,
                "predicted_cross_feeding_reference_flag": False,
                "predicted_enzyme_step1_promote_support_score": 0.0,
                "predicted_enzyme_step1_inhibit_risk_score": 0.8,
                "predicted_mechanism_support_score": 0.5,
                "predicted_parent_depletion_fraction": -0.5,
                "applicability_flag": True,
                "relation_weight": 1.0,
            },
        ]
    )


def test_infer_microbe_trait_priors_uses_species_override() -> None:
    frame = _build_base_frame()
    enriched = infer_microbe_trait_priors(frame)

    ben = enriched[enriched["nt_code"] == "NT_BEN"].iloc[0]
    risk = enriched[enriched["nt_code"] == "NT_RISK"].iloc[0]

    assert float(ben["anti_inflammatory_weight"]) >= 0.95
    assert float(ben["butyrate_weight"]) >= 0.95
    assert float(risk["pro_inflammatory_weight"]) >= 0.95
    assert float(risk["pathobiont_weight"]) >= 0.95


def test_compute_mechanism_layer_rewards_beneficial_profile() -> None:
    frame = _build_base_frame()
    frame = infer_microbe_trait_priors(frame)
    frame = attach_action_signals(
        frame,
        score_column="display_step1_predicted_effect_score",
        inhibit_probability_column="display_step1_predicted_inhibit_probability",
        promote_probability_column="predicted_promote_probability_refined",
    )

    result = compute_mechanism_layer(frame, relation_weight_column="relation_weight")
    scores = result["scores"]

    assert scores["anti_inflammatory_score"] > scores["pro_inflammatory_score"]
    assert scores["butyrate_support_score"] > 0.2
    assert scores["mechanism_benefit_score"] > scores["mechanism_risk_score"]
    assert 0.0 <= scores["competition_vs_crossfeeding_proxy"] <= 1.0


def test_compute_mechanism_layer_captures_risk_shift() -> None:
    frame = _build_base_frame().copy()
    frame.loc[frame["nt_code"] == "NT_BEN", "display_step1_predicted_effect_score"] = -0.8
    frame.loc[frame["nt_code"] == "NT_BEN", "display_step1_predicted_inhibit_probability"] = 0.9
    frame.loc[frame["nt_code"] == "NT_BEN", "predicted_promote_probability_refined"] = 0.1
    frame.loc[frame["nt_code"] == "NT_RISK", "display_step1_predicted_effect_score"] = 0.9
    frame.loc[frame["nt_code"] == "NT_RISK", "display_step1_predicted_inhibit_probability"] = 0.1
    frame.loc[frame["nt_code"] == "NT_RISK", "predicted_promote_probability_refined"] = 0.95

    frame = infer_microbe_trait_priors(frame)
    frame = attach_action_signals(
        frame,
        score_column="display_step1_predicted_effect_score",
        inhibit_probability_column="display_step1_predicted_inhibit_probability",
        promote_probability_column="predicted_promote_probability_refined",
    )
    result = compute_mechanism_layer(frame, relation_weight_column="relation_weight")
    scores = result["scores"]

    assert scores["pro_inflammatory_score"] > scores["anti_inflammatory_score"]
    assert scores["pathobiont_load"] > 0.2
    assert scores["mechanism_risk_score"] > scores["mechanism_benefit_score"]


def test_combine_disease_scores_keeps_ablation_interface_stable() -> None:
    combined = combine_disease_scores(raw_microbe_score=0.42, mechanism_balance_score=-0.10, raw_weight=0.65)
    assert round(combined, 4) == 0.238


def test_fusion_modes_match_expected_weights_and_score() -> None:
    assert fusion_weights("raw_only") == (1.0, 0.0)
    assert fusion_weights("mechanism_only") == (0.0, 1.0)
    assert fusion_weights("weighted_0.3_0.7") == (0.3, 0.7)

    raw = 0.8
    mechanism = 0.2
    assert round(fuse_disease_scores(raw, mechanism, "raw_only"), 4) == 0.8
    assert round(fuse_disease_scores(raw, mechanism, "mechanism_only"), 4) == 0.2
    assert round(fuse_disease_scores(raw, mechanism, "weighted_0.3_0.7"), 4) == 0.38
