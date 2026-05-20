from __future__ import annotations

import pandas as pd

from gut_drug_microbiome.web.service import evaluate_prediction_confidence


def _panel(rows: list[dict[str, object]]) -> pd.DataFrame:
    return pd.DataFrame(rows)


def test_unreasonable_prediction_triggers_low_confidence() -> None:
    frame = _panel(
        [
            {
                "microbe_label": "Faecalibacterium prausnitzii",
                "species_label": "Faecalibacterium prausnitzii",
                "predicted_effect_label": "inhibit",
                "predicted_inhibit_probability": 0.95,
                "predicted_effect_score": -0.72,
            },
            {
                "microbe_label": "Roseburia intestinalis",
                "species_label": "Roseburia intestinalis",
                "predicted_effect_label": "inhibit",
                "predicted_inhibit_probability": 0.88,
                "predicted_effect_score": -0.58,
            },
            {
                "microbe_label": "Eubacterium rectale",
                "species_label": "Eubacterium rectale",
                "predicted_effect_label": "inhibit",
                "predicted_inhibit_probability": 0.84,
                "predicted_effect_score": -0.52,
            },
            {
                "microbe_label": "Bifidobacterium longum",
                "species_label": "Bifidobacterium longum",
                "predicted_effect_label": "inhibit",
                "predicted_inhibit_probability": 0.74,
                "predicted_effect_score": -0.31,
            },
        ]
    )
    payload = evaluate_prediction_confidence(
        effect_frame=frame,
        step1_label_column="predicted_effect_label",
        step1_probability_column="predicted_inhibit_probability",
        step1_score_column="predicted_effect_score",
        drug_profile="eubiotic_modulator",
        molecular_weight=420.0,
        xlogp=3.2,
        mw_bounds=(100.0, 650.0),
        xlogp_bounds=(-2.0, 5.0),
    )
    assert payload["confidence_tier"] == "low"
    assert float(payload["confidence_score"]) < 0.45
    assert "主要风险来自" in str(payload["confidence_explanation"])


def test_ood_molecule_triggers_ood_warning() -> None:
    frame = _panel(
        [
            {
                "microbe_label": "Bifidobacterium longum",
                "species_label": "Bifidobacterium longum",
                "predicted_effect_label": "no_effect",
                "predicted_inhibit_probability": 0.12,
                "predicted_effect_score": 0.03,
            },
            {
                "microbe_label": "Blautia obeum",
                "species_label": "Blautia obeum",
                "predicted_effect_label": "no_effect",
                "predicted_inhibit_probability": 0.11,
                "predicted_effect_score": 0.02,
            },
        ]
    )
    payload = evaluate_prediction_confidence(
        effect_frame=frame,
        step1_label_column="predicted_effect_label",
        step1_probability_column="predicted_inhibit_probability",
        step1_score_column="predicted_effect_score",
        drug_profile="unknown",
        molecular_weight=40.0,
        xlogp=0.1,
        mw_bounds=(100.0, 650.0),
        xlogp_bounds=(-2.0, 5.0),
    )
    assert "OOD-molecule" in set(payload["warning_flags"])


def test_strong_core_butyrate_suppression_triggers_ecology_risk() -> None:
    frame = _panel(
        [
            {
                "microbe_label": "Faecalibacterium prausnitzii",
                "species_label": "Faecalibacterium prausnitzii",
                "predicted_effect_label": "inhibit",
                "predicted_inhibit_probability": 0.78,
                "predicted_effect_score": -0.44,
            },
            {
                "microbe_label": "Bifidobacterium longum",
                "species_label": "Bifidobacterium longum",
                "predicted_effect_label": "no_effect",
                "predicted_inhibit_probability": 0.15,
                "predicted_effect_score": 0.01,
            },
            {
                "microbe_label": "Akkermansia muciniphila",
                "species_label": "Akkermansia muciniphila",
                "predicted_effect_label": "no_effect",
                "predicted_inhibit_probability": 0.18,
                "predicted_effect_score": 0.00,
            },
        ]
    )
    payload = evaluate_prediction_confidence(
        effect_frame=frame,
        step1_label_column="predicted_effect_label",
        step1_probability_column="predicted_inhibit_probability",
        step1_score_column="predicted_effect_score",
        drug_profile="unknown",
        molecular_weight=320.0,
        xlogp=2.4,
        mw_bounds=(100.0, 650.0),
        xlogp_bounds=(-2.0, 5.0),
    )
    assert "ecology-risk" in set(payload["warning_flags"])


def test_reasonable_host_pathway_agent_not_mislabeled_over_suppression() -> None:
    frame = _panel(
        [
            {
                "microbe_label": "Bifidobacterium longum",
                "species_label": "Bifidobacterium longum",
                "predicted_effect_label": "no_effect",
                "predicted_inhibit_probability": 0.18,
                "predicted_effect_score": 0.03,
            },
            {
                "microbe_label": "Blautia obeum",
                "species_label": "Blautia obeum",
                "predicted_effect_label": "no_effect",
                "predicted_inhibit_probability": 0.12,
                "predicted_effect_score": 0.01,
            },
            {
                "microbe_label": "Akkermansia muciniphila",
                "species_label": "Akkermansia muciniphila",
                "predicted_effect_label": "promote",
                "predicted_inhibit_probability": 0.08,
                "predicted_effect_score": 0.17,
            },
            {
                "microbe_label": "Escherichia coli",
                "species_label": "Escherichia coli",
                "predicted_effect_label": "inhibit",
                "predicted_inhibit_probability": 0.56,
                "predicted_effect_score": -0.12,
            },
        ]
    )
    payload = evaluate_prediction_confidence(
        effect_frame=frame,
        step1_label_column="predicted_effect_label",
        step1_probability_column="predicted_inhibit_probability",
        step1_score_column="predicted_effect_score",
        drug_profile="host_pathway_agent",
        molecular_weight=300.0,
        xlogp=2.0,
        mw_bounds=(100.0, 650.0),
        xlogp_bounds=(-2.0, 5.0),
    )
    flags = set(payload["warning_flags"])
    assert "over-suppression" not in flags
    assert "drug-profile-conflict" not in flags
