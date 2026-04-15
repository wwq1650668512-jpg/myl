from __future__ import annotations

import pandas as pd

from gut_drug_microbiome.web.service import evaluate_prediction_confidence


def _panel(rows: list[dict[str, object]]) -> pd.DataFrame:
    return pd.DataFrame(rows)


def test_confidence_flags_over_suppression_and_profile_conflict() -> None:
    frame = _panel(
        [
            {
                "microbe_label": "Faecalibacterium prausnitzii",
                "species_label": "Faecalibacterium prausnitzii",
                "predicted_effect_label": "inhibit",
                "predicted_inhibit_probability": 0.92,
                "predicted_effect_score": -0.7,
            },
            {
                "microbe_label": "Roseburia intestinalis",
                "species_label": "Roseburia intestinalis",
                "predicted_effect_label": "inhibit",
                "predicted_inhibit_probability": 0.88,
                "predicted_effect_score": -0.6,
            },
            {
                "microbe_label": "Eubacterium rectale",
                "species_label": "Eubacterium rectale",
                "predicted_effect_label": "inhibit",
                "predicted_inhibit_probability": 0.86,
                "predicted_effect_score": -0.55,
            },
            {
                "microbe_label": "Escherichia coli",
                "species_label": "Escherichia coli",
                "predicted_effect_label": "inhibit",
                "predicted_inhibit_probability": 0.75,
                "predicted_effect_score": -0.35,
            },
        ]
    )
    payload = evaluate_prediction_confidence(
        effect_frame=frame,
        step1_label_column="predicted_effect_label",
        step1_probability_column="predicted_inhibit_probability",
        step1_score_column="predicted_effect_score",
        drug_profile="eubiotic_modulator",
        molecular_weight=700.0,
        xlogp=6.0,
        mw_bounds=(100.0, 650.0),
        xlogp_bounds=(-2.0, 5.0),
    )
    flags = set(payload["warning_flags"])
    assert "over-suppression" in flags
    assert "core-butyrate-suppression" in flags
    assert "drug-profile-conflict" in flags
    assert "OOD-molecule" in flags
    assert float(payload["confidence_score"]) < 0.45


def test_confidence_flags_ood_only() -> None:
    frame = _panel(
        [
            {
                "microbe_label": "Bifidobacterium longum",
                "species_label": "Bifidobacterium longum",
                "predicted_effect_label": "no_effect",
                "predicted_inhibit_probability": 0.12,
                "predicted_effect_score": 0.05,
            },
            {
                "microbe_label": "Blautia obeum",
                "species_label": "Blautia obeum",
                "predicted_effect_label": "no_effect",
                "predicted_inhibit_probability": 0.18,
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
        molecular_weight=50.0,
        xlogp=0.1,
        mw_bounds=(100.0, 650.0),
        xlogp_bounds=(-2.0, 5.0),
    )
    assert payload["warning_flags"] == ["OOD-molecule"]
    assert 0.6 <= float(payload["confidence_score"]) < 0.9


def test_confidence_high_when_no_risk_flags() -> None:
    frame = _panel(
        [
            {
                "microbe_label": "Bifidobacterium longum",
                "species_label": "Bifidobacterium longum",
                "predicted_effect_label": "no_effect",
                "predicted_inhibit_probability": 0.15,
                "predicted_effect_score": 0.04,
            },
            {
                "microbe_label": "Faecalibacterium prausnitzii",
                "species_label": "Faecalibacterium prausnitzii",
                "predicted_effect_label": "promote",
                "predicted_inhibit_probability": 0.05,
                "predicted_effect_score": 0.22,
            },
            {
                "microbe_label": "Roseburia intestinalis",
                "species_label": "Roseburia intestinalis",
                "predicted_effect_label": "no_effect",
                "predicted_inhibit_probability": 0.10,
                "predicted_effect_score": 0.01,
            },
        ]
    )
    payload = evaluate_prediction_confidence(
        effect_frame=frame,
        step1_label_column="predicted_effect_label",
        step1_probability_column="predicted_inhibit_probability",
        step1_score_column="predicted_effect_score",
        drug_profile="eubiotic_modulator",
        molecular_weight=350.0,
        xlogp=2.0,
        mw_bounds=(100.0, 650.0),
        xlogp_bounds=(-2.0, 5.0),
    )
    assert payload["warning_flags"] == []
    assert float(payload["confidence_score"]) >= 0.85
    assert payload["confidence_tier"] == "high"


def test_sulfonamide_profile_flags_antifolate_mismatch_when_core_not_suppressed() -> None:
    frame = _panel(
        [
            {
                "microbe_label": "Faecalibacterium prausnitzii",
                "species_label": "Faecalibacterium prausnitzii",
                "predicted_effect_label": "no_effect",
                "predicted_inhibit_probability": 0.12,
                "predicted_effect_score": -0.03,
            },
            {
                "microbe_label": "Roseburia intestinalis",
                "species_label": "Roseburia intestinalis",
                "predicted_effect_label": "no_effect",
                "predicted_inhibit_probability": 0.11,
                "predicted_effect_score": -0.02,
            },
            {
                "microbe_label": "Eubacterium rectale",
                "species_label": "Eubacterium rectale",
                "predicted_effect_label": "no_effect",
                "predicted_inhibit_probability": 0.10,
                "predicted_effect_score": -0.01,
            },
        ]
    )
    payload = evaluate_prediction_confidence(
        effect_frame=frame,
        step1_label_column="predicted_effect_label",
        step1_probability_column="predicted_inhibit_probability",
        step1_score_column="predicted_effect_score",
        drug_profile="sulfonamide_antifolate",
        molecular_weight=380.0,
        xlogp=2.8,
        mw_bounds=(100.0, 650.0),
        xlogp_bounds=(-2.0, 5.0),
    )
    assert "antifolate-mismatch" in set(payload["warning_flags"])


def test_sulfonamide_profile_passes_when_core_butyrate_is_suppressed() -> None:
    frame = _panel(
        [
            {
                "microbe_label": "Faecalibacterium prausnitzii",
                "species_label": "Faecalibacterium prausnitzii",
                "predicted_effect_label": "inhibit",
                "predicted_inhibit_probability": 0.81,
                "predicted_effect_score": -0.43,
            },
            {
                "microbe_label": "Roseburia intestinalis",
                "species_label": "Roseburia intestinalis",
                "predicted_effect_label": "inhibit",
                "predicted_inhibit_probability": 0.77,
                "predicted_effect_score": -0.38,
            },
            {
                "microbe_label": "Eubacterium rectale",
                "species_label": "Eubacterium rectale",
                "predicted_effect_label": "inhibit",
                "predicted_inhibit_probability": 0.73,
                "predicted_effect_score": -0.31,
            },
            {
                "microbe_label": "Escherichia coli",
                "species_label": "Escherichia coli",
                "predicted_effect_label": "no_effect",
                "predicted_inhibit_probability": 0.19,
                "predicted_effect_score": -0.04,
            },
        ]
    )
    payload = evaluate_prediction_confidence(
        effect_frame=frame,
        step1_label_column="predicted_effect_label",
        step1_probability_column="predicted_inhibit_probability",
        step1_score_column="predicted_effect_score",
        drug_profile="sulfonamide_antifolate",
        molecular_weight=380.0,
        xlogp=2.8,
        mw_bounds=(100.0, 650.0),
        xlogp_bounds=(-2.0, 5.0),
    )
    assert "antifolate-mismatch" not in set(payload["warning_flags"])
