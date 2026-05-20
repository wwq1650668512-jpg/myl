from __future__ import annotations

import pandas as pd

from gut_drug_microbiome.step1.compound_semantics import annotate_compound_semantics
from gut_drug_microbiome.step1.hybrid import _apply_step1_drug_profile_constraints


def test_compound_semantics_recognizes_sulfonamide_antifolate_family() -> None:
    frame = pd.DataFrame(
        [
            {"chemical_name": "Sulfasalazine", "therapeutic_class": "anti-inflammatory", "therapeutic_effect": "IBD"},
            {"chemical_name": "Sulfapyridine", "therapeutic_class": "sulfonamide antibiotic", "therapeutic_effect": ""},
        ]
    )
    annotated = annotate_compound_semantics(frame)
    families = set(annotated["compound_semantic_family"].fillna("").astype(str))
    assert "azo_prodrug_sulfonamide" in families or "sulfonamide_antifolate" in families
    assert "sulfonamide_antifolate" in families


def test_compound_semantics_extracts_structure_keywords_for_custom_smiles() -> None:
    frame = pd.DataFrame(
        [
            {
                "chemical_name": "custom_test_input",
                "therapeutic_class": "",
                "therapeutic_effect": "",
                "smiles": "CC(=O)NC1=CC=CC=C1",
            }
        ]
    )
    annotated = annotate_compound_semantics(frame)
    keywords = str(annotated.loc[0, "compound_semantic_keywords"] or "")
    assert "amide" in keywords
    assert "lactam" in keywords


def test_antifolate_constraint_boosts_core_butyrate_inhibit_pressure() -> None:
    frame = pd.DataFrame(
        [
            {
                "chemical_name": "Sulfasalazine",
                "therapeutic_class": "sulfonamide",
                "therapeutic_effect": "antifolate antibacterial",
                "compound_semantic_family": "sulfonamide_antifolate",
                "species_label": "Faecalibacterium prausnitzii",
                "microbe_label": "Faecalibacterium prausnitzii",
                "genus": "Faecalibacterium",
                "medium_preference": "strict anaerobe",
                "predicted_inhibit_probability": 0.22,
                "predicted_effect_score": -0.05,
            },
            {
                "chemical_name": "Sulfasalazine",
                "therapeutic_class": "sulfonamide",
                "therapeutic_effect": "antifolate antibacterial",
                "compound_semantic_family": "sulfonamide_antifolate",
                "species_label": "Escherichia coli",
                "microbe_label": "Escherichia coli",
                "genus": "Escherichia",
                "medium_preference": "facultative anaerobe",
                "predicted_inhibit_probability": 0.20,
                "predicted_effect_score": -0.03,
            },
        ]
    )
    constrained, summary = _apply_step1_drug_profile_constraints(
        frame,
        inhibit_probability_threshold=0.5,
        promote_score_threshold=0.2,
    )

    assert constrained.loc[0, "step1_drug_profile"] == "sulfonamide_antifolate"
    assert constrained.loc[0, "predicted_inhibit_probability"] > 0.40
    assert constrained.loc[0, "step1_constraint_reason"] == "antifolate_core_folate_vulnerability_boost"
    assert constrained.loc[0, "predicted_folate_vulnerability_score"] >= 0.55

    assert constrained.loc[1, "predicted_inhibit_probability"] <= constrained.loc[0, "predicted_inhibit_probability"]
    assert int(summary["antifolate_sensitive_rows_boosted"]) >= 1
