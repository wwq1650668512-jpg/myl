from __future__ import annotations

import pandas as pd

from gut_drug_microbiome.web.service import GutPredictionService


def _stub_service(disease_microbe_reference: pd.DataFrame) -> GutPredictionService:
    service = GutPredictionService.__new__(GutPredictionService)
    service.disease_microbe_reference = disease_microbe_reference
    service.disease_drug_reference = pd.DataFrame()
    service.step1_score_column = "predicted_effect_score"
    service.step1_probability_column = "predicted_inhibit_probability"
    service.step1_label_column = "predicted_effect_label"
    return service


def _single_microbe_work(*, applicability_flag: bool, scaffold_seen: bool, jaccard: float) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "nt_code": "NT_PREVOTELLA",
                "microbe_label": "Prevotella copri",
                "species_label": "Prevotella copri",
                "genus": "Prevotella",
                "family": "Prevotellaceae",
                "phylum": "Bacteroidetes",
                "display_step1_predicted_effect_label": "inhibit",
                "display_step1_predicted_effect_score": -0.85,
                "display_step1_predicted_inhibit_probability": 0.91,
                "predicted_promote_probability_refined": 0.04,
                "applicability_flag": applicability_flag,
                "scaffold_seen_in_training": scaffold_seen,
                "drug_max_fingerprint_jaccard": jaccard,
            }
        ]
    )


def _find_disease(rows: list[dict[str, object]], disease_name: str) -> dict[str, object]:
    for row in rows:
        if row.get("disease_name") == disease_name:
            return row
    raise AssertionError(f"Disease {disease_name!r} not found in candidate rows.")


def test_candidate_disease_relations_are_deduplicated() -> None:
    disease_reference = pd.DataFrame(
        [
            {
                "disease_name": "痔疮",
                "microbe_name_raw": "Prevotella",
                "genus_hint": "Prevotella",
                "taxon_level": "genus",
                "desired_step1_effect": "inhibit",
                "source_sheet": "microbe_to_disease",
                "relation_confidence": "high",
            },
            {
                "disease_name": "痔疮",
                "microbe_name_raw": "Prevotella",
                "genus_hint": "Prevotella",
                "taxon_level": "genus",
                "desired_step1_effect": "inhibit",
                "source_sheet": "disease_to_microbe",
                "relation_confidence": "medium",
            },
        ]
    )
    service = _stub_service(disease_reference)
    work = _single_microbe_work(applicability_flag=True, scaffold_seen=True, jaccard=0.95)

    rows = service._candidate_diseases_from_frame(work)
    hemorrhoid = _find_disease(rows, "痔疮")

    assert hemorrhoid["matched_relation_count"] == 1
    assert len(hemorrhoid["evidence_examples"]) == 1


def test_unseen_smiles_reduce_candidate_support_score() -> None:
    disease_reference = pd.DataFrame(
        [
            {
                "disease_name": "痔疮",
                "microbe_name_raw": "Prevotella",
                "genus_hint": "Prevotella",
                "taxon_level": "genus",
                "desired_step1_effect": "inhibit",
                "source_sheet": "microbe_to_disease",
                "relation_confidence": "high",
            }
        ]
    )
    service = _stub_service(disease_reference)

    seen_rows = service._candidate_diseases_from_frame(
        _single_microbe_work(applicability_flag=True, scaffold_seen=True, jaccard=0.95)
    )
    unseen_rows = service._candidate_diseases_from_frame(
        _single_microbe_work(applicability_flag=False, scaffold_seen=False, jaccard=0.10)
    )
    seen = _find_disease(seen_rows, "痔疮")
    unseen = _find_disease(unseen_rows, "痔疮")

    assert float(unseen["support_score"]) < float(seen["support_score"])
    assert float(unseen["disease_score_raw_only"]) < float(seen["disease_score_raw_only"])
