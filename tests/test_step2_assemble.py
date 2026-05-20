from __future__ import annotations

import pandas as pd

from gut_drug_microbiome.step2 import build_step2_input_tables
from gut_drug_microbiome.step2.normalize import NORMALIZED_STEP2_COLUMNS


def test_build_step2_input_tables_creates_candidates_and_labels(tmp_path) -> None:
    step1_predictions = pd.DataFrame(
        [
            {
                "prestwick_id": "P1",
                "nt_code": "NT1",
                "smiles": "CCO",
                "effect_label": "inhibit",
                "binary_effect_label": "inhibit",
                "effect_score": -0.4,
                "predicted_inhibit_probability": 0.9,
                "predicted_binary_effect_label": "inhibit",
                "predicted_effect_score": -0.3,
                "predicted_effect_label_hybrid": "inhibit",
                "predicted_effect_magnitude": 0.3,
            },
            {
                "prestwick_id": "P2",
                "nt_code": "NT2",
                "smiles": "CCC",
                "effect_label": "no_effect",
                "binary_effect_label": "no_effect",
                "effect_score": 0.0,
                "predicted_inhibit_probability": 0.1,
                "predicted_binary_effect_label": "no_effect",
                "predicted_effect_score": 0.05,
                "predicted_effect_label_hybrid": "promote",
                "predicted_effect_magnitude": 0.05,
            },
        ]
    )
    step1_path = tmp_path / "step1_predictions.csv"
    step1_predictions.to_csv(step1_path, index=False)

    label_row = {
        "pair_id": "P1::NT1",
        "prestwick_id": "P1",
        "nt_code": "NT1",
        "drug_name": "Drug 1",
        "microbe_name": "Microbe 1",
        "metabolism_label": "metabolized",
        "reaction_class": "reduction",
        "parent_depletion_fraction": -0.5,
        "product_ids": "CHEBI:1",
        "evidence_gene_ids": "geneA",
        "source_dataset": "zimmermann",
        "label_tier": "gold",
        "source_scope": "isolate",
        "source_record_id": "rec-1",
        "raw_metabolism_label": "yes",
        "raw_reaction_class": "reduction",
    }
    label_frame = pd.DataFrame([[label_row[column] for column in NORMALIZED_STEP2_COLUMNS]], columns=NORMALIZED_STEP2_COLUMNS)
    label_path = tmp_path / "step2_labels.csv"
    label_frame.to_csv(label_path, index=False)

    summary = build_step2_input_tables(step1_path, tmp_path / "output", [label_path])

    assert summary["n_candidate_pairs"] == 2
    assert summary["n_labeled_modeling_rows"] == 1

    modeling = pd.read_csv(tmp_path / "output" / "step2_modeling_table.csv")
    assert set(modeling["pair_id"]) == {"P1::NT1", "P2::NT2"}
    first = modeling.loc[modeling["pair_id"] == "P1::NT1"].iloc[0]
    second = modeling.loc[modeling["pair_id"] == "P2::NT2"].iloc[0]
    assert bool(first["step2_label_available"]) is True
    assert first["step2_metabolism_label"] == "metabolized"
    assert bool(second["step2_label_available"]) is False
