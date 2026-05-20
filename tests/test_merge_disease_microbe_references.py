from __future__ import annotations

import json

import pandas as pd

from scripts.merge_disease_microbe_references import merge_disease_microbe_references


def test_merge_disease_microbe_references_deduplicates_aliases(tmp_path) -> None:
    primary_path = tmp_path / "primary.csv"
    supplement_path = tmp_path / "supp.csv"
    output_path = tmp_path / "merged.csv"
    summary_path = tmp_path / "merged.summary.json"

    pd.DataFrame(
        [
            {
                "source_sheet": "microbe_to_disease",
                "disease_name": "结肠癌",
                "microbe_name_raw": "Fusobacterium nucleatum",
                "taxon_level": "species",
                "desired_step1_effect": "inhibit",
                "relation_confidence": "high",
            }
        ]
    ).to_csv(primary_path, index=False)

    pd.DataFrame(
        [
            {
                "source_sheet": "gmrepo_health_vs_disease",
                "disease_name": "Colorectal Neoplasms",
                "microbe_name_raw": "Fusobacterium nucleatum",
                "taxon_level": "species",
                "desired_step1_effect": "inhibit",
                "relation_confidence": "medium",
                "marker_nr_projects": 5,
            }
        ]
    ).to_csv(supplement_path, index=False)

    summary = merge_disease_microbe_references(
        primary_path=primary_path,
        supplement_paths=[supplement_path],
        output_path=output_path,
        summary_path=summary_path,
    )

    merged = pd.read_csv(output_path)
    assert len(merged) == 1
    assert merged.loc[0, "disease_name"] == "结直肠癌（CRC）"
    assert merged.loc[0, "relation_confidence"] == "high"
    assert summary["n_rows"] == 1
    loaded_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert loaded_summary["n_rows"] == 1
