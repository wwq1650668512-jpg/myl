from __future__ import annotations

import pandas as pd

from gut_drug_microbiome.web.service import GutPredictionService


def _stub_service() -> GutPredictionService:
    service = GutPredictionService.__new__(GutPredictionService)
    service.disease_microbe_reference = pd.DataFrame()
    service.disease_drug_reference = pd.DataFrame()
    return service


def test_canonicalize_disease_name_aliases() -> None:
    service = _stub_service()

    assert service._canonicalize_disease_name("克罗恩病") == "克罗恩病（CD）"
    assert service._canonicalize_disease_name("溃疡性结肠炎") == "溃疡性结肠炎（UC）"
    assert service._canonicalize_disease_name("结肠癌") == "结直肠癌（CRC）"
    assert service._canonicalize_disease_name("便秘") == "便秘（Constipation）"
    assert service._canonicalize_disease_name("腹泻") == "腹泻（Diarrhea）"
    assert service._canonicalize_disease_name("腹泻型肠易激综合征") == "肠易激综合征-腹泻型（IBS-D）"
    assert service._canonicalize_disease_name("便秘型肠易激综合征") == "肠易激综合征-便秘型（IBS-C）"


def test_build_disease_catalog_merges_alias_names() -> None:
    service = _stub_service()
    service.disease_microbe_reference = pd.DataFrame(
        [
            {"disease_name": "克罗恩病", "microbe_name_raw": "A", "taxon_level": "genus", "desired_step1_effect": "inhibit"},
            {"disease_name": "克罗恩病（CD）", "microbe_name_raw": "B", "taxon_level": "genus", "desired_step1_effect": "inhibit"},
            {"disease_name": "便秘", "microbe_name_raw": "C", "taxon_level": "genus", "desired_step1_effect": "promote"},
            {"disease_name": "便秘（Constipation）", "microbe_name_raw": "D", "taxon_level": "genus", "desired_step1_effect": "promote"},
        ]
    )
    service.disease_drug_reference = pd.DataFrame(
        [
            {"disease_name": "结肠癌", "marketed_drug_name_raw": "X", "marketed_drug_key": "x"},
            {"disease_name": "结直肠癌（CRC）", "marketed_drug_name_raw": "Y", "marketed_drug_key": "y"},
        ]
    )

    service._normalize_and_expand_disease_references()
    catalog = service._build_disease_catalog()
    names = {item["disease_name"] for item in catalog}

    assert "克罗恩病（CD）" in names
    assert "克罗恩病" not in names
    assert "便秘（Constipation）" in names
    assert "便秘" not in names
    assert "结直肠癌（CRC）" in names
    assert "结肠癌" not in names


def test_load_disease_microbe_reference_bundle_includes_supplements(tmp_path) -> None:
    service = _stub_service()
    primary_path = tmp_path / "primary.csv"
    supplement_path = tmp_path / "supplement.csv"

    pd.DataFrame(
        [
            {
                "source_sheet": "microbe_to_disease",
                "disease_name": "克罗恩病（CD）",
                "microbe_name_raw": "Faecalibacterium prausnitzii",
                "taxon_level": "species",
                "desired_step1_effect": "promote",
                "relation_confidence": "high",
                "mechanism_note": "primary",
            }
        ]
    ).to_csv(primary_path, index=False)
    pd.DataFrame(
        [
            {
                "source_sheet": "gmrepo_health_vs_disease",
                "disease_name": "Crohn Disease",
                "microbe_name_raw": "Escherichia coli",
                "taxon_level": "species",
                "desired_step1_effect": "inhibit",
                "relation_confidence": "medium",
                "mechanism_note": "supplement",
            }
        ]
    ).to_csv(supplement_path, index=False)

    merged = service._load_disease_microbe_reference_bundle(
        primary_path=primary_path,
        supplement_paths=[supplement_path],
    )

    assert len(merged) == 2
    assert set(merged["source_sheet"]) == {"microbe_to_disease", "gmrepo_health_vs_disease"}
