from __future__ import annotations

import pandas as pd

from gut_drug_microbiome.web.service import GutPredictionService


def _stub_similarity_service() -> GutPredictionService:
    service = GutPredictionService.__new__(GutPredictionService)
    service.fingerprint_columns = ["morgan_fp_000", "morgan_fp_001", "morgan_fp_002", "morgan_fp_003"]
    service.drug_similarity_table = pd.DataFrame(
        [
            {
                "prestwick_id": "Prestw-A",
                "chemical_name": "Imatinib-like A",
                "murcko_scaffold": "A",
                "morgan_fp_000": 1.0,
                "morgan_fp_001": 1.0,
                "morgan_fp_002": 1.0,
                "morgan_fp_003": 0.0,
            },
            {
                "prestwick_id": "Prestw-B",
                "chemical_name": "CRC-like B",
                "murcko_scaffold": "B",
                "morgan_fp_000": 1.0,
                "morgan_fp_001": 0.0,
                "morgan_fp_002": 0.0,
                "morgan_fp_003": 1.0,
            },
            {
                "prestwick_id": "Prestw-C",
                "chemical_name": "Unrelated C",
                "murcko_scaffold": "C",
                "morgan_fp_000": 0.0,
                "morgan_fp_001": 0.0,
                "morgan_fp_002": 1.0,
                "morgan_fp_003": 1.0,
            },
        ]
    )
    service.drug_similarity_matrix = (
        service.drug_similarity_table.loc[:, service.fingerprint_columns].fillna(0).to_numpy(dtype=bool)
    )
    service.biotransform_reference_table = pd.DataFrame(
        [
            {
                "prestwick_id": "Prestw-A",
                "chemical_name": "Imatinib-like A",
                "prestwick_key": "prestwa",
                "chemical_name_key": "imatiniblikea",
                "experimental_biotransform_product_count": 4.0,
                "experimental_biotransform_fraction_in_gut": 0.62,
                "experimental_biotransform_product_ids": "P1;P2",
                "experimental_biotransform_ec_numbers": "1.1.1.1",
                "experimental_biotransform_reaction_centers": "azo_reduction",
                "experimental_biotransform_primary_product_name": "Metabolite A",
                "morgan_fp_000": 1.0,
                "morgan_fp_001": 1.0,
                "morgan_fp_002": 1.0,
                "morgan_fp_003": 0.0,
            },
            {
                "prestwick_id": "Prestw-B",
                "chemical_name": "CRC-like B",
                "prestwick_key": "prestwb",
                "chemical_name_key": "crclikeb",
                "experimental_biotransform_product_count": 2.0,
                "experimental_biotransform_fraction_in_gut": 0.20,
                "experimental_biotransform_product_ids": "Q1",
                "experimental_biotransform_ec_numbers": "2.2.2.2",
                "experimental_biotransform_reaction_centers": "hydroxylation",
                "experimental_biotransform_primary_product_name": "Metabolite B",
                "morgan_fp_000": 1.0,
                "morgan_fp_001": 0.0,
                "morgan_fp_002": 0.0,
                "morgan_fp_003": 1.0,
            },
        ]
    )
    service.biotransform_reference_matrix = (
        service.biotransform_reference_table.loc[:, service.fingerprint_columns].fillna(0).to_numpy(dtype=bool)
    )
    service.disease_drug_reference = pd.DataFrame(
        [
            {
                "disease_name": "胃肠道间质瘤（GIST）",
                "marketed_drug_name_raw": "Imatinib-like A",
                "marketed_drug_key": "imatiniblikea",
            },
            {
                "disease_name": "结直肠癌（CRC）",
                "marketed_drug_name_raw": "CRC-like B",
                "marketed_drug_key": "crclikeb",
            },
        ]
    )
    return service


def test_similarity_disease_context_prefers_most_similar_marketed_disease() -> None:
    service = _stub_similarity_service()
    query = pd.Series(
        {
            "prestwick_id": "Custom-1",
            "chemical_name": "Custom kinase inhibitor",
            "morgan_fp_000": 1.0,
            "morgan_fp_001": 1.0,
            "morgan_fp_002": 1.0,
            "morgan_fp_003": 1.0,
        }
    )

    similar_drugs = service._top_similar_library_drugs(query, top_k=3, min_similarity=0.2)
    assert similar_drugs[0]["chemical_name"] == "Imatinib-like A"
    assert float(similar_drugs[0]["tanimoto_similarity"]) > float(similar_drugs[1]["tanimoto_similarity"])

    disease_context = service._similarity_disease_context(query)
    assert disease_context[0]["disease_name"] == "胃肠道间质瘤（GIST）"
    assert "Imatinib-like A" in disease_context[0]["matched_market_drugs"]


def test_gist_reference_files_contain_expected_entries() -> None:
    marketed = pd.read_csv("data/reference/disease_marketed_drug_catalog.csv")
    gist_drugs = marketed.loc[marketed["disease_name"].astype(str).eq("胃肠道间质瘤（GIST）"), "marketed_drug_key"]
    assert {"imatinib", "avapritinib", "sunitinib", "regorafenib", "ripretinib"}.issubset(set(gist_drugs))

    gist_microbes = pd.read_csv("data/reference/disease_microbe_gist_supplement.csv")
    assert gist_microbes["disease_name"].astype(str).eq("胃肠道间质瘤（GIST）").all()
    assert {"Proteobacteria", "Prevotella"}.issubset(set(gist_microbes["microbe_name_raw"].astype(str)))


def test_biotransform_sidecar_uses_direct_reference_when_available() -> None:
    service = _stub_similarity_service()
    query = pd.Series(
        {
            "prestwick_id": "Prestw-A",
            "chemical_name": "Imatinib-like A",
            "morgan_fp_000": 1.0,
            "morgan_fp_001": 1.0,
            "morgan_fp_002": 1.0,
            "morgan_fp_003": 0.0,
        }
    )

    payload = service._biotransform_sidecar_from_row(query)
    assert payload["biotransform_sidecar_enabled"] is True
    assert payload["biotransform_sidecar_mode"] == "direct_reference"
    assert float(payload["biotransform_sidecar_weighted_product_count"]) == 4.0
    assert payload["biotransform_sidecar_reference_drugs"] == "Imatinib-like A"


def test_biotransform_sidecar_similarity_transfer_for_custom_query() -> None:
    service = _stub_similarity_service()
    query = pd.Series(
        {
            "prestwick_id": "Custom-1",
            "chemical_name": "Custom kinase inhibitor",
            "morgan_fp_000": 1.0,
            "morgan_fp_001": 1.0,
            "morgan_fp_002": 1.0,
            "morgan_fp_003": 1.0,
        }
    )

    payload = service._biotransform_sidecar_from_row(query, min_similarity=0.2)
    assert payload["biotransform_sidecar_enabled"] is True
    assert payload["biotransform_sidecar_mode"] == "similarity_transfer"
    assert float(payload["biotransform_sidecar_support_score"]) > 0.0
    assert "Imatinib-like A" in str(payload["biotransform_sidecar_reference_drugs"])
