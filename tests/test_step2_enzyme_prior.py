from __future__ import annotations

import pandas as pd

from gut_drug_microbiome.step2.enzyme_prior import annotate_step2_with_enzyme_priors
from gut_drug_microbiome.step2.enzyme_prior import build_step2_enzyme_curation_template
from gut_drug_microbiome.step2.enzyme_prior import build_step2_enzyme_reference_tables


def test_build_step2_enzyme_reference_tables_creates_seeded_outputs(tmp_path) -> None:
    microbe_table = pd.DataFrame(
        [
            {"nt_code": "NT1", "microbe_label": "Bacteroides uniformis", "species_label": "Bacteroides uniformis", "genus": "Bacteroides"},
            {"nt_code": "NT2", "microbe_label": "Escherichia coli", "species_label": "Escherichia coli", "genus": "Escherichia"},
        ]
    )
    summary = build_step2_enzyme_reference_tables(
        microbe_table_path=microbe_table,
        enzyme_catalog_path=tmp_path / "catalog.csv",
        microbe_enzyme_long_path=tmp_path / "long.csv",
        microbe_enzyme_matrix_path=tmp_path / "matrix.csv",
        summary_path=tmp_path / "summary.json",
    )

    long_table = pd.read_csv(tmp_path / "long.csv")
    matrix = pd.read_csv(tmp_path / "matrix.csv")

    assert summary["n_microbes"] == 2
    assert {"NT1", "NT2"} == set(long_table["nt_code"])
    assert "ENZ001" in set(long_table["enzyme_id"])
    assert "ENZ003" in set(long_table.loc[long_table["nt_code"] == "NT2", "enzyme_id"])
    assert "ENZ001" in matrix.columns


def test_build_step2_enzyme_curation_template_creates_full_microbe_enzyme_grid(tmp_path) -> None:
    microbe_table = pd.DataFrame(
        [
            {"nt_code": "NT1", "microbe_label": "Bacteroides uniformis", "species_label": "Bacteroides uniformis", "genus": "Bacteroides"},
            {"nt_code": "NT2", "microbe_label": "Escherichia coli", "species_label": "Escherichia coli", "genus": "Escherichia"},
        ]
    )
    summary = build_step2_enzyme_curation_template(microbe_table, output_path=tmp_path / "template.csv")
    template = pd.read_csv(tmp_path / "template.csv")

    assert summary["n_microbes"] == 2
    assert summary["n_enzymes"] == 16
    assert len(template) == 32
    assert {"starter_presence_call", "curated_presence_call", "curated_pmid"} <= set(template.columns)
    assert template["review_priority"].isin(["high", "medium"]).all()


def test_species_literature_evidence_overrides_genus_prior(tmp_path) -> None:
    microbe_table = pd.DataFrame(
        [
            {
                "nt_code": "NT1",
                "microbe_label": "Bacteroides uniformis",
                "species_label": "Bacteroides uniformis",
                "species": "Bacteroides uniformis",
                "strain": "ATCC 8492",
                "genus": "Bacteroides",
                "family": "Bacteroidaceae",
                "phylum": "Bacteroidetes",
            }
        ]
    )
    curated = pd.DataFrame(
        [
            {
                "nt_code": "NT1",
                "enzyme_id": "ENZ001",
                "presence_call": "absent",
                "evidence_scope": "species_literature",
                "evidence_source": "manual_literature_curation",
                "literature_citation": "Example et al. 2024",
                "pmid": "12345678",
                "strain_match_level": "same_species",
                "evidence_note": "Species-level review found no detectable beta-glucuronidase for this strain set.",
                "curation_status": "reviewed_accepted",
            }
        ]
    )
    build_step2_enzyme_reference_tables(
        microbe_table_path=microbe_table,
        enzyme_catalog_path=tmp_path / "catalog.csv",
        microbe_enzyme_long_path=tmp_path / "long.csv",
        microbe_enzyme_matrix_path=tmp_path / "matrix.csv",
        evidence_ledger_path=tmp_path / "ledger.csv",
        curation_template_path=tmp_path / "template.csv",
        literature_evidence_path=curated,
        summary_path=None,
    )

    long_table = pd.read_csv(tmp_path / "long.csv")
    row = long_table.loc[(long_table["nt_code"] == "NT1") & (long_table["enzyme_id"] == "ENZ001")].iloc[0]

    assert row["evidence_scope"] == "species_literature"
    assert row["evidence_source"] == "manual_literature_curation"
    assert row["presence_call"] == "absent"
    assert float(row["presence_weight"]) == 0.0


def test_annotate_step2_with_enzyme_priors_matches_rutin_like_pair(tmp_path) -> None:
    microbe_table = pd.DataFrame(
        [{"nt_code": "NT1", "microbe_label": "Bacteroides uniformis", "species_label": "Bacteroides uniformis", "genus": "Bacteroides"}]
    )
    build_step2_enzyme_reference_tables(
        microbe_table_path=microbe_table,
        enzyme_catalog_path=tmp_path / "catalog.csv",
        microbe_enzyme_long_path=tmp_path / "long.csv",
        microbe_enzyme_matrix_path=tmp_path / "matrix.csv",
        summary_path=None,
    )

    pair_table = pd.DataFrame(
        [
            {
                "prestwick_id": "P1",
                "nt_code": "NT1",
                "chemical_name": "Rutin",
                "therapeutic_class": "",
                "therapeutic_effect": "",
                "compound_semantic_family": "flavonoid_glycoside",
                "compound_semantic_aliases": "rutin;rutinoside",
                "compound_semantic_keywords": "glycoside;rutinoside",
            }
        ]
    )
    annotated = annotate_step2_with_enzyme_priors(
        pair_table,
        microbe_enzyme_panel_path=tmp_path / "long.csv",
        enzyme_catalog_path=tmp_path / "catalog.csv",
    )
    row = annotated.iloc[0]
    assert bool(row["predicted_enzyme_prior_flag"]) is True
    assert int(row["predicted_enzyme_match_count"]) >= 1
    assert "ENZ008" in row["predicted_enzyme_ids"] or "ENZ007" in row["predicted_enzyme_ids"]
    assert float(row["predicted_enzyme_support_score"]) > 0


def test_annotate_step2_with_enzyme_priors_uses_smiles_derived_keywords(tmp_path) -> None:
    microbe_table = pd.DataFrame(
        [{"nt_code": "NT1", "microbe_label": "Bacteroides uniformis", "species_label": "Bacteroides uniformis", "genus": "Bacteroides"}]
    )
    build_step2_enzyme_reference_tables(
        microbe_table_path=microbe_table,
        enzyme_catalog_path=tmp_path / "catalog.csv",
        microbe_enzyme_long_path=tmp_path / "long.csv",
        microbe_enzyme_matrix_path=tmp_path / "matrix.csv",
        summary_path=None,
    )

    pair_table = pd.DataFrame(
        [
            {
                "prestwick_id": "P2",
                "nt_code": "NT1",
                "chemical_name": "custom_test_input",
                "smiles": "CC(=O)NC1=CC=CC=C1",
            }
        ]
    )
    annotated = annotate_step2_with_enzyme_priors(
        pair_table,
        microbe_enzyme_panel_path=tmp_path / "long.csv",
        enzyme_catalog_path=tmp_path / "catalog.csv",
    )
    row = annotated.iloc[0]
    assert bool(row["predicted_enzyme_prior_flag"]) is True
    assert "ENZ006" in str(row["predicted_enzyme_ids"])
    assert "amidase" in str(row["predicted_enzyme_names"])
    assert "amide" in str(row["compound_semantic_keywords"])
