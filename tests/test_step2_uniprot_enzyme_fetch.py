from __future__ import annotations

import pandas as pd

from gut_drug_microbiome.step2.uniprot_enzyme_fetch import _match_uniprot_record_to_project_enzymes
from gut_drug_microbiome.step2.uniprot_enzyme_fetch import build_uniprot_enzyme_query
from gut_drug_microbiome.step2.uniprot_enzyme_fetch import build_uniprot_taxonomy_query
from gut_drug_microbiome.step2.uniprot_enzyme_fetch import candidate_taxon_names_for_microbe


def test_candidate_taxon_names_for_microbe_normalizes_species_strings() -> None:
    row = pd.Series(
        {
            "species_name": "Bacteroides vulgatus (S1) DSM No.: 1447",
            "species_label": "Bacteroides vulgatus",
            "species": "Bacteroides vulgatus",
            "microbe_label": "Bacteroides vulgatus",
            "genus": "Bacteroides",
        }
    )
    candidates = candidate_taxon_names_for_microbe(row)

    assert candidates[0] == "Bacteroides vulgatus"
    assert "Bacteroides vulgatus" in candidates


def test_build_uniprot_taxonomy_query_includes_requested_names() -> None:
    query = build_uniprot_taxonomy_query(["Bacteroides vulgatus", "Escherichia coli"])

    assert 'VALUES ?query { "Bacteroides vulgatus" "Escherichia coli" }' in query
    assert "up:scientificName" in query


def test_build_uniprot_enzyme_query_targets_taxon_and_enzyme_links() -> None:
    query = build_uniprot_enzyme_query("http://purl.uniprot.org/taxonomy/821")

    assert "http://purl.uniprot.org/taxonomy/821" in query
    assert "up:enzyme" in query
    assert "rdfs:subClassOf" in query


def test_match_uniprot_record_to_project_enzymes_uses_ec_and_name_keywords() -> None:
    exact_matches = _match_uniprot_record_to_project_enzymes("3.2.1.31", "Beta-glucuronidase")
    keyword_matches = _match_uniprot_record_to_project_enzymes("", "Choloylglycine hydrolase family protein")

    assert ("ENZ001", "exact_ec") in exact_matches
    assert ("ENZ010", "name_keyword") in keyword_matches
