from __future__ import annotations

import pandas as pd

from gut_drug_microbiome.web.service import _enrich_microbe_taxonomy


def test_enrich_microbe_taxonomy_fills_genus_and_taxonomy_from_local_genus_priors() -> None:
    frame = pd.DataFrame(
        [
            {
                "nt_code": "NT1",
                "microbe_label": "Bacteroides vulgatus",
                "species_label": "Bacteroides vulgatus",
                "genus": "Bacteroides",
                "family": "Bacteroidaceae",
                "order": "Bacteroidales",
                "class": "Bacteroidia",
                "phylum": "Bacteroidetes",
                "gram_stain": "negative",
                "medium_preference": "mGAM",
                "starting_od_96_well_screen": 0.08,
            },
            {
                "nt_code": "NT2",
                "microbe_label": "Bacteroides dorei",
                "species_label": "Bacteroides dorei",
                "genus": None,
                "family": None,
                "order": None,
                "class": None,
                "phylum": None,
                "gram_stain": None,
                "medium_preference": None,
                "starting_od_96_well_screen": None,
            },
        ]
    )

    enriched = _enrich_microbe_taxonomy(frame)
    row = enriched.loc[enriched["nt_code"].eq("NT2")].iloc[0]

    assert row["genus"] == "Bacteroides"
    assert row["family"] == "Bacteroidaceae"
    assert row["phylum"] == "Bacteroidetes"
    assert row["gram_stain"] == "negative"
    assert row["medium_preference"] == "mGAM"
    assert float(row["starting_od_96_well_screen"]) == 0.08


def test_enrich_microbe_taxonomy_applies_curated_genus_fallbacks() -> None:
    frame = pd.DataFrame(
        [
            {
                "nt_code": "NT3",
                "microbe_label": "Faecalibacterium prausnitzii",
                "species_label": "Faecalibacterium prausnitzii",
                "genus": None,
                "family": None,
                "order": None,
                "class": None,
                "phylum": None,
                "gram_stain": None,
                "medium_preference": None,
            }
        ]
    )

    enriched = _enrich_microbe_taxonomy(frame)
    row = enriched.iloc[0]

    assert row["genus"] == "Faecalibacterium"
    assert row["family"] == "Ruminococcaceae"
    assert row["phylum"] == "Firmicutes"
    assert row["gram_stain"] == "positive"
