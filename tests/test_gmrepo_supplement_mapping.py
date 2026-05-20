from __future__ import annotations

from scripts.build_gmrepo_disease_microbe_supplement import _canonical_disease_name
from scripts.build_gmrepo_disease_microbe_supplement import _map_detail_row


def test_canonical_disease_name_maps_known_alias() -> None:
    assert _canonical_disease_name("D015179", "Colorectal Neoplasms") == "结直肠癌（CRC）"


def test_map_detail_row_health_vs_disease_direction() -> None:
    row = {
        "scientific_name": "Fusobacterium nucleatum",
        "taxon_rank_level": "species",
        "LDA": 3.2,
        "nrproj": 5,
        "conflict": 0,
        "project_id": "PRJNA12345",
    }
    mapped = _map_detail_row(
        row,
        disease_mesh_id="D015179",
        disease_name="结直肠癌（CRC）",
        mesh1="D006262",  # Health
        mesh2="D015179",  # Disease
    )
    assert mapped is not None
    assert mapped["disease_effect_on_microbe"] == "increase"
    assert mapped["desired_step1_effect"] == "inhibit"
    assert mapped["microbe_role_in_disease"] == "risk"


def test_map_detail_row_disease_vs_health_direction() -> None:
    row = {
        "scientific_name": "Faecalibacterium prausnitzii",
        "taxon_rank_level": "species",
        "LDA": -2.8,
        "nrproj": 4,
        "conflict": 0,
        "project_id": "PRJNA99999",
    }
    mapped = _map_detail_row(
        row,
        disease_mesh_id="D006262",  # now disease is mesh1
        disease_name="Health-like test disease",
        mesh1="D006262",  # Disease in this synthetic case
        mesh2="D999999",
    )
    assert mapped is not None
    assert mapped["disease_effect_on_microbe"] == "increase"
    assert mapped["desired_step1_effect"] == "inhibit"
