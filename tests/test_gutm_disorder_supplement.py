from __future__ import annotations

from scripts.build_gutm_disorder_disease_microbe_supplement import _canonical_project_disease_name
from scripts.build_gutm_disorder_disease_microbe_supplement import _infer_disease_side
from scripts.build_gutm_disorder_disease_microbe_supplement import _is_control_like_condition
from scripts.build_gutm_disorder_disease_microbe_supplement import _map_literature_row
from scripts.build_gutm_disorder_disease_microbe_supplement import _map_raw_row


def test_control_like_condition_does_not_map_to_project_disease() -> None:
    assert _is_control_like_condition("Non-Irritable Bowel Syndrome")
    assert _canonical_project_disease_name("Non-Irritable Bowel Syndrome") == ""


def test_canonical_project_disease_name_maps_project_aliases() -> None:
    assert _canonical_project_disease_name("Active Crohn's Disease") == "克罗恩病（CD）"
    assert _canonical_project_disease_name("Irritable Bowel Syndrome;Constipationt") == "肠易激综合征-便秘型（IBS-C）"
    assert _canonical_project_disease_name("Irritable Bowel Syndrome;Diarrhea") == "肠易激综合征-腹泻型（IBS-D）"


def test_infer_disease_side_skips_same_disease_stage_comparison() -> None:
    assert _infer_disease_side("Active Crohn's Disease", "Inactive Crohn's disease") is None


def test_map_literature_row_maps_direction_against_health_control() -> None:
    row = {
        "PMID": "17897884",
        "Condition1": "Crohn Disease",
        "Condition2": "Health",
        "Condition1ID": "D003424",
        "Condition2ID": "D006262",
        "GutMicrobe": "Enterococcus",
        "Classification": "genus",
        "Alteration": "increase",
        "SequencingTechnology": "16s rRNA gene sequencing",
    }
    mapped = _map_literature_row(row, node_label="Crohn Disease/Health", node_id="2")
    assert mapped is not None
    assert mapped["disease_name"] == "克罗恩病（CD）"
    assert mapped["disease_effect_on_microbe"] == "increase"
    assert mapped["desired_step1_effect"] == "inhibit"
    assert mapped["condition_comparator"] == "Health"


def test_map_literature_row_skips_cross_disease_comparison() -> None:
    row = {
        "PMID": "12345678",
        "Condition1": "Crohn Disease",
        "Condition2": "Colitis, Ulcerative",
        "GutMicrobe": "Faecalibacterium prausnitzii",
        "Classification": "species",
        "Alteration": "decrease",
    }
    assert _map_literature_row(row, node_label="Crohn Disease/Colitis, Ulcerative", node_id="999") is None


def test_map_raw_row_uses_lda_sign_for_direction() -> None:
    row = {
        "ProjectNumber": "PRJDB4176",
        "RunData": "PRJDB4176.txt",
        "Condition1": "Colorectal Neoplasms",
        "Condition2": "Health",
        "Condition1ID": "D015179",
        "Condition2ID": "D006262",
        "GutMicrobe": "Fusobacterium nucleatum",
        "Classification": "species",
        "LDAscore": "2.07939",
        "SequencingTechnology": "Whole metagenomic sequencing",
    }
    mapped = _map_raw_row(row, node_label="Colorectal Neoplasms/Health", node_id="31")
    assert mapped is not None
    assert mapped["disease_name"] == "结直肠癌（CRC）"
    assert mapped["disease_effect_on_microbe"] == "increase"
    assert mapped["desired_step1_effect"] == "inhibit"
    assert mapped["source_project_number"] == "PRJDB4176"
