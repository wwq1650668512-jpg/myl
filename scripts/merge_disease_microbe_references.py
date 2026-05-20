from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.utils.text import canonicalize_key as _canonicalize_key
from gut_drug_microbiome.utils.text import normalize_whitespace as _clean_text

CONFIDENCE_RANK = {"low": 1, "medium": 2, "high": 3}
DISEASE_ALIAS_MAP = {
    _canonicalize_key("克罗恩病", keep_cjk=True): "克罗恩病（CD）",
    _canonicalize_key("Crohn Disease", keep_cjk=True): "克罗恩病（CD）",
    _canonicalize_key("溃疡性结肠炎", keep_cjk=True): "溃疡性结肠炎（UC）",
    _canonicalize_key("Ulcerative Colitis", keep_cjk=True): "溃疡性结肠炎（UC）",
    _canonicalize_key("结肠癌", keep_cjk=True): "结直肠癌（CRC）",
    _canonicalize_key("Colorectal Neoplasms", keep_cjk=True): "结直肠癌（CRC）",
    _canonicalize_key("CRC", keep_cjk=True): "结直肠癌（CRC）",
    _canonicalize_key("便秘", keep_cjk=True): "便秘（Constipation）",
    _canonicalize_key("Constipation", keep_cjk=True): "便秘（Constipation）",
    _canonicalize_key("腹泻", keep_cjk=True): "腹泻（Diarrhea）",
    _canonicalize_key("Diarrhea", keep_cjk=True): "腹泻（Diarrhea）",
    _canonicalize_key("炎症性肠病", keep_cjk=True): "炎症性肠病（IBD）",
    _canonicalize_key("Inflammatory Bowel Diseases", keep_cjk=True): "炎症性肠病（IBD）",
    _canonicalize_key("肠易激综合征", keep_cjk=True): "肠易激综合征（IBS）",
    _canonicalize_key("Irritable Bowel Syndrome", keep_cjk=True): "肠易激综合征（IBS）",
}


def _canonicalize_disease_name(value: object) -> str:
    text = _clean_text(value)
    if not text:
        return ""
    mapped = DISEASE_ALIAS_MAP.get(_canonicalize_key(text, keep_cjk=True))
    return mapped or text


def _infer_taxon_level(value: object, microbe_name: object) -> str:
    text = _clean_text(value).lower()
    if text in {"species", "genus", "family", "order", "class", "phylum"}:
        return text
    microbe = _clean_text(microbe_name)
    if " " in microbe:
        return "species"
    return "genus" if microbe else "unknown"


def _normalize_desired_effect(value: object, disease_effect: object) -> str:
    effect = _clean_text(value).lower()
    if effect in {"promote", "inhibit", "unknown"}:
        return effect
    direction = _clean_text(disease_effect).lower()
    if direction == "increase":
        return "inhibit"
    if direction == "decrease":
        return "promote"
    return "unknown"


def _load_frame(path: Path, *, source_priority: int) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    frame = pd.read_csv(path, low_memory=False)
    if frame.empty:
        return frame
    work = frame.copy()
    work["disease_name"] = work.get("disease_name", pd.Series("", index=work.index)).map(_canonicalize_disease_name)
    work["microbe_name_raw"] = work.get("microbe_name_raw", pd.Series("", index=work.index)).map(_clean_text)
    work["taxon_level"] = [
        _infer_taxon_level(level, microbe)
        for level, microbe in zip(
            work.get("taxon_level", pd.Series("", index=work.index)),
            work.get("microbe_name_raw", pd.Series("", index=work.index)),
        )
    ]
    work["microbe_key"] = work.get("microbe_key", pd.Series("", index=work.index)).map(_clean_text)
    missing_key_mask = work["microbe_key"].eq("")
    work.loc[missing_key_mask, "microbe_key"] = work.loc[missing_key_mask, "microbe_name_raw"].map(
        lambda item: _canonicalize_key(item, keep_cjk=True)
    )
    work["desired_step1_effect"] = [
        _normalize_desired_effect(effect, direction)
        for effect, direction in zip(
            work.get("desired_step1_effect", pd.Series("", index=work.index)),
            work.get("disease_effect_on_microbe", pd.Series("", index=work.index)),
        )
    ]
    work["relation_confidence"] = work.get("relation_confidence", pd.Series("low", index=work.index)).map(_clean_text).str.lower()
    work.loc[~work["relation_confidence"].isin(CONFIDENCE_RANK), "relation_confidence"] = "low"
    if "source_database" not in work.columns:
        work["source_database"] = path.stem
    work["_merge_source_path"] = str(path)
    work["_merge_source_priority"] = int(source_priority)
    return work


def merge_disease_microbe_references(
    *,
    primary_path: Path,
    supplement_paths: list[Path],
    output_path: Path,
    summary_path: Path,
) -> dict[str, object]:
    paths = [primary_path] + supplement_paths
    frames: list[pd.DataFrame] = []
    priority = len(paths)
    for path in paths:
        frame = _load_frame(path, source_priority=priority)
        if not frame.empty:
            frames.append(frame)
        priority -= 1

    if frames:
        merged = pd.concat(frames, ignore_index=True, sort=False)
    else:
        merged = pd.DataFrame()

    if not merged.empty:
        merged["_confidence_rank"] = merged["relation_confidence"].map(CONFIDENCE_RANK).fillna(0).astype(int)
        merged["_abs_lda_score"] = pd.to_numeric(merged.get("lda_score", pd.Series(0.0, index=merged.index)), errors="coerce").abs().fillna(0.0)
        merged["_marker_nr_projects"] = pd.to_numeric(
            merged.get("marker_nr_projects", pd.Series(0, index=merged.index)),
            errors="coerce",
        ).fillna(0.0)
        merged["_has_mechanism_note"] = merged.get("mechanism_note", pd.Series("", index=merged.index)).map(_clean_text).ne("")
        merged.sort_values(
            [
                "_confidence_rank",
                "_merge_source_priority",
                "_marker_nr_projects",
                "_abs_lda_score",
                "_has_mechanism_note",
            ],
            ascending=[False, False, False, False, False],
            inplace=True,
        )
        dedup_keys = ["disease_name", "microbe_key", "taxon_level", "desired_step1_effect"]
        merged = merged.drop_duplicates(subset=dedup_keys, keep="first").reset_index(drop=True)
        merged.drop(
            columns=[
                "_confidence_rank",
                "_abs_lda_score",
                "_marker_nr_projects",
                "_has_mechanism_note",
                "_merge_source_priority",
            ],
            inplace=True,
            errors="ignore",
        )
        if "reference_id" in merged.columns:
            merged = merged.drop(columns=["reference_id"])
        merged.insert(0, "reference_id", [f"MDR_{index + 1:05d}" for index in range(len(merged))])
    else:
        merged = pd.DataFrame()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_path, index=False)

    summary = {
        "primary_path": str(primary_path),
        "supplement_paths": [str(path) for path in supplement_paths],
        "output_path": str(output_path),
        "n_rows": int(len(merged)),
        "n_diseases": int(merged["disease_name"].nunique()) if not merged.empty and "disease_name" in merged.columns else 0,
        "desired_step1_effect_counts": (
            {str(key): int(value) for key, value in merged["desired_step1_effect"].value_counts(dropna=False).to_dict().items()}
            if not merged.empty and "desired_step1_effect" in merged.columns
            else {}
        ),
        "source_database_counts": (
            {str(key): int(value) for key, value in merged["source_database"].value_counts(dropna=False).to_dict().items()}
            if not merged.empty and "source_database" in merged.columns
            else {}
        ),
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge primary and supplemental disease-microbe references into one de-duplicated table.")
    parser.add_argument(
        "--primary-path",
        default=ROOT / "data/reference/disease_microbe_dictionary.csv",
        type=Path,
    )
    parser.add_argument(
        "--supplement-path",
        action="append",
        default=[],
        type=Path,
        help="Supplemental disease-microbe CSV path. Can be repeated.",
    )
    parser.add_argument(
        "--output-path",
        default=ROOT / "data/reference/disease_microbe_dictionary_merged.csv",
        type=Path,
    )
    parser.add_argument(
        "--summary-path",
        default=ROOT / "data/reference/disease_microbe_dictionary_merged.summary.json",
        type=Path,
    )
    args = parser.parse_args()

    summary = merge_disease_microbe_references(
        primary_path=args.primary_path,
        supplement_paths=[path for path in args.supplement_path if path.exists()],
        output_path=args.output_path,
        summary_path=args.summary_path,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
