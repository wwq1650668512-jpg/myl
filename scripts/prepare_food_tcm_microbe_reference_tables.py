from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MICROBE_TABLE = ROOT / "data/processed/step1/step1_microbe_table.csv"
DEFAULT_OUTPUT_DIR = ROOT / "data/processed/food_tcm"

MICROBE_FOOD_COLUMNS = [
    "nt_code",
    "species_label",
    "species_name",
    "genus",
    "phylum",
    "food_name",
    "food_category",
    "food_component",
    "relation_direction",
    "relation_scope",
    "evidence_type",
    "host_context",
    "disease_context",
    "source_name",
    "source_url",
    "source_record_id",
    "pmid_or_doi",
    "curation_status",
    "notes",
]

MICROBE_TCM_COLUMNS = [
    "nt_code",
    "species_label",
    "species_name",
    "genus",
    "phylum",
    "tcm_name",
    "tcm_type",
    "tcm_component",
    "relation_direction",
    "relation_scope",
    "evidence_type",
    "host_context",
    "disease_context",
    "source_name",
    "source_url",
    "source_record_id",
    "pmid_or_doi",
    "curation_status",
    "notes",
]

SOURCE_REGISTRY_COLUMNS = [
    "domain",
    "source_name",
    "source_url",
    "citation_short",
    "coverage_note",
    "curation_status",
    "notes",
]


def _ensure_writable(path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} already exists. Re-run with --overwrite to replace it.")


def _seed_microbe_base(microbe_table: pd.DataFrame) -> pd.DataFrame:
    required_columns = ["nt_code", "species_label", "species_name", "genus", "phylum"]
    missing = [column for column in required_columns if column not in microbe_table.columns]
    if missing:
        raise ValueError(f"Microbe table is missing required columns: {missing}")
    return (
        microbe_table.loc[:, required_columns]
        .drop_duplicates(subset=["nt_code"])
        .sort_values(["species_label", "nt_code"], na_position="last")
        .reset_index(drop=True)
    )


def _seed_food_table(base: pd.DataFrame) -> pd.DataFrame:
    frame = base.copy()
    frame["food_name"] = ""
    frame["food_category"] = ""
    frame["food_component"] = ""
    frame["relation_direction"] = ""
    frame["relation_scope"] = ""
    frame["evidence_type"] = ""
    frame["host_context"] = ""
    frame["disease_context"] = ""
    frame["source_name"] = "FGMDI / gutMDisorder / literature"
    frame["source_url"] = ""
    frame["source_record_id"] = ""
    frame["pmid_or_doi"] = ""
    frame["curation_status"] = "pending_lookup"
    frame["notes"] = ""
    return frame[MICROBE_FOOD_COLUMNS]


def _seed_tcm_table(base: pd.DataFrame) -> pd.DataFrame:
    frame = base.copy()
    frame["tcm_name"] = ""
    frame["tcm_type"] = ""
    frame["tcm_component"] = ""
    frame["relation_direction"] = ""
    frame["relation_scope"] = ""
    frame["evidence_type"] = ""
    frame["host_context"] = ""
    frame["disease_context"] = ""
    frame["source_name"] = "MicrobeTCM / literature"
    frame["source_url"] = ""
    frame["source_record_id"] = ""
    frame["pmid_or_doi"] = ""
    frame["curation_status"] = "pending_lookup"
    frame["notes"] = ""
    return frame[MICROBE_TCM_COLUMNS]


def _seed_source_registry() -> pd.DataFrame:
    rows = [
        {
            "domain": "food",
            "source_name": "FGMDI",
            "source_url": "https://doi.org/10.1016/j.fbio.2024.104091",
            "citation_short": "Food Bioscience, 2024",
            "coverage_note": "1806 curated food-gut microbe associations; 495 gut microbes; 313 foods.",
            "curation_status": "recommended_primary_source",
            "notes": "Best starting point for food -> microbe links; also includes dietary patterns, dietary substances, plant compounds, herbal medicine, probiotics and prebiotics.",
        },
        {
            "domain": "food",
            "source_name": "gutMDisorder",
            "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6943049/",
            "citation_short": "NAR, 2020",
            "coverage_note": "2263 curated human microbe-disorder/intervention associations; 579 gut microbes; 77 intervention measures.",
            "curation_status": "recommended_secondary_source",
            "notes": "Useful for intervention measures including foods and broader host contexts.",
        },
        {
            "domain": "tcm",
            "source_name": "MicrobeTCM",
            "source_url": "https://www.microbetcm.com",
            "citation_short": "Pharmacol Res, 2024",
            "coverage_note": "725 microbes; 1032 herbs; 1468 herb-formulas; 15780 chemical compositions; 77 acupoints.",
            "curation_status": "recommended_primary_source",
            "notes": "Best starting point for herb/herb-formula/acupoint -> microbe links.",
        },
        {
            "domain": "microbe_host_axis",
            "source_name": "gutMGene v2.0",
            "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11701569/",
            "citation_short": "NAR, 2025",
            "coverage_note": "Microbe-metabolite-gene associations classified as causal or correlational.",
            "curation_status": "recommended_context_source",
            "notes": "Useful when a food or herb effect is mediated via microbial metabolites or host genes rather than direct abundance shifts.",
        },
    ]
    return pd.DataFrame(rows, columns=SOURCE_REGISTRY_COLUMNS)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Bootstrap food- and TCM-microbe relation templates for the current 83-microbe panel."
    )
    parser.add_argument(
        "--microbe-table",
        type=Path,
        default=DEFAULT_MICROBE_TABLE,
        help="Step 1 microbe table used to seed the lookup templates.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where food/TCM relation templates will be written.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing files if they already exist.",
    )
    args = parser.parse_args()

    microbe_table = pd.read_csv(args.microbe_table, low_memory=False)
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    food_table_path = output_dir / "microbe_food_relation_table.csv"
    tcm_table_path = output_dir / "microbe_tcm_relation_table.csv"
    source_registry_path = output_dir / "source_registry.csv"
    summary_path = output_dir / "template_summary.json"

    for path in [food_table_path, tcm_table_path, source_registry_path, summary_path]:
        _ensure_writable(path, overwrite=args.overwrite)

    base = _seed_microbe_base(microbe_table)
    food_table = _seed_food_table(base)
    tcm_table = _seed_tcm_table(base)
    source_registry = _seed_source_registry()

    food_table.to_csv(food_table_path, index=False)
    tcm_table.to_csv(tcm_table_path, index=False)
    source_registry.to_csv(source_registry_path, index=False)

    summary = {
        "output_dir": str(output_dir),
        "microbe_table": str(args.microbe_table),
        "food_table_path": str(food_table_path),
        "tcm_table_path": str(tcm_table_path),
        "source_registry_path": str(source_registry_path),
        "n_microbes": int(len(base)),
        "n_unique_species_labels": int(base["species_label"].nunique()),
        "n_source_registry_rows": int(len(source_registry)),
        "food_columns": MICROBE_FOOD_COLUMNS,
        "tcm_columns": MICROBE_TCM_COLUMNS,
    }
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
