from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MICROBE_TABLE = ROOT / "data/processed/step1/step1_microbe_table.csv"
DEFAULT_OUTPUT_DIR = ROOT / "data/processed/health_signature"

REFERENCE_COLUMNS = [
    "nt_code",
    "microbe_label",
    "species_label",
    "species_name",
    "canonical_species_name",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "gram_stain",
    "strain_hint",
    "culture_collection_hint",
    "ncbi_taxid",
    "ncbi_assembly_accession",
    "refseq_genome_accession",
    "gtdb_genome_id",
    "reference_genome_label",
    "reference_genome_source",
    "reference_fasta_path",
    "mapping_level",
    "mapping_confidence",
    "sequence_set_status",
    "curation_status",
    "source_name",
    "source_url",
    "notes",
]

TCG_PROXY_COLUMNS = [
    "nt_code",
    "microbe_label",
    "species_label",
    "canonical_species_name",
    "reference_genome_label",
    "tcg_membership",
    "tcg_mapping_level",
    "tcg_confidence",
    "tcg_target_genome_id",
    "tcg_target_label",
    "tcg_support_evidence",
    "tcg_source_name",
    "tcg_source_url",
    "ready_for_step3",
    "notes",
]

SEQUENCE_SET_COLUMNS = [
    "sequence_set_id",
    "sequence_role",
    "nt_code",
    "species_label",
    "sequence_accession",
    "sequence_source",
    "sequence_type",
    "file_path",
    "assembly_accession",
    "taxid",
    "is_representative",
    "quality_note",
    "notes",
]


def _ensure_writable(path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} already exists. Re-run with --overwrite to replace it.")


def _clean_text(value: object) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value).strip())


def _infer_genus(row: pd.Series) -> str:
    genus = _clean_text(row.get("genus"))
    if genus:
        return genus
    for candidate in [row.get("species_label"), row.get("microbe_label"), row.get("species_name")]:
        text = _clean_text(candidate)
        if text:
            return text.split()[0]
    return ""


def _infer_canonical_species_name(row: pd.Series) -> str:
    for candidate in [row.get("species_label"), row.get("microbe_label")]:
        text = _clean_text(candidate)
        if text:
            return re.sub(r"\s+\([^)]*\)$", "", text)
    species_name = _clean_text(row.get("species_name"))
    if species_name:
        return re.sub(r"\s+\([^)]*\).*", "", species_name)
    return ""


def _extract_culture_collection_hint(text: str) -> str:
    if not text:
        return ""
    matches = re.findall(r"\b(?:DSMZ?|ATCC|JCM|NCTC|CCUG|CIP|VPI)\s*(?:No\.:|No:|No\.|#)?\s*([A-Za-z0-9_\-\/]+)", text)
    if not matches:
        return ""
    return ";".join(dict.fromkeys(match.strip() for match in matches if match.strip()))


def _extract_strain_hint(text: str) -> str:
    if not text:
        return ""
    normalized = _clean_text(text)
    normalized = re.sub(r"\b(?:DSMZ?|ATCC|JCM|NCTC|CCUG|CIP|VPI)\b.*", "", normalized).strip(" ;,")
    return normalized


def _seed_reference_table(microbe_table: pd.DataFrame) -> pd.DataFrame:
    base_columns = [
        "nt_code",
        "microbe_label",
        "species_label",
        "species_name",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "gram_stain",
    ]
    missing = [column for column in base_columns if column not in microbe_table.columns]
    if missing:
        raise ValueError(f"Microbe table is missing required columns: {missing}")

    frame = (
        microbe_table[base_columns]
        .drop_duplicates(subset=["nt_code"])
        .sort_values(["species_label", "nt_code"], na_position="last")
        .reset_index(drop=True)
        .copy()
    )
    frame["canonical_species_name"] = frame.apply(_infer_canonical_species_name, axis=1)
    frame["genus"] = frame.apply(_infer_genus, axis=1)
    frame["strain_hint"] = frame["species_name"].map(_extract_strain_hint)
    frame["culture_collection_hint"] = frame["species_name"].map(_extract_culture_collection_hint)
    frame["ncbi_taxid"] = ""
    frame["ncbi_assembly_accession"] = ""
    frame["refseq_genome_accession"] = ""
    frame["gtdb_genome_id"] = ""
    frame["reference_genome_label"] = ""
    frame["reference_genome_source"] = ""
    frame["reference_fasta_path"] = ""
    frame["mapping_level"] = "unmapped"
    frame["mapping_confidence"] = ""
    frame["sequence_set_status"] = "pending"
    frame["curation_status"] = "pending"
    frame["source_name"] = ""
    frame["source_url"] = ""
    frame["notes"] = ""
    return frame[REFERENCE_COLUMNS]


def _seed_tcg_proxy_table(reference_table: pd.DataFrame) -> pd.DataFrame:
    frame = reference_table.loc[:, ["nt_code", "microbe_label", "species_label", "canonical_species_name"]].copy()
    frame["reference_genome_label"] = ""
    frame["tcg_membership"] = "unmapped"
    frame["tcg_mapping_level"] = "unmapped"
    frame["tcg_confidence"] = ""
    frame["tcg_target_genome_id"] = ""
    frame["tcg_target_label"] = ""
    frame["tcg_support_evidence"] = ""
    frame["tcg_source_name"] = "A core microbiome signature as an indicator of health"
    frame["tcg_source_url"] = "https://pubmed.ncbi.nlm.nih.gov/39378879/"
    frame["ready_for_step3"] = "no"
    frame["notes"] = ""
    return frame[TCG_PROXY_COLUMNS]


def _seed_sequence_set_table(reference_table: pd.DataFrame) -> pd.DataFrame:
    frame = reference_table.loc[:, ["nt_code", "species_label"]].copy()
    frame["sequence_set_id"] = frame["nt_code"].map(lambda value: f"{value}_reference")
    frame["sequence_role"] = "representative_reference"
    frame["sequence_accession"] = ""
    frame["sequence_source"] = ""
    frame["sequence_type"] = "genome"
    frame["file_path"] = ""
    frame["assembly_accession"] = ""
    frame["taxid"] = ""
    frame["is_representative"] = "yes"
    frame["quality_note"] = ""
    frame["notes"] = ""
    return frame[SEQUENCE_SET_COLUMNS]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Bootstrap reference-genome and core-health-signature mapping tables for the 83-microbe panel."
    )
    parser.add_argument(
        "--microbe-table",
        type=Path,
        default=DEFAULT_MICROBE_TABLE,
        help="Step 1 microbe table used to seed the mapping templates.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where the health-signature template CSVs will be written.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing template files if they already exist.",
    )
    args = parser.parse_args()

    microbe_table = pd.read_csv(args.microbe_table, low_memory=False)
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    reference_path = output_dir / "microbe_reference_genome_mapping.csv"
    tcg_proxy_path = output_dir / "microbe_tcg_proxy_mapping.csv"
    sequence_set_path = output_dir / "microbe_reference_sequence_sets.csv"
    summary_path = output_dir / "template_summary.json"

    for path in [reference_path, tcg_proxy_path, sequence_set_path, summary_path]:
        _ensure_writable(path, overwrite=args.overwrite)

    reference_table = _seed_reference_table(microbe_table)
    tcg_proxy_table = _seed_tcg_proxy_table(reference_table)
    sequence_set_table = _seed_sequence_set_table(reference_table)

    reference_table.to_csv(reference_path, index=False)
    tcg_proxy_table.to_csv(tcg_proxy_path, index=False)
    sequence_set_table.to_csv(sequence_set_path, index=False)

    summary = {
        "microbe_table": str(args.microbe_table),
        "output_dir": str(output_dir),
        "reference_path": str(reference_path),
        "tcg_proxy_path": str(tcg_proxy_path),
        "sequence_set_path": str(sequence_set_path),
        "n_microbe_rows": int(len(reference_table)),
        "n_unique_species_labels": int(reference_table["species_label"].nunique()),
        "n_unique_canonical_species_names": int(reference_table["canonical_species_name"].nunique()),
        "reference_columns": REFERENCE_COLUMNS,
        "tcg_proxy_columns": TCG_PROXY_COLUMNS,
        "sequence_set_columns": SEQUENCE_SET_COLUMNS,
    }
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
