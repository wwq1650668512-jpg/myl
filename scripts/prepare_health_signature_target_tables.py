from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "data/processed/health_signature"

TARGET_GENOME_COLUMNS = [
    "target_genome_id",
    "signature_name",
    "guild_membership",
    "guild_role",
    "species_name",
    "canonical_species_name",
    "strain_name",
    "ncbi_taxid",
    "ncbi_assembly_accession",
    "refseq_genome_accession",
    "gtdb_genome_id",
    "genome_label",
    "source_name",
    "source_url",
    "source_version",
    "evidence_level",
    "ready_for_matching",
    "notes",
]

SOURCE_REGISTRY_COLUMNS = [
    "signature_name",
    "source_name",
    "source_url",
    "source_version",
    "citation_short",
    "curation_status",
    "notes",
]


def _ensure_writable(path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} already exists. Re-run with --overwrite to replace it.")


def _seed_target_genomes() -> pd.DataFrame:
    return pd.DataFrame(columns=TARGET_GENOME_COLUMNS)


def _seed_source_registry() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "signature_name": "core_microbiome_signature_health",
                "source_name": "A core microbiome signature as an indicator of health",
                "source_url": "https://pubmed.ncbi.nlm.nih.gov/39378879/",
                "source_version": "accessed 2026-04-06",
                "citation_short": "Cell, 2024",
                "curation_status": "pending_member_extraction",
                "notes": "Fill health_signature_target_genomes.csv with the paper's genome/guild members before matching onto the 83-microbe panel.",
            }
        ],
        columns=SOURCE_REGISTRY_COLUMNS,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Bootstrap target genome tables for the core microbiome health-signature integration workflow."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where the target-genome template CSVs will be written.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing template files if they already exist.",
    )
    args = parser.parse_args()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    target_genome_path = output_dir / "health_signature_target_genomes.csv"
    source_registry_path = output_dir / "health_signature_source_registry.csv"
    summary_path = output_dir / "target_template_summary.json"

    for path in [target_genome_path, source_registry_path, summary_path]:
        _ensure_writable(path, overwrite=args.overwrite)

    target_genomes = _seed_target_genomes()
    source_registry = _seed_source_registry()

    target_genomes.to_csv(target_genome_path, index=False)
    source_registry.to_csv(source_registry_path, index=False)

    summary = {
        "output_dir": str(output_dir),
        "target_genome_path": str(target_genome_path),
        "source_registry_path": str(source_registry_path),
        "n_target_rows": int(len(target_genomes)),
        "n_source_rows": int(len(source_registry)),
        "target_genome_columns": TARGET_GENOME_COLUMNS,
        "source_registry_columns": SOURCE_REGISTRY_COLUMNS,
    }
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
