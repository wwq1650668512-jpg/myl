from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MICROBE_TABLE = ROOT / "data/processed/step1/step1_microbe_table.csv"
DEFAULT_OUTPUT_DIR = ROOT / "data/processed/amr"

MICROBE_REFERENCE_COLUMNS = [
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
    "reference_genome_id",
    "reference_genome_source",
    "has_beta_lactamase",
    "beta_lactamase_family",
    "has_efflux_amr",
    "has_target_alteration",
    "has_intrinsic_amr_evidence",
    "expected_beta_lactam_resistant",
    "amr_gene_summary",
    "amr_gene_source",
    "curation_status",
    "notes",
]
RULE_COLUMNS = [
    "rule_id",
    "drug_class",
    "drug_name",
    "species_label",
    "genus",
    "expected_phenotype",
    "rule_level",
    "mechanism_hint",
    "rule_strength",
    "source_name",
    "source_url",
    "source_version",
    "notes",
]
PHENOTYPE_PRIOR_COLUMNS = [
    "species_label",
    "genus",
    "drug_class",
    "drug_name",
    "n_tested",
    "resistant_fraction",
    "intermediate_fraction",
    "susceptible_fraction",
    "mic50",
    "mic90",
    "source_name",
    "source_url",
    "source_version",
    "evidence_note",
]


def _seed_microbe_reference(microbe_table: pd.DataFrame) -> pd.DataFrame:
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
        .sort_values(["genus", "species_label", "nt_code"], na_position="last")
        .reset_index(drop=True)
        .copy()
    )
    frame["reference_genome_id"] = ""
    frame["reference_genome_source"] = ""
    frame["has_beta_lactamase"] = ""
    frame["beta_lactamase_family"] = ""
    frame["has_efflux_amr"] = ""
    frame["has_target_alteration"] = ""
    frame["has_intrinsic_amr_evidence"] = ""
    frame["expected_beta_lactam_resistant"] = ""
    frame["amr_gene_summary"] = ""
    frame["amr_gene_source"] = ""
    frame["curation_status"] = "pending"
    frame["notes"] = ""
    return frame[MICROBE_REFERENCE_COLUMNS]


def _ensure_writable(path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} already exists. Re-run with --overwrite to replace it.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Bootstrap AMR reference tables for the current Step 1 / Step 3 microbe panel."
    )
    parser.add_argument(
        "--microbe-table",
        type=Path,
        default=DEFAULT_MICROBE_TABLE,
        help="Microbe panel table used to seed the AMR reference file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where the AMR template CSVs will be written.",
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

    microbe_reference_path = output_dir / "microbe_amr_reference.csv"
    rules_path = output_dir / "drug_resistance_rules.csv"
    phenotype_prior_path = output_dir / "species_drugclass_phenotype_prior.csv"
    summary_path = output_dir / "template_summary.json"

    for path in [microbe_reference_path, rules_path, phenotype_prior_path, summary_path]:
        _ensure_writable(path, overwrite=args.overwrite)

    microbe_reference = _seed_microbe_reference(microbe_table)
    rules = pd.DataFrame(columns=RULE_COLUMNS)
    phenotype_priors = pd.DataFrame(columns=PHENOTYPE_PRIOR_COLUMNS)

    microbe_reference.to_csv(microbe_reference_path, index=False)
    rules.to_csv(rules_path, index=False)
    phenotype_priors.to_csv(phenotype_prior_path, index=False)

    summary = {
        "microbe_table": str(args.microbe_table),
        "output_dir": str(output_dir),
        "microbe_reference_path": str(microbe_reference_path),
        "rules_path": str(rules_path),
        "phenotype_prior_path": str(phenotype_prior_path),
        "n_microbe_rows": int(len(microbe_reference)),
        "n_rule_rows": int(len(rules)),
        "n_phenotype_prior_rows": int(len(phenotype_priors)),
        "microbe_reference_columns": MICROBE_REFERENCE_COLUMNS,
        "rule_columns": RULE_COLUMNS,
        "phenotype_prior_columns": PHENOTYPE_PRIOR_COLUMNS,
    }
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
