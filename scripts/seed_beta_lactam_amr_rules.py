from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REFERENCE_PATH = ROOT / "data/processed/amr/microbe_amr_reference.csv"
DEFAULT_RULES_PATH = ROOT / "data/processed/amr/drug_resistance_rules.csv"


def _append_missing_rules(existing: pd.DataFrame, new_rows: list[dict[str, object]]) -> pd.DataFrame:
    existing = existing.copy()
    if "rule_id" not in existing.columns:
        existing["rule_id"] = ""
    existing_ids = set(existing["rule_id"].fillna("").astype(str))
    rows_to_add = [row for row in new_rows if str(row["rule_id"]) not in existing_ids]
    if not rows_to_add:
        return existing
    return pd.concat([existing, pd.DataFrame(rows_to_add)], ignore_index=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Seed the first batch of beta-lactam / penicillin AMR priors.")
    parser.add_argument("--reference-path", type=Path, default=DEFAULT_REFERENCE_PATH)
    parser.add_argument("--rules-path", type=Path, default=DEFAULT_RULES_PATH)
    args = parser.parse_args()

    reference = pd.read_csv(args.reference_path, low_memory=False)
    rules = pd.read_csv(args.rules_path, low_memory=False)

    for column in [
        "has_beta_lactamase",
        "beta_lactamase_family",
        "has_intrinsic_amr_evidence",
        "expected_beta_lactam_resistant",
        "amr_gene_summary",
        "amr_gene_source",
        "curation_status",
    ]:
        if column not in reference.columns:
            reference[column] = ""
        reference[column] = reference[column].astype("object")

    bacteroides_mask = reference["genus"].fillna("").astype(str).str.strip().eq("Bacteroides")
    reference.loc[bacteroides_mask, "has_beta_lactamase"] = "likely"
    reference.loc[bacteroides_mask, "beta_lactamase_family"] = "genus-level prior"
    reference.loc[bacteroides_mask, "has_intrinsic_amr_evidence"] = "reported"
    reference.loc[bacteroides_mask, "expected_beta_lactam_resistant"] = "penicillin-like beta-lactams"
    reference.loc[bacteroides_mask, "amr_gene_summary"] = (
        "Bacteroides spp. often carry beta-lactamase-associated resistance; validate per strain for class-specific exceptions."
    )
    reference.loc[bacteroides_mask, "amr_gene_source"] = "NCBI Bookshelf / PMC beta-lactamase literature"
    reference.loc[bacteroides_mask, "curation_status"] = "seeded_beta_lactam_prior"

    seeded_rules = [
        {
            "rule_id": "amr_bacteroides_beta_lactam_supporting",
            "drug_class": "beta_lactam",
            "drug_name": "",
            "species_label": "",
            "genus": "Bacteroides",
            "expected_phenotype": "resistant",
            "rule_level": "genus",
            "mechanism_hint": "beta-lactamase prior; use as conservative downweight only",
            "rule_strength": "supporting",
            "source_name": "NCBI Bookshelf Bacteroides fragilis overview",
            "source_url": "https://www.ncbi.nlm.nih.gov/sites/books/NBK553032/",
            "source_version": "accessed 2026-04-03",
            "notes": "Designed for service-layer correction when only a broad beta-lactam signal is available.",
        },
        {
            "rule_id": "amr_bacteroides_penicillin_strong",
            "drug_class": "penicillin",
            "drug_name": "",
            "species_label": "",
            "genus": "Bacteroides",
            "expected_phenotype": "resistant",
            "rule_level": "genus",
            "mechanism_hint": "beta-lactamase-mediated penicillin resistance prior",
            "rule_strength": "strong",
            "source_name": "Bacteroides beta-lactamase literature",
            "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC352644/",
            "source_version": "accessed 2026-04-03",
            "notes": "Use this stronger prior when the drug name/class indicates a penicillin-like agent.",
        },
        {
            "rule_id": "amr_bvulgatus_penicillin_strong",
            "drug_class": "penicillin",
            "drug_name": "",
            "species_label": "Bacteroides vulgatus",
            "genus": "Bacteroides",
            "expected_phenotype": "resistant",
            "rule_level": "species",
            "mechanism_hint": "reported beta-lactam resistance in B. vulgatus",
            "rule_strength": "strong",
            "source_name": "PMC Bacteroides vulgatus susceptibility paper",
            "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC187887/",
            "source_version": "accessed 2026-04-03",
            "notes": "Species-level prior used to outrank the broader genus rule.",
        },
        {
            "rule_id": "amr_buniformis_penicillin_moderate",
            "drug_class": "penicillin",
            "drug_name": "",
            "species_label": "Bacteroides uniformis",
            "genus": "Bacteroides",
            "expected_phenotype": "resistant",
            "rule_level": "species",
            "mechanism_hint": "genus-consistent beta-lactamase prior",
            "rule_strength": "moderate",
            "source_name": "NCBI Bookshelf Bacteroides fragilis overview",
            "source_url": "https://www.ncbi.nlm.nih.gov/sites/books/NBK553032/",
            "source_version": "accessed 2026-04-03",
            "notes": "Conservative species prior pending more direct phenotype data.",
        },
        {
            "rule_id": "amr_btheta_penicillin_moderate",
            "drug_class": "penicillin",
            "drug_name": "",
            "species_label": "Bacteroides thetaiotaomicron",
            "genus": "Bacteroides",
            "expected_phenotype": "resistant",
            "rule_level": "species",
            "mechanism_hint": "genus-consistent beta-lactamase prior",
            "rule_strength": "moderate",
            "source_name": "NCBI Bookshelf Bacteroides fragilis overview",
            "source_url": "https://www.ncbi.nlm.nih.gov/sites/books/NBK553032/",
            "source_version": "accessed 2026-04-03",
            "notes": "Conservative species prior pending more direct phenotype data.",
        },
        {
            "rule_id": "amr_bfragilis_penicillin_strong",
            "drug_class": "penicillin",
            "drug_name": "",
            "species_label": "Bacteroides fragilis",
            "genus": "Bacteroides",
            "expected_phenotype": "resistant",
            "rule_level": "species",
            "mechanism_hint": "beta-lactamase-mediated resistance prior",
            "rule_strength": "strong",
            "source_name": "NCBI Bookshelf Bacteroides fragilis overview",
            "source_url": "https://www.ncbi.nlm.nih.gov/sites/books/NBK553032/",
            "source_version": "accessed 2026-04-03",
            "notes": "Intended to catch obvious penicillin false positives in the current panel.",
        },
    ]

    updated_rules = _append_missing_rules(rules, seeded_rules)
    reference.to_csv(args.reference_path, index=False)
    updated_rules.to_csv(args.rules_path, index=False)

    print(
        {
            "reference_path": str(args.reference_path),
            "rules_path": str(args.rules_path),
            "n_bacteroides_rows_updated": int(bacteroides_mask.sum()),
            "n_rules_total": int(len(updated_rules)),
        }
    )


if __name__ == "__main__":
    main()
