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


def _rule_id(stem: str) -> str:
    return "".join(character if character.isalnum() else "_" for character in stem.strip().lower()).strip("_")


def main() -> None:
    parser = argparse.ArgumentParser(description="Seed broader high-risk false-positive AMR rules for antibiotics.")
    parser.add_argument("--reference-path", type=Path, default=DEFAULT_REFERENCE_PATH)
    parser.add_argument("--rules-path", type=Path, default=DEFAULT_RULES_PATH)
    args = parser.parse_args()

    reference = pd.read_csv(args.reference_path, low_memory=False)
    rules = pd.read_csv(args.rules_path, low_memory=False)

    gram_negative_genera = sorted(
        {
            str(genus).strip()
            for genus, gram_stain in zip(reference.get("genus", []), reference.get("gram_stain", []))
            if str(genus).strip() and str(gram_stain).strip().lower() == "negative"
        }
    )

    seeded_rules: list[dict[str, object]] = []
    for genus in gram_negative_genera:
        genus_key = _rule_id(genus)
        seeded_rules.extend(
            [
                {
                    "rule_id": f"amr_{genus_key}_vancomycin_strong",
                    "drug_class": "vancomycin",
                    "drug_name": "",
                    "species_label": "",
                    "genus": genus,
                    "expected_phenotype": "resistant",
                    "rule_level": "genus",
                    "mechanism_hint": "vancomycin has poor gram-negative activity due to the outer membrane barrier",
                    "rule_strength": "strong",
                    "source_name": "Vancomycin - StatPearls",
                    "source_url": "https://www.ncbi.nlm.nih.gov/books/NBK459263/",
                    "source_version": "accessed 2026-04-04",
                    "notes": "Use to suppress high-risk gram-negative false positives for explicit vancomycin-like inputs.",
                },
                {
                    "rule_id": f"amr_{genus_key}_glycopeptide_strong",
                    "drug_class": "glycopeptide",
                    "drug_name": "",
                    "species_label": "",
                    "genus": genus,
                    "expected_phenotype": "resistant",
                    "rule_level": "genus",
                    "mechanism_hint": "glycopeptide-class activity is largely limited to gram-positive bacteria",
                    "rule_strength": "strong",
                    "source_name": "Vancomycin - StatPearls",
                    "source_url": "https://www.ncbi.nlm.nih.gov/books/NBK459263/",
                    "source_version": "accessed 2026-04-04",
                    "notes": "Class-level prior for gram-negative genera currently present in the 83-microbe panel.",
                },
                {
                    "rule_id": f"amr_{genus_key}_daptomycin_strong",
                    "drug_class": "daptomycin",
                    "drug_name": "",
                    "species_label": "",
                    "genus": genus,
                    "expected_phenotype": "resistant",
                    "rule_level": "genus",
                    "mechanism_hint": "daptomycin activity is restricted to gram-positive organisms",
                    "rule_strength": "strong",
                    "source_name": "Daptomycin review (PMC4846043)",
                    "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4846043/",
                    "source_version": "accessed 2026-04-04",
                    "notes": "Use to suppress high-risk gram-negative false positives for explicit daptomycin inputs.",
                },
                {
                    "rule_id": f"amr_{genus_key}_lipopeptide_moderate",
                    "drug_class": "lipopeptide",
                    "drug_name": "",
                    "species_label": "",
                    "genus": genus,
                    "expected_phenotype": "resistant",
                    "rule_level": "genus",
                    "mechanism_hint": "lipopeptide antibiotic activity is expected to be gram-positive biased",
                    "rule_strength": "moderate",
                    "source_name": "Daptomycin review (PMC4846043)",
                    "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4846043/",
                    "source_version": "accessed 2026-04-04",
                    "notes": "Class-level prior kept one notch weaker than explicit daptomycin matching.",
                },
            ]
        )

    seeded_rules.extend(
        [
            {
                "rule_id": "amr_lactobacillus_vancomycin_strong",
                "drug_class": "vancomycin",
                "drug_name": "",
                "species_label": "",
                "genus": "Lactobacillus",
                "expected_phenotype": "resistant",
                "rule_level": "genus",
                "mechanism_hint": "intrinsic vancomycin resistance is well described in lactobacilli",
                "rule_strength": "strong",
                "source_name": "Gram-positive organisms intrinsically resistant to vancomycin",
                "source_url": "https://pubmed.ncbi.nlm.nih.gov/10467540/",
                "source_version": "accessed 2026-04-04",
                "notes": "Targets another common high-risk false-positive direction outside gram-negative genera.",
            },
            {
                "rule_id": "amr_lactobacillus_glycopeptide_moderate",
                "drug_class": "glycopeptide",
                "drug_name": "",
                "species_label": "",
                "genus": "Lactobacillus",
                "expected_phenotype": "resistant",
                "rule_level": "genus",
                "mechanism_hint": "intrinsic glycopeptide resistance prior in lactobacilli",
                "rule_strength": "moderate",
                "source_name": "Gram-positive organisms intrinsically resistant to vancomycin",
                "source_url": "https://pubmed.ncbi.nlm.nih.gov/10467540/",
                "source_version": "accessed 2026-04-04",
                "notes": "Broader glycopeptide prior kept weaker than explicit vancomycin matching.",
            },
        ]
    )

    updated_rules = _append_missing_rules(rules, seeded_rules)
    updated_rules.to_csv(args.rules_path, index=False)

    print(
        {
            "rules_path": str(args.rules_path),
            "n_gram_negative_genera_seeded": int(len(gram_negative_genera)),
            "n_rules_total": int(len(updated_rules)),
        }
    )


if __name__ == "__main__":
    main()
