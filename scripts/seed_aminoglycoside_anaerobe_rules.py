from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RULES_PATH = ROOT / "data/processed/amr/drug_resistance_rules.csv"

STRICT_ANAEROBE_GENERA = [
    "Akkermansia",
    "Bacteroides",
    "Bilophila",
    "Blautia",
    "Clostridium",
    "Collinsella",
    "Coprococcus",
    "Dorea",
    "Eggerthella",
    "Eubacterium",
    "Fusobacterium",
    "Odoribacter",
    "Parabacteroides",
    "Peptoclostridium",
    "Prevotella",
    "Roseburia",
    "Ruminococcus",
    "Veillonella",
]


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
    parser = argparse.ArgumentParser(description="Seed aminoglycoside x strict-anaerobe AMR priors.")
    parser.add_argument("--rules-path", type=Path, default=DEFAULT_RULES_PATH)
    args = parser.parse_args()

    rules = pd.read_csv(args.rules_path, low_memory=False)

    seeded_rules: list[dict[str, object]] = []
    for genus in STRICT_ANAEROBE_GENERA:
        genus_key = _rule_id(genus)
        seeded_rules.extend(
            [
                {
                    "rule_id": f"amr_{genus_key}_gentamicin_strong",
                    "drug_class": "gentamicin",
                    "drug_name": "",
                    "species_label": "",
                    "genus": genus,
                    "expected_phenotype": "resistant",
                    "rule_level": "genus",
                    "mechanism_hint": "aminoglycoside uptake is oxygen-dependent and strict anaerobes are expected to be intrinsically inactive",
                    "rule_strength": "strong",
                    "source_name": "Gentamicin - StatPearls",
                    "source_url": "https://www.ncbi.nlm.nih.gov/books/NBK557550/",
                    "source_version": "accessed 2026-04-05",
                    "notes": "Explicit gentamicin prior for strict anaerobe genera in the current panel.",
                },
                {
                    "rule_id": f"amr_{genus_key}_amikacin_strong",
                    "drug_class": "amikacin",
                    "drug_name": "",
                    "species_label": "",
                    "genus": genus,
                    "expected_phenotype": "resistant",
                    "rule_level": "genus",
                    "mechanism_hint": "aminoglycoside uptake is oxygen-dependent and strict anaerobes are expected to be intrinsically inactive",
                    "rule_strength": "strong",
                    "source_name": "Aminoglycosides - StatPearls",
                    "source_url": "https://www.ncbi.nlm.nih.gov/books/NBK541105/",
                    "source_version": "accessed 2026-04-05",
                    "notes": "Explicit amikacin prior for strict anaerobe genera in the current panel.",
                },
                {
                    "rule_id": f"amr_{genus_key}_tobramycin_strong",
                    "drug_class": "tobramycin",
                    "drug_name": "",
                    "species_label": "",
                    "genus": genus,
                    "expected_phenotype": "resistant",
                    "rule_level": "genus",
                    "mechanism_hint": "aminoglycoside uptake is oxygen-dependent and strict anaerobes are expected to be intrinsically inactive",
                    "rule_strength": "strong",
                    "source_name": "Aminoglycosides - StatPearls",
                    "source_url": "https://www.ncbi.nlm.nih.gov/books/NBK541105/",
                    "source_version": "accessed 2026-04-05",
                    "notes": "Explicit tobramycin prior for strict anaerobe genera in the current panel.",
                },
                {
                    "rule_id": f"amr_{genus_key}_streptomycin_moderate",
                    "drug_class": "streptomycin",
                    "drug_name": "",
                    "species_label": "",
                    "genus": genus,
                    "expected_phenotype": "resistant",
                    "rule_level": "genus",
                    "mechanism_hint": "aminoglycoside uptake is oxygen-dependent and strict anaerobes are expected to be intrinsically inactive",
                    "rule_strength": "moderate",
                    "source_name": "Aminoglycosides - StatPearls",
                    "source_url": "https://www.ncbi.nlm.nih.gov/books/NBK541105/",
                    "source_version": "accessed 2026-04-05",
                    "notes": "Kept one notch weaker than gentamicin / amikacin / tobramycin to stay conservative at the class edge.",
                },
                {
                    "rule_id": f"amr_{genus_key}_aminoglycoside_moderate",
                    "drug_class": "aminoglycoside",
                    "drug_name": "",
                    "species_label": "",
                    "genus": genus,
                    "expected_phenotype": "resistant",
                    "rule_level": "genus",
                    "mechanism_hint": "strict anaerobes are expected to lack the oxygen-dependent transport needed for aminoglycoside entry",
                    "rule_strength": "moderate",
                    "source_name": "Gentamicin - StatPearls",
                    "source_url": "https://www.ncbi.nlm.nih.gov/books/NBK557550/",
                    "source_version": "accessed 2026-04-05",
                    "notes": "Class-level prior for strict anaerobe genera in the current panel.",
                },
            ]
        )

    updated_rules = _append_missing_rules(rules, seeded_rules)
    updated_rules.to_csv(args.rules_path, index=False)

    print(
        {
            "rules_path": str(args.rules_path),
            "n_strict_anaerobe_genera_seeded": int(len(STRICT_ANAEROBE_GENERA)),
            "n_rules_total": int(len(updated_rules)),
        }
    )


if __name__ == "__main__":
    main()
