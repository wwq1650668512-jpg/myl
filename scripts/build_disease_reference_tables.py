from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.disease_knowledge import build_disease_reference_tables


def main() -> None:
    parser = argparse.ArgumentParser(description="Normalize curated disease-microbe and disease-drug Excel sheets.")
    parser.add_argument(
        "--disease-microbe-xlsx",
        default=ROOT / "data/肠道疾病与微生物的相互影响(1).xlsx",
        type=Path,
        help="Curated disease-microbe workbook path.",
    )
    parser.add_argument(
        "--disease-drug-xlsx",
        default=ROOT / "data/疾病与上市药物(1).xlsx",
        type=Path,
        help="Curated disease-drug workbook path.",
    )
    parser.add_argument(
        "--disease-microbe-output",
        default=ROOT / "data/reference/disease_microbe_dictionary.csv",
        type=Path,
        help="Normalized disease-microbe CSV output.",
    )
    parser.add_argument(
        "--disease-drug-output",
        default=ROOT / "data/reference/disease_marketed_drug_catalog.csv",
        type=Path,
        help="Normalized disease-drug CSV output.",
    )
    parser.add_argument(
        "--summary",
        default=ROOT / "data/reference/disease_reference_summary.json",
        type=Path,
        help="JSON summary output path.",
    )
    args = parser.parse_args()

    summary = build_disease_reference_tables(
        disease_microbe_xlsx_path=args.disease_microbe_xlsx,
        disease_drug_xlsx_path=args.disease_drug_xlsx,
        disease_microbe_output_path=args.disease_microbe_output,
        disease_drug_output_path=args.disease_drug_output,
        summary_path=args.summary,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
