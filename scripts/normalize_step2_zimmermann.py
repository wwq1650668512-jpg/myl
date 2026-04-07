from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step2.zimmermann_2019 import normalize_zimmermann_2019


def main() -> None:
    parser = argparse.ArgumentParser(description="Normalize Zimmermann 2019 Step 2 supplementary workbook.")
    parser.add_argument(
        "--input-path",
        default=ROOT / "data/raw/step2/zimmermann_2019/NIHMS1530152-supplement-Supplementary_Tables_1-21.xlsx",
        type=Path,
        help="Zimmermann 2019 supplementary workbook.",
    )
    parser.add_argument(
        "--output-dir",
        default=ROOT / "data/processed/step2/zimmermann_2019",
        type=Path,
        help="Directory for normalized Zimmermann 2019 outputs.",
    )
    args = parser.parse_args()

    summary = normalize_zimmermann_2019(
        input_path=args.input_path,
        output_dir=args.output_dir,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
