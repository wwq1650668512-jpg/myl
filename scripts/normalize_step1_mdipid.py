from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step1 import build_mdipid_silver_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Normalize MDIPID DEIM records into a Step 1 silver table.")
    parser.add_argument(
        "--raw-dir",
        default=ROOT / "data/raw/step1/mdipid",
        type=Path,
        help="Directory containing MDIPID raw files.",
    )
    parser.add_argument(
        "--output-dir",
        default=ROOT / "data/processed/step1",
        type=Path,
        help="Directory for the processed silver table.",
    )
    parser.add_argument(
        "--include-non-gut",
        action="store_true",
        help="Keep non-gut records instead of filtering to gut microbiota only.",
    )
    args = parser.parse_args()

    summary = build_mdipid_silver_table(
        raw_dir=args.raw_dir,
        output_dir=args.output_dir,
        gut_only=not args.include_non_gut,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
