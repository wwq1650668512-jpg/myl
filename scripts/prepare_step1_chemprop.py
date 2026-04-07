from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step1 import prepare_step1_chemprop_inputs


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare Chemprop-ready Step 1 datasets and descriptors.")
    parser.add_argument(
        "--modeling-table",
        default=ROOT / "data/processed/step1/step1_modeling_table.csv",
        type=Path,
        help="Normalized gold modeling table.",
    )
    parser.add_argument(
        "--silver-table",
        default=ROOT / "data/processed/step1/step1_silver_mdipid.csv",
        type=Path,
        help="Optional silver table appended to classification training only.",
    )
    parser.add_argument(
        "--output-dir",
        default=ROOT / "data/processed/step1/chemprop_scaffold",
        type=Path,
        help="Directory for Chemprop-ready CSV/NPZ files.",
    )
    parser.add_argument(
        "--split-mode",
        choices=["random", "drug", "scaffold", "microbe"],
        default="scaffold",
        help="Split strategy used to create train/val/test assignments.",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random seed.",
    )
    parser.add_argument(
        "--test-size",
        type=float,
        default=0.2,
        help="Fraction of gold data held out as test.",
    )
    parser.add_argument(
        "--val-size",
        type=float,
        default=0.1,
        help="Fraction of remaining train/val gold data used as validation.",
    )
    parser.add_argument(
        "--positive-label",
        default="inhibit",
        help="Positive class used for binary Chemprop classification.",
    )
    args = parser.parse_args()

    silver_table = args.silver_table if args.silver_table.exists() else None
    summary = prepare_step1_chemprop_inputs(
        modeling_table_path=args.modeling_table,
        output_dir=args.output_dir,
        silver_table_path=silver_table,
        split_mode=args.split_mode,
        random_state=args.random_state,
        test_size=args.test_size,
        val_size=args.val_size,
        positive_label=args.positive_label,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
