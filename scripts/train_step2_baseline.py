from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step2 import train_step2_baseline


def main() -> None:
    parser = argparse.ArgumentParser(description="Train the Step 2 baseline models.")
    parser.add_argument(
        "--modeling-table",
        default=ROOT / "data/processed/step2/zimmermann_2019/zimmermann_2019_modeling_table.csv",
        type=Path,
        help="Normalized Zimmermann 2019 Step 2 modeling table.",
    )
    parser.add_argument(
        "--output-dir",
        default=ROOT / "models/step2/zimmermann_scaffold_split",
        type=Path,
        help="Directory for trained models and metrics.",
    )
    parser.add_argument(
        "--split-mode",
        choices=["random", "drug", "scaffold", "microbe"],
        default="scaffold",
        help="Data split strategy for evaluation.",
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
        help="Fraction of held-out samples.",
    )
    parser.add_argument(
        "--n-estimators",
        type=int,
        default=300,
        help="Number of trees for both classifier and regressor.",
    )
    parser.add_argument(
        "--probability-threshold",
        type=float,
        default=0.3,
        help="Probability threshold for the predicted metabolism label.",
    )
    parser.add_argument(
        "--verbose",
        type=int,
        default=1,
        help="Verbosity level for coarse training progress logs.",
    )
    args = parser.parse_args()

    summary = train_step2_baseline(
        modeling_table_path=args.modeling_table,
        output_dir=args.output_dir,
        split_mode=args.split_mode,
        random_state=args.random_state,
        test_size=args.test_size,
        n_estimators=args.n_estimators,
        probability_threshold=args.probability_threshold,
        verbose=args.verbose,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
