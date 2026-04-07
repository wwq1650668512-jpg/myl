from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step1 import train_step1_baseline


def main() -> None:
    def _parse_source_weight(value: str) -> tuple[str, float]:
        if "=" not in value:
            raise argparse.ArgumentTypeError("Source weights must use the form source_dataset=weight")
        key, raw_weight = value.split("=", 1)
        key = key.strip()
        if not key:
            raise argparse.ArgumentTypeError("Source weight key cannot be empty")
        try:
            weight = float(raw_weight)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(f"Invalid weight: {raw_weight}") from exc
        return key, weight

    parser = argparse.ArgumentParser(description="Train the Step 1 baseline models.")
    parser.add_argument(
        "--modeling-table",
        default=ROOT / "data/processed/step1/step1_modeling_table.csv",
        type=Path,
        help="Normalized modeling table.",
    )
    parser.add_argument(
        "--output-dir",
        default=ROOT / "models/step1/baseline_drug_split",
        type=Path,
        help="Directory for trained models and metrics.",
    )
    parser.add_argument(
        "--silver-table",
        default=None,
        type=Path,
        help="Optional weak-supervision table to append to classifier training only.",
    )
    parser.add_argument(
        "--split-mode",
        choices=["random", "drug", "scaffold", "microbe"],
        default="drug",
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
        "--source-weight",
        action="append",
        default=[],
        type=_parse_source_weight,
        help="Optional per-source classifier sample weight in the form source_dataset=weight. Repeatable.",
    )
    parser.add_argument(
        "--default-gold-weight",
        type=float,
        default=1.0,
        help="Fallback sample weight applied to gold rows not explicitly covered by --source-weight.",
    )
    parser.add_argument(
        "--default-silver-weight",
        type=float,
        default=1.0,
        help="Fallback sample weight applied to silver rows not explicitly covered by --source-weight.",
    )
    parser.add_argument(
        "--verbose",
        type=int,
        default=1,
        help="Verbosity level for coarse training progress logs.",
    )
    args = parser.parse_args()
    source_weight_map = {key: value for key, value in args.source_weight}

    summary = train_step1_baseline(
        modeling_table_path=args.modeling_table,
        output_dir=args.output_dir,
        silver_table_path=args.silver_table,
        split_mode=args.split_mode,
        random_state=args.random_state,
        test_size=args.test_size,
        n_estimators=args.n_estimators,
        source_weight_map=source_weight_map,
        default_gold_weight=args.default_gold_weight,
        default_silver_weight=args.default_silver_weight,
        verbose=args.verbose,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
