from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step1 import train_step1_chemprop


def main() -> None:
    parser = argparse.ArgumentParser(description="Train a Chemprop Step 1 model from prepared scaffold files.")
    parser.add_argument(
        "--dataset-csv",
        required=True,
        type=Path,
        help="Prepared Chemprop CSV with smiles, target, and split columns.",
    )
    parser.add_argument(
        "--descriptors-path",
        required=True,
        type=Path,
        help="Prepared descriptor matrix .npz aligned with the dataset CSV.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Chemprop output directory.",
    )
    parser.add_argument(
        "--task-type",
        choices=["classification", "regression"],
        required=True,
        help="Chemprop task type.",
    )
    parser.add_argument(
        "--target-column",
        default="target",
        help="Target column name inside the prepared CSV.",
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=30,
        help="Number of Chemprop training epochs.",
    )
    parser.add_argument(
        "--extra-arg",
        action="append",
        default=[],
        help="Extra raw CLI argument forwarded to Chemprop. Repeatable.",
    )
    parser.add_argument(
        "--print-command-only",
        action="store_true",
        help="Print the resolved Chemprop CLI command without running it.",
    )
    args = parser.parse_args()

    summary = train_step1_chemprop(
        dataset_csv=args.dataset_csv,
        descriptors_path=args.descriptors_path,
        output_dir=args.output_dir,
        task_type=args.task_type,
        target_column=args.target_column,
        epochs=args.epochs,
        extra_args=args.extra_arg,
        dry_run=args.print_command_only,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
