from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step1 import summarize_step1_chemprop_run


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize Chemprop Step 1 test outputs into a compact metrics JSON.")
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Chemprop training output directory containing test.csv and model_0/test_predictions.csv.",
    )
    parser.add_argument(
        "--task-type",
        choices=["classification", "regression"],
        required=True,
        help="Chemprop task type used for the run.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.5,
        help="Probability threshold for binary classification summaries.",
    )
    args = parser.parse_args()

    summary = summarize_step1_chemprop_run(
        output_dir=args.output_dir,
        task_type=args.task_type,
        threshold=args.threshold,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
