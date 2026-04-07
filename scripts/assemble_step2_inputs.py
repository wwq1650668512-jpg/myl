from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step2 import build_step2_input_tables


def main() -> None:
    parser = argparse.ArgumentParser(description="Build Step 2 candidate and modeling tables from Step 1 hybrid outputs.")
    parser.add_argument(
        "--step1-predictions",
        default=ROOT / "predictions/step1/hybrid_scaffold_v1/predictions.csv",
        type=Path,
        help="Step 1 hybrid predictions table.",
    )
    parser.add_argument(
        "--output-dir",
        default=ROOT / "data/processed/step2",
        type=Path,
        help="Directory for Step 2 candidate/modeling tables.",
    )
    parser.add_argument(
        "--step2-label-table",
        action="append",
        default=[],
        type=Path,
        help="Optional normalized Step 2 label table. Repeatable.",
    )
    args = parser.parse_args()

    label_tables = [path for path in args.step2_label_table if path.exists()]
    summary = build_step2_input_tables(
        step1_predictions_path=args.step1_predictions,
        output_dir=args.output_dir,
        step2_label_table_paths=label_tables or None,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
