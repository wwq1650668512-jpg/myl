from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step1.weak_supervision import build_promote_literature_silver_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Normalize manually curated promote literature records into a Step 1 silver table.")
    parser.add_argument(
        "--input",
        default=ROOT / "data/reference/promote_literature_seed_table.csv",
        type=Path,
        help="Normalized literature seed table input path.",
    )
    parser.add_argument(
        "--output-dir",
        default=ROOT / "data/processed/step1",
        type=Path,
        help="Directory where the promote literature silver table is written.",
    )
    args = parser.parse_args()

    summary = build_promote_literature_silver_table(
        input_path=args.input,
        output_dir=args.output_dir,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
