from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step1 import build_step1_tables


def main() -> None:
    parser = argparse.ArgumentParser(description="Normalize Step 1 raw data into modeling tables.")
    parser.add_argument(
        "--raw-dir",
        default=ROOT / "data/raw/step1/maier_2018",
        type=Path,
        help="Directory containing raw Step 1 files.",
    )
    parser.add_argument(
        "--processed-dir",
        default=ROOT / "data/processed/step1",
        type=Path,
        help="Directory for standardized outputs.",
    )
    parser.add_argument(
        "--labels-config",
        default=ROOT / "configs/labeling_rules.yaml",
        type=Path,
        help="Path to the labeling configuration file.",
    )
    parser.add_argument(
        "--include-non-human-use",
        action="store_true",
        help="Keep all compounds instead of filtering to human-use drugs.",
    )
    parser.add_argument(
        "--include-non-primary-panel",
        action="store_true",
        help="Keep non-primary strains if present in future raw tables.",
    )
    args = parser.parse_args()

    summary = build_step1_tables(
        raw_dir=args.raw_dir,
        processed_dir=args.processed_dir,
        labels_config_path=args.labels_config,
        human_use_only=not args.include_non_human_use,
        primary_panel_only=not args.include_non_primary_panel,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
