from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step1 import build_masi_silver_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Normalize MASI records into a Step 1 silver table.")
    parser.add_argument(
        "--raw-dir",
        default=ROOT / "data/raw/step1/masi",
        type=Path,
        help="Directory containing MASI raw files.",
    )
    parser.add_argument(
        "--output-dir",
        default=ROOT / "data/processed/step1",
        type=Path,
        help="Directory for the processed silver table.",
    )
    parser.add_argument(
        "--drug-like-only",
        action="store_true",
        help="Keep only rows with drug-like scope signals inferred from MASI metadata.",
    )
    args = parser.parse_args()

    try:
        summary = build_masi_silver_table(
            raw_dir=args.raw_dir,
            output_dir=args.output_dir,
            drug_like_only=args.drug_like_only,
        )
    except FileNotFoundError as exc:
        raise SystemExit(str(exc)) from exc
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
