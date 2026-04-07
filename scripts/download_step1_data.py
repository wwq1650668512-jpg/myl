from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step1 import download_step1_data


def main() -> None:
    parser = argparse.ArgumentParser(description="Download Step 1 raw data.")
    parser.add_argument(
        "--raw-dir",
        default=ROOT / "data/raw/step1/maier_2018",
        type=Path,
        help="Destination directory for raw files.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Redownload files even if they already exist locally.",
    )
    args = parser.parse_args()

    manifest = download_step1_data(raw_dir=args.raw_dir, overwrite=args.overwrite)
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
