from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step1 import download_mdipid_data
from gut_drug_microbiome.step1 import download_masi_data


def main() -> None:
    parser = argparse.ArgumentParser(description="Download Step 1 weak-supervision files.")
    parser.add_argument(
        "--dataset",
        choices=["mdipid", "masi", "all"],
        default="all",
        help="Which weak-supervision dataset to download.",
    )
    parser.add_argument(
        "--raw-root",
        default=ROOT / "data/raw/step1",
        type=Path,
        help="Root directory containing weak-supervision raw datasets.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Redownload files even if they already exist.",
    )
    args = parser.parse_args()

    summary: dict[str, object] = {}
    if args.dataset in {"mdipid", "all"}:
        summary["mdipid"] = download_mdipid_data(args.raw_root / "mdipid", overwrite=args.overwrite)
    if args.dataset in {"masi", "all"}:
        summary["masi"] = download_masi_data(args.raw_root / "masi", overwrite=args.overwrite, allow_partial=True)

    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
