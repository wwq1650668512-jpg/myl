from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step1 import build_cross_feeding_reference_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a curated cross-feeding reference table from promote literature seeds.")
    parser.add_argument(
        "--input",
        default=ROOT / "data/reference/promote_literature_seed_table.csv",
        type=Path,
        help="Normalized promote literature seed table.",
    )
    parser.add_argument(
        "--output",
        default=ROOT / "data/reference/cross_feeding_edges.csv",
        type=Path,
        help="Cross-feeding reference CSV output path.",
    )
    parser.add_argument(
        "--extra-input",
        default=ROOT / "data/reference/cross_feeding_seed_template.csv",
        type=Path,
        help="Optional manual seed CSV with curated producer-consumer edges.",
    )
    args = parser.parse_args()
    summary = build_cross_feeding_reference_table(args.input, args.output, args.extra_input)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
