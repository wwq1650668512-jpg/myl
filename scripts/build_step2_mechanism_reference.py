from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step2 import build_step2_mechanism_reference


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a semi-mechanistic Step 2 reference from the labeled modeling table.")
    parser.add_argument(
        "--modeling-table",
        default=ROOT / "data/processed/step2/zimmermann_2019/zimmermann_2019_modeling_table.csv",
        type=Path,
        help="Step 2 modeling table with reaction class / product / gene evidence columns.",
    )
    parser.add_argument(
        "--output-path",
        default=ROOT / "models/step2/zimmermann_scaffold_split/mechanism_reference.joblib",
        type=Path,
        help="Joblib output path.",
    )
    args = parser.parse_args()

    reference = build_step2_mechanism_reference(args.modeling_table, output_path=args.output_path)
    print(
        json.dumps(
            {
                "output_path": str(args.output_path),
                "n_metabolized_rows": int(reference["n_metabolized_rows"]),
                "n_group_specs": int(len(reference["group_specs"])),
            },
            indent=2,
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
