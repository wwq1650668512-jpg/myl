from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step2 import normalize_step2_label_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Normalize a Step 2 metabolism label table into the project schema.")
    parser.add_argument("--input-path", required=True, type=Path, help="Source CSV/TSV/XLSX table.")
    parser.add_argument("--output-path", required=True, type=Path, help="Normalized CSV output path.")
    parser.add_argument("--source-dataset", required=True, help="Stable source identifier, e.g. zimmermann_2019.")
    parser.add_argument("--label-tier", default="gold", help="Label tier, e.g. gold/silver.")
    parser.add_argument("--source-scope", default="isolate", help="Scope such as isolate/community/mechanistic.")
    parser.add_argument("--sheet-name", default=None, help="Optional Excel sheet name or index.")
    args = parser.parse_args()

    sheet_name = args.sheet_name
    if sheet_name is not None and str(sheet_name).isdigit():
        sheet_name = int(sheet_name)

    summary = normalize_step2_label_table(
        input_path=args.input_path,
        output_path=args.output_path,
        source_dataset=args.source_dataset,
        label_tier=args.label_tier,
        source_scope=args.source_scope,
        sheet_name=sheet_name,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
