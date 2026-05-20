from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def merge_silver_tables(inputs: list[Path], output: Path) -> dict[str, object]:
    frames = [pd.read_csv(path, low_memory=False) for path in inputs]
    union_columns: list[str] = []
    for frame in frames:
        for column in frame.columns:
            if column not in union_columns:
                union_columns.append(column)

    aligned = [frame.reindex(columns=union_columns) for frame in frames]
    merged = pd.concat(aligned, ignore_index=True, sort=False)
    if "pair_id" in merged.columns:
        merged = merged.drop_duplicates(subset=["pair_id"], keep="first").reset_index(drop=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output, index=False)

    summary = {
        "inputs": [str(path) for path in inputs],
        "output": str(output),
        "n_rows": int(len(merged)),
        "label_counts": {
            str(key): int(value)
            for key, value in merged["effect_label"].fillna("missing").value_counts().to_dict().items()
        }
        if "effect_label" in merged.columns
        else {},
        "source_counts": {
            str(key): int(value)
            for key, value in merged["source_dataset"].fillna("missing").value_counts().to_dict().items()
        }
        if "source_dataset" in merged.columns
        else {},
    }
    output.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge multiple Step 1 silver tables into one CSV.")
    parser.add_argument(
        "--input",
        action="append",
        required=True,
        type=Path,
        help="Input silver table path. Repeatable.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output merged CSV path.",
    )
    args = parser.parse_args()
    summary = merge_silver_tables(args.input, args.output)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
