from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step2 import build_step2_enzyme_reference_tables


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build starter enzyme-prior tables linking the 83-microbe panel to enzyme functions."
    )
    parser.add_argument(
        "--microbe-table",
        default=ROOT / "data/processed/step1/step1_microbe_table.csv",
        type=Path,
        help="Expanded Step 1 microbe feature table used as the 83-microbe panel source.",
    )
    parser.add_argument(
        "--enzyme-catalog-path",
        default=ROOT / "data/reference/step2_enzyme_function_catalog.csv",
        type=Path,
        help="Output CSV path for the enzyme function catalog.",
    )
    parser.add_argument(
        "--microbe-enzyme-long-path",
        default=ROOT / "data/reference/step2_microbe_enzyme_prior_long.csv",
        type=Path,
        help="Output CSV path for the long microbe-enzyme prior table.",
    )
    parser.add_argument(
        "--microbe-enzyme-matrix-path",
        default=ROOT / "data/reference/step2_microbe_enzyme_prior_matrix.csv",
        type=Path,
        help="Output CSV path for the wide microbe-enzyme matrix.",
    )
    parser.add_argument(
        "--evidence-ledger-path",
        default=ROOT / "data/reference/step2_microbe_enzyme_evidence_ledger.csv",
        type=Path,
        help="Output CSV path for the combined starter + curated evidence ledger.",
    )
    parser.add_argument(
        "--curation-template-path",
        default=ROOT / "data/reference/step2_microbe_enzyme_curation_template.csv",
        type=Path,
        help="Output CSV path for the editable species/strain enzyme curation worksheet.",
    )
    parser.add_argument(
        "--literature-evidence-path",
        default=None,
        type=Path,
        help="Optional CSV path containing curated species/strain evidence rows to override starter priors.",
    )
    parser.add_argument(
        "--summary-path",
        default=ROOT / "data/reference/step2_enzyme_prior_summary.json",
        type=Path,
        help="Optional JSON summary path.",
    )
    args = parser.parse_args()

    summary = build_step2_enzyme_reference_tables(
        microbe_table_path=args.microbe_table,
        enzyme_catalog_path=args.enzyme_catalog_path,
        microbe_enzyme_long_path=args.microbe_enzyme_long_path,
        microbe_enzyme_matrix_path=args.microbe_enzyme_matrix_path,
        evidence_ledger_path=args.evidence_ledger_path,
        curation_template_path=args.curation_template_path,
        literature_evidence_path=args.literature_evidence_path,
        summary_path=args.summary_path,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
