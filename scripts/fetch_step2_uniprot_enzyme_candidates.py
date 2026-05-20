from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step2 import fetch_uniprot_enzyme_candidates


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch candidate microbe-enzyme evidence from UniProt SPARQL for the 83-microbe panel."
    )
    parser.add_argument(
        "--microbe-table",
        default=ROOT / "data/processed/step1/step1_microbe_table.csv",
        type=Path,
    )
    parser.add_argument(
        "--enzyme-catalog-path",
        default=ROOT / "data/reference/step2_enzyme_function_catalog.csv",
        type=Path,
    )
    parser.add_argument(
        "--output-path",
        default=ROOT / "data/reference/step2_uniprot_enzyme_candidate_evidence.csv",
        type=Path,
        help="Aggregated evidence rows aligned to the project enzyme table schema.",
    )
    parser.add_argument(
        "--raw-output-path",
        default=ROOT / "data/reference/step2_uniprot_protein_enzyme_hits.csv",
        type=Path,
        help="Raw UniProt protein-level enzyme hits.",
    )
    parser.add_argument(
        "--unresolved-output-path",
        default=ROOT / "data/reference/step2_uniprot_unresolved_taxa.csv",
        type=Path,
        help="Microbes that could not be mapped cleanly to UniProt taxonomy/protein hits.",
    )
    parser.add_argument(
        "--summary-output-path",
        default=ROOT / "data/reference/step2_uniprot_enzyme_fetch_summary.json",
        type=Path,
    )
    parser.add_argument(
        "--endpoint",
        default="https://sparql.uniprot.org/sparql",
    )
    parser.add_argument("--reviewed-only", action="store_true")
    parser.add_argument("--limit-microbes", type=int, default=None)
    parser.add_argument("--sleep-seconds", type=float, default=0.2)
    parser.add_argument("--timeout-seconds", type=int, default=90)
    parser.add_argument("--checkpoint-every", type=int, default=10)
    args = parser.parse_args()

    summary = fetch_uniprot_enzyme_candidates(
        microbe_table_path=args.microbe_table,
        enzyme_catalog_path=args.enzyme_catalog_path,
        evidence_output_path=args.output_path,
        raw_output_path=args.raw_output_path,
        unresolved_output_path=args.unresolved_output_path,
        summary_output_path=args.summary_output_path,
        endpoint=args.endpoint,
        reviewed_only=args.reviewed_only,
        limit_microbes=args.limit_microbes,
        sleep_seconds=args.sleep_seconds,
        timeout_seconds=args.timeout_seconds,
        checkpoint_every=args.checkpoint_every,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
