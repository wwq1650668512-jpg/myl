from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step2 import predict_step2_baseline


def main() -> None:
    parser = argparse.ArgumentParser(description="Predict Step 2 metabolism outcomes with the baseline models.")
    parser.add_argument(
        "--input-table",
        default=ROOT / "data/processed/step2/step2_candidate_pairs_full.csv",
        type=Path,
        help="Input candidate table to score.",
    )
    parser.add_argument(
        "--output-dir",
        default=ROOT / "predictions/step2/baseline_scaffold_v1",
        type=Path,
        help="Directory for prediction outputs.",
    )
    parser.add_argument(
        "--classifier-path",
        default=ROOT / "models/step2/zimmermann_scaffold_split/classifier_full.joblib",
        type=Path,
        help="Path to the fitted Step 2 classifier.",
    )
    parser.add_argument(
        "--regressor-path",
        default=ROOT / "models/step2/zimmermann_scaffold_split/regressor_full.joblib",
        type=Path,
        help="Path to the fitted Step 2 regressor.",
    )
    parser.add_argument(
        "--metrics-path",
        default=ROOT / "models/step2/zimmermann_scaffold_split/metrics.json",
        type=Path,
        help="Metrics JSON emitted during training.",
    )
    parser.add_argument(
        "--applicability-reference-path",
        default=ROOT / "models/step2/zimmermann_scaffold_split/applicability_reference.joblib",
        type=Path,
        help="Applicability reference emitted during training.",
    )
    parser.add_argument(
        "--mechanism-reference-path",
        default=ROOT / "models/step2/zimmermann_scaffold_split/mechanism_reference.joblib",
        type=Path,
        help="Semi-mechanistic Step 2 reference used to project reaction class / product / gene evidence.",
    )
    parser.add_argument(
        "--enzyme-microbe-panel-path",
        default=ROOT / "data/reference/step2_microbe_enzyme_prior_long.csv",
        type=Path,
        help="Optional curated microbe-enzyme prior table for enzyme-derived mechanism support.",
    )
    parser.add_argument(
        "--enzyme-function-catalog-path",
        default=ROOT / "data/reference/step2_enzyme_function_catalog.csv",
        type=Path,
        help="Optional enzyme function catalog matched against compound semantics.",
    )
    parser.add_argument(
        "--probability-threshold",
        type=float,
        default=None,
        help="Optional metabolism probability threshold. Defaults to the training config.",
    )
    parser.add_argument(
        "--similarity-threshold",
        type=float,
        default=0.25,
        help="Drug fingerprint Jaccard threshold used by the applicability flag.",
    )
    args = parser.parse_args()

    summary = predict_step2_baseline(
        input_table_path=args.input_table,
        output_dir=args.output_dir,
        classifier_path=args.classifier_path,
        regressor_path=args.regressor_path,
        metrics_path=args.metrics_path,
        applicability_reference_path=args.applicability_reference_path,
        mechanism_reference_path=args.mechanism_reference_path,
        enzyme_microbe_panel_path=args.enzyme_microbe_panel_path,
        enzyme_function_catalog_path=args.enzyme_function_catalog_path,
        probability_threshold=args.probability_threshold,
        similarity_threshold=args.similarity_threshold,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
