from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step1 import predict_step1_hybrid


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the Step 1 hybrid predictor on a standardized pair table.")
    parser.add_argument(
        "--input-table",
        default=ROOT / "data/processed/step1/step1_modeling_table.csv",
        type=Path,
        help="Standardized pair table with drug and microbe features.",
    )
    parser.add_argument(
        "--output-dir",
        default=ROOT / "predictions/step1/hybrid_scaffold_v1",
        type=Path,
        help="Directory for hybrid prediction outputs.",
    )
    parser.add_argument(
        "--classification-prepare-dir",
        default=ROOT / "data/processed/step1/chemprop_scaffold/classification",
        type=Path,
        help="Chemprop classification prepare directory containing descriptor assets.",
    )
    parser.add_argument(
        "--chemprop-model-path",
        default=ROOT / "models/step1/chemprop_scaffold_classification_v1/model_0/best.pt",
        type=Path,
        help="Trained Chemprop classification checkpoint.",
    )
    parser.add_argument(
        "--regressor-path",
        default=ROOT / "models/step1/gold_scaffold_split_rdkit_40/regressor.joblib",
        type=Path,
        help="RDKit baseline regressor joblib.",
    )
    parser.add_argument(
        "--regressor-metrics-path",
        default=ROOT / "models/step1/gold_scaffold_split_rdkit_40/metrics.json",
        type=Path,
        help="Metrics JSON that records the regressor feature schema.",
    )
    parser.add_argument(
        "--inhibit-probability-threshold",
        type=float,
        default=0.5,
        help="Threshold for mapping Chemprop inhibit probability to the binary inhibit label.",
    )
    parser.add_argument(
        "--promote-score-threshold",
        type=float,
        default=0.2,
        help="Threshold for mapping positive effect_score to the promote label when Chemprop does not call inhibit.",
    )
    parser.add_argument(
        "--accelerator",
        default="cpu",
        help="Accelerator forwarded to Chemprop predict.",
    )
    parser.add_argument(
        "--devices",
        default="1",
        help="Devices forwarded to Chemprop predict.",
    )
    args = parser.parse_args()

    summary = predict_step1_hybrid(
        input_table_path=args.input_table,
        output_dir=args.output_dir,
        classification_prepare_dir=args.classification_prepare_dir,
        chemprop_model_path=args.chemprop_model_path,
        regressor_path=args.regressor_path,
        regressor_metrics_path=args.regressor_metrics_path,
        inhibit_probability_threshold=args.inhibit_probability_threshold,
        promote_score_threshold=args.promote_score_threshold,
        accelerator=args.accelerator,
        devices=args.devices,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
