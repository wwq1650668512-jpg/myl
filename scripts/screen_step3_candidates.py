from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step3 import BUILTIN_SCENARIOS
from gut_drug_microbiome.step3 import run_step3_simulation


def _slugify(value: str) -> str:
    text = re.sub(r"[^a-zA-Z0-9]+", "_", value.strip().lower()).strip("_")
    return text or "drug"


def _load_default_drug_queries(integrated_predictions: Path, limit: int) -> list[str]:
    frame = pd.read_csv(integrated_predictions, low_memory=False, usecols=["prestwick_id", "chemical_name"])
    frame = frame.drop_duplicates().sort_values(["chemical_name", "prestwick_id"]).reset_index(drop=True)
    return frame["chemical_name"].head(limit).astype(str).tolist()


def main() -> None:
    parser = argparse.ArgumentParser(description="Batch screen candidate drugs with the Step 3 simulator.")
    parser.add_argument(
        "--integrated-predictions",
        default=ROOT / "predictions/step2/baseline_scaffold_v1/predictions.csv",
        type=Path,
        help="Integrated Step 1 + Step 2 prediction table.",
    )
    parser.add_argument(
        "--output-dir",
        default=ROOT / "predictions/step3/candidate_screen",
        type=Path,
        help="Directory for candidate-screen outputs.",
    )
    parser.add_argument(
        "--drug-query",
        action="append",
        default=[],
        help="Drug query to simulate. Repeatable.",
    )
    parser.add_argument(
        "--drug-list-file",
        default=None,
        type=Path,
        help="Optional text file with one drug query per line.",
    )
    parser.add_argument(
        "--scenario",
        choices=sorted(BUILTIN_SCENARIOS),
        default="healthy_reference",
        help="Builtin scenario used for candidate screening.",
    )
    parser.add_argument(
        "--community-table",
        default=None,
        type=Path,
        help="Optional real cohort community table. When provided, it overrides the builtin scenario.",
    )
    parser.add_argument("--n-steps", type=int, default=14, help="Number of discrete time steps.")
    parser.add_argument("--initial-dose", type=float, default=1.0, help="Initial parent drug amount.")
    parser.add_argument("--repeat-dose", type=float, default=1.0, help="Repeated dose amount.")
    parser.add_argument("--dosing-interval", type=int, default=1, help="Dose interval in time steps.")
    parser.add_argument("--drug-clearance-rate", type=float, default=0.12, help="Per-step parent drug clearance.")
    parser.add_argument("--product-clearance-rate", type=float, default=0.18, help="Per-step product clearance.")
    parser.add_argument("--metabolism-scale", type=float, default=0.85, help="Scaling factor for microbial metabolism.")
    parser.add_argument("--effect-scale", type=float, default=0.55, help="Scaling factor for Step 1 drug pressure.")
    parser.add_argument("--ecology-strength", type=float, default=0.20, help="Pull back toward the baseline community.")
    parser.add_argument(
        "--default-limit",
        type=int,
        default=0,
        help="If no drug query is provided, optionally auto-screen the first N drugs in the integrated table.",
    )
    args = parser.parse_args()

    queries = [item.strip() for item in args.drug_query if item and item.strip()]
    if args.drug_list_file is not None:
        queries.extend(
            [
                line.strip()
                for line in args.drug_list_file.read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]
        )
    if not queries and args.default_limit > 0:
        queries = _load_default_drug_queries(args.integrated_predictions, limit=args.default_limit)
    if not queries:
        raise SystemExit("Provide --drug-query, --drug-list-file, or --default-limit.")

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    summaries = []
    for drug_query in queries:
        drug_output_dir = output_dir / _slugify(drug_query)
        summary = run_step3_simulation(
            integrated_predictions_path=args.integrated_predictions,
            output_dir=drug_output_dir,
            drug_query=drug_query,
            scenario_name=args.scenario,
            community_table_path=args.community_table,
            n_steps=args.n_steps,
            initial_dose=args.initial_dose,
            repeat_dose=args.repeat_dose,
            dosing_interval=args.dosing_interval,
            drug_clearance_rate=args.drug_clearance_rate,
            product_clearance_rate=args.product_clearance_rate,
            metabolism_scale=args.metabolism_scale,
            effect_scale=args.effect_scale,
            ecology_strength=args.ecology_strength,
        )
        summaries.append(summary)

    ranking = pd.DataFrame(summaries).sort_values(
        ["development_score", "final_health_index", "final_parent_retention_ratio"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    if not ranking.empty:
        ranking["development_rank"] = np.arange(1, len(ranking) + 1)
        if len(ranking) == 1:
            ranking["development_percentile"] = 100.0
        else:
            ranking["development_percentile"] = 100.0 * (len(ranking) - ranking["development_rank"]) / (len(ranking) - 1)
    ranking.to_csv(output_dir / "candidate_ranking.csv", index=False)
    (output_dir / "candidate_ranking.json").write_text(
        json.dumps(ranking.to_dict(orient="records"), indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(json.dumps({"output_dir": str(output_dir), "n_candidates": int(len(ranking))}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
