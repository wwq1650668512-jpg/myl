from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step3 import BUILTIN_SCENARIOS
from gut_drug_microbiome.step3 import run_step3_simulation


def _slugify(value: str) -> str:
    text = re.sub(r"[^a-zA-Z0-9]+", "_", value.strip().lower()).strip("_")
    return text or "simulation"


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Step 3 gut community simulations.")
    parser.add_argument(
        "--integrated-predictions",
        default=ROOT / "predictions/step2/baseline_scaffold_v1/predictions.csv",
        type=Path,
        help="Integrated Step 1 + Step 2 prediction table.",
    )
    parser.add_argument(
        "--drug-query",
        required=True,
        help="Prestwick ID or chemical name query.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        type=Path,
        help="Directory for simulation outputs. Defaults to predictions/step3/<drug_slug>.",
    )
    parser.add_argument(
        "--scenario",
        choices=sorted(BUILTIN_SCENARIOS),
        default="healthy_reference",
        help="Builtin community scenario to simulate.",
    )
    parser.add_argument(
        "--all-scenarios",
        action="store_true",
        help="Run every builtin scenario and write an aggregated summary.",
    )
    parser.add_argument(
        "--community-table",
        default=None,
        type=Path,
        help="Optional custom community table with microbe ids/names and abundances.",
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
        "--experimental-multi-product-enabled",
        action="store_true",
        help="Enable an additive experimental metabolite-branching score without changing the base score.",
    )
    parser.add_argument(
        "--experimental-branching-scale",
        type=float,
        default=0.35,
        help="Additional direct product burden scale under the experimental mode.",
    )
    parser.add_argument(
        "--experimental-secondary-metabolism-rate",
        type=float,
        default=0.10,
        help="Recursive downstream-burden rate for the experimental metabolite pool.",
    )
    args = parser.parse_args()

    drug_slug = _slugify(args.drug_query)
    root_output_dir = args.output_dir or (ROOT / "predictions/step3" / drug_slug)

    run_kwargs = {
        "integrated_predictions_path": args.integrated_predictions,
        "drug_query": args.drug_query,
        "community_table_path": args.community_table,
        "n_steps": args.n_steps,
        "initial_dose": args.initial_dose,
        "repeat_dose": args.repeat_dose,
        "dosing_interval": args.dosing_interval,
        "drug_clearance_rate": args.drug_clearance_rate,
        "product_clearance_rate": args.product_clearance_rate,
        "metabolism_scale": args.metabolism_scale,
        "effect_scale": args.effect_scale,
        "ecology_strength": args.ecology_strength,
        "experimental_multi_product_enabled": args.experimental_multi_product_enabled,
        "experimental_branching_scale": args.experimental_branching_scale,
        "experimental_secondary_metabolism_rate": args.experimental_secondary_metabolism_rate,
    }

    if args.all_scenarios and args.community_table is None:
        summaries = []
        for scenario_name in sorted(BUILTIN_SCENARIOS):
            scenario_output_dir = root_output_dir / scenario_name
            summary = run_step3_simulation(
                output_dir=scenario_output_dir,
                scenario_name=scenario_name,
                **run_kwargs,
            )
            summaries.append(summary)
        summary_frame = pd.DataFrame(summaries)
        summary_frame.to_csv(root_output_dir / "scenario_grid_summary.csv", index=False)
        (root_output_dir / "scenario_grid_summary.json").write_text(
            json.dumps(summaries, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        print(json.dumps({"output_dir": str(root_output_dir), "n_scenarios": len(summaries)}, indent=2, ensure_ascii=False))
        return

    summary = run_step3_simulation(
        output_dir=root_output_dir,
        scenario_name=args.scenario,
        **run_kwargs,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
