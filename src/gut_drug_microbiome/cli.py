from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Sequence

import pandas as pd

from gut_drug_microbiome.step1 import build_step1_tables
from gut_drug_microbiome.step1 import train_step1_baseline
from gut_drug_microbiome.step2 import build_step2_input_tables
from gut_drug_microbiome.step2 import build_step2_enzyme_reference_tables
from gut_drug_microbiome.step2 import fetch_uniprot_enzyme_candidates
from gut_drug_microbiome.step3 import BUILTIN_SCENARIOS
from gut_drug_microbiome.step3 import run_step3_simulation
from gut_drug_microbiome.web import serve_web_app


ROOT = Path(__file__).resolve().parents[2]

# 用于命令行参数解析与命令处理的辅助函数。这些函数仅限当前模块内部使用，放置于此是为了避免与核心逻辑模块产生循环导入问题。
def _parse_source_weight(value: str) -> tuple[str, float]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("Source weights must use the form source_dataset=weight")
    key, raw_weight = value.split("=", 1)
    key = key.strip()
    if not key:
        raise argparse.ArgumentTypeError("Source weight key cannot be empty")
    try:
        weight = float(raw_weight)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Invalid weight: {raw_weight}") from exc
    return key, weight

# 简单的字符串清理函数，用于生成文件夹名称等。保留字母数字字符，其他字符替换为下划线，连续下划线合并，首尾下划线去除。如果结果为空，则返回 "simulation" 作为默认名称。
def _slugify(value: str) -> str:
    text = re.sub(r"[^a-zA-Z0-9]+", "_", value.strip().lower()).strip("_")
    return text or "simulation"

# 以下是每个命令的处理函数。这些函数负责调用核心逻辑模块中的函数来执行具体的任务，并将结果以适当的格式返回。每个函数都接受一个 argparse.Namespace 对象作为参数，该对象包含了从命令行解析得到的参数值。
def _command_step1_normalize(args: argparse.Namespace) -> dict[str, object]:
    return build_step1_tables(
        raw_dir=args.raw_dir,
        processed_dir=args.processed_dir,
        labels_config_path=args.labels_config,
        human_use_only=not args.include_non_human_use,
        primary_panel_only=not args.include_non_primary_panel,
    )

# 训练 Step 1 基线模型的命令处理函数。它将从命令行参数中提取必要的信息，并调用 train_step1_baseline 函数来执行训练过程。训练完成后，它会返回一个包含输出目录和其他相关信息的字典。
def _command_step1_train_baseline(args: argparse.Namespace) -> dict[str, object]:
    source_weight_map = {key: value for key, value in args.source_weight}
    return train_step1_baseline(
        modeling_table_path=args.modeling_table,
        output_dir=args.output_dir,
        silver_table_path=args.silver_table,
        promote_feature_table_path=args.promote_feature_table,
        split_mode=args.split_mode,
        random_state=args.random_state,
        test_size=args.test_size,
        n_estimators=args.n_estimators,
        source_weight_map=source_weight_map,
        default_gold_weight=args.default_gold_weight,
        default_silver_weight=args.default_silver_weight,
        enable_promote_head=args.enable_promote_head,
        verbose=args.verbose,
    )


def _command_step2_assemble(args: argparse.Namespace) -> dict[str, object]:
    label_tables = [path for path in args.step2_label_table if path.exists()]
    return build_step2_input_tables(
        step1_predictions_path=args.step1_predictions,
        output_dir=args.output_dir,
        step2_label_table_paths=label_tables or None,
    )


def _command_step2_build_enzyme_priors(args: argparse.Namespace) -> dict[str, object]:
    return build_step2_enzyme_reference_tables(
        microbe_table_path=args.microbe_table,
        enzyme_catalog_path=args.enzyme_catalog_path,
        microbe_enzyme_long_path=args.microbe_enzyme_long_path,
        microbe_enzyme_matrix_path=args.microbe_enzyme_matrix_path,
        evidence_ledger_path=args.evidence_ledger_path,
        curation_template_path=args.curation_template_path,
        literature_evidence_path=args.literature_evidence_path,
        summary_path=args.summary_path,
    )


def _command_step2_fetch_uniprot_enzymes(args: argparse.Namespace) -> dict[str, object]:
    return fetch_uniprot_enzyme_candidates(
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


def _command_step3_simulate(args: argparse.Namespace) -> dict[str, object]:
    root_output_dir = args.output_dir or (ROOT / "predictions/step3" / _slugify(args.drug_query))
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
        return {"output_dir": str(root_output_dir), "n_scenarios": len(summaries)}

    return run_step3_simulation(
        output_dir=root_output_dir,
        scenario_name=args.scenario,
        **run_kwargs,
    )


def _command_web_serve(args: argparse.Namespace) -> None:
    demo_ranking = args.demo_ranking if args.demo_ranking.exists() else None
    supplement_paths = [path for path in args.disease_microbe_supplement if path.exists()]
    serve_web_app(
        host=args.host,
        port=args.port,
        integrated_predictions_path=args.integrated_predictions,
        demo_ranking_path=demo_ranking,
        disease_microbe_reference_path=args.disease_microbe_reference,
        disease_microbe_supplement_paths=tuple(supplement_paths),
        disease_drug_reference_path=args.disease_drug_reference,
        static_dir=args.static_dir,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Unified command line interface for the gut drug-microbiome project.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    step1_parser = subparsers.add_parser("step1", help="Step 1 data processing and training commands.")
    step1_subparsers = step1_parser.add_subparsers(dest="step1_command", required=True)

    step1_normalize = step1_subparsers.add_parser("normalize", help="Normalize Step 1 raw data.")
    step1_normalize.add_argument("--raw-dir", default=ROOT / "data/raw/step1/maier_2018", type=Path)
    step1_normalize.add_argument("--processed-dir", default=ROOT / "data/processed/step1", type=Path)
    step1_normalize.add_argument("--labels-config", default=ROOT / "configs/labeling_rules.yaml", type=Path)
    step1_normalize.add_argument("--include-non-human-use", action="store_true")
    step1_normalize.add_argument("--include-non-primary-panel", action="store_true")
    step1_normalize.set_defaults(handler=_command_step1_normalize)

    step1_train = step1_subparsers.add_parser("train-baseline", help="Train the Step 1 baseline models.")
    step1_train.add_argument("--modeling-table", default=ROOT / "data/processed/step1/step1_modeling_table.csv", type=Path)
    step1_train.add_argument("--output-dir", default=ROOT / "models/step1/baseline_drug_split", type=Path)
    step1_train.add_argument("--silver-table", default=None, type=Path)
    step1_train.add_argument("--promote-feature-table", default=None, type=Path)
    step1_train.add_argument("--split-mode", choices=["random", "drug", "scaffold", "microbe"], default="drug")
    step1_train.add_argument("--random-state", type=int, default=42)
    step1_train.add_argument("--test-size", type=float, default=0.2)
    step1_train.add_argument("--n-estimators", type=int, default=300)
    step1_train.add_argument("--source-weight", action="append", default=[], type=_parse_source_weight)
    step1_train.add_argument("--default-gold-weight", type=float, default=1.0)
    step1_train.add_argument("--default-silver-weight", type=float, default=1.0)
    step1_train.add_argument("--enable-promote-head", action="store_true")
    step1_train.add_argument("--verbose", type=int, default=1)
    step1_train.set_defaults(handler=_command_step1_train_baseline)

    step2_parser = subparsers.add_parser("step2", help="Step 2 assembly commands.")
    step2_subparsers = step2_parser.add_subparsers(dest="step2_command", required=True)
    step2_assemble = step2_subparsers.add_parser("assemble", help="Build Step 2 candidate and modeling tables.")
    step2_assemble.add_argument(
        "--step1-predictions",
        default=ROOT / "predictions/step1/hybrid_scaffold_v1/predictions.csv",
        type=Path,
    )
    step2_assemble.add_argument("--output-dir", default=ROOT / "data/processed/step2", type=Path)
    step2_assemble.add_argument("--step2-label-table", action="append", default=[], type=Path)
    step2_assemble.set_defaults(handler=_command_step2_assemble)

    step2_enzyme = step2_subparsers.add_parser(
        "build-enzyme-priors",
        help="Build starter microbe-enzyme and enzyme-function reference tables for Step 2.",
    )
    step2_enzyme.add_argument("--microbe-table", default=ROOT / "data/processed/step1/step1_microbe_table.csv", type=Path)
    step2_enzyme.add_argument(
        "--enzyme-catalog-path",
        default=ROOT / "data/reference/step2_enzyme_function_catalog.csv",
        type=Path,
    )
    step2_enzyme.add_argument(
        "--microbe-enzyme-long-path",
        default=ROOT / "data/reference/step2_microbe_enzyme_prior_long.csv",
        type=Path,
    )
    step2_enzyme.add_argument(
        "--microbe-enzyme-matrix-path",
        default=ROOT / "data/reference/step2_microbe_enzyme_prior_matrix.csv",
        type=Path,
    )
    step2_enzyme.add_argument(
        "--evidence-ledger-path",
        default=ROOT / "data/reference/step2_microbe_enzyme_evidence_ledger.csv",
        type=Path,
    )
    step2_enzyme.add_argument(
        "--curation-template-path",
        default=ROOT / "data/reference/step2_microbe_enzyme_curation_template.csv",
        type=Path,
    )
    step2_enzyme.add_argument(
        "--literature-evidence-path",
        default=None,
        type=Path,
    )
    step2_enzyme.add_argument(
        "--summary-path",
        default=ROOT / "data/reference/step2_enzyme_prior_summary.json",
        type=Path,
    )
    step2_enzyme.set_defaults(handler=_command_step2_build_enzyme_priors)

    step2_uniprot = step2_subparsers.add_parser(
        "fetch-uniprot-enzymes",
        help="Fetch candidate species-level enzyme evidence for the microbe panel from UniProt SPARQL.",
    )
    step2_uniprot.add_argument("--microbe-table", default=ROOT / "data/processed/step1/step1_microbe_table.csv", type=Path)
    step2_uniprot.add_argument(
        "--enzyme-catalog-path",
        default=ROOT / "data/reference/step2_enzyme_function_catalog.csv",
        type=Path,
    )
    step2_uniprot.add_argument(
        "--output-path",
        default=ROOT / "data/reference/step2_uniprot_enzyme_candidate_evidence.csv",
        type=Path,
    )
    step2_uniprot.add_argument(
        "--raw-output-path",
        default=ROOT / "data/reference/step2_uniprot_protein_enzyme_hits.csv",
        type=Path,
    )
    step2_uniprot.add_argument(
        "--unresolved-output-path",
        default=ROOT / "data/reference/step2_uniprot_unresolved_taxa.csv",
        type=Path,
    )
    step2_uniprot.add_argument(
        "--summary-output-path",
        default=ROOT / "data/reference/step2_uniprot_enzyme_fetch_summary.json",
        type=Path,
    )
    step2_uniprot.add_argument("--endpoint", default="https://sparql.uniprot.org/sparql")
    step2_uniprot.add_argument("--reviewed-only", action="store_true")
    step2_uniprot.add_argument("--limit-microbes", default=None, type=int)
    step2_uniprot.add_argument("--sleep-seconds", default=0.2, type=float)
    step2_uniprot.add_argument("--timeout-seconds", default=90, type=int)
    step2_uniprot.add_argument("--checkpoint-every", default=10, type=int)
    step2_uniprot.set_defaults(handler=_command_step2_fetch_uniprot_enzymes)

    step3_parser = subparsers.add_parser("step3", help="Step 3 simulation commands.")
    step3_subparsers = step3_parser.add_subparsers(dest="step3_command", required=True)
    step3_simulate = step3_subparsers.add_parser("simulate", help="Run the Step 3 simulator.")
    step3_simulate.add_argument(
        "--integrated-predictions",
        default=ROOT / "predictions/step2/baseline_scaffold_v1/predictions.csv",
        type=Path,
    )
    step3_simulate.add_argument("--drug-query", required=True)
    step3_simulate.add_argument("--output-dir", default=None, type=Path)
    step3_simulate.add_argument("--scenario", choices=sorted(BUILTIN_SCENARIOS), default="healthy_reference")
    step3_simulate.add_argument("--all-scenarios", action="store_true")
    step3_simulate.add_argument("--community-table", default=None, type=Path)
    step3_simulate.add_argument("--n-steps", type=int, default=14)
    step3_simulate.add_argument("--initial-dose", type=float, default=1.0)
    step3_simulate.add_argument("--repeat-dose", type=float, default=1.0)
    step3_simulate.add_argument("--dosing-interval", type=int, default=1)
    step3_simulate.add_argument("--drug-clearance-rate", type=float, default=0.12)
    step3_simulate.add_argument("--product-clearance-rate", type=float, default=0.18)
    step3_simulate.add_argument("--metabolism-scale", type=float, default=0.85)
    step3_simulate.add_argument("--effect-scale", type=float, default=0.55)
    step3_simulate.add_argument("--ecology-strength", type=float, default=0.20)
    step3_simulate.add_argument("--experimental-multi-product-enabled", action="store_true")
    step3_simulate.add_argument("--experimental-branching-scale", type=float, default=0.35)
    step3_simulate.add_argument("--experimental-secondary-metabolism-rate", type=float, default=0.10)
    step3_simulate.set_defaults(handler=_command_step3_simulate)

    web_parser = subparsers.add_parser("web", help="Local web app commands.")
    web_subparsers = web_parser.add_subparsers(dest="web_command", required=True)
    web_serve = web_subparsers.add_parser("serve", help="Run the local web app.")
    web_serve.add_argument("--host", default="127.0.0.1")
    web_serve.add_argument("--port", type=int, default=8080)
    web_serve.add_argument(
        "--integrated-predictions",
        default=ROOT / "predictions/step2/baseline_scaffold_v1_83/predictions.csv",
        type=Path,
    )
    web_serve.add_argument(
        "--demo-ranking",
        default=ROOT / "predictions/step3/candidate_screen_demo/candidate_ranking.csv",
        type=Path,
    )
    web_serve.add_argument(
        "--disease-microbe-reference",
        default=ROOT / "data/reference/disease_microbe_dictionary.csv",
        type=Path,
    )
    web_serve.add_argument(
        "--disease-microbe-supplement",
        action="append",
        default=[
            ROOT / "data/reference/disease_microbe_gmrepo_supplement.csv",
            ROOT / "data/reference/disease_microbe_gutm_disorder_supplement.csv",
        ],
        type=Path,
        help="Optional supplemental disease-microbe CSV paths. Can be repeated.",
    )
    web_serve.add_argument(
        "--disease-drug-reference",
        default=ROOT / "data/reference/disease_marketed_drug_catalog.csv",
        type=Path,
    )
    web_serve.add_argument("--static-dir", default=ROOT / "webapp/static", type=Path)
    web_serve.set_defaults(handler=_command_web_serve)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    result = args.handler(args)
    if result is not None:
        print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0
