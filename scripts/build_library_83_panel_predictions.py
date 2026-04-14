from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step1 import predict_step1_hybrid
from gut_drug_microbiome.step1 import annotate_compound_semantics
from gut_drug_microbiome.step1 import refine_step1_promote_with_step2
from gut_drug_microbiome.step2 import build_step2_input_tables
from gut_drug_microbiome.step2 import predict_step2_baseline


DEFAULT_SOURCE_PREDICTIONS = ROOT / "predictions/step2/baseline_scaffold_v1/predictions.csv"
DEFAULT_MICROBE_TABLE = ROOT / "data/processed/step1/step1_microbe_table.csv"
DEFAULT_OUTPUT_DIR = ROOT / "predictions/step2/baseline_scaffold_v1_83"
DEFAULT_STEP1_CHEMPROP_PREPARE_DIR = ROOT / "data/processed/step1/chemprop_scaffold/classification"
DEFAULT_STEP1_CHEMPROP_MODEL_PATH = ROOT / "models/step1/chemprop_scaffold_classification_v1/model_0/best.pt"
DEFAULT_STEP1_REGRESSOR_PATH = ROOT / "models/step1/gold_scaffold_split_rdkit_40/regressor.joblib"
DEFAULT_STEP1_REGRESSOR_METRICS_PATH = ROOT / "models/step1/gold_scaffold_split_rdkit_40/metrics.json"
DEFAULT_STEP1_PROMOTE_CLASSIFIER_PATH = (
    ROOT / "models/step1/promote_aux_scaffold_mdipid_plus_promote_literature_v2_40/promote_classifier.joblib"
)
DEFAULT_STEP1_PROMOTE_METRICS_PATH = (
    ROOT / "models/step1/promote_aux_scaffold_mdipid_plus_promote_literature_v2_40/metrics.json"
)
DEFAULT_CROSS_FEEDING_REFERENCE_PATH = ROOT / "data/reference/cross_feeding_edges.csv"
DEFAULT_STEP2_CLASSIFIER_PATH = ROOT / "models/step2/zimmermann_scaffold_split/classifier_full.joblib"
DEFAULT_STEP2_REGRESSOR_PATH = ROOT / "models/step2/zimmermann_scaffold_split/regressor_full.joblib"
DEFAULT_STEP2_METRICS_PATH = ROOT / "models/step2/zimmermann_scaffold_split/metrics.json"
DEFAULT_STEP2_APPLICABILITY_REFERENCE_PATH = ROOT / "models/step2/zimmermann_scaffold_split/applicability_reference.joblib"
DEFAULT_ENZYME_MICROBE_PANEL_PATH = ROOT / "data/reference/step2_microbe_enzyme_prior_long.csv"
DEFAULT_ENZYME_FUNCTION_CATALOG_PATH = ROOT / "data/reference/step2_enzyme_function_catalog.csv"

DRUG_COLUMNS = [
    "prestwick_id",
    "chemical_name",
    "cid_flat",
    "cid_active",
    "cid_main",
    "main_component_smiles",
    "smiles",
    "molecular_formula",
    "molecular_weight",
    "xlogp",
    "tpsa",
    "complexity",
    "volume3d",
    "therapeutic_class",
    "therapeutic_effect",
    "atc_codes",
    "atc_primary_l1",
    "atc_primary_l3",
    "atc_primary_l4",
    "target_species",
    "human_use",
    "veterinary",
    "dose_umol",
    "estimated_intestine_concentration_um",
    "plasma_concentration_um",
    "fraction_excreted_in_feces",
    "fraction_excreted_in_urine",
    "estimated_colon_concentration_um",
    "screen_conc_20_um_as_ug_ml",
    "smiles_length",
    "smiles_uppercase_count",
    "smiles_ring_index_count",
    "smiles_branch_count",
    "smiles_double_bond_count",
    "smiles_halogen_count",
    "canonical_smiles_rdkit",
    "inchikey",
    "murcko_scaffold",
    "rdkit_formula",
    "rdkit_valid_smiles",
    "rdkit_exact_mol_wt",
    "rdkit_logp",
    "rdkit_tpsa",
    "rdkit_molar_refractivity",
    "rdkit_formal_charge",
    "rdkit_heavy_atom_count",
    "rdkit_hbond_donor_count",
    "rdkit_hbond_acceptor_count",
    "rdkit_rotatable_bond_count",
    "rdkit_ring_count",
    "rdkit_aromatic_ring_count",
    "rdkit_aliphatic_ring_count",
    "rdkit_hetero_atom_count",
    "rdkit_fraction_csp3",
]
DRUG_COLUMNS += [f"morgan_fp_{index:03d}" for index in range(256)]

MICROBE_COLUMNS = [
    "nt_code",
    "microbe_label",
    "species_label",
    "species_name",
    "species",
    "strain",
    "speci_cluster",
    "is_unique",
    "biosafety",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "gram_stain",
    "medium_preference",
    "starting_od_96_well_screen",
    "starting_od_384_well_screen",
]


def build_library_pair_table(source_predictions_path: Path, microbe_table_path: Path) -> pd.DataFrame:
    source = pd.read_csv(source_predictions_path, low_memory=False)
    source = source.rename(
        columns={
            "cid_flat_drug": "cid_flat",
            "target_species_drug": "target_species",
            "human_use_drug": "human_use",
            "veterinary_drug": "veterinary",
        }
    )

    existing_drug_columns = [column for column in DRUG_COLUMNS if column in source.columns]
    unique_drug_columns = []
    for column in existing_drug_columns:
        if column not in unique_drug_columns:
            unique_drug_columns.append(column)

    drug_table = source.loc[:, unique_drug_columns].drop_duplicates(subset=["prestwick_id"]).reset_index(drop=True)
    drug_table = annotate_compound_semantics(drug_table)
    microbe_table = pd.read_csv(microbe_table_path, low_memory=False)
    existing_microbe_columns = [column for column in MICROBE_COLUMNS if column in microbe_table.columns]
    microbe_table = microbe_table.loc[:, existing_microbe_columns].drop_duplicates(subset=["nt_code"]).reset_index(drop=True)

    drug_table["_merge_key"] = 1
    microbe_table["_merge_key"] = 1
    pair_table = microbe_table.merge(drug_table, on="_merge_key", how="inner").drop(columns=["_merge_key"])
    pair_table["pair_id"] = pair_table["prestwick_id"].astype(str) + "::" + pair_table["nt_code"].astype(str)
    pair_table["effect_label"] = np.nan
    pair_table["binary_effect_label"] = np.nan
    pair_table["effect_score"] = np.nan
    pair_table["source_dataset"] = "library_panel_83"
    pair_table["label_tier"] = "inference"
    return pair_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Build 83-microbe library predictions for the web app.")
    parser.add_argument(
        "--source-predictions",
        default=DEFAULT_SOURCE_PREDICTIONS,
        type=Path,
        help="Legacy integrated library predictions used as the drug metadata source.",
    )
    parser.add_argument(
        "--microbe-table",
        default=DEFAULT_MICROBE_TABLE,
        type=Path,
        help="Expanded Step 1 microbe feature table.",
    )
    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        type=Path,
        help="Directory for the expanded 83-microbe library predictions.",
    )
    parser.add_argument(
        "--step1-chemprop-prepare-dir",
        default=DEFAULT_STEP1_CHEMPROP_PREPARE_DIR,
        type=Path,
        help="Prepared Chemprop descriptor directory for Step 1.",
    )
    parser.add_argument(
        "--step1-chemprop-model-path",
        default=DEFAULT_STEP1_CHEMPROP_MODEL_PATH,
        type=Path,
        help="Chemprop model path for Step 1 classification.",
    )
    parser.add_argument(
        "--step1-regressor-path",
        default=DEFAULT_STEP1_REGRESSOR_PATH,
        type=Path,
        help="RDKit regressor path for Step 1 regression.",
    )
    parser.add_argument(
        "--step1-regressor-metrics-path",
        default=DEFAULT_STEP1_REGRESSOR_METRICS_PATH,
        type=Path,
        help="Metrics JSON for the Step 1 regressor.",
    )
    parser.add_argument(
        "--step1-promote-classifier-path",
        default=DEFAULT_STEP1_PROMOTE_CLASSIFIER_PATH,
        type=Path,
        help="Optional Step 1 promote auxiliary classifier used during Step 2-aware rescoring.",
    )
    parser.add_argument(
        "--step1-promote-metrics-path",
        default=DEFAULT_STEP1_PROMOTE_METRICS_PATH,
        type=Path,
        help="Metrics JSON for the optional Step 1 promote auxiliary classifier.",
    )
    parser.add_argument(
        "--cross-feeding-reference-path",
        default=DEFAULT_CROSS_FEEDING_REFERENCE_PATH,
        type=Path,
        help="Curated cross-feeding reference table used for conservative promote rescoring.",
    )
    parser.add_argument(
        "--step2-classifier-path",
        default=DEFAULT_STEP2_CLASSIFIER_PATH,
        type=Path,
        help="Step 2 classifier path.",
    )
    parser.add_argument(
        "--step2-regressor-path",
        default=DEFAULT_STEP2_REGRESSOR_PATH,
        type=Path,
        help="Step 2 regressor path.",
    )
    parser.add_argument(
        "--step2-metrics-path",
        default=DEFAULT_STEP2_METRICS_PATH,
        type=Path,
        help="Metrics JSON for Step 2 baseline inference.",
    )
    parser.add_argument(
        "--step2-applicability-reference-path",
        default=DEFAULT_STEP2_APPLICABILITY_REFERENCE_PATH,
        type=Path,
        help="Applicability reference for Step 2 baseline inference.",
    )
    parser.add_argument(
        "--enzyme-microbe-panel-path",
        default=DEFAULT_ENZYME_MICROBE_PANEL_PATH,
        type=Path,
        help="Optional curated microbe-enzyme prior table for the expanded panel.",
    )
    parser.add_argument(
        "--enzyme-function-catalog-path",
        default=DEFAULT_ENZYME_FUNCTION_CATALOG_PATH,
        type=Path,
        help="Optional enzyme function catalog used for enzyme-derived support.",
    )
    args = parser.parse_args()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    pair_table = build_library_pair_table(
        source_predictions_path=args.source_predictions,
        microbe_table_path=args.microbe_table,
    )
    pair_table_path = output_dir / "library_83_pair_table.csv"
    pair_table.to_csv(pair_table_path, index=False)

    step1_output_dir = output_dir / "step1"
    predict_step1_hybrid(
        input_table_path=pair_table_path,
        output_dir=step1_output_dir,
        classification_prepare_dir=args.step1_chemprop_prepare_dir,
        chemprop_model_path=args.step1_chemprop_model_path,
        regressor_path=args.step1_regressor_path,
        regressor_metrics_path=args.step1_regressor_metrics_path,
    )

    step2_input_dir = output_dir / "step2_inputs"
    step2_input_summary = build_step2_input_tables(
        step1_predictions_path=step1_output_dir / "predictions.csv",
        output_dir=step2_input_dir,
        step2_label_table_paths=None,
    )

    step2_output_dir = output_dir
    step2_summary = predict_step2_baseline(
        input_table_path=step2_input_dir / "step2_candidate_pairs_full.csv",
        output_dir=step2_output_dir,
        classifier_path=args.step2_classifier_path,
        regressor_path=args.step2_regressor_path,
        metrics_path=args.step2_metrics_path,
        applicability_reference_path=args.step2_applicability_reference_path,
        enzyme_microbe_panel_path=args.enzyme_microbe_panel_path,
        enzyme_function_catalog_path=args.enzyme_function_catalog_path,
    )

    predictions = pd.read_csv(output_dir / "predictions.csv", low_memory=False)
    predictions = refine_step1_promote_with_step2(
        predictions,
        promote_classifier_path=args.step1_promote_classifier_path,
        promote_metrics_path=args.step1_promote_metrics_path,
        cross_feeding_reference_path=args.cross_feeding_reference_path,
    )
    predictions.to_csv(output_dir / "predictions.csv", index=False)
    summary = {
        "source_predictions": str(args.source_predictions),
        "microbe_table": str(args.microbe_table),
        "output_dir": str(output_dir),
        "pair_table_path": str(pair_table_path),
        "step1_predictions_path": str(step1_output_dir / "predictions.csv"),
        "step2_predictions_path": str(output_dir / "predictions.csv"),
        "n_pairs": int(len(predictions)),
        "n_drugs": int(predictions["prestwick_id"].nunique()),
        "n_microbes": int(predictions["nt_code"].nunique()),
        "step2_input_summary_path": str(step2_input_dir / "step2_summary.json"),
        "step2_summary_path": str(output_dir / "summary.json"),
        "step2_input_summary": step2_input_summary,
        "step2_summary": step2_summary,
    }
    (output_dir / "build_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"output_dir": str(output_dir), "n_pairs": summary["n_pairs"], "n_microbes": summary["n_microbes"]}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
