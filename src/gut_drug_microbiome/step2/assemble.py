from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .normalize import NORMALIZED_STEP2_COLUMNS
from .normalize import read_table_auto


STEP1_RENAME_MAP = {
    "effect_label": "step1_observed_effect_label",
    "binary_effect_label": "step1_observed_binary_effect_label",
    "effect_score": "step1_observed_effect_score",
    "predicted_inhibit_probability": "step1_predicted_inhibit_probability",
    "predicted_binary_effect_label": "step1_predicted_binary_effect_label",
    "predicted_effect_score": "step1_predicted_effect_score",
    "predicted_effect_label_hybrid": "step1_predicted_effect_label_hybrid",
    "predicted_effect_magnitude": "step1_predicted_effect_magnitude",
}

STEP2_LABEL_RENAME_MAP = {
    "metabolism_label": "step2_metabolism_label",
    "reaction_class": "step2_reaction_class",
    "parent_depletion_fraction": "step2_parent_depletion_fraction",
    "product_ids": "step2_product_ids",
    "evidence_gene_ids": "step2_evidence_gene_ids",
    "source_dataset": "step2_source_dataset",
    "label_tier": "step2_label_tier",
    "source_scope": "step2_source_scope",
    "source_record_id": "step2_source_record_id",
    "raw_metabolism_label": "step2_raw_metabolism_label",
    "raw_reaction_class": "step2_raw_reaction_class",
}


def _ensure_pair_id(frame: pd.DataFrame) -> pd.DataFrame:
    """Ensure a pair_id column exists for drug-microbe join operations."""
    result = frame.copy()
    if "pair_id" not in result.columns and {"prestwick_id", "nt_code"}.issubset(result.columns):
        result["pair_id"] = result["prestwick_id"].astype(str) + "::" + result["nt_code"].astype(str)
    return result


def _prepare_step1_candidate_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Rename and augment Step 1 outputs so they can seed Step 2 candidate pairs."""
    candidate = _ensure_pair_id(frame)
    candidate = candidate.rename(columns=STEP1_RENAME_MAP)
    candidate["step1_has_smiles"] = candidate.get("smiles", pd.Series(np.nan, index=candidate.index)).notna()
    candidate["step1_predicted_inhibit_flag"] = candidate.get(
        "step1_predicted_binary_effect_label",
        pd.Series(np.nan, index=candidate.index),
    ).eq("inhibit")
    candidate["step1_predicted_promote_flag"] = candidate.get(
        "step1_predicted_effect_label_hybrid",
        pd.Series(np.nan, index=candidate.index),
    ).eq("promote")
    candidate["step1_predicted_no_effect_flag"] = candidate.get(
        "step1_predicted_effect_label_hybrid",
        pd.Series(np.nan, index=candidate.index),
    ).eq("no_effect")
    return candidate


def _load_step2_label_tables(paths: list[str | Path] | None) -> pd.DataFrame | None:
    """Load and concatenate normalized Step 2 label tables after schema validation."""
    if not paths:
        return None
    frames = []
    for path in paths:
        table = read_table_auto(path)
        missing = [column for column in NORMALIZED_STEP2_COLUMNS if column not in table.columns]
        if missing:
            raise ValueError(f"Step 2 normalized label table missing columns {missing}: {path}")
        frames.append(table.loc[:, NORMALIZED_STEP2_COLUMNS].copy())
    if not frames:
        return None
    return pd.concat(frames, ignore_index=True, sort=False)


def build_step2_input_tables(
    step1_predictions_path: str | Path,
    output_dir: str | Path,
    step2_label_table_paths: list[str | Path] | None = None,
) -> dict[str, object]:
    """Build Step 2 candidate and modeling tables by combining Step 1 outputs with labels.

    Args:
        step1_predictions_path: Step 1 prediction or modeling table used as candidate pairs.
        output_dir: Directory where Step 2 CSV outputs and summary JSON are written.
        step2_label_table_paths: Optional normalized Step 2 label tables to merge in.

    Returns:
        A summary dictionary with counts and generated file paths.
    """
    step1_predictions_path = Path(step1_predictions_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    step1_raw = read_table_auto(step1_predictions_path)
    candidate_full = _prepare_step1_candidate_frame(step1_raw)
    candidate_full = candidate_full.drop_duplicates(subset=["pair_id"]).reset_index(drop=True)

    candidate_full_path = output_dir / "step2_candidate_pairs_full.csv"
    candidate_full.to_csv(candidate_full_path, index=False)

    candidate_slim_columns = [
        "pair_id",
        "prestwick_id",
        "nt_code",
        "smiles",
        "step1_observed_effect_label",
        "step1_observed_binary_effect_label",
        "step1_observed_effect_score",
        "step1_predicted_inhibit_probability",
        "step1_predicted_binary_effect_label",
        "step1_predicted_effect_score",
        "step1_predicted_effect_label_hybrid",
        "step1_predicted_effect_magnitude",
        "step1_has_smiles",
        "step1_predicted_inhibit_flag",
        "step1_predicted_promote_flag",
        "step1_predicted_no_effect_flag",
    ]
    existing_candidate_slim_columns = [column for column in candidate_slim_columns if column in candidate_full.columns]
    candidate_slim = candidate_full.loc[:, existing_candidate_slim_columns].copy()
    candidate_slim_path = output_dir / "step2_candidate_pairs_slim.csv"
    candidate_slim.to_csv(candidate_slim_path, index=False)

    labels_long = _load_step2_label_tables(step2_label_table_paths)
    if labels_long is not None:
        labels_long_path = output_dir / "step2_label_table_long.csv"
        labels_long.to_csv(labels_long_path, index=False)
        modeling = candidate_full.merge(
            labels_long.rename(columns=STEP2_LABEL_RENAME_MAP),
            on=["pair_id", "prestwick_id", "nt_code"],
            how="left",
        )
        modeling["step2_label_available"] = modeling["step2_metabolism_label"].notna()
        modeling["step2_has_resolved_products"] = modeling["step2_product_ids"].fillna("").astype(str).str.strip().ne("")
    else:
        labels_long_path = None
        modeling = candidate_full.copy()
        for column in STEP2_LABEL_RENAME_MAP.values():
            modeling[column] = np.nan
        modeling["step2_label_available"] = False
        modeling["step2_has_resolved_products"] = False

    modeling_path = output_dir / "step2_modeling_table.csv"
    modeling.to_csv(modeling_path, index=False)

    summary = {
        "step1_predictions_path": str(step1_predictions_path),
        "output_dir": str(output_dir),
        "step2_label_table_paths": [] if step2_label_table_paths is None else [str(path) for path in step2_label_table_paths],
        "n_candidate_pairs": int(len(candidate_full)),
        "n_pairs_with_smiles": int(candidate_full["step1_has_smiles"].sum()),
        "candidate_step1_hybrid_counts": {
            str(key): int(value)
            for key, value in candidate_full["step1_predicted_effect_label_hybrid"].fillna("missing").value_counts().to_dict().items()
        },
        "n_step2_label_rows": int(0 if labels_long is None else len(labels_long)),
        "n_labeled_modeling_rows": int(modeling["step2_label_available"].sum()),
        "candidate_full_path": str(candidate_full_path),
        "candidate_slim_path": str(candidate_slim_path),
        "labels_long_path": None if labels_long_path is None else str(labels_long_path),
        "modeling_path": str(modeling_path),
    }
    if labels_long is not None:
        summary["step2_metabolism_label_counts"] = {
            str(key): int(value)
            for key, value in labels_long["metabolism_label"].fillna("missing").value_counts().to_dict().items()
        }

    (output_dir / "step2_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    return summary
