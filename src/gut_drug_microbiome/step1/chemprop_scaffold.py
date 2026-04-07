from __future__ import annotations

import json
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import dump
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.metrics import accuracy_score
from sklearn.metrics import auc
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import f1_score
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
from sklearn.metrics import precision_recall_curve
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
from sklearn.metrics import r2_score
from sklearn.metrics import roc_auc_score

from .train_baseline import _make_split


CHEMPROP_NUMERIC_DESCRIPTOR_COLUMNS = [
    "molecular_weight",
    "xlogp",
    "tpsa",
    "complexity",
    "volume3d",
    "dose_umol",
    "estimated_intestine_concentration_um",
    "plasma_concentration_um",
    "fraction_excreted_in_feces",
    "fraction_excreted_in_urine",
    "estimated_colon_concentration_um",
    "screen_conc_20_um_as_ug_ml",
    "starting_od_96_well_screen",
    "starting_od_384_well_screen",
]

CHEMPROP_CATEGORICAL_DESCRIPTOR_COLUMNS = [
    "target_species",
    "therapeutic_class",
    "therapeutic_effect",
    "atc_primary_l1",
    "atc_primary_l3",
    "atc_primary_l4",
    "human_use",
    "veterinary",
    "species_label",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "gram_stain",
    "medium_preference",
    "biosafety",
    "is_unique",
]

SMILES_PRIORITY = [
    "canonical_smiles_rdkit",
    "main_component_smiles",
    "canonical_smiles",
    "iso_smiles",
    "smiles",
]


def _resolve_chemprop_bin(require_installed: bool = True) -> str:
    """Locate the Chemprop CLI executable in the current environment."""
    candidate_paths = []

    # Prefer the CLI installed alongside the current Python interpreter so a
    # dedicated venv works even when PATH is not activated.
    candidate_paths.append(Path(sys.executable).parent / "chemprop")

    # Some environments expose only the resolved interpreter path, so keep a
    # second candidate based on the symlink target parent as a fallback.
    try:
        candidate_paths.append(Path(sys.executable).resolve().parent / "chemprop")
    except OSError:
        pass

    for candidate in candidate_paths:
        if candidate.exists():
            return str(candidate)

    chemprop_bin = shutil.which("chemprop")
    if chemprop_bin is not None:
        return chemprop_bin

    if require_installed:
        raise RuntimeError(
            "Chemprop CLI is not installed in the current environment. "
            "Prepare data with prepare_step1_chemprop_inputs first, then install Chemprop in a dedicated environment."
        )
    return "chemprop"


def _pick_smiles(frame: pd.DataFrame) -> pd.Series:
    """Pick one preferred SMILES string per row from the available chemistry columns."""
    result = pd.Series(np.nan, index=frame.index, dtype=object)
    for column in SMILES_PRIORITY:
        if column not in frame.columns:
            continue
        values = frame[column]
        if values.dtype != object:
            values = values.astype(object)
        mask = result.isna() & values.notna() & values.astype(str).str.strip().ne("")
        result.loc[mask] = values.loc[mask]
    return result


def _available_columns(frame: pd.DataFrame, desired: list[str]) -> list[str]:
    """Return requested descriptor columns that are present in the frame."""
    return [column for column in desired if column in frame.columns]


def _prepare_descriptor_frame(frame: pd.DataFrame, numeric_features: list[str], categorical_features: list[str]) -> pd.DataFrame:
    """Cast descriptor columns into numeric/string form before preprocessing."""
    prepared = frame.copy()
    for column in numeric_features:
        if column in prepared.columns:
            prepared[column] = pd.to_numeric(prepared[column], errors="coerce")
    for column in categorical_features:
        if column in prepared.columns:
            prepared[column] = prepared[column].map(lambda value: np.nan if pd.isna(value) else str(value))
    return prepared


def _descriptor_features(train_frame: pd.DataFrame, full_frame: pd.DataFrame) -> tuple[list[str], list[str]]:
    """Choose descriptor columns that are both available and informative in training data."""
    numeric_features = [
        column
        for column in _available_columns(full_frame, CHEMPROP_NUMERIC_DESCRIPTOR_COLUMNS)
        if pd.to_numeric(train_frame[column], errors="coerce").notna().any()
    ]
    categorical_features = _available_columns(full_frame, CHEMPROP_CATEGORICAL_DESCRIPTOR_COLUMNS)
    return numeric_features, categorical_features


def _fit_descriptor_preprocessor(
    train_frame: pd.DataFrame,
    full_frame: pd.DataFrame,
) -> tuple[ColumnTransformer, list[str], list[str]]:
    """Fit the descriptor preprocessing pipeline used by Chemprop side features."""
    numeric_features, categorical_features = _descriptor_features(train_frame, full_frame)
    numeric_pipeline = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
        ]
    )
    categorical_pipeline = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="constant", fill_value="missing")),
            ("encoder", OneHotEncoder(handle_unknown="ignore", sparse_output=False)),
        ]
    )
    preprocessor = ColumnTransformer(
        transformers=[
            ("numeric", numeric_pipeline, numeric_features),
            ("categorical", categorical_pipeline, categorical_features),
        ],
        remainder="drop",
        sparse_threshold=0.0,
    )
    preprocessor.fit(train_frame.loc[:, numeric_features + categorical_features])
    return preprocessor, numeric_features, categorical_features


def _split_gold_frame(
    gold_frame: pd.DataFrame,
    split_mode: str,
    random_state: int,
    test_size: float,
    val_size: float,
) -> pd.Series:
    """Assign train/val/test split labels to gold data using grouped splitting rules."""
    target = gold_frame["binary_effect_label"]
    train_idx, test_idx = _make_split(
        frame=gold_frame,
        split_mode=split_mode,
        target=target,
        random_state=random_state,
        test_size=test_size,
    )

    trainval_frame = gold_frame.loc[train_idx].reset_index(drop=False).rename(columns={"index": "original_index"})
    inner_target = trainval_frame["binary_effect_label"]
    inner_train_idx, inner_val_idx = _make_split(
        frame=trainval_frame,
        split_mode=split_mode,
        target=inner_target,
        random_state=random_state,
        test_size=val_size,
    )

    split = pd.Series("unused", index=gold_frame.index, dtype=object)
    split.loc[trainval_frame.loc[inner_train_idx, "original_index"]] = "train"
    split.loc[trainval_frame.loc[inner_val_idx, "original_index"]] = "val"
    split.loc[gold_frame.index[test_idx]] = "test"
    return split


def _finalize_dataset_frame(frame: pd.DataFrame, target_column: str) -> pd.DataFrame:
    """Export the subset of columns expected by Chemprop training and evaluation."""
    columns_to_keep = [
        "pair_id",
        "smiles",
        target_column,
        "split",
        "source_dataset",
        "label_tier",
        "prestwick_id",
        "nt_code",
        "effect_label",
        "binary_effect_label",
        "effect_score",
    ]
    existing_columns = [column for column in columns_to_keep if column in frame.columns]
    result = frame[existing_columns].copy()
    result = result.rename(columns={target_column: "target"})
    result.insert(0, "row_id", np.arange(len(result)))
    return result


def _save_descriptors(
    dataset_frame: pd.DataFrame,
    train_mask: pd.Series,
    output_dir: Path,
) -> dict[str, object]:
    """Fit descriptor preprocessing, save side-feature arrays, and return metadata."""
    train_frame = dataset_frame.loc[train_mask].copy()
    preprocessor, numeric_features, categorical_features = _fit_descriptor_preprocessor(train_frame, dataset_frame)
    transformed = preprocessor.transform(dataset_frame.loc[:, numeric_features + categorical_features])
    feature_names = [str(value) for value in preprocessor.get_feature_names_out()]

    descriptors_path = output_dir / "descriptors.npz"
    np.savez_compressed(descriptors_path, transformed)
    (output_dir / "descriptor_feature_names.json").write_text(
        json.dumps(feature_names, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    preprocessor_path = output_dir / "descriptor_preprocessor.joblib"
    dump(preprocessor, preprocessor_path)
    descriptor_schema = {
        "numeric_descriptor_features": numeric_features,
        "categorical_descriptor_features": categorical_features,
        "feature_names": feature_names,
    }
    (output_dir / "descriptor_schema.json").write_text(
        json.dumps(descriptor_schema, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return {
        "descriptors_path": str(descriptors_path),
        "descriptor_preprocessor_path": str(preprocessor_path),
        "descriptor_schema_path": str(output_dir / "descriptor_schema.json"),
        "n_descriptor_features": int(transformed.shape[1]),
        "numeric_descriptor_features": numeric_features,
        "categorical_descriptor_features": categorical_features,
    }


def prepare_step1_chemprop_inputs(
    modeling_table_path: str | Path,
    output_dir: str | Path,
    silver_table_path: str | Path | None = None,
    split_mode: str = "scaffold",
    random_state: int = 42,
    test_size: float = 0.2,
    val_size: float = 0.1,
    positive_label: str = "inhibit",
) -> dict[str, object]:
    """Prepare Chemprop-ready Step 1 datasets, splits, and descriptor artifacts.

    Args:
        modeling_table_path: Normalized Step 1 modeling table.
        output_dir: Root directory where classification/regression datasets are written.
        silver_table_path: Optional silver table appended only to classification training.
        split_mode: Split strategy shared with the baseline pipeline.
        random_state: Random seed for split reproducibility.
        test_size: Fraction reserved for the held-out test split.
        val_size: Fraction of the training pool reserved for validation.
        positive_label: Step 1 class mapped to the positive binary target.

    Returns:
        A summary dictionary describing exported datasets and descriptor files.
    """
    modeling_table_path = Path(modeling_table_path)
    output_dir = Path(output_dir)
    classification_dir = output_dir / "classification"
    regression_dir = output_dir / "regression"
    classification_dir.mkdir(parents=True, exist_ok=True)
    regression_dir.mkdir(parents=True, exist_ok=True)

    gold = pd.read_csv(modeling_table_path, low_memory=False)
    gold = gold.dropna(subset=["effect_label", "binary_effect_label", "effect_score"]).reset_index(drop=True).copy()
    gold["smiles"] = _pick_smiles(gold)
    gold = gold.dropna(subset=["smiles"]).reset_index(drop=True)
    gold["source_dataset"] = gold.get("source_dataset", "maier_2018")
    gold["label_tier"] = gold.get("label_tier", "gold")
    gold["split"] = _split_gold_frame(
        gold_frame=gold,
        split_mode=split_mode,
        random_state=random_state,
        test_size=test_size,
        val_size=val_size,
    )

    gold["classification_target"] = gold["effect_label"].eq(positive_label).astype(int)
    gold["regression_target"] = gold["effect_score"].astype(float)

    silver = None
    if silver_table_path is not None and Path(silver_table_path).exists():
        silver = pd.read_csv(silver_table_path, low_memory=False)
        silver = silver.dropna(subset=["effect_label"]).reset_index(drop=True).copy()
        silver["smiles"] = _pick_smiles(silver)
        silver = silver.dropna(subset=["smiles"]).reset_index(drop=True)
        silver["classification_target"] = silver["effect_label"].eq(positive_label).astype(int)
        silver["split"] = "train"
        silver["source_dataset"] = silver.get("source_dataset", "silver")
        silver["label_tier"] = silver.get("label_tier", "silver")

    classification_frame = gold.copy()
    if silver is not None:
        classification_frame = pd.concat([classification_frame, silver], ignore_index=True, sort=False)
    classification_frame = _prepare_descriptor_frame(
        classification_frame,
        numeric_features=_available_columns(classification_frame, CHEMPROP_NUMERIC_DESCRIPTOR_COLUMNS),
        categorical_features=_available_columns(classification_frame, CHEMPROP_CATEGORICAL_DESCRIPTOR_COLUMNS),
    )

    regression_frame = _prepare_descriptor_frame(
        gold.copy(),
        numeric_features=_available_columns(gold, CHEMPROP_NUMERIC_DESCRIPTOR_COLUMNS),
        categorical_features=_available_columns(gold, CHEMPROP_CATEGORICAL_DESCRIPTOR_COLUMNS),
    )

    classification_csv = classification_dir / "dataset.csv"
    classification_export = _finalize_dataset_frame(classification_frame, "classification_target")
    classification_export.to_csv(classification_csv, index=False)
    classification_descriptor_summary = _save_descriptors(
        classification_frame,
        train_mask=classification_frame["split"].eq("train"),
        output_dir=classification_dir,
    )

    regression_csv = regression_dir / "dataset.csv"
    regression_export = _finalize_dataset_frame(regression_frame, "regression_target")
    regression_export.to_csv(regression_csv, index=False)
    regression_descriptor_summary = _save_descriptors(
        regression_frame,
        train_mask=regression_frame["split"].eq("train"),
        output_dir=regression_dir,
    )

    summary = {
        "modeling_table_path": str(modeling_table_path),
        "silver_table_path": None if silver_table_path is None else str(silver_table_path),
        "output_dir": str(output_dir),
        "split_mode": split_mode,
        "random_state": random_state,
        "test_size": test_size,
        "val_size": val_size,
        "positive_label": positive_label,
        "classification": {
            "dataset_csv": str(classification_csv),
            **classification_descriptor_summary,
            "n_rows": int(len(classification_export)),
            "split_counts": {key: int(value) for key, value in classification_export["split"].value_counts().to_dict().items()},
            "target_counts": {str(key): int(value) for key, value in classification_export["target"].value_counts().to_dict().items()},
        },
        "regression": {
            "dataset_csv": str(regression_csv),
            **regression_descriptor_summary,
            "n_rows": int(len(regression_export)),
            "split_counts": {key: int(value) for key, value in regression_export["split"].value_counts().to_dict().items()},
        },
    }
    (output_dir / "chemprop_prepare_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return summary


def build_chemprop_train_command(
    dataset_csv: str | Path,
    descriptors_path: str | Path,
    output_dir: str | Path,
    task_type: str,
    target_column: str = "target",
    epochs: int = 30,
    extra_args: list[str] | None = None,
    require_installed: bool = True,
) -> list[str]:
    """Build the Chemprop CLI command for one training run."""
    dataset_csv = Path(dataset_csv)
    descriptors_path = Path(descriptors_path)
    output_dir = Path(output_dir)

    chemprop_bin = _resolve_chemprop_bin(require_installed=require_installed)

    command = [
        chemprop_bin,
        "train",
        "--data-path",
        str(dataset_csv),
        "--task-type",
        task_type,
        "--output-dir",
        str(output_dir),
        "--smiles-columns",
        "smiles",
        "--target-columns",
        target_column,
        "--splits-column",
        "split",
        "--descriptors-path",
        str(descriptors_path),
        "--epochs",
        str(epochs),
        "--num-workers",
        "0",
    ]
    if extra_args:
        command.extend(extra_args)
    return command


def train_step1_chemprop(
    dataset_csv: str | Path,
    descriptors_path: str | Path,
    output_dir: str | Path,
    task_type: str,
    target_column: str = "target",
    epochs: int = 30,
    extra_args: list[str] | None = None,
    dry_run: bool = False,
) -> dict[str, object]:
    """Launch a Chemprop training job and optionally summarize the resulting run."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    command = build_chemprop_train_command(
        dataset_csv=dataset_csv,
        descriptors_path=descriptors_path,
        output_dir=output_dir,
        task_type=task_type,
        target_column=target_column,
        epochs=epochs,
        extra_args=extra_args,
        require_installed=not dry_run,
    )
    request_summary = {
        "command": command,
        "command_string": shlex.join(command),
        "dry_run": dry_run,
        "dataset_csv": str(dataset_csv),
        "descriptors_path": str(descriptors_path),
        "output_dir": str(output_dir),
        "task_type": task_type,
        "target_column": target_column,
        "epochs": epochs,
        "extra_args": [] if extra_args is None else list(extra_args),
    }
    (output_dir / "chemprop_train_request.json").write_text(
        json.dumps(request_summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    if dry_run:
        return request_summary

    subprocess.run(command, check=True)
    summary = dict(request_summary)
    try:
        summary["metrics_summary"] = summarize_step1_chemprop_run(
            output_dir=output_dir,
            task_type=task_type,
        )
    except (FileNotFoundError, ValueError):
        # Keep training usable even if a future Chemprop version changes its
        # prediction export layout; the dedicated summary script can still be
        # rerun manually after inspection.
        pass
    return summary


def summarize_step1_chemprop_run(
    output_dir: str | Path,
    task_type: str,
    threshold: float = 0.5,
) -> dict[str, object]:
    """Summarize Chemprop test predictions into standard classification/regression metrics."""
    output_dir = Path(output_dir)
    test_csv = output_dir / "test.csv"
    test_predictions_csv = output_dir / "model_0" / "test_predictions.csv"

    if not test_csv.exists():
        raise FileNotFoundError(f"Expected Chemprop test split at {test_csv}")
    if not test_predictions_csv.exists():
        raise FileNotFoundError(f"Expected Chemprop test predictions at {test_predictions_csv}")

    test_frame = pd.read_csv(test_csv, low_memory=False)
    prediction_frame = pd.read_csv(test_predictions_csv, low_memory=False)

    if len(test_frame) != len(prediction_frame):
        raise ValueError(
            "Chemprop test rows and prediction rows do not align: "
            f"{len(test_frame)} != {len(prediction_frame)}"
        )

    observed = pd.to_numeric(test_frame["target"], errors="coerce")
    predicted = pd.to_numeric(prediction_frame["target"], errors="coerce")
    if observed.isna().any() or predicted.isna().any():
        raise ValueError("Chemprop evaluation inputs contain non-numeric target values.")

    summary: dict[str, object] = {
        "output_dir": str(output_dir),
        "task_type": task_type,
        "n_test": int(len(test_frame)),
        "test_csv": str(test_csv),
        "test_predictions_csv": str(test_predictions_csv),
    }

    if task_type == "classification":
        y_true = observed.astype(int)
        y_score = predicted.astype(float)
        y_pred = (y_score >= threshold).astype(int)
        precision, recall, _ = precision_recall_curve(y_true, y_score)
        summary.update(
            {
                "threshold": float(threshold),
                "positive_rate": float(y_true.mean()),
                "roc_auc": float(roc_auc_score(y_true, y_score)),
                "pr_auc": float(auc(recall, precision)),
                "accuracy": float(accuracy_score(y_true, y_pred)),
                "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
                "f1": float(f1_score(y_true, y_pred)),
            }
        )
    elif task_type == "regression":
        y_true = observed.astype(float)
        y_pred = predicted.astype(float)
        summary.update(
            {
                "target_mean": float(y_true.mean()),
                "target_std": float(y_true.std(ddof=0)),
                "rmse": float(mean_squared_error(y_true, y_pred) ** 0.5),
                "mae": float(mean_absolute_error(y_true, y_pred)),
                "r2": float(r2_score(y_true, y_pred)),
                "spearman": float(y_true.corr(y_pred, method="spearman")),
            }
        )
    else:
        raise ValueError(f"Unsupported task_type: {task_type}")

    (output_dir / "metrics_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return summary
