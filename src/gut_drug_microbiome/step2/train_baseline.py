from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import dump
from scipy.stats import spearmanr
from sklearn.base import clone
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.impute import SimpleImputer
from sklearn.metrics import accuracy_score
from sklearn.metrics import average_precision_score
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import f1_score
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupShuffleSplit
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

from .mechanism import build_step2_mechanism_reference


NUMERIC_PREFIXES = ("morgan_fp_",)
CATEGORICAL_PREFIXES = ("murcko_scaffold",)

DEFAULT_NUMERIC_FEATURES = [
    "molecular_weight",
    "xlogp",
    "tpsa",
    "complexity",
    "volume3d",
    "estimated_colon_concentration_um",
    "smiles_length",
    "smiles_uppercase_count",
    "smiles_ring_index_count",
    "smiles_branch_count",
    "smiles_double_bond_count",
    "smiles_halogen_count",
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

DEFAULT_CATEGORICAL_FEATURES = [
    "therapeutic_class",
    "therapeutic_effect",
    "therapeutic_indication",
    "atc_primary_l1",
    "atc_primary_l3",
    "atc_primary_l4",
    "human_use",
    "veterinary",
    "species_label",
    "species_name",
    "phylum",
    "phylum_or_description",
    "genus",
    "species_epithet",
    "rdkit_formula",
]


def _available_columns(frame: pd.DataFrame, desired: list[str]) -> list[str]:
    """Return the subset of requested feature columns that exists in the input table."""
    return [column for column in desired if column in frame.columns]


def _prefixed_columns(frame: pd.DataFrame, prefixes: tuple[str, ...]) -> list[str]:
    """Collect all columns whose names begin with one of the requested prefixes."""
    columns: list[str] = []
    for column in frame.columns:
        if any(column.startswith(prefix) for prefix in prefixes):
            columns.append(column)
    return sorted(columns)


def _prepare_feature_frame(
    frame: pd.DataFrame,
    numeric_features: list[str],
    categorical_features: list[str],
) -> pd.DataFrame:
    """Cast training/inference features into numeric or string form."""
    prepared = frame.copy()
    for column in numeric_features:
        if column in prepared.columns:
            prepared[column] = pd.to_numeric(prepared[column], errors="coerce")
    for column in categorical_features:
        if column in prepared.columns:
            prepared[column] = prepared[column].map(lambda value: np.nan if pd.isna(value) else str(value))
    return prepared


def _build_preprocessor(frame: pd.DataFrame) -> tuple[ColumnTransformer, list[str], list[str]]:
    """Construct the sklearn preprocessing pipeline and selected Step 2 feature lists."""
    numeric_features = _available_columns(frame, DEFAULT_NUMERIC_FEATURES) + _prefixed_columns(frame, NUMERIC_PREFIXES)
    categorical_features = _available_columns(frame, DEFAULT_CATEGORICAL_FEATURES) + _prefixed_columns(
        frame, CATEGORICAL_PREFIXES
    )

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
    return preprocessor, numeric_features, categorical_features


def _pick_group_column(frame: pd.DataFrame, choices: list[str]) -> pd.Series:
    """Choose the first available grouping column, or fall back to row index strings."""
    for choice in choices:
        if choice in frame.columns:
            return frame[choice]
    return frame.index.astype(str)


def _make_split(
    frame: pd.DataFrame,
    split_mode: str,
    target: pd.Series,
    random_state: int,
    test_size: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Create a train/test split, optionally grouped by drug, scaffold, or microbe."""
    index = np.arange(len(frame))
    if split_mode == "random":
        train_idx, test_idx = train_test_split(
            index,
            test_size=test_size,
            random_state=random_state,
            stratify=target,
        )
        return train_idx, test_idx

    if split_mode == "drug":
        groups = _pick_group_column(frame, ["prestwick_id", "drug_id", "drug_name"])
    elif split_mode == "scaffold":
        groups = frame["murcko_scaffold"].astype(object) if "murcko_scaffold" in frame.columns else pd.Series(np.nan, index=frame.index)
        fallback = _pick_group_column(frame, ["prestwick_id", "drug_id", "drug_name"])
        groups = groups.where(groups.notna(), fallback.astype(str))
    elif split_mode == "microbe":
        groups = _pick_group_column(frame, ["nt_code", "microbe_id", "microbe_name"])
    else:
        raise ValueError(f"Unsupported split_mode: {split_mode}")

    splitter = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=random_state)
    train_idx, test_idx = next(splitter.split(frame, y=target, groups=groups))
    return train_idx, test_idx


def _safe_binary_metric(
    metric_fn,
    y_true: pd.Series,
    y_score: np.ndarray,
) -> float | None:
    """Compute a binary metric only when both classes are present."""
    if y_true.nunique() < 2:
        return None
    return float(metric_fn(y_true, y_score))


def _build_applicability_reference(
    frame: pd.DataFrame,
    fingerprint_columns: list[str],
) -> dict[str, object]:
    """Package training-domain metadata used later for applicability scoring."""
    drug_key_column = "prestwick_id" if "prestwick_id" in frame.columns else "drug_id"
    microbe_key_column = "nt_code" if "nt_code" in frame.columns else "microbe_id"

    reference = {
        "drug_key_column": drug_key_column,
        "microbe_key_column": microbe_key_column,
        "fingerprint_columns": fingerprint_columns,
        "seen_scaffolds": sorted(
            {
                str(value)
                for value in frame.get("murcko_scaffold", pd.Series(dtype=object)).dropna().astype(str).tolist()
                if value.strip()
            }
        ),
        "seen_genera": sorted(
            {
                str(value)
                for value in frame.get("genus", pd.Series(dtype=object)).dropna().astype(str).tolist()
                if value.strip()
            }
        ),
        "seen_phyla": sorted(
            {
                str(value)
                for value in frame.get("phylum", pd.Series(dtype=object)).dropna().astype(str).tolist()
                if value.strip()
            }
        ),
    }

    drug_reference_columns = [drug_key_column]
    for optional in ["drug_name", "murcko_scaffold"]:
        if optional in frame.columns and optional not in drug_reference_columns:
            drug_reference_columns.append(optional)
    drug_reference_columns.extend(fingerprint_columns)

    drug_reference = frame.loc[:, [column for column in drug_reference_columns if column in frame.columns]].copy()
    drug_reference = drug_reference.drop_duplicates(subset=[drug_key_column]).reset_index(drop=True)
    reference["drug_reference"] = drug_reference
    return reference


def train_step2_baseline(
    modeling_table_path: str | Path,
    output_dir: str | Path,
    split_mode: str = "scaffold",
    random_state: int = 42,
    test_size: float = 0.2,
    n_estimators: int = 300,
    probability_threshold: float = 0.3,
    verbose: int = 1,
) -> dict[str, object]:
    """Train Step 2 ExtraTrees models and save metrics, predictions, and references.

    Args:
        modeling_table_path: Path to the normalized Step 2 modeling table.
        output_dir: Directory where model artifacts and summaries are written.
        split_mode: Evaluation split strategy such as random, drug, scaffold, or microbe.
        random_state: Random seed for reproducibility.
        test_size: Fraction of rows held out for evaluation.
        n_estimators: Tree count shared by classifier and regressor.
        probability_threshold: Probability cutoff used when producing predicted labels.
        verbose: Whether to print coarse training progress logs.

    Returns:
        A summary dictionary with metrics, feature lists, and artifact paths.
    """
    def _log(message: str) -> None:
        if verbose:
            print(message, flush=True)

    modeling_table_path = Path(modeling_table_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    _log(f"[step2] loading modeling table from {modeling_table_path}")
    raw_frame = pd.read_csv(modeling_table_path, low_memory=False)

    label_column = "metabolism_label"
    regression_column = "parent_depletion_fraction"
    valid_labels = {"metabolized", "not_metabolized"}
    frame = raw_frame[raw_frame[label_column].isin(valid_labels)].copy()
    frame[regression_column] = pd.to_numeric(frame[regression_column], errors="coerce")
    frame = frame.dropna(subset=[regression_column]).reset_index(drop=True)
    if frame.empty:
        raise RuntimeError("No usable Step 2 labeled rows were found after filtering.")

    _log("[step2] preparing feature lists and preprocessing tables")
    preprocessor, numeric_features, categorical_features = _build_preprocessor(frame)
    frame = _prepare_feature_frame(frame, numeric_features=numeric_features, categorical_features=categorical_features)
    feature_columns = numeric_features + categorical_features

    classifier = Pipeline(
        steps=[
            ("preprocessor", clone(preprocessor)),
            (
                "model",
                ExtraTreesClassifier(
                    n_estimators=n_estimators,
                    random_state=random_state,
                    n_jobs=-1,
                    class_weight="balanced_subsample",
                    min_samples_leaf=2,
                ),
            ),
        ]
    )
    regressor = Pipeline(
        steps=[
            ("preprocessor", clone(preprocessor)),
            (
                "model",
                ExtraTreesRegressor(
                    n_estimators=n_estimators,
                    random_state=random_state,
                    n_jobs=-1,
                    min_samples_leaf=2,
                ),
            ),
        ]
    )

    train_idx, test_idx = _make_split(
        frame=frame,
        split_mode=split_mode,
        target=frame[label_column],
        random_state=random_state,
        test_size=test_size,
    )
    train_frame = frame.loc[train_idx].copy()
    test_frame = frame.loc[test_idx].copy()

    x_train = train_frame.loc[:, feature_columns]
    x_test = test_frame.loc[:, feature_columns]
    y_train_cls = train_frame[label_column]
    y_test_cls = test_frame[label_column]
    y_train_reg = train_frame[regression_column]
    y_test_reg = test_frame[regression_column]

    _log(
        "[step2] fitting classifier and regressor "
        f"(split={split_mode}, n_estimators={n_estimators}, train={len(train_frame)}, test={len(test_frame)})"
    )
    classifier.fit(x_train, y_train_cls)
    regressor.fit(x_train, y_train_reg)

    _log("[step2] generating predictions and metrics")
    if hasattr(classifier.named_steps["model"], "predict_proba"):
        proba_matrix = classifier.predict_proba(x_test)
        class_labels = classifier.named_steps["model"].classes_.tolist()
        if "metabolized" in class_labels:
            metabolized_probability = proba_matrix[:, class_labels.index("metabolized")]
        else:
            metabolized_probability = np.full(len(x_test), np.nan)
    else:
        class_labels = []
        metabolized_probability = np.full(len(x_test), np.nan)
    if class_labels and "metabolized" in class_labels:
        cls_pred = np.where(metabolized_probability >= probability_threshold, "metabolized", "not_metabolized")
    else:
        cls_pred = classifier.predict(x_test)
    cls_labels = sorted(set(y_test_cls.astype(str).tolist()) | set(pd.Series(cls_pred).astype(str).tolist()))
    reg_pred = regressor.predict(x_test)

    cls_metrics = {
        "label_counts_train": {str(key): int(value) for key, value in y_train_cls.value_counts().to_dict().items()},
        "label_counts_test": {str(key): int(value) for key, value in y_test_cls.value_counts().to_dict().items()},
        "accuracy": float(accuracy_score(y_test_cls, cls_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_test_cls, cls_pred)),
        "macro_f1": float(f1_score(y_test_cls, cls_pred, average="macro")),
        "classification_report": classification_report(y_test_cls, cls_pred, output_dict=True),
        "confusion_matrix": confusion_matrix(y_test_cls, cls_pred, labels=cls_labels).tolist(),
        "labels": cls_labels,
        "model_classes": class_labels,
        "roc_auc": _safe_binary_metric(
            roc_auc_score,
            y_test_cls.eq("metabolized"),
            metabolized_probability,
        ),
        "average_precision": _safe_binary_metric(
            average_precision_score,
            y_test_cls.eq("metabolized"),
            metabolized_probability,
        ),
    }

    spearman_value = spearmanr(y_test_reg, reg_pred).statistic
    reg_metrics = {
        "rmse": float(np.sqrt(mean_squared_error(y_test_reg, reg_pred))),
        "mae": float(mean_absolute_error(y_test_reg, reg_pred)),
        "r2": float(r2_score(y_test_reg, reg_pred)),
        "spearman": float(0.0 if np.isnan(spearman_value) else spearman_value),
    }

    predictions = test_frame.loc[
        :,
        [
            column
            for column in [
                "pair_id",
                "prestwick_id",
                "nt_code",
                "drug_name",
                "microbe_name",
                label_column,
                regression_column,
            ]
            if column in test_frame.columns
        ],
    ].copy()
    predictions["predicted_metabolized_probability"] = metabolized_probability
    predictions["predicted_metabolism_label"] = np.where(
        predictions["predicted_metabolized_probability"] >= probability_threshold,
        "metabolized",
        "not_metabolized",
    )
    predictions["predicted_parent_depletion_fraction"] = reg_pred
    predictions["predicted_parent_depletion_magnitude"] = np.abs(reg_pred)
    predictions.to_csv(output_dir / "predictions.csv", index=False)

    _log("[step2] fitting full-data deployment models")
    classifier_full = clone(classifier)
    regressor_full = clone(regressor)
    classifier_full.fit(frame.loc[:, feature_columns], frame[label_column])
    regressor_full.fit(frame.loc[:, feature_columns], frame[regression_column])

    dump(classifier, output_dir / "classifier.joblib")
    dump(regressor, output_dir / "regressor.joblib")
    dump(classifier_full, output_dir / "classifier_full.joblib")
    dump(regressor_full, output_dir / "regressor_full.joblib")

    fingerprint_columns = _prefixed_columns(frame, NUMERIC_PREFIXES)
    applicability_reference = _build_applicability_reference(frame, fingerprint_columns=fingerprint_columns)
    dump(applicability_reference, output_dir / "applicability_reference.joblib")
    build_step2_mechanism_reference(frame, output_path=output_dir / "mechanism_reference.joblib")

    summary = {
        "model": "ExtraTrees",
        "split_mode": split_mode,
        "random_state": random_state,
        "test_size": test_size,
        "n_estimators": n_estimators,
        "probability_threshold": float(probability_threshold),
        "n_rows": int(len(frame)),
        "n_train": int(len(train_idx)),
        "n_test": int(len(test_idx)),
        "numeric_features": numeric_features,
        "categorical_features": categorical_features,
        "fingerprint_columns": fingerprint_columns,
        "classification": cls_metrics,
        "regression": reg_metrics,
        "artifacts": {
            "classifier_path": str(output_dir / "classifier.joblib"),
            "regressor_path": str(output_dir / "regressor.joblib"),
            "classifier_full_path": str(output_dir / "classifier_full.joblib"),
            "regressor_full_path": str(output_dir / "regressor_full.joblib"),
            "applicability_reference_path": str(output_dir / "applicability_reference.joblib"),
            "mechanism_reference_path": str(output_dir / "mechanism_reference.joblib"),
            "predictions_path": str(output_dir / "predictions.csv"),
        },
    }
    (output_dir / "metrics.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    _log(f"[step2] done: metrics saved to {output_dir / 'metrics.json'}")
    return summary
