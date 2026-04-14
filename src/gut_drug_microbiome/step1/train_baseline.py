# train_baseline是一个用于训练Step 1的基线模型的模块，包含一个主要函数train_step1_baseline，该函数接受预处理后的建模表和可选的弱监督表，构建特征列表和预处理管道，进行数据分割，训练ExtraTrees分类器和回归器，并生成预测结果和评估指标。训练完成后，它会将模型、预测结果和评估指标保存到指定的输出目录中，以供后续分析和比较使用。
from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping

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
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import f1_score
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from sklearn.model_selection import GroupShuffleSplit
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder


NUMERIC_PREFIXES = ("morgan_fp_",)
CATEGORICAL_PREFIXES = ("murcko_scaffold",)
PROMOTE_STEP2_NUMERIC_FEATURES = [
    "predicted_metabolized_probability",
    "predicted_parent_depletion_fraction",
    "drug_max_fingerprint_jaccard",
    "predicted_reaction_confidence",
    "predicted_mechanism_support_score",
    "predicted_candidate_product_count",
    "predicted_evidence_gene_count",
    "predicted_enzyme_match_count",
    "predicted_enzyme_presence_score",
    "predicted_enzyme_support_score",
    "predicted_enzyme_step1_promote_support_score",
    "predicted_enzyme_step1_inhibit_risk_score",
    "applicability_flag",
    "predicted_mechanism_projection_flag",
    "scaffold_seen_in_training",
    "microbe_genus_seen_in_training",
    "microbe_phylum_seen_in_training",
]

DEFAULT_NUMERIC_FEATURES = [
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
    "smiles_length",
    "smiles_uppercase_count",
    "smiles_ring_index_count",
    "smiles_branch_count",
    "smiles_double_bond_count",
    "smiles_halogen_count",
    "starting_od_96_well_screen",
    "starting_od_384_well_screen",
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
    "rdkit_formula",
]


def _available_columns(frame: pd.DataFrame, desired: list[str]) -> list[str]:
    """Return the subset of desired feature columns that actually exist in a table."""
    return [column for column in desired if column in frame.columns]


def _numeric_columns_with_observations(frame: pd.DataFrame, columns: list[str]) -> list[str]:
    """Keep only numeric feature columns that have at least one observed value."""
    kept: list[str] = []
    for column in columns:
        values = pd.to_numeric(frame[column], errors="coerce")
        if values.notna().any():
            kept.append(column)
    return kept


def _prefixed_columns(frame: pd.DataFrame, prefixes: tuple[str, ...]) -> list[str]:
    """Collect all columns whose names start with any requested feature prefix."""
    columns: list[str] = []
    for column in frame.columns:
        if any(column.startswith(prefix) for prefix in prefixes):
            columns.append(column)
    return sorted(columns)


def _make_split(
    frame: pd.DataFrame,
    split_mode: str,
    target: pd.Series,
    random_state: int,
    test_size: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Build a train/test split, optionally grouped by drug, scaffold, or microbe."""
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
        groups = frame["prestwick_id"]
    elif split_mode == "scaffold":
        groups = frame["murcko_scaffold"].astype(object)
        fallback = frame["prestwick_id"] if "prestwick_id" in frame.columns else frame.index.astype(str)
        groups = groups.where(groups.notna(), fallback.astype(str))
    elif split_mode == "microbe":
        groups = frame["nt_code"]
    else:
        raise ValueError(f"Unsupported split_mode: {split_mode}")

    splitter = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=random_state)
    train_idx, test_idx = next(splitter.split(frame, y=target, groups=groups))
    return train_idx, test_idx


def _build_preprocessor(
    frame: pd.DataFrame,
    extra_numeric_features: list[str] | None = None,
    extra_categorical_features: list[str] | None = None,
) -> tuple[ColumnTransformer, list[str], list[str]]:
    """Construct the sklearn preprocessing pipeline and list the selected features."""
    numeric_candidates = DEFAULT_NUMERIC_FEATURES + ([] if extra_numeric_features is None else extra_numeric_features)
    categorical_candidates = DEFAULT_CATEGORICAL_FEATURES + (
        [] if extra_categorical_features is None else extra_categorical_features
    )
    numeric_features = _available_columns(frame, numeric_candidates) + _prefixed_columns(frame, NUMERIC_PREFIXES)
    categorical_features = _available_columns(frame, categorical_candidates) + _prefixed_columns(
        frame, CATEGORICAL_PREFIXES
    )
    numeric_features = list(dict.fromkeys(numeric_features))
    categorical_features = list(dict.fromkeys(categorical_features))
    numeric_features = _numeric_columns_with_observations(frame, numeric_features)

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


def _prepare_feature_frame(frame: pd.DataFrame, numeric_features: list[str], categorical_features: list[str]) -> pd.DataFrame:
    """Cast feature columns into numeric or string form before model inference/training."""
    prepared = frame.copy()
    for column in numeric_features:
        if column in prepared.columns:
            prepared[column] = pd.to_numeric(prepared[column], errors="coerce")
    for column in categorical_features:
        if column in prepared.columns:
            prepared[column] = prepared[column].map(lambda value: np.nan if pd.isna(value) else str(value))
    return prepared


def _normalize_source_weight_map(source_weight_map: Mapping[str, float] | None) -> dict[str, float]:
    """Convert optional source weights into a normalized string-to-float mapping."""
    if source_weight_map is None:
        return {}
    normalized: dict[str, float] = {}
    for key, value in source_weight_map.items():
        normalized[str(key)] = float(value)
    return normalized


def _classifier_sample_weights(
    frame: pd.DataFrame,
    source_weight_map: Mapping[str, float] | None,
    default_gold_weight: float,
    default_silver_weight: float,
) -> np.ndarray:
    """Compute classifier sample weights from dataset source and gold/silver tier."""
    source_weight_map = _normalize_source_weight_map(source_weight_map)
    source_series = frame["source_dataset"] if "source_dataset" in frame.columns else pd.Series("unknown", index=frame.index)
    tier_series = frame["label_tier"] if "label_tier" in frame.columns else pd.Series(np.nan, index=frame.index)

    weights: list[float] = []
    for source_value, tier_value in zip(source_series, tier_series):
        source_name = "unknown" if pd.isna(source_value) else str(source_value)
        tier_name = "" if pd.isna(tier_value) else str(tier_value).lower()
        if source_name in source_weight_map:
            weights.append(float(source_weight_map[source_name]))
        elif tier_name == "gold":
            weights.append(float(default_gold_weight))
        else:
            weights.append(float(default_silver_weight))
    return np.asarray(weights, dtype=float)


def _merge_promote_feature_table(frame: pd.DataFrame, promote_feature_table_path: str | Path | None) -> pd.DataFrame:
    """Left-join optional Step 2-derived promote features onto a Step 1 table."""
    if promote_feature_table_path is None:
        return frame
    promote_feature_table_path = Path(promote_feature_table_path)
    if not promote_feature_table_path.exists():
        raise FileNotFoundError(f"Promote feature table not found: {promote_feature_table_path}")

    external = pd.read_csv(promote_feature_table_path, low_memory=False)
    join_columns = [column for column in ["pair_id", "prestwick_id", "nt_code"] if column in frame.columns and column in external.columns]
    if not join_columns:
        raise ValueError("Promote feature table must share at least one of pair_id/prestwick_id/nt_code")

    feature_columns = [column for column in PROMOTE_STEP2_NUMERIC_FEATURES if column in external.columns]
    if not feature_columns:
        return frame

    external = external.loc[:, join_columns + feature_columns].drop_duplicates(subset=join_columns).copy()
    merged = frame.merge(external, on=join_columns, how="left", suffixes=("", "_promote_aux"))
    for column in feature_columns:
        aux_column = f"{column}_promote_aux"
        if aux_column in merged.columns:
            if column in frame.columns:
                merged[column] = merged[column].combine_first(merged[aux_column])
            else:
                merged[column] = merged[aux_column]
            merged = merged.drop(columns=[aux_column])
    return merged


def train_step1_baseline(
    modeling_table_path: str | Path,
    output_dir: str | Path,
    silver_table_path: str | Path | None = None,
    promote_feature_table_path: str | Path | None = None,
    split_mode: str = "drug",
    random_state: int = 42,
    test_size: float = 0.2,
    n_estimators: int = 300,
    source_weight_map: Mapping[str, float] | None = None,
    default_gold_weight: float = 1.0,
    default_silver_weight: float = 1.0,
    enable_promote_head: bool = False,
    verbose: int = 1,
) -> dict:
    """Train the Step 1 ExtraTrees classifier and regressor and save their artifacts.

    Args:
        modeling_table_path: Path to the normalized Step 1 modeling table.
        output_dir: Directory where models, predictions, and metrics are written.
        silver_table_path: Optional weak-supervision table appended to classifier training.
        promote_feature_table_path: Optional Step 2-derived feature table used only by the promote auxiliary head.
        split_mode: Evaluation split strategy such as random, drug, scaffold, or microbe.
        random_state: Random seed for reproducibility.
        test_size: Fraction of gold rows held out for evaluation.
        n_estimators: Tree count shared by classifier and regressor.
        source_weight_map: Optional source-specific classifier weights.
        default_gold_weight: Fallback classifier weight for gold rows.
        default_silver_weight: Fallback classifier weight for silver rows.
        enable_promote_head: Whether to train an auxiliary promote-vs-not-promote classifier.
        verbose: Whether to print coarse progress logs.

    Returns:
        A summary dictionary with feature lists, metrics, and artifact metadata.
    """
    def _log(message: str) -> None:
        if verbose:
            print(message, flush=True)

    modeling_table_path = Path(modeling_table_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    _log(f"[step1] loading gold table from {modeling_table_path}")
    gold_frame = pd.read_csv(modeling_table_path, low_memory=False)
    gold_frame = gold_frame.dropna(subset=["effect_label", "binary_effect_label", "effect_score"]).reset_index(drop=True)
    gold_frame = _merge_promote_feature_table(gold_frame, promote_feature_table_path)

    silver_frame = None
    if silver_table_path is not None:
        silver_table_path = Path(silver_table_path)
        if silver_table_path.exists():
            _log(f"[step1] loading silver table from {silver_table_path}")
            silver_frame = pd.read_csv(silver_table_path, low_memory=False)
            silver_frame = silver_frame.dropna(subset=["effect_label"]).reset_index(drop=True)
            if "binary_effect_label" not in silver_frame.columns:
                silver_frame["binary_effect_label"] = silver_frame["effect_label"]
            if "effect_score" not in silver_frame.columns:
                silver_frame["effect_score"] = np.nan
            silver_frame = _merge_promote_feature_table(silver_frame, promote_feature_table_path)

    classification_target = "effect_label"
    target_counts = gold_frame[classification_target].value_counts()
    if target_counts.shape[0] < 2:
        classification_target = "binary_effect_label"
        target_counts = gold_frame[classification_target].value_counts()
    if target_counts.shape[0] < 2:
        raise RuntimeError("Not enough label diversity to train the classifier.")

    _log("[step1] preparing feature lists and preprocessing tables")
    feature_frame = gold_frame if silver_frame is None else pd.concat([gold_frame, silver_frame], ignore_index=True, sort=False)
    preprocessor, numeric_features, categorical_features = _build_preprocessor(feature_frame)
    gold_frame = _prepare_feature_frame(gold_frame, numeric_features, categorical_features)
    if silver_frame is not None:
        silver_frame = _prepare_feature_frame(silver_frame, numeric_features, categorical_features)

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
        frame=gold_frame,
        split_mode=split_mode,
        target=gold_frame[classification_target],
        random_state=random_state,
        test_size=test_size,
    )

    feature_columns = numeric_features + categorical_features
    gold_train = gold_frame.loc[train_idx].copy()
    gold_test = gold_frame.loc[test_idx].copy()
    if silver_frame is not None:
        classifier_train = pd.concat([gold_train, silver_frame], ignore_index=True, sort=False)
    else:
        classifier_train = gold_train

    x_train_cls = classifier_train.loc[:, feature_columns]
    y_train_cls = classifier_train[classification_target]
    x_test = gold_test.loc[:, feature_columns]
    y_test_cls = gold_test[classification_target]
    x_train_reg = gold_train.loc[:, feature_columns]
    y_train_reg = gold_train["effect_score"]
    y_test_reg = gold_test["effect_score"]
    classifier_sample_weights = _classifier_sample_weights(
        classifier_train,
        source_weight_map=source_weight_map,
        default_gold_weight=default_gold_weight,
        default_silver_weight=default_silver_weight,
    )

    _log(
        "[step1] fitting classifier "
        f"(split={split_mode}, n_estimators={n_estimators}, gold_train={len(gold_train)}, "
        f"silver_train={0 if silver_frame is None else len(silver_frame)}, gold_test={len(gold_test)})"
    )
    classifier.fit(x_train_cls, y_train_cls, model__sample_weight=classifier_sample_weights)
    _log("[step1] fitting regressor")
    regressor.fit(x_train_reg, y_train_reg)

    promote_classifier = None
    promote_metrics: dict[str, object] | None = None
    if enable_promote_head:
        promote_preprocessor, promote_numeric_features, promote_categorical_features = _build_preprocessor(
            feature_frame,
            extra_numeric_features=PROMOTE_STEP2_NUMERIC_FEATURES,
        )
        promote_feature_columns = promote_numeric_features + promote_categorical_features
        classifier_train_promote = _prepare_feature_frame(
            classifier_train,
            numeric_features=promote_numeric_features,
            categorical_features=promote_categorical_features,
        )
        gold_test_promote = _prepare_feature_frame(
            gold_test,
            numeric_features=promote_numeric_features,
            categorical_features=promote_categorical_features,
        )
        y_train_promote = np.where(classifier_train_promote["effect_label"].eq("promote"), "promote", "not_promote")
        y_test_promote = np.where(gold_test_promote["effect_label"].eq("promote"), "promote", "not_promote")
        if pd.Series(y_train_promote).nunique() >= 2:
            _log("[step1] fitting promote auxiliary classifier")
            promote_classifier = Pipeline(
                steps=[
                    ("preprocessor", clone(promote_preprocessor)),
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
            promote_classifier.fit(
                classifier_train_promote.loc[:, promote_feature_columns],
                y_train_promote,
                model__sample_weight=classifier_sample_weights,
            )
            promote_pred = promote_classifier.predict(gold_test_promote.loc[:, promote_feature_columns])
            if hasattr(promote_classifier.named_steps["model"], "predict_proba"):
                promote_proba_matrix = promote_classifier.predict_proba(gold_test_promote.loc[:, promote_feature_columns])
                promote_class_labels = promote_classifier.named_steps["model"].classes_.tolist()
                if "promote" in promote_class_labels:
                    promote_prob = promote_proba_matrix[:, promote_class_labels.index("promote")]
                else:
                    promote_prob = np.zeros(len(gold_test_promote), dtype=float)
            else:
                promote_prob = np.where(pd.Series(promote_pred).eq("promote"), 1.0, 0.0)
            predictions_promote_probability = promote_prob
            promote_metrics = {
                "trained": True,
                "numeric_features": promote_numeric_features,
                "categorical_features": promote_categorical_features,
                "label_counts_train": pd.Series(y_train_promote).value_counts().to_dict(),
                "label_counts_test": pd.Series(y_test_promote).value_counts().to_dict(),
                "accuracy": float(accuracy_score(y_test_promote, promote_pred)),
                "balanced_accuracy": float(balanced_accuracy_score(y_test_promote, promote_pred)),
                "macro_f1": float(f1_score(y_test_promote, promote_pred, average="macro")),
            }
        else:
            predictions_promote_probability = np.full(len(gold_test), np.nan)
            promote_metrics = {
                "trained": False,
                "reason": "not_enough_promote_label_diversity",
                "numeric_features": promote_numeric_features,
                "categorical_features": promote_categorical_features,
                "label_counts_train": pd.Series(y_train_promote).value_counts().to_dict(),
            }
    else:
        predictions_promote_probability = np.full(len(gold_test), np.nan)

    _log("[step1] generating predictions and metrics")
    cls_pred = classifier.predict(x_test)
    reg_pred = regressor.predict(x_test)
    cls_labels = sorted(set(y_test_cls.astype(str).tolist()) | set(pd.Series(cls_pred).astype(str).tolist()))

    cls_metrics = {
        "classification_target": classification_target,
        "label_counts_train": y_train_cls.value_counts().to_dict(),
        "label_counts_test": y_test_cls.value_counts().to_dict(),
        "accuracy": float(accuracy_score(y_test_cls, cls_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_test_cls, cls_pred)),
        "macro_f1": float(f1_score(y_test_cls, cls_pred, average="macro")),
        "classification_report": classification_report(y_test_cls, cls_pred, output_dict=True),
        "confusion_matrix": confusion_matrix(
            y_test_cls,
            cls_pred,
            labels=cls_labels,
        ).tolist(),
        "labels": cls_labels,
    }

    spearman_value = spearmanr(y_test_reg, reg_pred).statistic
    reg_metrics = {
        "rmse": float(np.sqrt(mean_squared_error(y_test_reg, reg_pred))),
        "mae": float(mean_absolute_error(y_test_reg, reg_pred)),
        "r2": float(r2_score(y_test_reg, reg_pred)),
        "spearman": float(0.0 if np.isnan(spearman_value) else spearman_value),
    }

    predictions = gold_test.loc[:, ["pair_id", "prestwick_id", "nt_code", "effect_label", "binary_effect_label", "effect_score"]].copy()
    predictions["predicted_effect_label"] = cls_pred
    predictions["predicted_effect_score"] = reg_pred
    if enable_promote_head:
        predictions["predicted_promote_probability"] = predictions_promote_probability
    predictions.to_csv(output_dir / "predictions.csv", index=False)

    _log("[step1] saving trained artifacts")
    dump(classifier, output_dir / "classifier.joblib")
    dump(regressor, output_dir / "regressor.joblib")
    if promote_classifier is not None:
        dump(promote_classifier, output_dir / "promote_classifier.joblib")

    summary = {
        "model": "ExtraTrees",
        "split_mode": split_mode,
        "random_state": random_state,
        "test_size": test_size,
        "n_estimators": n_estimators,
        "n_train_gold": int(len(train_idx)),
        "n_train_silver": int(0 if silver_frame is None else len(silver_frame)),
        "n_test_gold": int(len(test_idx)),
        "silver_table_path": None if silver_table_path is None else str(silver_table_path),
        "promote_feature_table_path": None if promote_feature_table_path is None else str(promote_feature_table_path),
        "source_weight_map": _normalize_source_weight_map(source_weight_map),
        "default_gold_weight": float(default_gold_weight),
        "default_silver_weight": float(default_silver_weight),
        "enable_promote_head": bool(enable_promote_head),
        "sample_weight_summary": {
            "min": float(classifier_sample_weights.min()),
            "max": float(classifier_sample_weights.max()),
            "mean": float(classifier_sample_weights.mean()),
        },
        "numeric_features": numeric_features,
        "categorical_features": categorical_features,
        "classification": cls_metrics,
        "regression": reg_metrics,
    }
    if promote_metrics is not None:
        summary["promote_auxiliary"] = promote_metrics

    (output_dir / "metrics.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    _log(f"[step1] done: metrics saved to {output_dir / 'metrics.json'}")
    return summary
