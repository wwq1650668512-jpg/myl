from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import load

from .mechanism import Step2MechanismProjector


def _ensure_identifier_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Ensure inference inputs contain a stable pair identifier."""
    result = frame.copy()
    if "pair_id" not in result.columns:
        if {"prestwick_id", "nt_code"}.issubset(result.columns):
            result["pair_id"] = (
                result["prestwick_id"].fillna("unknown_drug").astype(str)
                + "::"
                + result["nt_code"].fillna("unknown_microbe").astype(str)
            )
        else:
            result["pair_id"] = [f"row_{index}" for index in range(len(result))]
    return result


def _ensure_columns(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """Add missing model feature columns as NaN placeholders."""
    result = frame.copy()
    for column in columns:
        if column not in result.columns:
            result[column] = np.nan
    return result


def _prepare_feature_frame(
    frame: pd.DataFrame,
    numeric_features: list[str],
    categorical_features: list[str],
) -> pd.DataFrame:
    """Cast inference features into the formats expected by the trained pipelines."""
    prepared = frame.copy()
    for column in numeric_features:
        if column in prepared.columns:
            prepared[column] = pd.to_numeric(prepared[column], errors="coerce")
    for column in categorical_features:
        if column in prepared.columns:
            prepared[column] = prepared[column].map(lambda value: np.nan if pd.isna(value) else str(value))
    return prepared


def _max_jaccard_similarity(query: np.ndarray, reference: np.ndarray) -> float | None:
    """Compute the maximum Jaccard similarity between one query fingerprint and references."""
    if query.size == 0 or reference.size == 0:
        return None
    intersection = np.logical_and(reference, query).sum(axis=1)
    union = np.logical_or(reference, query).sum(axis=1)
    valid = union > 0
    if not valid.any():
        return None
    similarity = intersection[valid] / union[valid]
    return float(similarity.max())


def _build_drug_similarity_map(
    frame: pd.DataFrame,
    drug_key_column: str,
    fingerprint_columns: list[str],
    reference_drug_table: pd.DataFrame,
) -> dict[str, float | None]:
    """Compute nearest-training-drug fingerprint similarity for each unique query drug."""
    if not fingerprint_columns or reference_drug_table.empty or drug_key_column not in frame.columns:
        return {}

    query_drugs = frame.loc[:, [drug_key_column] + fingerprint_columns].copy()
    query_drugs = query_drugs.drop_duplicates(subset=[drug_key_column]).reset_index(drop=True)
    query_drugs[fingerprint_columns] = query_drugs[fingerprint_columns].fillna(0)

    reference = reference_drug_table.loc[:, [drug_key_column] + fingerprint_columns].copy()
    reference = reference.drop_duplicates(subset=[drug_key_column]).reset_index(drop=True)
    reference[fingerprint_columns] = reference[fingerprint_columns].fillna(0)
    reference_matrix = reference.loc[:, fingerprint_columns].to_numpy(dtype=bool)

    similarity_map: dict[str, float | None] = {}
    for _, row in query_drugs.iterrows():
        drug_key = row[drug_key_column]
        if pd.isna(drug_key):
            continue
        similarity_map[str(drug_key)] = _max_jaccard_similarity(
            row.loc[fingerprint_columns].to_numpy(dtype=bool),
            reference_matrix,
        )
    return similarity_map


def predict_step2_baseline(
    input_table_path: str | Path,
    output_dir: str | Path,
    classifier_path: str | Path,
    regressor_path: str | Path,
    metrics_path: str | Path,
    applicability_reference_path: str | Path | None = None,
    mechanism_reference_path: str | Path | None = None,
    probability_threshold: float | None = None,
    similarity_threshold: float = 0.25,
) -> dict[str, object]:
    """Run the Step 2 baseline models and export predictions plus applicability metadata.

    Args:
        input_table_path: Candidate input table for inference.
        output_dir: Directory where prediction CSVs and summary JSON are written.
        classifier_path: Trained Step 2 classifier pipeline.
        regressor_path: Trained Step 2 regressor pipeline.
        metrics_path: Metrics JSON containing feature schemas and thresholds.
        applicability_reference_path: Optional reference for applicability calculations.
        mechanism_reference_path: Optional mechanism reference for projection annotations.
        probability_threshold: Optional override for the metabolized probability cutoff.
        similarity_threshold: Minimum fingerprint similarity used in applicability checks.

    Returns:
        A summary dictionary with counts and generated artifact paths.
    """
    input_table_path = Path(input_table_path)
    output_dir = Path(output_dir)
    classifier_path = Path(classifier_path)
    regressor_path = Path(regressor_path)
    metrics_path = Path(metrics_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(input_table_path, low_memory=False)
    raw = _ensure_identifier_columns(raw)

    classifier = load(classifier_path)
    regressor = load(regressor_path)
    metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    if applicability_reference_path is not None:
        applicability_reference = load(applicability_reference_path)
    else:
        applicability_reference = {}
    mechanism_projector = Step2MechanismProjector.from_joblib(mechanism_reference_path)

    numeric_features = [str(value) for value in metrics.get("numeric_features", [])]
    categorical_features = [str(value) for value in metrics.get("categorical_features", [])]
    threshold = metrics.get("probability_threshold", 0.5) if probability_threshold is None else probability_threshold
    feature_columns = numeric_features + categorical_features

    inference = _ensure_columns(raw, feature_columns)
    inference = _prepare_feature_frame(inference, numeric_features=numeric_features, categorical_features=categorical_features)

    predictions = raw.copy()
    if hasattr(classifier.named_steps["model"], "predict_proba"):
        proba_matrix = classifier.predict_proba(inference.loc[:, feature_columns])
        class_labels = classifier.named_steps["model"].classes_.tolist()
        if "metabolized" in class_labels:
            predictions["predicted_metabolized_probability"] = proba_matrix[:, class_labels.index("metabolized")]
        else:
            predictions["predicted_metabolized_probability"] = np.nan
    else:
        class_labels = []
        predictions["predicted_metabolized_probability"] = np.nan
    predictions["predicted_metabolism_label"] = np.where(
        predictions["predicted_metabolized_probability"] >= threshold,
        "metabolized",
        "not_metabolized",
    )
    predictions["predicted_parent_depletion_fraction"] = regressor.predict(inference.loc[:, feature_columns])
    predictions["predicted_parent_depletion_magnitude"] = predictions["predicted_parent_depletion_fraction"].abs()

    drug_key_column = str(applicability_reference.get("drug_key_column", "prestwick_id"))
    fingerprint_columns = [str(value) for value in applicability_reference.get("fingerprint_columns", [])]
    reference_drug_table = applicability_reference.get("drug_reference", pd.DataFrame())
    if isinstance(reference_drug_table, pd.DataFrame):
        drug_similarity_map = _build_drug_similarity_map(
            frame=raw,
            drug_key_column=drug_key_column,
            fingerprint_columns=fingerprint_columns,
            reference_drug_table=reference_drug_table,
        )
    else:
        drug_similarity_map = {}

    predictions["drug_max_fingerprint_jaccard"] = predictions.get(
        drug_key_column,
        pd.Series(np.nan, index=predictions.index),
    ).map(lambda value: drug_similarity_map.get(str(value)) if pd.notna(value) else np.nan)
    seen_scaffolds = {str(value) for value in applicability_reference.get("seen_scaffolds", [])}
    predictions["scaffold_seen_in_training"] = predictions.get(
        "murcko_scaffold",
        pd.Series(np.nan, index=predictions.index),
    ).map(lambda value: False if pd.isna(value) else str(value) in seen_scaffolds)

    seen_genera = {str(value) for value in applicability_reference.get("seen_genera", [])}
    seen_phyla = {str(value) for value in applicability_reference.get("seen_phyla", [])}
    predictions["microbe_genus_seen_in_training"] = predictions.get(
        "genus",
        pd.Series(np.nan, index=predictions.index),
    ).map(lambda value: False if pd.isna(value) else str(value) in seen_genera)
    predictions["microbe_phylum_seen_in_training"] = predictions.get(
        "phylum",
        pd.Series(np.nan, index=predictions.index),
    ).map(lambda value: False if pd.isna(value) else str(value) in seen_phyla)
    predictions["applicability_flag"] = (
        predictions["microbe_genus_seen_in_training"] | predictions["microbe_phylum_seen_in_training"]
    ) & (
        predictions["scaffold_seen_in_training"]
        | (predictions["drug_max_fingerprint_jaccard"].fillna(-1) >= similarity_threshold)
    )
    predictions = mechanism_projector.annotate_frame(
        predictions,
        predicted_probability_column="predicted_metabolized_probability",
        predicted_label_column="predicted_metabolism_label",
    )

    predictions_path = output_dir / "predictions.csv"
    predictions.to_csv(predictions_path, index=False)

    slim_columns = [
        "pair_id",
        "prestwick_id",
        "nt_code",
        "smiles",
        "predicted_metabolized_probability",
        "predicted_metabolism_label",
        "predicted_parent_depletion_fraction",
        "predicted_parent_depletion_magnitude",
        "drug_max_fingerprint_jaccard",
        "scaffold_seen_in_training",
        "microbe_genus_seen_in_training",
        "microbe_phylum_seen_in_training",
        "applicability_flag",
        "predicted_mechanism_projection_flag",
        "predicted_reaction_class",
        "predicted_reaction_confidence",
        "predicted_candidate_product_count",
        "predicted_evidence_gene_count",
    ]
    existing_slim_columns = [column for column in slim_columns if column in predictions.columns]
    predictions_slim_path = output_dir / "predictions_slim.csv"
    predictions.loc[:, existing_slim_columns].to_csv(predictions_slim_path, index=False)

    summary = {
        "input_table_path": str(input_table_path),
        "output_dir": str(output_dir),
        "classifier_path": str(classifier_path),
        "regressor_path": str(regressor_path),
        "metrics_path": str(metrics_path),
        "applicability_reference_path": None
        if applicability_reference_path is None
        else str(applicability_reference_path),
        "mechanism_reference_path": None if mechanism_reference_path is None else str(mechanism_reference_path),
        "probability_threshold": float(threshold),
        "similarity_threshold": float(similarity_threshold),
        "n_rows": int(len(predictions)),
        "predicted_metabolism_label_counts": {
            str(key): int(value)
            for key, value in predictions["predicted_metabolism_label"].value_counts(dropna=False).to_dict().items()
        },
        "n_applicable_rows": int(predictions["applicability_flag"].sum()),
        "n_mechanism_projected_rows": int(predictions["predicted_mechanism_projection_flag"].fillna(False).sum()),
        "predictions_path": str(predictions_path),
        "predictions_slim_path": str(predictions_slim_path),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    return summary
