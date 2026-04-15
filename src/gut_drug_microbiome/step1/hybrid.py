# hybrid是一个集成了Chemprop分类器和基于RDKit特征的回归模型的预测模块，用于生成Step 1的混合效果标签和相关预测结果。它首先确保输入数据具有必要的标识符列，然后分别调用Chemprop分类器和回归模型进行预测，最后将这些预测结果融合成一个综合的效果标签，并输出详细的预测文件和摘要信息。
from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import load
from rdkit import Chem

from .chemprop_scaffold import _pick_smiles
from .chemprop_scaffold import _prepare_descriptor_frame
from .chemprop_scaffold import _resolve_chemprop_bin
from .train_baseline import _prepare_feature_frame

STEP1_PROFILE_EUBIOTIC = "eubiotic_modulator"
STEP1_PROFILE_HOST = "host_pathway_agent"
STEP1_PROFILE_ANTIFOLATE = "sulfonamide_antifolate"


def _pick_first_series(frame: pd.DataFrame, candidates: list[str], default_value: object = np.nan) -> pd.Series:
    """Return the first existing candidate column or a default-filled fallback series."""
    for candidate in candidates:
        if candidate in frame.columns:
            return frame[candidate]
    return pd.Series(default_value, index=frame.index)


def _normalize_lookup_value(value: object) -> str:
    """Normalize free-text lookup keys so literature references can match inference rows."""
    if pd.isna(value):
        return ""
    text = str(value).strip().lower()
    return "".join(character for character in text if character.isalnum())


def _normalize_lookup_series(values: pd.Series) -> pd.Series:
    """Vectorized wrapper for normalization of lookup keys."""
    return values.map(_normalize_lookup_value)


def _split_lookup_terms(value: object) -> list[str]:
    """Split a pipe/comma/semicolon separated field into normalized lookup terms."""
    if pd.isna(value):
        return []
    raw = str(value)
    terms: list[str] = []
    for piece in raw.replace(";", "|").replace(",", "|").split("|"):
        normalized = _normalize_lookup_value(piece)
        if normalized:
            terms.append(normalized)
    return list(dict.fromkeys(terms))


def _coerce_bool_series(values: pd.Series) -> pd.Series:
    """Map mixed boolean-like values onto numeric support scores in [0, 1]."""
    mapped = values.map(
        lambda value: np.nan
        if pd.isna(value)
        else (
            float(value)
            if isinstance(value, (bool, np.bool_, int, float, np.integer, np.floating))
            else {"true": 1.0, "false": 0.0, "yes": 1.0, "no": 0.0, "1": 1.0, "0": 0.0}.get(
                str(value).strip().lower(),
                np.nan,
            )
        )
    )
    return pd.to_numeric(mapped, errors="coerce")


def _is_valid_smiles_series(frame: pd.DataFrame, smiles_column: str) -> pd.Series:
    """Detect parseable SMILES rows so Chemprop inference can skip invalid strings."""
    if "rdkit_valid_smiles" in frame.columns:
        return pd.to_numeric(frame["rdkit_valid_smiles"], errors="coerce").fillna(0.0) > 0
    return frame[smiles_column].map(lambda value: isinstance(value, str) and Chem.MolFromSmiles(value) is not None)


def _sigmoid(values: pd.Series | np.ndarray, scale: float) -> np.ndarray:
    """Convert signed margins into soft probabilities with a stable sigmoid."""
    safe_scale = max(float(scale), 1e-6)
    numeric = np.asarray(values, dtype=float) / safe_scale
    numeric = np.clip(numeric, -30.0, 30.0)
    return 1.0 / (1.0 + np.exp(-numeric))


def _predict_promote_probability_from_auxiliary_head(
    frame: pd.DataFrame,
    promote_classifier_path: str | Path | None,
    promote_metrics_path: str | Path | None,
) -> pd.Series:
    """Run an optional promote auxiliary head on a frame that already contains Step 2 features."""
    if promote_classifier_path is None or promote_metrics_path is None:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    promote_classifier_path = Path(promote_classifier_path)
    promote_metrics_path = Path(promote_metrics_path)
    if not promote_classifier_path.exists() or not promote_metrics_path.exists():
        return pd.Series(np.nan, index=frame.index, dtype=float)

    metrics = json.loads(promote_metrics_path.read_text(encoding="utf-8"))
    promote_summary = metrics.get("promote_auxiliary", {})
    if not isinstance(promote_summary, dict) or not promote_summary.get("trained", False):
        return pd.Series(np.nan, index=frame.index, dtype=float)

    numeric_features = [str(value) for value in promote_summary.get("numeric_features", [])]
    categorical_features = [str(value) for value in promote_summary.get("categorical_features", [])]
    feature_columns = numeric_features + categorical_features
    if not feature_columns:
        return pd.Series(np.nan, index=frame.index, dtype=float)

    classifier = load(promote_classifier_path)
    inference = _ensure_columns(frame, feature_columns)
    inference = _prepare_feature_frame(
        inference,
        numeric_features=numeric_features,
        categorical_features=categorical_features,
    )
    if hasattr(classifier.named_steps["model"], "predict_proba"):
        proba_matrix = classifier.predict_proba(inference.loc[:, feature_columns])
        class_labels = classifier.named_steps["model"].classes_.tolist()
        if "promote" in class_labels:
            return pd.Series(proba_matrix[:, class_labels.index("promote")], index=frame.index, dtype=float)
    return pd.Series(np.nan, index=frame.index, dtype=float)


def _load_cross_feeding_reference(path: str | Path | None) -> pd.DataFrame:
    """Load a curated cross-feeding reference table if one is available."""
    expected_columns = [
        "reference_id",
        "compound_name_normalized",
        "compound_aliases",
        "compound_family",
        "match_keywords",
        "consumer_microbe_label",
        "producer_microbe_label",
        "evidence_type",
        "source_pmid",
        "source_title",
        "evidence_level",
        "notes",
    ]
    if path is None:
        return pd.DataFrame(columns=expected_columns)
    path = Path(path)
    if not path.exists():
        return pd.DataFrame(columns=expected_columns)

    reference = pd.read_csv(path, low_memory=False)
    for column in expected_columns:
        if column not in reference.columns:
            reference[column] = np.nan
    reference["compound_lookup_key"] = _normalize_lookup_series(reference["compound_name_normalized"])
    reference["consumer_lookup_key"] = _normalize_lookup_series(reference["consumer_microbe_label"])
    reference["compound_alias_keys"] = reference["compound_aliases"].map(_split_lookup_terms)
    reference["compound_family_key"] = _normalize_lookup_series(reference["compound_family"])
    reference["match_keyword_keys"] = reference["match_keywords"].map(_split_lookup_terms)
    reference = reference.drop_duplicates(subset=["compound_lookup_key", "consumer_lookup_key", "producer_microbe_label"], keep="first")
    return reference.loc[
        :,
        expected_columns
        + [
            "compound_lookup_key",
            "consumer_lookup_key",
            "compound_alias_keys",
            "compound_family_key",
            "match_keyword_keys",
        ],
    ].copy()


def _match_cross_feeding_reference(
    frame: pd.DataFrame,
    reference: pd.DataFrame,
) -> pd.DataFrame:
    """Match inference rows to curated cross-feeding edges using exact, alias, or keyword rules."""
    empty = pd.DataFrame(
        {
            "reference_id": pd.Series(np.nan, index=frame.index, dtype=object),
            "producer_microbe_label": pd.Series(np.nan, index=frame.index, dtype=object),
            "source_pmid": pd.Series(np.nan, index=frame.index, dtype=object),
            "evidence_level": pd.Series(np.nan, index=frame.index, dtype=object),
            "match_mode": pd.Series(np.nan, index=frame.index, dtype=object),
            "matched_term": pd.Series(np.nan, index=frame.index, dtype=object),
        }
    )
    if reference.empty or frame.empty:
        return empty

    compound_text = (
        _pick_first_series(frame, ["compound_name_normalized", "chemical_name"], default_value="").fillna("").astype(str)
        + " "
        + _pick_first_series(frame, ["compound_semantic_family"], default_value="").fillna("").astype(str)
        + " "
        + _pick_first_series(frame, ["compound_semantic_aliases"], default_value="").fillna("").astype(str)
        + " "
        + _pick_first_series(frame, ["compound_semantic_keywords"], default_value="").fillna("").astype(str)
        + " "
        + _pick_first_series(frame, ["therapeutic_class"], default_value="").fillna("").astype(str)
        + " "
        + _pick_first_series(frame, ["therapeutic_effect"], default_value="").fillna("").astype(str)
    )
    compound_text_key = _normalize_lookup_series(compound_text)
    consumer_lookup_key = _normalize_lookup_series(
        _pick_first_series(frame, ["species_label", "microbe_label", "species_name"], default_value="")
    )

    output = empty.copy()
    reference_records = reference.to_dict(orient="records")
    for index in frame.index:
        consumer_key = consumer_lookup_key.loc[index]
        if not consumer_key:
            continue
        row_compound_key = compound_text_key.loc[index]
        best_record: dict[str, object] | None = None
        best_priority = -1
        best_mode = ""
        best_term = ""

        for record in reference_records:
            if consumer_key != str(record.get("consumer_lookup_key", "")):
                continue
            exact_key = str(record.get("compound_lookup_key", ""))
            alias_keys = [str(term) for term in record.get("compound_alias_keys", [])]
            keyword_keys = [str(term) for term in record.get("match_keyword_keys", [])]
            family_key = str(record.get("compound_family_key", ""))

            priority = -1
            mode = ""
            matched_term = ""
            if exact_key and row_compound_key == exact_key:
                priority = 4
                mode = "exact_compound"
                matched_term = exact_key
            elif alias_keys and row_compound_key in alias_keys:
                priority = 3
                mode = "compound_alias"
                matched_term = row_compound_key
            else:
                for keyword in keyword_keys:
                    if keyword and keyword in row_compound_key:
                        priority = 2
                        mode = "keyword_family"
                        matched_term = keyword
                        break
                if priority < 0 and family_key and family_key in row_compound_key:
                    priority = 1
                    mode = "compound_family"
                    matched_term = family_key

            if priority > best_priority:
                best_priority = priority
                best_record = record
                best_mode = mode
                best_term = matched_term

        if best_record is None:
            continue
        output.at[index, "reference_id"] = best_record.get("reference_id")
        output.at[index, "producer_microbe_label"] = best_record.get("producer_microbe_label")
        output.at[index, "source_pmid"] = best_record.get("source_pmid")
        output.at[index, "evidence_level"] = best_record.get("evidence_level")
        output.at[index, "match_mode"] = best_mode
        output.at[index, "matched_term"] = best_term
    return output


def refine_step1_promote_with_step2(
    frame: pd.DataFrame,
    promote_classifier_path: str | Path | None = None,
    promote_metrics_path: str | Path | None = None,
    cross_feeding_reference_path: str | Path | None = None,
    inhibit_probability_threshold: float = 0.5,
    promote_score_threshold: float = 0.2,
    promote_probability_threshold: float = 0.35,
) -> pd.DataFrame:
    """Rescore Step 1 promote probability using downstream Step 2 metabolism evidence."""
    result = frame.copy()

    inhibit_probability = pd.to_numeric(
        _pick_first_series(result, ["step1_predicted_inhibit_probability", "predicted_inhibit_probability"]),
        errors="coerce",
    )
    effect_score = pd.to_numeric(
        _pick_first_series(result, ["step1_predicted_effect_score", "predicted_effect_score"]),
        errors="coerce",
    )
    base_label = _pick_first_series(
        result,
        ["step1_predicted_effect_label_hybrid", "predicted_effect_label_hybrid"],
        default_value="no_effect",
    ).astype(str)

    score_margin = effect_score.fillna(promote_score_threshold) - float(promote_score_threshold)
    promote_from_score = _sigmoid(score_margin, scale=0.05)
    anti_inhibit = 1.0 - inhibit_probability.fillna(0.5).clip(0.0, 1.0)
    heuristic_promote_probability = np.clip(promote_from_score * (0.5 + 0.5 * anti_inhibit), 0.0, 1.0)

    auxiliary_promote_probability = _predict_promote_probability_from_auxiliary_head(
        result,
        promote_classifier_path=promote_classifier_path,
        promote_metrics_path=promote_metrics_path,
    )
    base_promote_probability = pd.Series(heuristic_promote_probability, index=result.index, dtype=float)
    head_available = auxiliary_promote_probability.notna()
    base_promote_probability.loc[head_available] = (
        0.35 * base_promote_probability.loc[head_available]
        + 0.65 * auxiliary_promote_probability.loc[head_available]
    )

    cross_feeding_reference = _load_cross_feeding_reference(cross_feeding_reference_path)
    if cross_feeding_reference.empty:
        cross_feeding_match = pd.Series(False, index=result.index)
        cross_feeding_support_microbe = pd.Series(np.nan, index=result.index, dtype=object)
        cross_feeding_reference_pmid = pd.Series(np.nan, index=result.index, dtype=object)
        cross_feeding_evidence_level = pd.Series(np.nan, index=result.index, dtype=object)
        cross_feeding_match_mode = pd.Series(np.nan, index=result.index, dtype=object)
        cross_feeding_matched_term = pd.Series(np.nan, index=result.index, dtype=object)
    else:
        merged_reference = _match_cross_feeding_reference(result, cross_feeding_reference).reindex(result.index)
        cross_feeding_match = merged_reference["reference_id"].notna()
        cross_feeding_support_microbe = merged_reference["producer_microbe_label"]
        cross_feeding_reference_pmid = merged_reference["source_pmid"]
        cross_feeding_evidence_level = merged_reference["evidence_level"]
        cross_feeding_match_mode = merged_reference["match_mode"]
        cross_feeding_matched_term = merged_reference["matched_term"]

    metabolism_component = pd.to_numeric(
        _pick_first_series(result, ["predicted_metabolized_probability"]),
        errors="coerce",
    ).clip(0.0, 1.0)
    depletion_component = pd.to_numeric(
        _pick_first_series(result, ["predicted_parent_depletion_fraction"]),
        errors="coerce",
    ).clip(0.0, 1.0)
    mechanism_component = pd.concat(
        [
            _coerce_bool_series(_pick_first_series(result, ["predicted_mechanism_projection_flag"])),
            pd.to_numeric(_pick_first_series(result, ["predicted_reaction_confidence"]), errors="coerce"),
            pd.to_numeric(_pick_first_series(result, ["predicted_mechanism_support_score"]), errors="coerce"),
        ],
        axis=1,
    ).max(axis=1, skipna=True).clip(0.0, 1.0)
    enzyme_support_component = pd.concat(
        [
            _coerce_bool_series(_pick_first_series(result, ["predicted_enzyme_prior_flag"])),
            pd.to_numeric(_pick_first_series(result, ["predicted_enzyme_support_score"]), errors="coerce"),
            pd.to_numeric(_pick_first_series(result, ["predicted_enzyme_step1_promote_support_score"]), errors="coerce"),
        ],
        axis=1,
    ).max(axis=1, skipna=True).clip(0.0, 1.0)
    enzyme_risk_component = pd.to_numeric(
        _pick_first_series(result, ["predicted_enzyme_step1_inhibit_risk_score"]),
        errors="coerce",
    ).clip(0.0, 1.0)
    context_component = pd.concat(
        [
            pd.to_numeric(_pick_first_series(result, ["drug_max_fingerprint_jaccard"]), errors="coerce"),
            _coerce_bool_series(_pick_first_series(result, ["applicability_flag"])),
            _coerce_bool_series(_pick_first_series(result, ["scaffold_seen_in_training"])),
            _coerce_bool_series(_pick_first_series(result, ["microbe_genus_seen_in_training"])),
            _coerce_bool_series(_pick_first_series(result, ["microbe_phylum_seen_in_training"])),
        ],
        axis=1,
    ).max(axis=1, skipna=True).clip(0.0, 1.0)

    support_components = pd.DataFrame(
        {
            "metabolism": metabolism_component,
            "depletion": depletion_component,
            "mechanism": mechanism_component,
            "enzyme": enzyme_support_component,
            "cross_feeding": cross_feeding_match.astype(float),
            "context": context_component,
        },
        index=result.index,
    )
    weights = pd.Series(
        {"metabolism": 0.34, "depletion": 0.16, "mechanism": 0.14, "enzyme": 0.14, "cross_feeding": 0.14, "context": 0.08}
    )
    weighted_numerator = support_components.fillna(0.0).mul(weights, axis=1).sum(axis=1)
    weighted_denominator = support_components.notna().mul(weights, axis=1).sum(axis=1)
    support_score = weighted_numerator.div(weighted_denominator.where(weighted_denominator > 0, np.nan)).fillna(0.0)

    metabolizer_gate = metabolism_component.fillna(0.0)
    strong_metabolizer = metabolizer_gate >= 0.6
    weak_metabolizer = (metabolizer_gate >= 0.3) & ~strong_metabolizer
    mechanism_support = (
        pd.concat([depletion_component, mechanism_component, enzyme_support_component], axis=1).max(axis=1, skipna=True).fillna(0.0)
        >= 0.35
    )

    uplift_factor = pd.Series(0.85, index=result.index, dtype=float)
    uplift_factor.loc[weak_metabolizer] = 0.95 + 0.15 * support_score.loc[weak_metabolizer]
    uplift_factor.loc[strong_metabolizer] = 1.0 + 0.35 * support_score.loc[strong_metabolizer]
    uplift_factor.loc[strong_metabolizer & mechanism_support] = (
        1.05 + 0.45 * support_score.loc[strong_metabolizer & mechanism_support]
    )
    uplift_factor.loc[(metabolizer_gate < 0.3) & cross_feeding_match] = np.maximum(
        uplift_factor.loc[(metabolizer_gate < 0.3) & cross_feeding_match],
        1.02 + 0.18 * support_score.loc[(metabolizer_gate < 0.3) & cross_feeding_match],
    )
    enzyme_supported = enzyme_support_component.fillna(0.0) >= 0.35
    uplift_factor.loc[strong_metabolizer & enzyme_supported] = np.maximum(
        uplift_factor.loc[strong_metabolizer & enzyme_supported],
        1.02 + 0.24 * support_score.loc[strong_metabolizer & enzyme_supported],
    )

    refined_promote_probability = (
        base_promote_probability * uplift_factor * (1.0 - 0.12 * enzyme_risk_component.fillna(0.0))
    ).clip(0.0, 1.0)
    non_metabolizer_mask = metabolizer_gate < 0.3
    constrained_non_metabolizer_mask = non_metabolizer_mask & ~cross_feeding_match
    refined_promote_probability.loc[constrained_non_metabolizer_mask] = np.minimum(
        refined_promote_probability.loc[constrained_non_metabolizer_mask],
        base_promote_probability.loc[constrained_non_metabolizer_mask],
    )
    cross_feeding_promote_candidate = (
        cross_feeding_match
        & non_metabolizer_mask
        & (base_promote_probability >= float(promote_probability_threshold) * 0.8)
    )
    promote_mask = (inhibit_probability.fillna(0.0) < float(inhibit_probability_threshold)) & (
        ((effect_score.fillna(-np.inf) >= float(promote_score_threshold)) & (refined_promote_probability >= float(promote_probability_threshold)))
        | (refined_promote_probability >= 0.6)
        | cross_feeding_promote_candidate
    )
    refined_label = np.where(
        inhibit_probability.fillna(0.0) >= float(inhibit_probability_threshold),
        "inhibit",
        np.where(promote_mask, "promote", "no_effect"),
    )
    support_type = np.where(
        strong_metabolizer & mechanism_support,
        "self_metabolism_supported",
        np.where(
            cross_feeding_match,
            "cross_feeding_reference",
            np.where(
                strong_metabolizer & enzyme_supported,
                "enzyme_prior_supported",
                np.where(weak_metabolizer | strong_metabolizer, "self_metabolism_consistent", "weak_or_unspecified"),
            ),
        ),
    )
    evidence_type = np.where(
        cross_feeding_match,
        "cross_feeding_supported_promote",
        np.where(
            strong_metabolizer & mechanism_support,
            "self_metabolism_supported_promote",
            np.where(
                strong_metabolizer & enzyme_supported,
                "enzyme_prior_supported_promote",
                np.where(base_promote_probability >= float(promote_probability_threshold), "direct_or_host_like_promote", "unspecified_promote"),
            ),
        ),
    )

    result["predicted_promote_probability_base"] = base_promote_probability
    result["predicted_promote_probability_refined"] = refined_promote_probability
    result["predicted_promote_support_score"] = support_score
    result["predicted_promote_support_type"] = support_type
    result["predicted_promote_evidence_type"] = evidence_type
    result["predicted_cross_feeding_reference_flag"] = cross_feeding_match
    result["predicted_cross_feeding_support_microbe"] = cross_feeding_support_microbe
    result["predicted_cross_feeding_reference_pmid"] = cross_feeding_reference_pmid
    result["predicted_cross_feeding_evidence_level"] = cross_feeding_evidence_level
    result["predicted_cross_feeding_match_mode"] = cross_feeding_match_mode
    result["predicted_cross_feeding_matched_term"] = cross_feeding_matched_term
    result["predicted_effect_label_step2_refined"] = refined_label
    result["predicted_effect_label_step2_refined_changed"] = base_label.ne(result["predicted_effect_label_step2_refined"])
    return result

# _ensure_identifier_columns：确保输入数据框中存在一个稳定的pair_id列，如果没有，则根据prestwick_id和nt_code生成一个组合列，或者使用行索引生成唯一标识符，以便后续的预测结果能够正确地与输入数据对应。
def _ensure_identifier_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Guarantee that inference rows have a stable pair_id column."""
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

# _ensure_columns：检查数据框中是否存在指定的列，如果缺失则添加这些列并填充NaN值，以确保后续的模型预测步骤能够顺利进行而不会因为缺少必要的特征列而出错。
def _ensure_columns(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """Add any missing columns needed by a downstream model as NaN placeholders."""
    missing_columns = [column for column in columns if column not in frame.columns]
    if not missing_columns:
        return frame.copy()
    padding = pd.DataFrame({column: np.nan for column in missing_columns}, index=frame.index)
    return pd.concat([frame.copy(), padding], axis=1)

# _hybrid_effect_label：根据Chemprop的抑制概率、辅助促进概率和回归效果分数，结合预设的阈值，生成一个综合的效果标签。优先考虑抑制标签，如果满足抑制概率阈值则返回"inhibit"，否则如果满足促进概率或效果分数的阈值则返回"promote"，否则返回"no_effect"。
def _hybrid_effect_label(
    inhibit_probability: float | None,
    effect_score: float | None,
    inhibit_probability_threshold: float,
    promote_score_threshold: float,
) -> str:
    """Fuse classification probability and regression score into one hybrid label."""
    if inhibit_probability is not None and not np.isnan(inhibit_probability):
        if inhibit_probability >= inhibit_probability_threshold:
            return "inhibit"
    if effect_score is not None and not np.isnan(effect_score):
        if effect_score >= promote_score_threshold:
            return "promote"
    return "no_effect"


def _infer_step1_drug_profile_series(frame: pd.DataFrame) -> pd.Series:
    """Infer lightweight drug profile tags used by Step1 realism constraints."""
    name_text = _pick_first_series(frame, ["chemical_name"], default_value="").fillna("").astype(str)
    class_text = _pick_first_series(frame, ["therapeutic_class"], default_value="").fillna("").astype(str)
    effect_text = _pick_first_series(frame, ["therapeutic_effect"], default_value="").fillna("").astype(str)
    semantic_family = _pick_first_series(frame, ["compound_semantic_family"], default_value="").fillna("").astype(str)
    merged = (name_text + " " + class_text + " " + effect_text + " " + semantic_family).map(_normalize_lookup_value)

    profile = pd.Series("unknown", index=frame.index, dtype=object)
    eubiotic_mask = merged.str.contains("rifaximin", regex=False)
    host_mask = (
        merged.str.contains("lubiprostone", regex=False)
        | merged.str.contains("secretagogue", regex=False)
        | merged.str.contains("chloridechannelactivator", regex=False)
    )
    antifolate_mask = (
        merged.str.contains("sulfasalazine", regex=False)
        | merged.str.contains("sasp", regex=False)
        | merged.str.contains("sulfapyridine", regex=False)
        | merged.str.contains("sulfonamide", regex=False)
        | merged.str.contains("sulfamethoxazole", regex=False)
        | merged.str.contains("sulfadiazine", regex=False)
        | merged.str.contains("sulfisoxazole", regex=False)
        | merged.str.contains("cotrimoxazole", regex=False)
        | merged.str.contains("antifolate", regex=False)
    )
    profile.loc[eubiotic_mask] = STEP1_PROFILE_EUBIOTIC
    profile.loc[host_mask] = STEP1_PROFILE_HOST
    profile.loc[antifolate_mask] = STEP1_PROFILE_ANTIFOLATE
    return profile


def _core_butyrate_mask(frame: pd.DataFrame) -> pd.Series:
    """Identify core butyrate-producer taxa for Rifaximin constraint checks."""
    species = _pick_first_series(frame, ["species_label", "microbe_label", "species_name"], default_value="").fillna("").astype(str)
    genus = _pick_first_series(frame, ["genus"], default_value="").fillna("").astype(str)
    species_key = species.map(_normalize_lookup_value)
    genus_key = genus.map(_normalize_lookup_value)
    return (
        species_key.str.contains("faecalibacteriumprausnitzii", regex=False)
        | species_key.str.contains("eubacteriumrectale", regex=False)
        | genus_key.str.contains("roseburia", regex=False)
    )


def _antifolate_vulnerability_score(frame: pd.DataFrame) -> pd.Series:
    """Estimate microbe susceptibility for sulfonamide-like antifolate pressure."""
    species = _pick_first_series(frame, ["species_label", "microbe_label", "species_name"], default_value="").fillna("").astype(str)
    genus = _pick_first_series(frame, ["genus"], default_value="").fillna("").astype(str)
    medium = _pick_first_series(frame, ["medium_preference"], default_value="").fillna("").astype(str)

    species_key = species.map(_normalize_lookup_value)
    genus_key = genus.map(_normalize_lookup_value)
    medium_key = medium.map(_normalize_lookup_value)
    core_mask = _core_butyrate_mask(frame)

    strict_anaerobe_mask = (
        medium_key.str.contains("anaerob", regex=False)
        | genus_key.isin({"faecalibacterium", "roseburia", "eubacterium", "anaerostipes", "coprococcus", "butyrivibrio"})
    )
    likely_folate_denovo_mask = (
        core_mask
        | genus_key.isin({"faecalibacterium", "roseburia", "eubacterium", "anaerostipes", "coprococcus", "butyrivibrio"})
        | species_key.str.contains("faecalibacteriumprausnitzii", regex=False)
    )

    score = pd.Series(0.18, index=frame.index, dtype=float)
    score = score + core_mask.astype(float) * 0.44
    score = score + strict_anaerobe_mask.astype(float) * 0.20
    score = score + likely_folate_denovo_mask.astype(float) * 0.18
    return score.clip(0.0, 1.0)


def _apply_step1_drug_profile_constraints(
    predictions: pd.DataFrame,
    inhibit_probability_threshold: float,
    promote_score_threshold: float,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Apply profile-aware realism constraints before exporting final Step1 effect calls."""
    result = predictions.copy()
    inhibit_prob = pd.to_numeric(result["predicted_inhibit_probability"], errors="coerce")
    effect_score = pd.to_numeric(result["predicted_effect_score"], errors="coerce")
    profile = _infer_step1_drug_profile_series(result)
    result["step1_drug_profile"] = profile
    result["step1_constraint_applied"] = False
    result["step1_constraint_reason"] = "none"

    summary = {
        "eubiotic_core_rows_adjusted": 0,
        "host_rows_adjusted": 0,
        "host_high_evidence_rows_retained": 0,
        "antifolate_rows_adjusted": 0,
        "antifolate_sensitive_rows_boosted": 0,
    }

    # Constraint 1: eubiotic_modulator should not indiscriminately strongly inhibit core butyrate producers.
    core_mask = _core_butyrate_mask(result)
    eubiotic_mask = profile.eq(STEP1_PROFILE_EUBIOTIC)
    strong_inhibit_mask = (inhibit_prob.fillna(0.0) >= 0.65) | effect_score.fillna(0.0).lt(-0.06)
    adjust_eubiotic = eubiotic_mask & core_mask & strong_inhibit_mask
    if adjust_eubiotic.any():
        summary["eubiotic_core_rows_adjusted"] = int(adjust_eubiotic.sum())
        inhibit_prob.loc[adjust_eubiotic] = inhibit_prob.loc[adjust_eubiotic].clip(upper=min(0.45, inhibit_probability_threshold - 0.01))
        effect_score.loc[adjust_eubiotic] = effect_score.loc[adjust_eubiotic].clip(lower=-0.02)
        result.loc[adjust_eubiotic, "step1_constraint_applied"] = True
        result.loc[adjust_eubiotic, "step1_constraint_reason"] = "eubiotic_core_butyrate_inhibit_clipped"

    # Constraint 2: host_pathway_agent defaults to weak direct microbe effects; keep only high-confidence tails.
    host_mask = profile.eq(STEP1_PROFILE_HOST)
    if host_mask.any():
        host_prob = inhibit_prob.loc[host_mask].fillna(0.0)
        host_score = effect_score.loc[host_mask].fillna(0.0)
        high_evidence = (host_prob >= 0.88) | host_score.abs().ge(0.28)
        keep_high = host_mask.copy()
        keep_high.loc[host_mask] = high_evidence.values
        soften_mask = host_mask & ~keep_high

        # Global downscale for host-pathway agents.
        inhibit_prob.loc[host_mask] = (inhibit_prob.loc[host_mask].fillna(0.0) * 0.35).clip(0.0, 1.0)
        effect_score.loc[host_mask] = effect_score.loc[host_mask].fillna(0.0) * 0.30

        # Default to no_effect unless high-evidence tail is retained.
        if soften_mask.any():
            inhibit_prob.loc[soften_mask] = inhibit_prob.loc[soften_mask].clip(upper=max(0.0, inhibit_probability_threshold - 0.01))
            effect_score.loc[soften_mask] = effect_score.loc[soften_mask].clip(lower=-0.08, upper=0.08)
            result.loc[soften_mask, "step1_constraint_applied"] = True
            result.loc[soften_mask, "step1_constraint_reason"] = "host_pathway_global_soften_to_no_effect"
        if keep_high.any():
            result.loc[keep_high, "step1_constraint_applied"] = True
            result.loc[keep_high, "step1_constraint_reason"] = "host_pathway_high_evidence_retained"

        summary["host_rows_adjusted"] = int(host_mask.sum())
        summary["host_high_evidence_rows_retained"] = int(keep_high.sum())

    # Constraint 3: sulfonamide-antifolate agents should preferentially inhibit folate-vulnerable anaerobes.
    antifolate_mask = profile.eq(STEP1_PROFILE_ANTIFOLATE)
    if antifolate_mask.any():
        vulnerability_score = _antifolate_vulnerability_score(result)
        result["predicted_folate_vulnerability_score"] = vulnerability_score
        sensitive_mask = antifolate_mask & (vulnerability_score >= 0.55)
        background_mask = antifolate_mask & ~sensitive_mask

        if sensitive_mask.any():
            inhibit_prob.loc[sensitive_mask] = (
                inhibit_prob.loc[sensitive_mask].fillna(0.0) + 0.18 + 0.22 * vulnerability_score.loc[sensitive_mask]
            ).clip(0.0, 1.0)
            effect_score.loc[sensitive_mask] = (
                effect_score.loc[sensitive_mask].fillna(0.0) - (0.06 + 0.18 * vulnerability_score.loc[sensitive_mask])
            )
            result.loc[sensitive_mask, "step1_constraint_applied"] = True
            result.loc[sensitive_mask, "step1_constraint_reason"] = "antifolate_core_folate_vulnerability_boost"

        if background_mask.any():
            inhibit_prob.loc[background_mask] = (
                inhibit_prob.loc[background_mask].fillna(0.0) * (0.90 + 0.08 * vulnerability_score.loc[background_mask])
            ).clip(0.0, 1.0)

        summary["antifolate_rows_adjusted"] = int(antifolate_mask.sum())
        summary["antifolate_sensitive_rows_boosted"] = int(sensitive_mask.sum())

    result["predicted_inhibit_probability"] = inhibit_prob
    result["predicted_effect_score"] = effect_score
    result["predicted_binary_effect_label"] = np.where(
        result["predicted_inhibit_probability"].fillna(0.0) >= inhibit_probability_threshold,
        "inhibit",
        "no_effect",
    )
    result["predicted_effect_label_hybrid"] = [
        _hybrid_effect_label(
            inhibit_probability=float(inhibit_probability) if pd.notna(inhibit_probability) else None,
            effect_score=float(score) if pd.notna(score) else None,
            inhibit_probability_threshold=inhibit_probability_threshold,
            promote_score_threshold=promote_score_threshold,
        )
        for inhibit_probability, score in zip(result["predicted_inhibit_probability"], result["predicted_effect_score"])
    ]
    result["predicted_effect_magnitude"] = pd.to_numeric(result["predicted_effect_score"], errors="coerce").abs()
    return result, summary

# _predict_chemprop_inhibit_probability：使用Chemprop分类器对输入数据进行预测，返回每行的抑制概率。它首先准备输入数据，包括选择合适的SMILES列、处理缺失值、计算化学描述符等，然后调用Chemprop的命令行接口进行预测，并将结果与输入数据合并返回。
def _predict_chemprop_inhibit_probability(
    frame: pd.DataFrame,
    classification_prepare_dir: str | Path,
    chemprop_model_path: str | Path,
    output_dir: str | Path,
    accelerator: str = "cpu",
    devices: str = "1",
) -> pd.DataFrame:
    """Run the Chemprop classifier and return inhibit probabilities for each row."""
    classification_prepare_dir = Path(classification_prepare_dir)
    chemprop_model_path = Path(chemprop_model_path)
    output_dir = Path(output_dir)

    descriptor_preprocessor_path = classification_prepare_dir / "descriptor_preprocessor.joblib"
    descriptor_schema_path = classification_prepare_dir / "descriptor_schema.json"
    if not descriptor_preprocessor_path.exists():
        raise FileNotFoundError(
            f"Expected Chemprop descriptor preprocessor at {descriptor_preprocessor_path}. "
            "Re-run scripts/prepare_step1_chemprop.py after updating the repository."
        )
    if not descriptor_schema_path.exists():
        raise FileNotFoundError(
            f"Expected Chemprop descriptor schema at {descriptor_schema_path}. "
            "Re-run scripts/prepare_step1_chemprop.py after updating the repository."
        )

    descriptor_preprocessor = load(descriptor_preprocessor_path)
    descriptor_schema = json.loads(descriptor_schema_path.read_text(encoding="utf-8"))
    numeric_features = [str(value) for value in descriptor_schema.get("numeric_descriptor_features", [])]
    categorical_features = [str(value) for value in descriptor_schema.get("categorical_descriptor_features", [])]

    inference = frame.copy()
    inference["smiles"] = _pick_smiles(inference)
    non_empty_mask = inference["smiles"].notna() & inference["smiles"].astype(str).str.strip().ne("")
    valid_mask = non_empty_mask & _is_valid_smiles_series(inference, "smiles")

    result = inference.loc[:, ["hybrid_row_id", "pair_id", "prestwick_id", "nt_code", "smiles"]].copy()
    result["predicted_inhibit_probability"] = np.nan

    if not valid_mask.any():
        return result

    valid = inference.loc[valid_mask].copy()
    valid = _ensure_columns(valid, numeric_features + categorical_features)
    valid = _prepare_descriptor_frame(valid, numeric_features=numeric_features, categorical_features=categorical_features)

    descriptors = descriptor_preprocessor.transform(valid.loc[:, numeric_features + categorical_features])
    predict_input = valid.loc[:, ["hybrid_row_id", "pair_id", "prestwick_id", "nt_code", "smiles"]].copy()

    with tempfile.TemporaryDirectory(prefix="chemprop_predict_", dir=output_dir) as tmp_dir_name:
        tmp_dir = Path(tmp_dir_name)
        predict_csv = tmp_dir / "predict_input.csv"
        descriptors_path = tmp_dir / "predict_descriptors.npz"
        predictions_path = tmp_dir / "predictions.csv"
        predict_input.to_csv(predict_csv, index=False)
        np.savez_compressed(descriptors_path, descriptors)

        command = [
            _resolve_chemprop_bin(require_installed=True),
            "predict",
            "--test-path",
            str(predict_csv),
            "--smiles-columns",
            "smiles",
            "--descriptors-path",
            str(descriptors_path),
            "--model-paths",
            str(chemprop_model_path),
            "--output",
            str(predictions_path),
            "--accelerator",
            accelerator,
            "--devices",
            devices,
            "--num-workers",
            "0",
        ]
        subprocess.run(command, check=True)
        predicted = pd.read_csv(predictions_path, low_memory=False)

    predicted = predicted.rename(columns={"target": "predicted_inhibit_probability"})
    predicted = predicted.loc[:, ["hybrid_row_id", "predicted_inhibit_probability"]].copy()
    result = result.merge(predicted, on="hybrid_row_id", how="left", suffixes=("", "_new"))
    if "predicted_inhibit_probability_new" in result.columns:
        result["predicted_inhibit_probability"] = result["predicted_inhibit_probability_new"].combine_first(
            result["predicted_inhibit_probability"]
        )
        result = result.drop(columns=["predicted_inhibit_probability_new"])
    return result

# _predict_regression_effect_score：使用基于RDKit特征的回归模型对输入数据进行预测，返回每行的效果分数。它首先加载训练好的回归模型和相关的特征列表，然后确保输入数据包含这些特征列，并进行必要的预处理，最后调用模型进行预测并将结果与输入数据合并返回。
def _predict_regression_effect_score(
    frame: pd.DataFrame,
    regressor_path: str | Path,
    regressor_metrics_path: str | Path,
) -> pd.DataFrame:
    """Run the RDKit baseline regressor and return predicted effect scores."""
    regressor_path = Path(regressor_path)
    regressor_metrics_path = Path(regressor_metrics_path)

    regressor = load(regressor_path)
    regressor_metrics = json.loads(regressor_metrics_path.read_text(encoding="utf-8"))
    numeric_features = [str(value) for value in regressor_metrics.get("numeric_features", [])]
    categorical_features = [str(value) for value in regressor_metrics.get("categorical_features", [])]
    feature_columns = numeric_features + categorical_features

    inference = _ensure_columns(frame, feature_columns)
    inference = _prepare_feature_frame(inference, numeric_features=numeric_features, categorical_features=categorical_features)
    predicted_score = regressor.predict(inference.loc[:, feature_columns])

    result = frame.loc[:, ["hybrid_row_id", "pair_id", "prestwick_id", "nt_code"]].copy()
    result["predicted_effect_score"] = predicted_score
    return result

# _predict_promote_probability：运行可选的辅助促进分类器并返回促进概率。它首先加载促进分类器和相关的特征列表，然后确保输入数据包含这些特征列，并进行必要的预处理，最后调用分类器进行预测并将结果与输入数据合并返回。
# predict_step1_hybrid：生成Step 1的混合预测结果，通过结合Chemprop的抑制概率和基于RDKit特征的回归效果分数，来生成一个综合的效果标签。它还会输出详细的预测文件和一个包含统计摘要的JSON文件。
def predict_step1_hybrid(
    input_table_path: str | Path,
    output_dir: str | Path,
    classification_prepare_dir: str | Path,
    chemprop_model_path: str | Path,
    regressor_path: str | Path,
    regressor_metrics_path: str | Path,
    inhibit_probability_threshold: float = 0.5,
    promote_score_threshold: float = 0.2,
    accelerator: str = "cpu",
    devices: str = "1",
) -> dict[str, object]:
    """Generate Step 1 hybrid predictions by combining Chemprop and baseline regression.

    Args:
        input_table_path: Input candidate pair table for inference.
        output_dir: Directory for prediction CSVs and summary JSON.
        classification_prepare_dir: Chemprop preparation directory with descriptor metadata.
        chemprop_model_path: Trained Chemprop classification checkpoint.
        regressor_path: Trained baseline regression model.
        regressor_metrics_path: Metrics JSON containing regressor feature lists.
        inhibit_probability_threshold: Probability cutoff for an inhibit call.
        promote_score_threshold: Regression cutoff for a promote call.
        accelerator: Chemprop inference accelerator name.
        devices: Device specification passed to Chemprop.

    Returns:
        A summary dictionary describing the generated prediction files.
    """
    input_table_path = Path(input_table_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(input_table_path, low_memory=False)
    raw = _ensure_identifier_columns(raw)
    raw = raw.reset_index(drop=True).copy()
    raw["hybrid_row_id"] = np.arange(len(raw))

    chemprop_predictions = _predict_chemprop_inhibit_probability(
        frame=raw,
        classification_prepare_dir=classification_prepare_dir,
        chemprop_model_path=chemprop_model_path,
        output_dir=output_dir,
        accelerator=accelerator,
        devices=devices,
    )
    regression_predictions = _predict_regression_effect_score(
        frame=raw,
        regressor_path=regressor_path,
        regressor_metrics_path=regressor_metrics_path,
    )

    predictions = raw.merge(
        chemprop_predictions.loc[:, ["hybrid_row_id", "predicted_inhibit_probability"]],
        on="hybrid_row_id",
        how="left",
    )
    predictions = predictions.merge(
        regression_predictions.loc[:, ["hybrid_row_id", "predicted_effect_score"]],
        on="hybrid_row_id",
        how="left",
    )
    predictions["predicted_binary_effect_label"] = np.where(
        predictions["predicted_inhibit_probability"] >= inhibit_probability_threshold,
        "inhibit",
        "no_effect",
    )
    predictions["predicted_effect_label_hybrid"] = [
        _hybrid_effect_label(
            inhibit_probability=inhibit_probability,
            effect_score=effect_score,
            inhibit_probability_threshold=inhibit_probability_threshold,
            promote_score_threshold=promote_score_threshold,
        )
        for inhibit_probability, effect_score in zip(
            predictions["predicted_inhibit_probability"],
            predictions["predicted_effect_score"],
        )
    ]
    predictions["predicted_effect_magnitude"] = predictions["predicted_effect_score"].abs()
    predictions, constraint_summary = _apply_step1_drug_profile_constraints(
        predictions,
        inhibit_probability_threshold=float(inhibit_probability_threshold),
        promote_score_threshold=float(promote_score_threshold),
    )

    predictions_output = predictions.drop(columns=["hybrid_row_id"])
    predictions_output_path = output_dir / "predictions.csv"
    predictions_output.to_csv(predictions_output_path, index=False)

    slim_columns = [
        "pair_id",
        "prestwick_id",
        "nt_code",
        "smiles",
        "effect_label",
        "binary_effect_label",
        "effect_score",
        "predicted_inhibit_probability",
        "predicted_binary_effect_label",
        "predicted_effect_score",
        "predicted_effect_label_hybrid",
        "predicted_effect_magnitude",
    ]
    existing_slim_columns = [column for column in slim_columns if column in predictions_output.columns]
    predictions_slim_path = output_dir / "predictions_slim.csv"
    predictions_output.loc[:, existing_slim_columns].to_csv(predictions_slim_path, index=False)

    summary = {
        "input_table_path": str(input_table_path),
        "output_dir": str(output_dir),
        "classification_prepare_dir": str(classification_prepare_dir),
        "chemprop_model_path": str(chemprop_model_path),
        "regressor_path": str(regressor_path),
        "regressor_metrics_path": str(regressor_metrics_path),
        "inhibit_probability_threshold": float(inhibit_probability_threshold),
        "promote_score_threshold": float(promote_score_threshold),
        "n_rows": int(len(predictions_output)),
        "n_rows_with_smiles": int(predictions_output["predicted_inhibit_probability"].notna().sum()),
        "n_rows_with_regression_score": int(predictions_output["predicted_effect_score"].notna().sum()),
        "predicted_binary_effect_label_counts": {
            str(key): int(value)
            for key, value in predictions_output["predicted_binary_effect_label"].value_counts().to_dict().items()
        },
        "predicted_effect_label_hybrid_counts": {
            str(key): int(value)
            for key, value in predictions_output["predicted_effect_label_hybrid"].value_counts().to_dict().items()
        },
        "step1_drug_profile_counts": {
            str(key): int(value)
            for key, value in predictions_output.get("step1_drug_profile", pd.Series(dtype=object)).value_counts().to_dict().items()
        },
        "step1_constraint_summary": constraint_summary,
        "predictions_path": str(predictions_output_path),
        "predictions_slim_path": str(predictions_slim_path),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    return summary
