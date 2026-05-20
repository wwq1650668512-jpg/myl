from __future__ import annotations

import json
import math
import re
import tempfile
import uuid
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from gut_drug_microbiome.amr import AmrRuleEngine
from gut_drug_microbiome.disease_knowledge import build_disease_adjusted_community
from gut_drug_microbiome.mechanism_layer import attach_action_signals
from gut_drug_microbiome.mechanism_layer import compute_mechanism_layer
from gut_drug_microbiome.mechanism_layer import fuse_disease_scores
from gut_drug_microbiome.mechanism_layer import infer_microbe_trait_priors
from gut_drug_microbiome.step1 import annotate_compound_semantics
from gut_drug_microbiome.step1 import predict_step1_hybrid
from gut_drug_microbiome.step1 import refine_step1_promote_with_step2
from gut_drug_microbiome.step1.chem_features import enrich_drug_table_with_rdkit
from gut_drug_microbiome.step2 import Step2MechanismProjector
from gut_drug_microbiome.step2 import build_step2_input_tables
from gut_drug_microbiome.step2 import predict_step2_baseline
from gut_drug_microbiome.step3 import BUILTIN_SCENARIOS
from gut_drug_microbiome.step3 import run_step3_simulation
from gut_drug_microbiome.utils.chem import compute_smiles_descriptors as _compute_smiles_descriptors
from gut_drug_microbiome.utils.text import canonicalize_key as _canonicalize_key


ROOT = Path(__file__).resolve().parents[3]
EXPANDED_INTEGRATED_PREDICTIONS_PATH = ROOT / "predictions/step2/baseline_scaffold_v1_83/predictions.csv"
LEGACY_INTEGRATED_PREDICTIONS_PATH = ROOT / "predictions/step2/baseline_scaffold_v1/predictions.csv"
DEFAULT_INTEGRATED_PREDICTIONS_PATH = (
    EXPANDED_INTEGRATED_PREDICTIONS_PATH
    if EXPANDED_INTEGRATED_PREDICTIONS_PATH.exists()
    else LEGACY_INTEGRATED_PREDICTIONS_PATH
)
DEFAULT_DEMO_RANKING_PATH = ROOT / "predictions/step3/candidate_screen_demo/candidate_ranking.csv"
DEFAULT_TEMP_ROOT = Path("/tmp/gut_drug_microbiome_webapp")
DEFAULT_STEP1_MICROBE_TABLE_PATH = ROOT / "data/processed/step1/step1_microbe_table.csv"
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
DEFAULT_STEP2_MECHANISM_REFERENCE_PATH = ROOT / "models/step2/zimmermann_scaffold_split/mechanism_reference.joblib"
DEFAULT_STEP2_ENZYME_MICROBE_PANEL_PATH = ROOT / "data/reference/step2_microbe_enzyme_prior_long.csv"
DEFAULT_STEP2_ENZYME_FUNCTION_CATALOG_PATH = ROOT / "data/reference/step2_enzyme_function_catalog.csv"
DEFAULT_STEP2_BIOTRANSFORM_REFERENCE_PATH = (
    ROOT / "predictions/step2/baseline_scaffold_v1_83/predictions_biotransform_experimental.csv"
)
DEFAULT_STEP3_COHORT_ROOT = ROOT / "data/processed/step3/cohorts"
DEFAULT_STEP3_HEALTH_SIGNATURE_PROXY_PATH = ROOT / "data/processed/health_signature/microbe_tcg_proxy_mapping.csv"
DEFAULT_DISEASE_MICROBE_REFERENCE_PATH = ROOT / "data/reference/disease_microbe_dictionary.csv"
DEFAULT_DISEASE_MICROBE_SUPPLEMENT_PATHS = (
    ROOT / "data/reference/disease_microbe_gmrepo_supplement.csv",
    ROOT / "data/reference/disease_microbe_gutm_disorder_supplement.csv",
    ROOT / "data/reference/disease_microbe_gist_supplement.csv",
)
DEFAULT_DISEASE_DRUG_REFERENCE_PATH = ROOT / "data/reference/disease_marketed_drug_catalog.csv"

REQUIRED_DISEASE_CATALOG_ENTRIES = [
    "肠易激综合征（IBS）",
    "肠易激综合征-腹泻型（IBS-D）",
    "肠易激综合征-便秘型（IBS-C）",
]
BENCHMARK_DISEASE_PANEL_ENTRIES = [
    "克罗恩病（CD）",
    "溃疡性结肠炎（UC）",
    "便秘（Constipation）",
    "肛周脓肿（Anorectal Abscess）",
]
IBS_STANDARD_NAME = "肠易激综合征（IBS）"
IBS_D_STANDARD_NAME = "肠易激综合征-腹泻型（IBS-D）"
IBS_C_STANDARD_NAME = "肠易激综合征-便秘型（IBS-C）"
CORE_BUTYRATE_REGEX = r"Faecalibacterium\s+prausnitzii|Roseburia|Eubacterium\s+rectale"
CD_STANDARD_NAME = "克罗恩病（CD）"
UC_STANDARD_NAME = "溃疡性结肠炎（UC）"
CRC_STANDARD_NAME = "结直肠癌（CRC）"
CONSTIPATION_STANDARD_NAME = "便秘（Constipation）"
DIARRHEA_STANDARD_NAME = "腹泻（Diarrhea）"
DISEASE_RELATION_LEVEL_WEIGHTS = {
    "species": 1.0,
    "genus": 0.75,
    "family": 0.45,
    "phylum": 0.30,
    "class": 0.25,
    "order": 0.25,
}
DISEASE_RELATION_SOURCE_WEIGHTS = {
    "microbe_to_disease": 1.0,
    "disease_to_microbe": 0.7,
    "gmrepo_health_vs_disease": 0.95,
}
DISEASE_RELATION_CONFIDENCE_WEIGHTS = {"high": 1.0, "medium": 0.9, "low": 0.6}
SIMILARITY_MIN_THRESHOLD = 0.35
SIMILARITY_DRUG_TOP_K = 8
SIMILARITY_DISEASE_TOP_K = 5
BIOTRANSFORM_MIN_THRESHOLD = 0.22
BIOTRANSFORM_TOP_K = 3
BIOTRANSFORM_DIRECT_MATCH_THRESHOLD = 0.999

DESIRED_COLUMNS = [
    "pair_id",
    "prestwick_id",
    "chemical_name",
    "therapeutic_class",
    "therapeutic_effect",
    "atc_primary_l1",
    "atc_primary_l3",
    "atc_primary_l4",
    "smiles",
    "canonical_smiles_rdkit",
    "molecular_formula",
    "molecular_weight",
    "xlogp",
    "tpsa",
    "murcko_scaffold",
    "nt_code",
    "microbe_label",
    "species_label",
    "species_name",
    "genus",
    "family",
    "phylum",
    "gram_stain",
    "medium_preference",
    "biosafety",
    "step1_observed_effect_label",
    "step1_observed_binary_effect_label",
    "step1_observed_effect_score",
    "step1_predicted_inhibit_probability",
    "step1_predicted_binary_effect_label",
    "step1_predicted_effect_score",
    "step1_predicted_effect_label_hybrid",
    "step1_predicted_effect_magnitude",
    "predicted_inhibit_probability",
    "predicted_binary_effect_label",
    "predicted_effect_score",
    "predicted_effect_label_hybrid",
    "predicted_effect_magnitude",
    "predicted_promote_probability_base",
    "predicted_promote_probability_refined",
    "predicted_promote_support_score",
    "predicted_promote_support_type",
    "predicted_promote_evidence_type",
    "predicted_cross_feeding_reference_flag",
    "predicted_cross_feeding_support_microbe",
    "predicted_cross_feeding_reference_pmid",
    "predicted_cross_feeding_evidence_level",
    "predicted_cross_feeding_match_mode",
    "predicted_cross_feeding_matched_term",
    "predicted_effect_label_step2_refined",
    "predicted_effect_label_step2_refined_changed",
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
    "predicted_reaction_support_pairs",
    "predicted_mechanism_support_score",
    "predicted_mechanism_support_scopes",
    "predicted_candidate_product_ids",
    "predicted_candidate_product_count",
    "predicted_evidence_gene_ids",
    "predicted_evidence_gene_count",
    "predicted_enzyme_prior_flag",
    "predicted_enzyme_match_count",
    "predicted_enzyme_ids",
    "predicted_enzyme_names",
    "predicted_enzyme_reaction_classes",
    "predicted_enzyme_bond_targets",
    "predicted_enzyme_presence_score",
    "predicted_enzyme_support_score",
    "predicted_enzyme_step1_promote_support_score",
    "predicted_enzyme_step1_inhibit_risk_score",
]

STEP1_LABEL_CANDIDATES = ["predicted_effect_label_step2_refined", "step1_predicted_effect_label_hybrid", "predicted_effect_label_hybrid"]
STEP1_BINARY_CANDIDATES = ["step1_predicted_binary_effect_label", "predicted_binary_effect_label"]
STEP1_PROBABILITY_CANDIDATES = ["step1_predicted_inhibit_probability", "predicted_inhibit_probability"]
STEP1_SCORE_CANDIDATES = ["step1_predicted_effect_score", "predicted_effect_score"]
STEP1_MAGNITUDE_CANDIDATES = ["step1_predicted_effect_magnitude", "predicted_effect_magnitude"]
STEP1_OBSERVED_LABEL_CANDIDATES = ["step1_observed_effect_label", "effect_label"]
STEP1_OBSERVED_BINARY_CANDIDATES = ["step1_observed_binary_effect_label", "binary_effect_label"]
STEP1_OBSERVED_SCORE_CANDIDATES = ["step1_observed_effect_score", "effect_score"]
MICROBE_TAXONOMY_COLUMNS = ["species_label", "microbe_label", "genus", "family", "order", "class", "phylum", "gram_stain", "medium_preference"]
MICROBE_CONTEXT_NUMERIC_COLUMNS = ["starting_od_96_well_screen", "starting_od_384_well_screen"]
GENUS_TAXONOMY_FALLBACK = {
    "Actinomyces": {
        "family": "Actinomycetaceae",
        "order": "Actinomycetales",
        "class": "Actinobacteria",
        "phylum": "Actinobacteria",
        "gram_stain": "positive",
        "medium_preference": "mGAM",
    },
    "Alistipes": {
        "family": "Rikenellaceae",
        "order": "Bacteroidales",
        "class": "Bacteroidia",
        "phylum": "Bacteroidetes",
        "gram_stain": "negative",
        "medium_preference": "mGAM",
    },
    "Butyrivibrio": {
        "family": "Lachnospiraceae",
        "order": "Clostridiales",
        "class": "Clostridia",
        "phylum": "Firmicutes",
        "gram_stain": "positive",
        "medium_preference": "mGAM",
    },
    "Desulfovibrio": {
        "family": "Desulfovibrionaceae",
        "order": "Desulfovibrionales",
        "class": "Deltaproteobacteria",
        "phylum": "Proteobacteria",
        "gram_stain": "negative",
        "medium_preference": "mGAM supplemented with 60 mM sodium formiate and 10 mM taurine",
    },
    "Dialister": {
        "family": "Veillonellaceae",
        "order": "Selenomonadales",
        "class": "Negativicutes",
        "phylum": "Firmicutes",
        "gram_stain": "negative",
        "medium_preference": "mGAM",
    },
    "Faecalibacterium": {
        "family": "Ruminococcaceae",
        "order": "Clostridiales",
        "class": "Clostridia",
        "phylum": "Firmicutes",
        "gram_stain": "positive",
        "medium_preference": "mGAM",
    },
    "Haemophilus": {
        "family": "Pasteurellaceae",
        "order": "Pasteurellales",
        "class": "Gammaproteobacteria",
        "phylum": "Proteobacteria",
        "gram_stain": "negative",
        "medium_preference": "Todd-Hewitt+0.6% sodium lactate",
    },
    "Parvimonas": {
        "family": "Peptostreptococcaceae",
        "order": "Clostridiales",
        "class": "Clostridia",
        "phylum": "Firmicutes",
        "gram_stain": "positive",
        "medium_preference": "mGAM",
    },
    "Pseudoflavonifractor": {
        "family": "Ruminococcaceae",
        "order": "Clostridiales",
        "class": "Clostridia",
        "phylum": "Firmicutes",
        "gram_stain": "positive",
        "medium_preference": "mGAM",
    },
}


def _pick_column(frame: pd.DataFrame, candidates: list[str]) -> str | None:
    for candidate in candidates:
        if candidate in frame.columns:
            return candidate
    return None


def _safe_float(value: object) -> float | None:
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(numeric):
        return None
    return float(numeric)


def _safe_bool(value: object) -> bool | None:
    if pd.isna(value):
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no"}:
        return False
    return None


def _fingerprint_array_from_row(row: pd.Series, columns: list[str]) -> np.ndarray:
    if not columns:
        return np.zeros((0,), dtype=bool)
    values = pd.to_numeric(row.reindex(columns), errors="coerce").fillna(0)
    return values.to_numpy(dtype=bool)


def _tanimoto_similarity_vector(query_bits: np.ndarray, reference_matrix: np.ndarray) -> np.ndarray:
    if query_bits.size == 0 or reference_matrix.size == 0:
        return np.zeros((reference_matrix.shape[0],), dtype=float)
    query = query_bits.astype(bool)
    reference = reference_matrix.astype(bool)
    intersection = np.logical_and(reference, query).sum(axis=1).astype(float)
    union = np.logical_or(reference, query).sum(axis=1).astype(float)
    similarities = np.zeros((reference.shape[0],), dtype=float)
    valid = union > 0
    similarities[valid] = intersection[valid] / union[valid]
    return similarities


def _string_list_from_semicolon(value: object, *, max_items: int = 8) -> list[str]:
    text = str(value or "").strip()
    if not text or text.lower() in {"nan", "none", "null", "n/a", "na"}:
        return []
    items: list[str] = []
    for raw in text.split(";"):
        item = raw.strip()
        if item and item.lower() not in {"nan", "none", "null", "n/a", "na"} and item not in items:
            items.append(item)
        if len(items) >= max_items:
            break
    return items


def _clean_json(value: object) -> object:
    if isinstance(value, dict):
        return {str(key): _clean_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean_json(item) for item in value]
    if isinstance(value, pd.DataFrame):
        return [_clean_json(record) for record in value.to_dict(orient="records")]
    if isinstance(value, pd.Series):
        return {str(key): _clean_json(item) for key, item in value.to_dict().items()}
    if value is None:
        return None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        numeric = float(value)
        if math.isnan(numeric) or math.isinf(numeric):
            return None
        return numeric
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if pd.isna(value):
        return None
    return value


def _existing_usecols(path: Path, desired: list[str]) -> list[str]:
    columns = pd.read_csv(path, nrows=0).columns.tolist()
    return [column for column in desired if column in columns]


def _normalize_null_like_strings(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """Normalize textual null markers ('nan', 'null', etc.) to NaN for selected columns."""
    result = frame.copy()
    for column in columns:
        if column not in result.columns:
            continue
        series = result[column]
        if not pd.api.types.is_object_dtype(series) and not pd.api.types.is_string_dtype(series):
            continue
        normalized = series.astype(str).str.strip()
        lower = normalized.str.lower()
        mask = lower.isin({"", "nan", "none", "null", "n/a", "na"})
        result.loc[mask, column] = np.nan
    return result


def _mode_non_empty(series: pd.Series) -> str | None:
    normalized = series.dropna().astype(str).str.strip()
    normalized = normalized[~normalized.str.lower().isin({"", "nan", "none", "null", "n/a", "na"})]
    if normalized.empty:
        return None
    return str(normalized.value_counts().index[0])


def _extract_genus_name(raw_name: object) -> str | None:
    text = str(raw_name or "").strip()
    if not text:
        return None
    ascii_tokens = re.findall(r"[A-Za-z]+", text)
    if not ascii_tokens:
        return None
    return ascii_tokens[0]


def _enrich_microbe_taxonomy(frame: pd.DataFrame) -> pd.DataFrame:
    """Fill missing taxonomy/context fields for the 83-microbe panel using table-local priors."""
    result = frame.copy()
    result = _normalize_null_like_strings(result, MICROBE_TAXONOMY_COLUMNS + MICROBE_CONTEXT_NUMERIC_COLUMNS)

    for column in MICROBE_TAXONOMY_COLUMNS + MICROBE_CONTEXT_NUMERIC_COLUMNS:
        if column not in result.columns:
            result[column] = np.nan

    inferred_species = result["species_label"].fillna(result["microbe_label"]).map(_extract_genus_name)
    inferred_microbe = result["microbe_label"].map(_extract_genus_name)
    result["genus"] = result["genus"].fillna(inferred_species).fillna(inferred_microbe)

    genus_clean = result["genus"].fillna("").astype(str).str.strip()
    genus_key = genus_clean.str.lower()
    taxonomy_targets = ["family", "order", "class", "phylum", "gram_stain", "medium_preference"]

    for column in taxonomy_targets:
        known = result.loc[genus_clean.ne("") & result[column].notna(), ["genus", column]].copy()
        if known.empty:
            continue
        known["genus_key"] = known["genus"].astype(str).str.strip().str.lower()
        lookup = known.groupby("genus_key")[column].agg(_mode_non_empty).dropna().to_dict()
        fill_values = genus_key.map(lookup)
        result[column] = result[column].fillna(fill_values)

    for genus_name, fallback in GENUS_TAXONOMY_FALLBACK.items():
        mask = genus_clean.eq(genus_name)
        if not mask.any():
            continue
        for column, value in fallback.items():
            if column in result.columns:
                result.loc[mask, column] = result.loc[mask, column].fillna(value)

    for column in MICROBE_CONTEXT_NUMERIC_COLUMNS:
        numeric = pd.to_numeric(result[column], errors="coerce")
        result[column] = numeric
        if column == "starting_od_384_well_screen":
            # Keep this field conservative to avoid introducing artificial assay precision.
            continue
        group_median = result.groupby(genus_key)[column].transform("median")
        result[column] = result[column].fillna(group_median)

    return result


def _is_ibs_like_disease(disease_name: object) -> bool:
    key = _canonicalize_key(disease_name)
    return any(token in key for token in ["ibs", "肠易激", "ibsd", "ibsc"])


def _species_base_key(raw_name: object) -> str:
    text = str(raw_name or "")
    ascii_tokens = re.findall(r"[A-Za-z]+", text)
    if len(ascii_tokens) >= 2:
        return _canonicalize_key(" ".join(ascii_tokens[:2]))
    return _canonicalize_key(text)


def _clip01(value: float) -> float:
    return float(np.clip(float(value), 0.0, 1.0))


def _mean_boolean_like(series: pd.Series, default: float = 1.0) -> float:
    if series.empty:
        return float(default)
    text = series.fillna("").astype(str).str.strip().str.lower()
    mapped = text.map({"true": 1.0, "false": 0.0, "1": 1.0, "0": 0.0, "yes": 1.0, "no": 0.0})
    numeric = pd.to_numeric(series, errors="coerce")
    values = mapped.where(mapped.notna(), numeric)
    values = values.dropna()
    if values.empty:
        return float(default)
    return _clip01(float(values.mean()))


def _numeric_bounds(series: pd.Series, low_quantile: float = 0.01, high_quantile: float = 0.99) -> tuple[float | None, float | None]:
    numeric = pd.to_numeric(series, errors="coerce").dropna()
    if numeric.empty:
        return None, None
    if len(numeric) < 20:
        return float(numeric.min()), float(numeric.max())
    return float(numeric.quantile(low_quantile)), float(numeric.quantile(high_quantile))


def _first_numeric_value(row: pd.Series, columns: list[str]) -> float | None:
    for column in columns:
        if column not in row.index:
            continue
        value = _safe_float(row.get(column))
        if value is not None:
            return float(value)
    return None


def _sigmoid_100(value: float, center: float = 0.0, scale: float = 8.0) -> float:
    safe_scale = max(float(scale), 1e-6)
    return float(100.0 / (1.0 + math.exp(-(float(value) - float(center)) / safe_scale)))


def _inject_placebo_deltas(
    active_metrics: pd.DataFrame,
    placebo_metrics: pd.DataFrame,
) -> pd.DataFrame:
    if active_metrics.empty or placebo_metrics.empty or "timepoint" not in active_metrics.columns or "timepoint" not in placebo_metrics.columns:
        return active_metrics
    placebo_renamed = placebo_metrics.copy()
    placebo_renamed = placebo_renamed.rename(
        columns={
            column: f"{column}_placebo"
            for column in placebo_renamed.columns
            if column != "timepoint"
        }
    )
    merged = active_metrics.merge(placebo_renamed, on="timepoint", how="left")
    delta_targets = [
        "health_index",
        "development_score",
        "experimental_development_score",
        "development_score_balance",
        "experimental_development_score_balance",
        "benefit_subscore",
        "risk_subscore",
        "experimental_risk_subscore",
        "parent_retention_ratio",
        "experimental_aggregate_metabolite_pool",
        "disease_target_alignment_score",
    ]
    for column in delta_targets:
        placebo_column = f"{column}_placebo"
        if column not in merged.columns or placebo_column not in merged.columns:
            continue
        merged[f"{column}_delta_vs_placebo"] = pd.to_numeric(merged[column], errors="coerce").fillna(0.0) - pd.to_numeric(
            merged[placebo_column],
            errors="coerce",
        ).fillna(0.0)
    if "development_score_delta_vs_placebo" in merged.columns:
        merged["development_score_normalized_vs_placebo"] = merged["development_score_delta_vs_placebo"].map(
            lambda item: _sigmoid_100(item, center=0.0, scale=8.0)
        )
    drop_columns = [
        column
        for column in placebo_renamed.columns
        if column != "timepoint" and column in merged.columns
    ]
    return merged.drop(columns=drop_columns)


def _attach_placebo_summary_deltas(
    active_summary: dict[str, object],
    placebo_summary: dict[str, object],
) -> dict[str, object]:
    result = dict(active_summary)
    delta_targets = [
        "development_score",
        "experimental_development_score",
        "development_score_balance",
        "experimental_development_score_balance",
        "health_index",
        "parent_retention_ratio",
        "experimental_aggregate_metabolite_pool",
        "benefit_subscore",
        "risk_subscore",
        "experimental_risk_subscore",
        "disease_target_alignment_score",
    ]
    for column in delta_targets:
        active_value = _safe_float(
            result.get(
                f"final_{column}",
                result.get(
                    f"{column}_final",
                    result.get(column),
                ),
            )
        )
        placebo_value = _safe_float(
            placebo_summary.get(
                f"final_{column}",
                placebo_summary.get(
                    f"{column}_final",
                    placebo_summary.get(column),
                ),
            )
        )
        if active_value is None or placebo_value is None:
            continue
        result[f"{column}_delta_vs_placebo"] = float(active_value - placebo_value)
    if "development_score_delta_vs_placebo" in result:
        result["development_score_normalized_vs_placebo"] = _sigmoid_100(
            float(result["development_score_delta_vs_placebo"]),
            center=0.0,
            scale=8.0,
        )
    result["placebo_baseline"] = {
        "final_development_score": _safe_float(placebo_summary.get("development_score")),
        "final_experimental_development_score": _safe_float(placebo_summary.get("experimental_development_score")),
        "final_development_score_balance": _safe_float(placebo_summary.get("development_score_balance")),
        "final_experimental_development_score_balance": _safe_float(
            placebo_summary.get("experimental_development_score_balance")
        ),
        "final_health_index": _safe_float(placebo_summary.get("final_health_index")),
        "final_parent_retention_ratio": _safe_float(placebo_summary.get("final_parent_retention_ratio")),
        "final_experimental_aggregate_metabolite_pool": _safe_float(
            placebo_summary.get("final_experimental_aggregate_metabolite_pool")
        ),
        "final_benefit_subscore": _safe_float(placebo_summary.get("benefit_subscore_final")),
        "final_risk_subscore": _safe_float(placebo_summary.get("risk_subscore_final")),
        "final_experimental_risk_subscore": _safe_float(placebo_summary.get("experimental_risk_subscore_final")),
        "final_disease_target_alignment_score": _safe_float(placebo_summary.get("disease_target_alignment_score_final")),
    }
    return result


def _warning_reason_text(flag: str, breakdown: dict[str, object]) -> str:
    inhibit_fraction = _safe_float(breakdown.get("inhibit_fraction"))
    strong_core = pd.to_numeric(pd.Series([breakdown.get("strong_core_butyrate_count")]), errors="coerce").iloc[0]
    core_total = pd.to_numeric(pd.Series([breakdown.get("core_butyrate_total")]), errors="coerce").iloc[0]
    molecular_weight = _safe_float(breakdown.get("molecular_weight"))
    xlogp = _safe_float(breakdown.get("xlogp"))
    profile = str(breakdown.get("drug_profile") or "unknown")

    if flag == "over-suppression":
        if inhibit_fraction is not None:
            return f"整体抑制比例偏高（inhibit_fraction={inhibit_fraction:.2f}）"
        return "整体抑制比例偏高"
    if flag == "core-butyrate-suppression":
        if not pd.isna(strong_core) and not pd.isna(core_total):
            return f"核心产丁酸菌出现强抑制（{int(strong_core)}/{max(int(core_total), 1)}）"
        return "核心产丁酸菌出现强抑制"
    if flag == "ecology-risk":
        if not pd.isna(strong_core) and not pd.isna(core_total):
            return f"存在菌群生态风险信号（核心产丁酸菌受抑 {int(strong_core)}/{max(int(core_total), 1)}）"
        return "存在菌群生态风险信号"
    if flag == "drug-profile-conflict":
        if profile and profile != "unknown":
            return f"预测方向与药物类型（{profile}）不一致"
        return "预测方向与药物类型不一致"
    if flag == "antifolate-mismatch":
        if not pd.isna(strong_core) and not pd.isna(core_total):
            return f"磺胺类预期靶向效应不足（核心产丁酸菌强抑制仅 {int(strong_core)}/{max(int(core_total), 1)}）"
        return "磺胺类预期靶向效应不足（未见核心产丁酸菌抑制）"
    if flag == "OOD-molecule":
        details = []
        if molecular_weight is not None:
            details.append(f"MW={molecular_weight:.1f}")
        if xlogp is not None:
            details.append(f"logP={xlogp:.2f}")
        if details:
            return f"分子性质超出训练分布（{'，'.join(details)}）"
        return "分子性质超出训练分布"
    return flag


def build_confidence_explanation(
    *,
    confidence_score: float,
    confidence_tier: str,
    warning_flags: list[str],
    confidence_breakdown: dict[str, object],
) -> str:
    tier_cn = {"high": "高", "medium": "中", "low": "低"}.get(str(confidence_tier), str(confidence_tier))
    unique_flags = sorted(set(str(flag) for flag in warning_flags if str(flag).strip()))
    if not unique_flags:
        return (
            f"当前预测置信度{tier_cn}（{confidence_score:.2f}），未触发主要风险警告。"
            "模型在当前分子与菌群组合上未发现明显异常信号。"
        )

    reasons = [_warning_reason_text(flag, confidence_breakdown) for flag in unique_flags]
    lead_reason = reasons[0]
    if len(reasons) == 1:
        return f"当前预测置信度{tier_cn}（{confidence_score:.2f}），主要风险来自：{lead_reason}。"

    return (
        f"当前预测置信度{tier_cn}（{confidence_score:.2f}），主要风险来自：{lead_reason}。"
        f"另检测到 {len(reasons) - 1} 项风险信号。"
    )


def evaluate_prediction_confidence(
    *,
    effect_frame: pd.DataFrame,
    step1_label_column: str,
    step1_probability_column: str,
    step1_score_column: str,
    drug_profile: str,
    molecular_weight: float | None,
    xlogp: float | None,
    mw_bounds: tuple[float | None, float | None],
    xlogp_bounds: tuple[float | None, float | None],
) -> dict[str, object]:
    """Estimate prediction confidence and warning flags from ecology/risk heuristics.

    Rule rationale:
    - Over-suppression, core-butyrate suppression, and profile conflicts indicate biologically risky outputs.
    - MW/logP OOD indicates lower model familiarity with this molecule.
    """
    work = effect_frame.copy()
    warnings: list[str] = []
    penalties: list[dict[str, object]] = []

    label = work.get(step1_label_column, pd.Series("", index=work.index)).fillna("").astype(str).str.lower()
    inhibit_prob = pd.to_numeric(work.get(step1_probability_column, pd.Series(np.nan, index=work.index)), errors="coerce").fillna(0.0)
    effect_score = pd.to_numeric(work.get(step1_score_column, pd.Series(np.nan, index=work.index)), errors="coerce").fillna(0.0)

    panel_size = int(len(work))
    inhibit_count = int(label.eq("inhibit").sum())
    inhibit_fraction = float(inhibit_count / max(panel_size, 1))
    if inhibit_fraction > 0.70:
        warnings.append("over-suppression")
        penalties.append({"rule": "global_strong_inhibition", "penalty": 0.30, "evidence": f"inhibit_fraction={inhibit_fraction:.3f}"})

    butyrate_key = (
        work.get("microbe_label", pd.Series("", index=work.index)).fillna("").astype(str)
        + " "
        + work.get("species_label", pd.Series("", index=work.index)).fillna("").astype(str)
    )
    core_mask = butyrate_key.str.contains(CORE_BUTYRATE_REGEX, case=False, regex=True, na=False)
    core_total = int(core_mask.sum())
    strong_core = int(
        (core_mask & label.eq("inhibit") & ((inhibit_prob >= 0.50) | (effect_score <= -0.20))).sum()
    )
    core_strong_fraction = float(strong_core / max(core_total, 1))
    if strong_core > 0:
        warnings.append("ecology-risk")
        if core_strong_fraction >= 0.40:
            warnings.append("core-butyrate-suppression")
            penalties.append(
                {
                    "rule": "strong_core_butyrate_suppression",
                    "penalty": 0.25,
                    "evidence": f"strong_core={strong_core}/{max(core_total,1)}",
                }
            )
        else:
            penalties.append(
                {
                    "rule": "mild_core_butyrate_suppression",
                    "penalty": 0.10,
                    "evidence": f"strong_core={strong_core}/{max(core_total,1)}",
                }
            )

    profile_conflict = False
    profile_key = str(drug_profile or "unknown").strip().lower()
    if profile_key == "eubiotic_modulator":
        profile_conflict = inhibit_fraction > 0.70 or core_strong_fraction >= 0.40
    elif profile_key in {"host_secretagogue", "host_pathway_agent"}:
        profile_conflict = inhibit_fraction > 0.35
    elif profile_key == "disruptive_antibiotic":
        profile_conflict = inhibit_fraction < 0.20
    elif profile_key == "contextual_antimicrobial":
        profile_conflict = inhibit_fraction < 0.20
    elif profile_key == "sulfonamide_antifolate":
        profile_conflict = False
        if core_total > 0 and core_strong_fraction < 0.30:
            warnings.append("antifolate-mismatch")
            penalties.append(
                {
                    "rule": "sulfonamide_antifolate_core_mismatch",
                    "penalty": 0.22,
                    "evidence": (
                        f"profile={profile_key}, strong_core={strong_core}/{max(core_total,1)}, "
                        f"core_strong_fraction={core_strong_fraction:.3f}"
                    ),
                }
            )
    if profile_conflict:
        warnings.append("drug-profile-conflict")
        penalties.append(
            {
                "rule": "drug_profile_conflict",
                "penalty": 0.20,
                "evidence": f"profile={profile_key}, inhibit_fraction={inhibit_fraction:.3f}, core_strong_fraction={core_strong_fraction:.3f}",
            }
        )

    mw_low, mw_high = mw_bounds
    xlogp_low, xlogp_high = xlogp_bounds
    mw_ood = molecular_weight is not None and mw_low is not None and mw_high is not None and (
        molecular_weight < mw_low or molecular_weight > mw_high
    )
    xlogp_ood = xlogp is not None and xlogp_low is not None and xlogp_high is not None and (
        xlogp < xlogp_low or xlogp > xlogp_high
    )
    if mw_ood or xlogp_ood:
        warnings.append("OOD-molecule")
        penalties.append(
            {
                "rule": "molecule_out_of_distribution",
                "penalty": 0.22,
                "evidence": f"mw={molecular_weight}, bounds=({mw_low},{mw_high}); xlogp={xlogp}, bounds=({xlogp_low},{xlogp_high})",
            }
        )

    confidence = 0.90
    confidence -= float(sum(float(item.get("penalty", 0.0)) for item in penalties))
    confidence = _clip01(confidence)
    confidence = max(0.05, confidence)

    if confidence >= 0.75:
        confidence_tier = "high"
    elif confidence >= 0.45:
        confidence_tier = "medium"
    else:
        confidence_tier = "low"

    confidence_score = round(float(confidence), 4)
    warning_flags = sorted(set(warnings))
    confidence_breakdown = {
        "panel_size": panel_size,
        "inhibit_fraction": round(inhibit_fraction, 4),
        "strong_core_butyrate_count": strong_core,
        "core_butyrate_total": core_total,
        "core_butyrate_strong_fraction": round(core_strong_fraction, 4),
        "drug_profile": profile_key,
        "molecular_weight": molecular_weight,
        "xlogp": xlogp,
        "mw_bounds": [mw_low, mw_high],
        "xlogp_bounds": [xlogp_low, xlogp_high],
        "penalties": penalties,
    }
    return {
        "confidence_score": confidence_score,
        "confidence_tier": confidence_tier,
        "warning_flags": warning_flags,
        "confidence_breakdown": confidence_breakdown,
        "confidence_explanation": build_confidence_explanation(
            confidence_score=confidence_score,
            confidence_tier=confidence_tier,
            warning_flags=warning_flags,
            confidence_breakdown=confidence_breakdown,
        ),
    }


class GutPredictionService:
    def __init__(
        self,
        integrated_predictions_path: str | Path = DEFAULT_INTEGRATED_PREDICTIONS_PATH,
        demo_ranking_path: str | Path | None = DEFAULT_DEMO_RANKING_PATH,
        temp_root: str | Path = DEFAULT_TEMP_ROOT,
        step1_microbe_table_path: str | Path = DEFAULT_STEP1_MICROBE_TABLE_PATH,
        step1_chemprop_prepare_dir: str | Path = DEFAULT_STEP1_CHEMPROP_PREPARE_DIR,
        step1_chemprop_model_path: str | Path = DEFAULT_STEP1_CHEMPROP_MODEL_PATH,
        step1_regressor_path: str | Path = DEFAULT_STEP1_REGRESSOR_PATH,
        step1_regressor_metrics_path: str | Path = DEFAULT_STEP1_REGRESSOR_METRICS_PATH,
        step1_promote_classifier_path: str | Path | None = DEFAULT_STEP1_PROMOTE_CLASSIFIER_PATH,
        step1_promote_metrics_path: str | Path | None = DEFAULT_STEP1_PROMOTE_METRICS_PATH,
        cross_feeding_reference_path: str | Path | None = DEFAULT_CROSS_FEEDING_REFERENCE_PATH,
        step2_classifier_path: str | Path = DEFAULT_STEP2_CLASSIFIER_PATH,
        step2_regressor_path: str | Path = DEFAULT_STEP2_REGRESSOR_PATH,
        step2_metrics_path: str | Path = DEFAULT_STEP2_METRICS_PATH,
        step2_applicability_reference_path: str | Path = DEFAULT_STEP2_APPLICABILITY_REFERENCE_PATH,
        step2_mechanism_reference_path: str | Path = DEFAULT_STEP2_MECHANISM_REFERENCE_PATH,
        step2_enzyme_microbe_panel_path: str | Path | None = DEFAULT_STEP2_ENZYME_MICROBE_PANEL_PATH,
        step2_enzyme_function_catalog_path: str | Path | None = DEFAULT_STEP2_ENZYME_FUNCTION_CATALOG_PATH,
        step2_biotransform_reference_path: str | Path | None = DEFAULT_STEP2_BIOTRANSFORM_REFERENCE_PATH,
        step3_cohort_root: str | Path = DEFAULT_STEP3_COHORT_ROOT,
        step3_health_signature_proxy_path: str | Path = DEFAULT_STEP3_HEALTH_SIGNATURE_PROXY_PATH,
        disease_microbe_reference_path: str | Path | None = DEFAULT_DISEASE_MICROBE_REFERENCE_PATH,
        disease_microbe_supplement_paths: Sequence[str | Path] | None = DEFAULT_DISEASE_MICROBE_SUPPLEMENT_PATHS,
        disease_drug_reference_path: str | Path | None = DEFAULT_DISEASE_DRUG_REFERENCE_PATH,
    ) -> None:
        self.integrated_predictions_path = Path(integrated_predictions_path)
        self.demo_ranking_path = None if demo_ranking_path is None else Path(demo_ranking_path)
        self.temp_root = Path(temp_root)
        self.temp_root.mkdir(parents=True, exist_ok=True)
        self.step1_microbe_table_path = Path(step1_microbe_table_path)
        self.step1_chemprop_prepare_dir = Path(step1_chemprop_prepare_dir)
        self.step1_chemprop_model_path = Path(step1_chemprop_model_path)
        self.step1_regressor_path = Path(step1_regressor_path)
        self.step1_regressor_metrics_path = Path(step1_regressor_metrics_path)
        self.step1_promote_classifier_path = (
            None if step1_promote_classifier_path is None else Path(step1_promote_classifier_path)
        )
        self.step1_promote_metrics_path = None if step1_promote_metrics_path is None else Path(step1_promote_metrics_path)
        self.cross_feeding_reference_path = (
            None if cross_feeding_reference_path is None else Path(cross_feeding_reference_path)
        )
        self.step2_classifier_path = Path(step2_classifier_path)
        self.step2_regressor_path = Path(step2_regressor_path)
        self.step2_metrics_path = Path(step2_metrics_path)
        self.step2_applicability_reference_path = Path(step2_applicability_reference_path)
        self.step2_mechanism_reference_path = Path(step2_mechanism_reference_path)
        self.step2_enzyme_microbe_panel_path = (
            None if step2_enzyme_microbe_panel_path is None else Path(step2_enzyme_microbe_panel_path)
        )
        self.step2_enzyme_function_catalog_path = (
            None if step2_enzyme_function_catalog_path is None else Path(step2_enzyme_function_catalog_path)
        )
        self.step2_biotransform_reference_path = (
            None if step2_biotransform_reference_path is None else Path(step2_biotransform_reference_path)
        )
        self.step3_cohort_root = Path(step3_cohort_root)
        self.step3_health_signature_proxy_path = Path(step3_health_signature_proxy_path)
        self.disease_microbe_reference_path = (
            None if disease_microbe_reference_path is None else Path(disease_microbe_reference_path)
        )
        self.disease_microbe_supplement_paths = [
            Path(path) for path in (disease_microbe_supplement_paths or []) if path is not None
        ]
        self.disease_drug_reference_path = None if disease_drug_reference_path is None else Path(disease_drug_reference_path)

        usecols = _existing_usecols(self.integrated_predictions_path, DESIRED_COLUMNS)
        self.frame = pd.read_csv(self.integrated_predictions_path, usecols=usecols, low_memory=False)
        self.frame = _normalize_null_like_strings(
            self.frame,
            columns=[
                "predicted_reaction_class",
                "predicted_enzyme_names",
                "predicted_enzyme_ids",
                "predicted_enzyme_reaction_classes",
                "predicted_enzyme_bond_targets",
                "predicted_candidate_product_ids",
                "predicted_evidence_gene_ids",
            ],
        )
        self.frame = refine_step1_promote_with_step2(
            self.frame,
            promote_classifier_path=self.step1_promote_classifier_path,
            promote_metrics_path=self.step1_promote_metrics_path,
            cross_feeding_reference_path=self.cross_feeding_reference_path,
        )
        self.frame["drug_search_key"] = self.frame["prestwick_id"].map(_canonicalize_key)
        self.frame["drug_name_key"] = self.frame["chemical_name"].map(_canonicalize_key)
        self.frame["microbe_search_key"] = self.frame["nt_code"].map(_canonicalize_key)
        self.frame["microbe_name_key"] = self.frame["microbe_label"].map(_canonicalize_key)
        self.frame["microbe_species_key"] = self.frame["species_label"].map(_canonicalize_key)
        self.mw_ood_bounds = _numeric_bounds(self.frame.get("molecular_weight", pd.Series(dtype=float)))
        self.xlogp_ood_bounds = _numeric_bounds(self.frame.get("xlogp", pd.Series(dtype=float)))

        self.step1_label_column = _pick_column(self.frame, STEP1_LABEL_CANDIDATES)
        self.step1_binary_column = _pick_column(self.frame, STEP1_BINARY_CANDIDATES)
        self.step1_probability_column = _pick_column(self.frame, STEP1_PROBABILITY_CANDIDATES)
        self.step1_score_column = _pick_column(self.frame, STEP1_SCORE_CANDIDATES)
        self.step1_magnitude_column = _pick_column(self.frame, STEP1_MAGNITUDE_CANDIDATES)
        self.step1_observed_label_column = _pick_column(self.frame, STEP1_OBSERVED_LABEL_CANDIDATES)
        self.step1_observed_binary_column = _pick_column(self.frame, STEP1_OBSERVED_BINARY_CANDIDATES)
        self.step1_observed_score_column = _pick_column(self.frame, STEP1_OBSERVED_SCORE_CANDIDATES)

        self.drug_table = (
            self.frame.loc[
                :,
                [
                    column
                    for column in [
                        "prestwick_id",
                        "chemical_name",
                        "therapeutic_class",
                        "therapeutic_effect",
                        "atc_primary_l1",
                        "atc_primary_l3",
                        "atc_primary_l4",
                        "molecular_formula",
                        "molecular_weight",
                        "xlogp",
                        "tpsa",
                        "murcko_scaffold",
                        "canonical_smiles_rdkit",
                        "smiles",
                    ]
                    if column in self.frame.columns
                ],
            ]
            .drop_duplicates(subset=["prestwick_id"])
            .sort_values(["chemical_name", "prestwick_id"])
            .reset_index(drop=True)
        )
        drug_similarity_features = enrich_drug_table_with_rdkit(
            self.drug_table.loc[
                :,
                [column for column in ["canonical_smiles_rdkit", "smiles"] if column in self.drug_table.columns],
            ].copy(),
            smiles_columns=["canonical_smiles_rdkit", "smiles"],
        )
        drug_similarity_features = drug_similarity_features.loc[:, ~drug_similarity_features.columns.duplicated()].copy()
        similarity_feature_columns = [column for column in drug_similarity_features.columns if str(column).startswith("morgan_fp_")]
        if similarity_feature_columns:
            self.drug_table = pd.concat(
                [self.drug_table.reset_index(drop=True), drug_similarity_features.loc[:, similarity_feature_columns].reset_index(drop=True)],
                axis=1,
            )
        if "murcko_scaffold" in drug_similarity_features.columns and "murcko_scaffold" not in self.drug_table.columns:
            self.drug_table["murcko_scaffold"] = drug_similarity_features["murcko_scaffold"]
        self.fingerprint_columns = sorted(column for column in self.drug_table.columns if str(column).startswith("morgan_fp_"))
        self.drug_table["search_key"] = self.drug_table["prestwick_id"].map(_canonicalize_key)
        self.drug_table["name_key"] = self.drug_table["chemical_name"].map(_canonicalize_key)
        similarity_columns = [
            column
            for column in [
                "prestwick_id",
                "chemical_name",
                "murcko_scaffold",
                "canonical_smiles_rdkit",
                "smiles",
            ]
            if column in self.drug_table.columns
        ] + self.fingerprint_columns
        self.drug_similarity_table = (
            self.drug_table.loc[:, similarity_columns]
            .drop_duplicates(subset=["prestwick_id"])
            .sort_values(["chemical_name", "prestwick_id"])
            .reset_index(drop=True)
        )
        self.drug_similarity_matrix = (
            self.drug_similarity_table.loc[:, self.fingerprint_columns].fillna(0).to_numpy(dtype=bool)
            if self.fingerprint_columns
            else np.zeros((len(self.drug_similarity_table), 0), dtype=bool)
        )
        self.biotransform_reference_table = self._load_biotransform_reference_table()
        self.biotransform_reference_matrix = (
            self.biotransform_reference_table.loc[:, self.fingerprint_columns].fillna(0).to_numpy(dtype=bool)
            if (not self.biotransform_reference_table.empty and self.fingerprint_columns)
            else np.zeros((len(self.biotransform_reference_table), 0), dtype=bool)
        )
        self.frame = self._annotate_biotransform_sidecar(self.frame)

        self.library_microbe_table = (
            self.frame.loc[
                :,
                [
                    column
                    for column in [
                        "nt_code",
                        "microbe_label",
                        "species_label",
                        "species_name",
                        "genus",
                        "family",
                        "phylum",
                        "gram_stain",
                        "medium_preference",
                        "biosafety",
                    ]
                    if column in self.frame.columns
                ],
            ]
            .drop_duplicates(subset=["nt_code"])
            .sort_values(["microbe_label", "nt_code"])
            .reset_index(drop=True)
        )
        self.library_microbe_table["search_key"] = self.library_microbe_table["nt_code"].map(_canonicalize_key)
        self.library_microbe_table["name_key"] = self.library_microbe_table["microbe_label"].map(_canonicalize_key)
        self.library_microbe_table["species_key"] = self.library_microbe_table["species_label"].map(_canonicalize_key)

        microbe_feature_table = pd.read_csv(self.step1_microbe_table_path, low_memory=False)
        microbe_feature_table = _enrich_microbe_taxonomy(microbe_feature_table)
        self.microbe_feature_table = (
            microbe_feature_table.drop_duplicates(subset=["nt_code"])
            .sort_values(["microbe_label", "nt_code"])
            .reset_index(drop=True)
        )
        self.microbe_table = self.microbe_feature_table.loc[
            :,
            [
                column
                for column in [
                    "nt_code",
                    "microbe_label",
                    "species_label",
                    "species_name",
                    "genus",
                    "family",
                    "phylum",
                    "gram_stain",
                    "medium_preference",
                    "biosafety",
                ]
                if column in self.microbe_feature_table.columns
            ],
        ].copy()
        self.microbe_table["search_key"] = self.microbe_table["nt_code"].map(_canonicalize_key)
        self.microbe_table["name_key"] = self.microbe_table["microbe_label"].map(_canonicalize_key)
        self.microbe_table["species_key"] = self.microbe_table["species_label"].map(_canonicalize_key)
        self.amr_engine = AmrRuleEngine()
        self.step2_mechanism_projector = Step2MechanismProjector.from_joblib(self.step2_mechanism_reference_path)
        self.cohort_communities = self._discover_cohort_communities()
        self.demo_ranking = self._load_demo_ranking()
        self.disease_microbe_reference = self._load_disease_microbe_reference_bundle(
            primary_path=self.disease_microbe_reference_path,
            supplement_paths=self.disease_microbe_supplement_paths,
        )
        self.disease_drug_reference = self._load_optional_reference(self.disease_drug_reference_path)
        self._normalize_and_expand_disease_references()
        self.disease_catalog = self._build_disease_catalog()
        self.custom_sessions: dict[str, dict[str, object]] = {}

    def _load_optional_reference(self, path: Path | None) -> pd.DataFrame:
        """Load an optional CSV reference table, returning an empty frame when unavailable."""
        if path is None or not path.exists():
            return pd.DataFrame()
        return pd.read_csv(path, low_memory=False)

    def _load_disease_microbe_reference_bundle(
        self,
        primary_path: Path | None,
        supplement_paths: Sequence[Path],
    ) -> pd.DataFrame:
        """Load primary + supplemental disease-microbe references and perform lightweight de-duplication."""
        candidate_paths: list[Path] = []
        if primary_path is not None:
            candidate_paths.append(primary_path)
        candidate_paths.extend(path for path in supplement_paths if path is not None)

        frames: list[pd.DataFrame] = []
        for path in candidate_paths:
            frame = self._load_optional_reference(path)
            if frame.empty:
                continue
            work = frame.copy()
            if "source_database" not in work.columns:
                work["source_database"] = path.stem
            frames.append(work)

        if not frames:
            return pd.DataFrame()

        merged = pd.concat(frames, ignore_index=True, sort=False)
        dedup_columns = [
            column
            for column in [
                "source_sheet",
                "disease_name",
                "microbe_name_raw",
                "taxon_level",
                "desired_step1_effect",
                "mechanism_note",
            ]
            if column in merged.columns
        ]
        if dedup_columns:
            merged = merged.drop_duplicates(subset=dedup_columns).reset_index(drop=True)
        return merged

    def _canonicalize_disease_name(self, disease_name: object) -> str:
        text = str(disease_name or "").strip()
        if not text:
            return ""
        text_lower = text.lower()
        key_ascii = _canonicalize_key(text)
        key_cjk = _canonicalize_key(text, keep_cjk=True)
        if any(token in key_cjk for token in ["腹泻型肠易激"]) or any(
            token in key_ascii for token in ["diarrheapredominantibs", "ibsd"]
        ):
            return IBS_D_STANDARD_NAME
        if any(token in key_cjk for token in ["便秘型肠易激"]) or any(
            token in key_ascii for token in ["constipationpredominantibs", "ibsc"]
        ):
            return IBS_C_STANDARD_NAME
        if "ibs" in key_ascii or "肠易激" in key_cjk:
            return IBS_STANDARD_NAME
        if "crohn" in key_ascii or "克罗恩" in key_cjk:
            return CD_STANDARD_NAME
        if "ulcerativecolitis" in key_ascii or "溃疡性结肠炎" in key_cjk:
            return UC_STANDARD_NAME
        if any(token in key_ascii for token in ["colorectalcancer", "crc"]) or any(
            token in key_cjk for token in ["结直肠癌", "结肠癌"]
        ):
            return CRC_STANDARD_NAME
        if "便秘" in key_cjk or "constipation" in text_lower:
            return CONSTIPATION_STANDARD_NAME
        if "腹泻" in key_cjk or "diarrhea" in text_lower:
            return DIARRHEA_STANDARD_NAME
        return text

    def _disease_lookup_key(self, disease_name: object) -> str:
        text = str(disease_name or "").strip()
        canonical = _canonicalize_key(text, keep_cjk=True)
        if canonical:
            return canonical
        return text.lower()

    def _normalize_and_expand_disease_references(self) -> None:
        """Normalize disease naming and ensure IBS subtype entries exist for candidate generation."""
        if not self.disease_microbe_reference.empty and "disease_name" in self.disease_microbe_reference.columns:
            reference = self.disease_microbe_reference.copy()
            reference["disease_name"] = reference["disease_name"].map(self._canonicalize_disease_name)
            ibs_rows = reference[reference["disease_name"].map(_canonicalize_key).eq(_canonicalize_key(IBS_STANDARD_NAME))].copy()
            for subtype in [IBS_D_STANDARD_NAME, IBS_C_STANDARD_NAME]:
                subtype_key = _canonicalize_key(subtype)
                has_subtype = reference["disease_name"].map(_canonicalize_key).eq(subtype_key).any()
                if not has_subtype and not ibs_rows.empty:
                    clone = ibs_rows.copy()
                    clone["disease_name"] = subtype
                    if "mechanism_note" in clone.columns:
                        clone["mechanism_note"] = (
                            "Auto-expanded from IBS baseline to preserve subtype candidate coverage."
                        )
                    if "relation_confidence" in clone.columns:
                        clone["relation_confidence"] = clone["relation_confidence"].fillna("medium")
                    reference = pd.concat([reference, clone], ignore_index=True)
            self.disease_microbe_reference = reference

        if not self.disease_drug_reference.empty and "disease_name" in self.disease_drug_reference.columns:
            reference = self.disease_drug_reference.copy()
            reference["disease_name"] = reference["disease_name"].map(self._canonicalize_disease_name)
            if "marketed_drug_name_raw" in reference.columns:
                reference["marketed_drug_key"] = reference["marketed_drug_name_raw"].map(_canonicalize_key)
            self.disease_drug_reference = reference

    def _build_disease_catalog(self) -> list[dict[str, object]]:
        """Summarize the loaded disease references for bootstrap and UI dropdowns."""
        disease_names: set[str] = set()
        if not self.disease_microbe_reference.empty and "disease_name" in self.disease_microbe_reference.columns:
            disease_names |= {
                str(value).strip()
                for value in self.disease_microbe_reference["disease_name"].dropna().astype(str).tolist()
                if str(value).strip()
            }
        if not self.disease_drug_reference.empty and "disease_name" in self.disease_drug_reference.columns:
            disease_names |= {
                str(value).strip()
                for value in self.disease_drug_reference["disease_name"].dropna().astype(str).tolist()
                if str(value).strip()
            }
        disease_names |= set(REQUIRED_DISEASE_CATALOG_ENTRIES)

        catalog: list[dict[str, object]] = []
        for disease_name in sorted(disease_names):
            disease_key = self._disease_lookup_key(disease_name)
            microbe_count = 0
            if not self.disease_microbe_reference.empty:
                microbe_count = int(
                    self.disease_microbe_reference["disease_name"].map(self._disease_lookup_key).eq(disease_key).sum()
                )
            marketed_count = 0
            if not self.disease_drug_reference.empty:
                marketed_count = int(
                    self.disease_drug_reference["disease_name"].map(self._disease_lookup_key).eq(disease_key).sum()
                )
            catalog.append(
                {
                    "disease_name": disease_name,
                    "disease_key": disease_key,
                    "microbe_relation_count": microbe_count,
                    "marketed_drug_count": marketed_count,
                }
            )
        return catalog

    def _load_demo_ranking(self) -> list[dict[str, object]]:
        if self.demo_ranking_path is None or not self.demo_ranking_path.exists():
            return []
        ranking = pd.read_csv(self.demo_ranking_path, low_memory=False)
        keep_columns = [
            column
            for column in [
                "chemical_name",
                "prestwick_id",
                "scenario_name",
                "development_score",
                "final_health_index",
                "final_parent_retention_ratio",
            ]
            if column in ranking.columns
        ]
        ranking = ranking.loc[:, keep_columns].head(8).copy()
        return _clean_json(ranking)

    def _discover_cohort_communities(self) -> list[dict[str, object]]:
        if not self.step3_cohort_root.exists():
            return []
        communities: list[dict[str, object]] = []
        for path in sorted(self.step3_cohort_root.rglob("*_community.csv")):
            try:
                relative_path = path.relative_to(ROOT)
            except ValueError:
                relative_path = path
            dataset_name = path.parent.name
            sample_id = path.stem.removesuffix("_community")
            communities.append(
                {
                    "community_id": _canonicalize_key(str(relative_path)),
                    "label": f"{dataset_name} / {sample_id}",
                    "dataset_name": dataset_name,
                    "sample_id": sample_id,
                    "community_table_path": str(relative_path),
                }
            )
        return communities

    def _resolve_community_table_path(self, community_table_path: str | None) -> Path | None:
        if community_table_path is None or not str(community_table_path).strip():
            return None
        raw_value = str(community_table_path).strip()
        for item in self.cohort_communities:
            if raw_value == item["community_table_path"] or _canonicalize_key(raw_value) == item["community_id"]:
                candidate = ROOT / str(item["community_table_path"])
                if candidate.exists():
                    return candidate
        candidate = Path(raw_value).expanduser()
        if not candidate.is_absolute():
            candidate = ROOT / candidate
        if not candidate.exists() or not candidate.is_file():
            raise ValueError(f"找不到 community_table: {raw_value}")
        return candidate

    def _annotate_step2_mechanism(self, frame: pd.DataFrame) -> pd.DataFrame:
        return self.step2_mechanism_projector.annotate_frame(
            frame,
            predicted_probability_column="predicted_metabolized_probability",
            predicted_label_column="predicted_metabolism_label",
        )

    def _resolve_drug_id(self, drug_query: str) -> str:
        query_key = _canonicalize_key(drug_query)
        if not query_key:
            raise ValueError("drug_query 不能为空。")

        exact = self.drug_table[
            self.drug_table["search_key"].eq(query_key) | self.drug_table["name_key"].eq(query_key)
        ]
        if not exact.empty:
            return str(exact.iloc[0]["prestwick_id"])

        contains = self.drug_table[
            self.drug_table["chemical_name"].astype(str).str.contains(str(drug_query), case=False, na=False)
        ]
        if len(contains) == 1:
            return str(contains.iloc[0]["prestwick_id"])
        if len(contains) > 1:
            options = contains.loc[:, ["prestwick_id", "chemical_name"]].head(8).to_dict(orient="records")
            raise ValueError(f"drug_query 匹配到多个药物，请更具体一些。候选示例: {options}")
        raise ValueError(f"找不到药物: {drug_query}")

    def _resolve_microbe_id(self, microbe_query: str) -> str:
        query_key = _canonicalize_key(microbe_query)
        if not query_key:
            raise ValueError("microbe_query 不能为空。")

        exact = self.microbe_table[
            self.microbe_table["search_key"].eq(query_key)
            | self.microbe_table["name_key"].eq(query_key)
            | self.microbe_table["species_key"].eq(query_key)
        ]
        if not exact.empty:
            return str(exact.iloc[0]["nt_code"])

        contains = self.microbe_table[
            self.microbe_table["microbe_label"].astype(str).str.contains(str(microbe_query), case=False, na=False)
        ]
        if len(contains) == 1:
            return str(contains.iloc[0]["nt_code"])
        if len(contains) > 1:
            options = contains.loc[:, ["nt_code", "microbe_label"]].head(8).to_dict(orient="records")
            raise ValueError(f"microbe_query 匹配到多个微生物，请更具体一些。候选示例: {options}")
        raise ValueError(f"找不到微生物: {microbe_query}")

    def _drug_frame(self, drug_query: str) -> pd.DataFrame:
        drug_id = self._resolve_drug_id(drug_query)
        frame = self.frame[self.frame["prestwick_id"].eq(drug_id)].copy()
        if frame.empty:
            raise ValueError(f"找不到 drug_id={drug_id} 对应的预测结果。")
        return frame.reset_index(drop=True)

    def _custom_frame(self, session_id: str) -> pd.DataFrame:
        if session_id not in self.custom_sessions:
            raise ValueError(f"找不到 session_id={session_id} 的新药预测结果。")
        return self.custom_sessions[session_id]["frame"].copy()  # type: ignore[return-value]

    def _custom_session(self, session_id: str) -> dict[str, object]:
        if session_id not in self.custom_sessions:
            raise ValueError(f"找不到 session_id={session_id} 的新药预测结果。")
        return self.custom_sessions[session_id]

    def _pair_row(self, drug_query: str, microbe_query: str) -> pd.Series:
        frame = self._annotate_step2_mechanism(self._annotate_amr(self._drug_frame(drug_query)))
        microbe_id = self._resolve_microbe_id(microbe_query)
        matched = frame[frame["nt_code"].eq(microbe_id)]
        if matched.empty:
            raise ValueError(f"当前数据中不存在 microbe={microbe_id} 这个 pair。")
        return matched.iloc[0]

    def _pair_row_from_frame(self, frame: pd.DataFrame, microbe_query: str) -> pd.Series:
        frame = self._annotate_step2_mechanism(self._annotate_amr(frame))
        microbe_id = self._resolve_microbe_id(microbe_query)
        matched = frame[frame["nt_code"].eq(microbe_id)]
        if matched.empty:
            raise ValueError(f"当前结果中不存在 microbe={microbe_id} 这个 pair。")
        return matched.iloc[0]

    def _annotate_amr(self, frame: pd.DataFrame) -> pd.DataFrame:
        return self.amr_engine.annotate_frame(
            frame,
            label_column=self.step1_label_column,
            probability_column=self.step1_probability_column,
            score_column=self.step1_score_column,
        )

    def _drug_metadata(self, row: pd.Series) -> dict[str, object]:
        return _clean_json(
            {
                "prestwick_id": row.get("prestwick_id"),
                "chemical_name": row.get("chemical_name"),
                "therapeutic_class": row.get("therapeutic_class"),
                "therapeutic_effect": row.get("therapeutic_effect"),
                "atc_primary_l1": row.get("atc_primary_l1"),
                "atc_primary_l3": row.get("atc_primary_l3"),
                "atc_primary_l4": row.get("atc_primary_l4"),
                "molecular_formula": row.get("molecular_formula"),
                "molecular_weight": _safe_float(row.get("molecular_weight")),
                "xlogp": _safe_float(row.get("xlogp")),
                "tpsa": _safe_float(row.get("tpsa")),
                "murcko_scaffold": row.get("murcko_scaffold"),
                "canonical_smiles_rdkit": row.get("canonical_smiles_rdkit"),
                "smiles": row.get("smiles"),
            }
        )

    def _microbe_metadata(self, row: pd.Series) -> dict[str, object]:
        return _clean_json(
            {
                "nt_code": row.get("nt_code"),
                "microbe_label": row.get("microbe_label"),
                "species_label": row.get("species_label"),
                "species_name": row.get("species_name"),
                "genus": row.get("genus"),
                "family": row.get("family"),
                "phylum": row.get("phylum"),
                "gram_stain": row.get("gram_stain"),
                "medium_preference": row.get("medium_preference"),
                "biosafety": row.get("biosafety"),
            }
        )

    def _load_biotransform_reference_table(self) -> pd.DataFrame:
        reference_path = self.step2_biotransform_reference_path
        if reference_path is None or not reference_path.exists():
            return pd.DataFrame()

        reference = pd.read_csv(reference_path, low_memory=False)
        if reference.empty:
            return pd.DataFrame()

        has_annotation = pd.to_numeric(
            reference.get("experimental_biotransform_product_count", pd.Series(0, index=reference.index)),
            errors="coerce",
        ).fillna(0).gt(0)
        has_annotation = has_annotation | pd.to_numeric(
            reference.get("experimental_biotransform_fraction_in_gut", pd.Series(0.0, index=reference.index)),
            errors="coerce",
        ).fillna(0.0).gt(0)
        has_annotation = has_annotation | reference.get(
            "experimental_biotransform_product_ids",
            pd.Series("", index=reference.index),
        ).fillna("").astype(str).str.strip().ne("")
        reference = reference.loc[has_annotation].copy()
        if reference.empty:
            return pd.DataFrame()

        keep_columns = [
            column
            for column in [
                "prestwick_id",
                "chemical_name",
                "canonical_smiles_rdkit",
                "smiles",
                "experimental_biotransform_drug_name",
                "experimental_biotransform_product_count",
                "experimental_biotransform_product_ids",
                "experimental_biotransform_ec_numbers",
                "experimental_biotransform_reaction_centers",
                "experimental_biotransform_support_rows",
                "experimental_biotransform_summary_path",
                "experimental_biotransform_fraction_in_gut",
                "experimental_biotransform_primary_product_id",
                "experimental_biotransform_primary_product_name",
                "experimental_biotransform_drugbank_id",
            ]
            if column in reference.columns
        ]
        reference = (
            reference.loc[:, keep_columns]
            .drop_duplicates(subset=["prestwick_id"] if "prestwick_id" in reference.columns else None)
            .reset_index(drop=True)
        )
        if reference.empty:
            return pd.DataFrame()

        feature_source = reference.loc[
            :,
            [column for column in ["canonical_smiles_rdkit", "smiles"] if column in reference.columns],
        ].copy()
        feature_frame = enrich_drug_table_with_rdkit(feature_source, smiles_columns=list(feature_source.columns))
        feature_frame = feature_frame.loc[:, ~feature_frame.columns.duplicated()].copy()
        fp_columns = [column for column in feature_frame.columns if str(column).startswith("morgan_fp_")]
        if fp_columns:
            reference = pd.concat(
                [reference.reset_index(drop=True), feature_frame.loc[:, fp_columns].reset_index(drop=True)],
                axis=1,
            )
        reference["prestwick_key"] = reference.get("prestwick_id", pd.Series("", index=reference.index)).map(_canonicalize_key)
        reference["chemical_name_key"] = reference.get("chemical_name", pd.Series("", index=reference.index)).map(_canonicalize_key)
        reference["experimental_biotransform_product_count"] = pd.to_numeric(
            reference.get("experimental_biotransform_product_count"),
            errors="coerce",
        ).fillna(0.0)
        reference["experimental_biotransform_fraction_in_gut"] = pd.to_numeric(
            reference.get("experimental_biotransform_fraction_in_gut"),
            errors="coerce",
        ).fillna(0.0)
        return reference

    def _biotransform_sidecar_from_row(
        self,
        row: pd.Series,
        *,
        top_k: int = BIOTRANSFORM_TOP_K,
        min_similarity: float = BIOTRANSFORM_MIN_THRESHOLD,
    ) -> dict[str, object]:
        if self.biotransform_reference_table.empty or not self.fingerprint_columns:
            return {
                "biotransform_sidecar_enabled": False,
                "biotransform_sidecar_mode": "unavailable",
                "biotransform_sidecar_reason": "no_reference_annotations",
            }

        query_bits = _fingerprint_array_from_row(row, self.fingerprint_columns)
        if query_bits.size == 0 or not query_bits.any():
            smiles_columns = [column for column in ["canonical_smiles_rdkit", "smiles"] if column in row.index]
            if smiles_columns:
                feature_source = pd.DataFrame([{column: row.get(column) for column in smiles_columns}])
                feature_frame = enrich_drug_table_with_rdkit(feature_source, smiles_columns=smiles_columns)
                feature_frame = feature_frame.loc[:, ~feature_frame.columns.duplicated()].copy()
                query_bits = _fingerprint_array_from_row(feature_frame.iloc[0], self.fingerprint_columns)
        if query_bits.size == 0 or not query_bits.any():
            return {
                "biotransform_sidecar_enabled": False,
                "biotransform_sidecar_mode": "unavailable",
                "biotransform_sidecar_reason": "query_fingerprint_missing",
            }

        ranked = self.biotransform_reference_table.copy()
        ranked["biotransform_similarity"] = _tanimoto_similarity_vector(query_bits, self.biotransform_reference_matrix)
        query_prestwick_key = _canonicalize_key(row.get("prestwick_id"))
        query_name_key = _canonicalize_key(row.get("chemical_name"))
        direct_ranked = ranked.copy()
        if query_prestwick_key:
            direct_ranked = direct_ranked[
                direct_ranked.get("prestwick_key", pd.Series("", index=direct_ranked.index)).astype(str).eq(query_prestwick_key)
            ].copy()
        elif query_name_key:
            direct_ranked = direct_ranked[
                direct_ranked.get("chemical_name_key", pd.Series("", index=direct_ranked.index)).astype(str).eq(query_name_key)
            ].copy()
        if not direct_ranked.empty:
            ranked = direct_ranked.sort_values(["biotransform_similarity", "chemical_name"], ascending=[False, True]).head(top_k)
            direct_mode = True
        else:
            direct_mode = False
        if not direct_mode:
            if query_prestwick_key:
                ranked = ranked[
                    ranked.get("prestwick_key", pd.Series("", index=ranked.index)).astype(str).ne(query_prestwick_key)
                ].copy()
            if query_name_key:
                ranked = ranked[
                    ranked.get("chemical_name_key", pd.Series("", index=ranked.index)).astype(str).ne(query_name_key)
                ].copy()
            ranked = ranked[
                pd.to_numeric(ranked["biotransform_similarity"], errors="coerce").fillna(0.0).ge(float(min_similarity))
            ].copy()
        if ranked.empty:
            return {
                "biotransform_sidecar_enabled": False,
                "biotransform_sidecar_mode": "similarity_transfer",
                "biotransform_sidecar_reason": "no_similar_annotated_drug",
                "biotransform_sidecar_support_score": 0.0,
            }

        ranked = ranked.sort_values(["biotransform_similarity", "chemical_name"], ascending=[False, True]).head(top_k)
        similarities = pd.to_numeric(ranked["biotransform_similarity"], errors="coerce").fillna(0.0).to_numpy(dtype=float)
        weights = similarities if float(similarities.sum()) > 0 else np.ones_like(similarities)
        weighted_product_count = float(
            np.average(
                pd.to_numeric(ranked["experimental_biotransform_product_count"], errors="coerce").fillna(0.0),
                weights=weights,
            )
        )
        weighted_fraction_in_gut = float(
            np.average(
                pd.to_numeric(ranked["experimental_biotransform_fraction_in_gut"], errors="coerce").fillna(0.0),
                weights=weights,
            )
        )
        max_similarity = float(similarities.max()) if similarities.size else 0.0
        support_score = float(
            np.clip(
                0.65 * max_similarity
                + 0.20 * min(weighted_fraction_in_gut, 1.0)
                + 0.15 * min(weighted_product_count / 5.0, 1.0),
                0.0,
                1.0,
            )
        )
        if max_similarity >= BIOTRANSFORM_DIRECT_MATCH_THRESHOLD:
            confidence = "high"
        elif max_similarity >= 0.45:
            confidence = "medium"
        else:
            confidence = "low"

        reference_drugs = ranked["chemical_name"].fillna("").astype(str).tolist()
        product_ids: list[str] = []
        ec_numbers: list[str] = []
        reaction_centers: list[str] = []
        primary_products: list[str] = []
        for _, reference_row in ranked.iterrows():
            product_ids.extend(_string_list_from_semicolon(reference_row.get("experimental_biotransform_product_ids")))
            ec_numbers.extend(_string_list_from_semicolon(reference_row.get("experimental_biotransform_ec_numbers")))
            reaction_centers.extend(_string_list_from_semicolon(reference_row.get("experimental_biotransform_reaction_centers")))
            primary_name = str(reference_row.get("experimental_biotransform_primary_product_name") or "").strip()
            if primary_name:
                primary_products.append(primary_name)

        def _dedupe(values: list[str], max_items: int = 8) -> list[str]:
            ordered: list[str] = []
            for value in values:
                if value and value not in ordered:
                    ordered.append(value)
                if len(ordered) >= max_items:
                    break
            return ordered

        return {
            "biotransform_sidecar_enabled": True,
            "biotransform_sidecar_mode": "direct_reference" if direct_mode else "similarity_transfer",
            "biotransform_sidecar_reason": "matched_direct_reference" if direct_mode else "matched_similar_annotated_drugs",
            "biotransform_sidecar_support_score": round(support_score, 4),
            "biotransform_sidecar_similarity": round(max_similarity, 4),
            "biotransform_sidecar_confidence": confidence,
            "biotransform_sidecar_reference_count": int(len(ranked)),
            "biotransform_sidecar_weighted_product_count": round(weighted_product_count, 3),
            "biotransform_sidecar_weighted_fraction_in_gut": round(weighted_fraction_in_gut, 4),
            "biotransform_sidecar_reference_drugs": ";".join(_dedupe(reference_drugs, max_items=top_k)),
            "biotransform_sidecar_reference_product_ids": ";".join(_dedupe(product_ids)),
            "biotransform_sidecar_reference_ec_numbers": ";".join(_dedupe(ec_numbers)),
            "biotransform_sidecar_reference_reaction_centers": ";".join(_dedupe(reaction_centers)),
            "biotransform_sidecar_primary_products": ";".join(_dedupe(primary_products)),
        }

    def _annotate_biotransform_sidecar(self, frame: pd.DataFrame) -> pd.DataFrame:
        if frame.empty:
            return frame.copy()

        key_columns = [column for column in ["prestwick_id", "chemical_name"] if column in frame.columns]
        if not key_columns:
            return frame.copy()

        unique_drugs = frame.loc[:, key_columns + [column for column in self.fingerprint_columns if column in frame.columns] + [
            column for column in ["canonical_smiles_rdkit", "smiles"] if column in frame.columns
        ]].drop_duplicates(subset=["prestwick_id"] if "prestwick_id" in frame.columns else ["chemical_name"])
        annotations: list[dict[str, object]] = []
        for _, drug_row in unique_drugs.iterrows():
            payload = self._biotransform_sidecar_from_row(drug_row)
            entry = {column: drug_row.get(column) for column in key_columns}
            entry.update(payload)
            annotations.append(entry)

        if not annotations:
            return frame.copy()

        annotation_frame = pd.DataFrame(annotations)
        merge_columns = ["prestwick_id"] if "prestwick_id" in frame.columns and "prestwick_id" in annotation_frame.columns else ["chemical_name"]
        return frame.merge(annotation_frame, on=merge_columns, how="left")

    def _top_similar_library_drugs(
        self,
        row: pd.Series,
        *,
        top_k: int = SIMILARITY_DRUG_TOP_K,
        min_similarity: float = SIMILARITY_MIN_THRESHOLD,
    ) -> list[dict[str, object]]:
        """Return the nearest library drugs for one query molecule using RDKit Morgan fingerprints."""
        if self.drug_similarity_table.empty or not self.fingerprint_columns:
            return []

        query_bits = _fingerprint_array_from_row(row, self.fingerprint_columns)
        if query_bits.size == 0 or not query_bits.any():
            return []

        similarities = _tanimoto_similarity_vector(query_bits, self.drug_similarity_matrix)
        ranked = self.drug_similarity_table.copy()
        ranked["tanimoto_similarity"] = similarities

        query_prestwick = str(row.get("prestwick_id", "")).strip()
        query_name_key = _canonicalize_key(row.get("chemical_name"))
        if query_prestwick:
            ranked = ranked[ranked["prestwick_id"].astype(str).ne(query_prestwick)].copy()
        if query_name_key:
            ranked = ranked[ranked["chemical_name"].map(_canonicalize_key).ne(query_name_key)].copy()

        ranked = ranked[
            pd.to_numeric(ranked["tanimoto_similarity"], errors="coerce").fillna(0.0).ge(float(min_similarity))
        ].copy()
        if ranked.empty:
            return []

        ranked = ranked.sort_values(["tanimoto_similarity", "chemical_name"], ascending=[False, True]).head(top_k)
        records: list[dict[str, object]] = []
        for _, similar_row in ranked.iterrows():
            drug_key = _canonicalize_key(similar_row.get("chemical_name"))
            matched_diseases = []
            if not self.disease_drug_reference.empty:
                matched = self.disease_drug_reference[
                    self.disease_drug_reference["marketed_drug_key"].map(str).eq(drug_key)
                ].copy()
                matched_diseases = (
                    matched["disease_name"].dropna().astype(str).drop_duplicates().sort_values().tolist()
                    if not matched.empty
                    else []
                )
            records.append(
                {
                    "prestwick_id": similar_row.get("prestwick_id"),
                    "chemical_name": similar_row.get("chemical_name"),
                    "tanimoto_similarity": round(float(similar_row.get("tanimoto_similarity", 0.0)), 4),
                    "murcko_scaffold": similar_row.get("murcko_scaffold"),
                    "matched_diseases": matched_diseases,
                }
            )
        return records

    def _similarity_disease_context(self, row: pd.Series) -> list[dict[str, object]]:
        """Aggregate disease hints from structurally similar marketed library drugs."""
        similar_drugs = self._top_similar_library_drugs(row)
        if not similar_drugs or self.disease_drug_reference.empty:
            return []

        disease_records: dict[str, dict[str, object]] = {}
        total_similarity = 0.0
        for similar_drug in similar_drugs:
            similarity = float(similar_drug.get("tanimoto_similarity", 0.0) or 0.0)
            if similarity <= 0:
                continue
            drug_key = _canonicalize_key(similar_drug.get("chemical_name"))
            matched = self.disease_drug_reference[
                self.disease_drug_reference["marketed_drug_key"].map(str).eq(drug_key)
            ].copy()
            if matched.empty:
                continue
            total_similarity += similarity
            for disease_name, group in matched.groupby("disease_name", dropna=False):
                disease_name_str = str(disease_name).strip()
                if not disease_name_str:
                    continue
                record = disease_records.setdefault(
                    disease_name_str,
                    {
                        "disease_name": disease_name_str,
                        "weighted_similarity_sum": 0.0,
                        "max_similarity": 0.0,
                        "matched_drug_count": 0,
                        "matched_market_drugs": [],
                    },
                )
                record["weighted_similarity_sum"] = float(record["weighted_similarity_sum"] + similarity)
                record["max_similarity"] = float(max(float(record["max_similarity"]), similarity))
                record["matched_drug_count"] = int(record["matched_drug_count"] + 1)
                record["matched_market_drugs"] = sorted(
                    {
                        *record["matched_market_drugs"],
                        *group["marketed_drug_name_raw"].dropna().astype(str).tolist(),
                    }
                )[:5]

        if not disease_records:
            return []

        denominator = total_similarity if total_similarity > 0 else 1.0
        ranked = []
        for record in disease_records.values():
            weighted_similarity = float(record["weighted_similarity_sum"])
            ranked.append(
                {
                    "disease_name": record["disease_name"],
                    "support_score": round(float(min(1.0, weighted_similarity / denominator)), 4),
                    "weighted_similarity_sum": round(weighted_similarity, 4),
                    "max_similarity": round(float(record["max_similarity"]), 4),
                    "matched_drug_count": int(record["matched_drug_count"]),
                    "matched_market_drugs": record["matched_market_drugs"],
                }
            )
        ranked.sort(
            key=lambda item: (
                float(item.get("support_score", 0.0)),
                float(item.get("max_similarity", 0.0)),
                int(item.get("matched_drug_count", 0)),
            ),
            reverse=True,
        )
        return ranked[:SIMILARITY_DISEASE_TOP_K]

    def _marketed_disease_context(self, drug_name: str | None) -> list[dict[str, object]]:
        """Return marketed-disease context rows matching the current drug name when possible."""
        if not drug_name or self.disease_drug_reference.empty:
            return []
        drug_key = _canonicalize_key(drug_name)
        if not drug_key:
            return []

        matched = self.disease_drug_reference[
            self.disease_drug_reference["marketed_drug_key"].map(str).eq(drug_key)
        ].copy()
        if matched.empty:
            matched = self.disease_drug_reference[
                self.disease_drug_reference["marketed_drug_name_raw"]
                .astype(str)
                .str.contains(str(drug_name), case=False, na=False, regex=False)
            ].copy()
        if matched.empty:
            return []

        grouped = (
            matched.groupby("disease_name", dropna=False)["marketed_drug_name_raw"]
            .apply(lambda series: sorted({str(value).strip() for value in series if str(value).strip()}))
            .reset_index()
        )
        return [
            {
                "disease_name": row["disease_name"],
                "matched_market_drugs": row["marketed_drug_name_raw"][:5],
            }
            for _, row in grouped.iterrows()
        ]

    def _match_disease_relation(self, work: pd.DataFrame, relation: pd.Series) -> pd.Series:
        """Match one disease relation row onto the currently available microbe panel."""
        taxon_level = str(relation.get("taxon_level", "unknown"))
        microbe_key = _canonicalize_key(relation.get("microbe_name_raw"))
        genus_hint = _canonicalize_key(relation.get("genus_hint"))
        if taxon_level == "species":
            species_keys = work["species_label"].map(_canonicalize_key).fillna("")
            microbe_label_keys = work["microbe_label"].map(_canonicalize_key).fillna("")
            base_species_key = _species_base_key(relation.get("microbe_name_raw"))
            exact_mask = species_keys.eq(microbe_key) | microbe_label_keys.eq(microbe_key)
            if not base_species_key:
                return exact_mask

            # Species aliases in curated references often carry suffixes like "(AIEC)" or "产毒株".
            # We tolerate these by matching the canonical binomial base key in both directions.
            tolerant_mask = (
                species_keys.eq(base_species_key)
                | microbe_label_keys.eq(base_species_key)
                | species_keys.str.contains(base_species_key, case=False, regex=False, na=False)
                | microbe_label_keys.str.contains(base_species_key, case=False, regex=False, na=False)
                | species_keys.map(lambda item: bool(item) and item in base_species_key)
                | microbe_label_keys.map(lambda item: bool(item) and item in base_species_key)
            )
            return exact_mask | tolerant_mask
        if taxon_level == "genus":
            return work["genus"].map(_canonicalize_key).eq(genus_hint or microbe_key)
        if taxon_level == "family":
            return work["family"].map(_canonicalize_key).eq(microbe_key)
        if taxon_level == "phylum":
            return work["phylum"].map(_canonicalize_key).eq(microbe_key)
        if taxon_level == "class":
            return work["class"].map(_canonicalize_key).eq(microbe_key) if "class" in work.columns else pd.Series(False, index=work.index)
        if taxon_level == "order":
            return work["order"].map(_canonicalize_key).eq(microbe_key) if "order" in work.columns else pd.Series(False, index=work.index)
        return work["genus"].map(_canonicalize_key).eq(genus_hint) if genus_hint else pd.Series(False, index=work.index)

    def _relation_dedup_key(self, relation: pd.Series, desired_effect: str) -> tuple[str, str, str]:
        taxon_level = str(relation.get("taxon_level", "unknown"))
        if taxon_level == "species":
            microbe_key = _species_base_key(relation.get("microbe_name_raw"))
        elif taxon_level == "genus":
            microbe_key = _canonicalize_key(relation.get("genus_hint")) or _canonicalize_key(relation.get("microbe_name_raw"))
        else:
            microbe_key = _canonicalize_key(relation.get("microbe_name_raw"))
        return desired_effect, taxon_level, microbe_key

    def _relation_weight(self, relation: pd.Series, *, ibs_like: bool, desired_effect: str) -> float:
        source_sheet = str(relation.get("source_sheet", ""))
        relation_confidence = _canonicalize_key(relation.get("relation_confidence"))
        weight = (
            DISEASE_RELATION_LEVEL_WEIGHTS.get(str(relation.get("taxon_level", "unknown")), 0.2)
            * DISEASE_RELATION_SOURCE_WEIGHTS.get(source_sheet, 0.5)
            * DISEASE_RELATION_CONFIDENCE_WEIGHTS.get(relation_confidence, 0.8)
        )
        if ibs_like and source_sheet == "disease_to_microbe":
            # For IBS-like disorders, curated symptom-context relations are predominantly disease_to_microbe.
            # Raise their contribution slightly so host-symptom evidence is not overwhelmed by inflammatory priors.
            weight *= 1.35
            if desired_effect == "promote":
                weight *= 1.10
        return float(weight)

    def _unseen_smiles_reliability_context(self, work: pd.DataFrame) -> dict[str, float]:
        applicability_rate = _mean_boolean_like(
            work.get("applicability_flag", pd.Series(np.nan, index=work.index)),
            default=1.0,
        )
        scaffold_seen_rate = _mean_boolean_like(
            work.get("scaffold_seen_in_training", pd.Series(np.nan, index=work.index)),
            default=1.0,
        )
        jaccard_series = pd.to_numeric(
            work.get("drug_max_fingerprint_jaccard", pd.Series(np.nan, index=work.index)),
            errors="coerce",
        ).dropna()
        jaccard_mean = float(jaccard_series.mean()) if not jaccard_series.empty else 1.0
        similarity_reliability = _clip01((jaccard_mean - 0.15) / 0.65)
        global_reliability = _clip01(0.50 * applicability_rate + 0.30 * similarity_reliability + 0.20 * scaffold_seen_rate)
        score_multiplier = float(0.55 + 0.45 * global_reliability)
        return {
            "applicability_rate": float(applicability_rate),
            "scaffold_seen_rate": float(scaffold_seen_rate),
            "fingerprint_jaccard_mean": float(jaccard_mean),
            "similarity_reliability": float(similarity_reliability),
            "global_reliability": float(global_reliability),
            "score_multiplier": score_multiplier,
        }

    def _candidate_diseases_from_frame(self, work: pd.DataFrame) -> list[dict[str, object]]:
        """Score diseases whose curated microbe patterns are directionally consistent with Step 1 outputs."""
        if self.disease_microbe_reference.empty:
            return []

        disease_records: list[dict[str, object]] = []
        unseen_smiles_context = self._unseen_smiles_reliability_context(work)
        score_series = pd.to_numeric(
            work.get("display_step1_predicted_effect_score", work.get(self.step1_score_column, pd.Series(np.nan, index=work.index))),
            errors="coerce",
        )
        inhibit_series = pd.to_numeric(
            work.get(
                "display_step1_predicted_inhibit_probability",
                work.get(self.step1_probability_column, pd.Series(np.nan, index=work.index)),
            ),
            errors="coerce",
        )
        promote_series = pd.to_numeric(
            work.get("predicted_promote_probability_refined", pd.Series(np.nan, index=work.index)),
            errors="coerce",
        )
        mechanism_work = attach_action_signals(
            infer_microbe_trait_priors(work),
            score_column="display_step1_predicted_effect_score"
            if "display_step1_predicted_effect_score" in work.columns
            else (self.step1_score_column or "predicted_effect_score"),
            inhibit_probability_column="display_step1_predicted_inhibit_probability"
            if "display_step1_predicted_inhibit_probability" in work.columns
            else (self.step1_probability_column or "predicted_inhibit_probability"),
            promote_probability_column="predicted_promote_probability_refined",
        )

        for disease_name, group in self.disease_microbe_reference.groupby("disease_name", dropna=False):
            disease_name_str = str(disease_name)
            ibs_like = _is_ibs_like_disease(disease_name_str)
            relation_scores: list[float] = []
            evidence_rows: list[dict[str, object]] = []
            mechanism_rows: list[pd.DataFrame] = []
            deduplicated_relations: dict[tuple[str, str, str], tuple[pd.Series, float]] = {}
            for _, relation in group.iterrows():
                desired_effect = str(relation.get("desired_step1_effect", "unknown"))
                if desired_effect not in {"promote", "inhibit"}:
                    continue
                dedup_key = self._relation_dedup_key(relation, desired_effect)
                if not dedup_key[2]:
                    continue
                relation_weight = self._relation_weight(relation, ibs_like=ibs_like, desired_effect=desired_effect)
                existing = deduplicated_relations.get(dedup_key)
                if existing is None or relation_weight > existing[1]:
                    deduplicated_relations[dedup_key] = (relation, relation_weight)

            for relation, weight in deduplicated_relations.values():
                desired_effect = str(relation.get("desired_step1_effect", "unknown"))
                mask = self._match_disease_relation(mechanism_work, relation)
                matched = mechanism_work.loc[mask].copy()
                if matched.empty:
                    continue
                matched_scores = score_series.loc[matched.index].fillna(0.0)
                matched_inhibit = inhibit_series.loc[matched.index].fillna(0.0)
                matched_promote = promote_series.loc[matched.index].fillna(0.0)

                if desired_effect == "promote":
                    local_support = matched_scores.clip(lower=0.0) + 0.75 * matched_promote
                else:
                    local_support = (-matched_scores).clip(lower=0.0) + 0.75 * matched_inhibit
                relation_score = float(local_support.mean())
                local_applicability = _mean_boolean_like(
                    matched.get("applicability_flag", pd.Series(np.nan, index=matched.index)),
                    default=unseen_smiles_context["applicability_rate"],
                )
                local_reliability = _clip01(0.70 * local_applicability + 0.30 * unseen_smiles_context["global_reliability"])
                reliability_multiplier = float(0.55 + 0.45 * local_reliability)
                weighted_score = relation_score * weight * reliability_multiplier
                relation_scores.append(weighted_score)
                matched = matched.copy()
                matched["relation_weight"] = float(weight)
                matched["desired_step1_effect"] = desired_effect
                mechanism_rows.append(matched)
                strongest_idx = local_support.idxmax()
                strongest = matched.loc[strongest_idx]
                evidence_rows.append(
                    {
                        "microbe_name_raw": relation.get("microbe_name_raw"),
                        "desired_step1_effect": desired_effect,
                        "matched_microbe": strongest.get("microbe_label"),
                        "matched_effect_label": strongest.get("display_step1_predicted_effect_label", strongest.get(self.step1_label_column)),
                        "matched_effect_score": _safe_float(strongest.get("display_step1_predicted_effect_score", strongest.get(self.step1_score_column))),
                        "matched_promote_probability": _safe_float(strongest.get("predicted_promote_probability_refined")),
                        "matched_inhibit_probability": _safe_float(strongest.get("display_step1_predicted_inhibit_probability", strongest.get(self.step1_probability_column))),
                        "relation_score": round(weighted_score, 4),
                        "reliability_multiplier": round(reliability_multiplier, 4),
                        "taxon_level": relation.get("taxon_level"),
                        "source_sheet": relation.get("source_sheet"),
                    }
                )

            if not relation_scores:
                continue

            evidence_deduplicated: dict[tuple[str, str, str, str], dict[str, object]] = {}
            for item in evidence_rows:
                evidence_key = (
                    _canonicalize_key(item.get("microbe_name_raw")),
                    _canonicalize_key(item.get("matched_microbe")),
                    str(item.get("desired_step1_effect")),
                    str(item.get("taxon_level")),
                )
                if evidence_key not in evidence_deduplicated or float(item.get("relation_score", 0.0)) > float(
                    evidence_deduplicated[evidence_key].get("relation_score", 0.0)
                ):
                    evidence_deduplicated[evidence_key] = item
            evidence_rows = sorted(
                evidence_deduplicated.values(),
                key=lambda item: float(item.get("relation_score", 0.0)),
                reverse=True,
            )

            if mechanism_rows:
                mechanism_result = compute_mechanism_layer(
                    pd.concat(mechanism_rows, ignore_index=True),
                    relation_weight_column="relation_weight",
                    top_n_contributors=4,
                )
            else:
                mechanism_result = compute_mechanism_layer(pd.DataFrame())

            raw_microbe_score = float(np.mean(relation_scores))
            mechanism_scores = mechanism_result["scores"]
            mechanism_balance = float(mechanism_scores.get("mechanism_balance_score", 0.0))
            if ibs_like and mechanism_balance > 0:
                # Mild mechanism emphasis for IBS-like ranking: prioritize barrier/butyrate-compatible signals.
                mechanism_balance = float(mechanism_balance * 1.15)
                mechanism_scores = dict(mechanism_scores)
                mechanism_scores["mechanism_balance_score"] = mechanism_balance
            default_fusion_mode = "weighted_0.65_0.35"
            disease_score_mechanism = fuse_disease_scores(
                raw_score=raw_microbe_score,
                mechanism_score=mechanism_balance,
                fusion_mode=default_fusion_mode,
            )
            marketed_examples = []
            if not self.disease_drug_reference.empty:
                disease_key = self._disease_lookup_key(disease_name)
                marketed_examples = (
                    self.disease_drug_reference.loc[
                        self.disease_drug_reference["disease_name"].map(self._disease_lookup_key).eq(disease_key),
                        "marketed_drug_name_raw",
                    ]
                    .dropna()
                    .astype(str)
                    .head(5)
                    .tolist()
                )
            disease_records.append(
                {
                    "disease_name": disease_name_str,
                    "support_score": round(disease_score_mechanism, 4),
                    "disease_score_mechanism": round(disease_score_mechanism, 4),
                    "disease_score_raw_only": round(raw_microbe_score, 4),
                    "mechanism_delta": round(float(disease_score_mechanism - raw_microbe_score), 4),
                    "mechanism_scores": mechanism_scores,
                    "mechanism_top_contributors": mechanism_result["top_contributors"],
                    "mechanism_model_version": "v1_literature_minimal_2026_04",
                    "fusion_mode": default_fusion_mode,
                    "matched_relation_count": int(len(relation_scores)),
                    "matched_microbe_count": int(
                        pd.concat(mechanism_rows, ignore_index=True)["nt_code"].nunique() if mechanism_rows else 0
                    ),
                    "unseen_smiles_reliability": unseen_smiles_context,
                    "marketed_drug_examples": marketed_examples,
                    "evidence_examples": evidence_rows[:5],
                }
            )

        deduplicated_records: dict[str, dict[str, object]] = {}
        for record in disease_records:
            disease_key = self._disease_lookup_key(record.get("disease_name"))
            existing = deduplicated_records.get(disease_key)
            if existing is None:
                deduplicated_records[disease_key] = record
                continue
            current_rank = (
                float(record.get("disease_score_mechanism", 0.0)),
                float(record.get("disease_score_raw_only", 0.0)),
                int(record.get("matched_relation_count", 0)),
            )
            existing_rank = (
                float(existing.get("disease_score_mechanism", 0.0)),
                float(existing.get("disease_score_raw_only", 0.0)),
                int(existing.get("matched_relation_count", 0)),
            )
            if current_rank > existing_rank:
                deduplicated_records[disease_key] = record
        disease_records = list(deduplicated_records.values())

        disease_records.sort(
            key=lambda item: (
                item["disease_score_mechanism"],
                item["disease_score_raw_only"],
                item["matched_relation_count"],
            ),
            reverse=True,
        )
        top_records = list(disease_records[:8])
        selected_keys = {_canonicalize_key(item["disease_name"]) for item in top_records}
        for required in REQUIRED_DISEASE_CATALOG_ENTRIES:
            required_key = _canonicalize_key(required)
            if required_key in selected_keys:
                continue
            matched = next((item for item in disease_records if _canonicalize_key(item["disease_name"]) == required_key), None)
            if matched is not None:
                top_records.append(matched)
                selected_keys.add(required_key)
        for required in BENCHMARK_DISEASE_PANEL_ENTRIES:
            required_key = _canonicalize_key(required)
            if required_key in selected_keys:
                continue
            matched = next((item for item in disease_records if _canonicalize_key(item["disease_name"]) == required_key), None)
            if matched is not None:
                top_records.append(matched)
                selected_keys.add(required_key)
                continue
            top_records.append(
                {
                    "disease_name": required,
                    "support_score": 0.0,
                    "disease_score_mechanism": 0.0,
                    "disease_score_raw_only": 0.0,
                    "mechanism_delta": 0.0,
                    "mechanism_scores": {
                        "anti_inflammatory_score": 0.0,
                        "pro_inflammatory_score": 0.0,
                        "butyrate_support_score": 0.0,
                        "barrier_protection_score": 0.0,
                        "toxin_risk_score": 0.0,
                        "mucus_degradation_score": 0.0,
                        "pathobiont_load": 0.0,
                        "competition_vs_crossfeeding_proxy": 0.0,
                        "mechanism_benefit_score": 0.0,
                        "mechanism_risk_score": 0.0,
                        "mechanism_balance_score": 0.0,
                    },
                    "mechanism_top_contributors": [],
                    "mechanism_model_version": "v1_literature_minimal_2026_04",
                    "fusion_mode": "weighted_0.65_0.35",
                    "matched_relation_count": 0,
                    "matched_microbe_count": 0,
                    "marketed_drug_examples": [],
                    "evidence_examples": [],
                    "coverage_note": "forced_benchmark_panel_entry_without_matched_relations",
                }
            )
            selected_keys.add(required_key)
        return top_records

    def _infer_drug_profile_from_frame(self, work: pd.DataFrame, row: pd.Series) -> str:
        if "step1_drug_profile" in work.columns:
            observed = work["step1_drug_profile"].dropna().astype(str)
            if not observed.empty:
                top = observed.value_counts().idxmax()
                if str(top).strip():
                    return str(top).strip()

        name_key = _canonicalize_key(row.get("chemical_name"))
        if "rifaximin" in name_key:
            return "eubiotic_modulator"
        if "vancomycin" in name_key:
            return "disruptive_antibiotic"
        if "lubiprostone" in name_key:
            return "host_secretagogue"
        if "metronidazole" in name_key:
            return "contextual_antimicrobial"
        if (
            "sulfasalazine" in name_key
            or "sasp" in name_key
            or "sulfapyridine" in name_key
            or "sulfonamide" in name_key
            or "sulfamethoxazole" in name_key
            or "sulfadiazine" in name_key
            or "sulfisoxazole" in name_key
            or "cotrimoxazole" in name_key
            or "trimethoprimsulfamethoxazole" in name_key
        ):
            return "sulfonamide_antifolate"

        cls_key = _canonicalize_key(row.get("therapeutic_class"))
        if "sulfonamide" in cls_key or "antifolate" in cls_key:
            return "sulfonamide_antifolate"
        if "antibiotic" in cls_key:
            return "disruptive_antibiotic"
        return "unknown"

    def _write_disease_adjusted_community(self, output_path: Path, disease_name: str) -> Path:
        """Create a temporary community table reflecting one curated disease profile."""
        if self.disease_microbe_reference.empty:
            raise ValueError("当前没有可用的疾病-微生物参考表。")
        canonical_disease_name = self._canonicalize_disease_name(disease_name)
        community = build_disease_adjusted_community(
            microbe_metadata=self.microbe_feature_table,
            disease_name=canonical_disease_name,
            disease_microbe_reference=self.disease_microbe_reference,
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        community.to_csv(output_path, index=False)
        return output_path

    def _build_disease_target_profile_from_frame(
        self,
        work: pd.DataFrame,
        disease_name: str | None,
    ) -> tuple[dict[str, float], dict[str, object]]:
        if self.disease_microbe_reference.empty or not disease_name or not str(disease_name).strip():
            return {}, {"enabled": False, "reason": "empty_reference_or_missing_disease"}
        canonical_disease_name = self._canonicalize_disease_name(disease_name)
        disease_key = self._disease_lookup_key(canonical_disease_name)
        reference = self.disease_microbe_reference[
            self.disease_microbe_reference["disease_name"].map(self._disease_lookup_key).eq(disease_key)
        ].copy()
        if reference.empty:
            return {}, {
                "enabled": False,
                "reason": "disease_not_found",
                "requested_disease_name": str(disease_name),
                "canonical_disease_name": canonical_disease_name,
            }
        ibs_like = _is_ibs_like_disease(canonical_disease_name)
        deduplicated_relations: dict[tuple[str, str, str], tuple[pd.Series, float]] = {}
        for _, relation in reference.iterrows():
            desired_effect = str(relation.get("desired_step1_effect", "unknown"))
            if desired_effect not in {"promote", "inhibit"}:
                continue
            dedup_key = self._relation_dedup_key(relation, desired_effect)
            if not dedup_key[2]:
                continue
            relation_weight = self._relation_weight(relation, ibs_like=ibs_like, desired_effect=desired_effect)
            existing = deduplicated_relations.get(dedup_key)
            if existing is None or relation_weight > existing[1]:
                deduplicated_relations[dedup_key] = (relation, relation_weight)

        target_weights: dict[str, float] = {}
        matched_relation_count = 0
        for relation, relation_weight in deduplicated_relations.values():
            desired_effect = str(relation.get("desired_step1_effect", "unknown"))
            desired_sign = 1.0 if desired_effect == "promote" else -1.0
            matched_nt_codes = (
                work.loc[self._match_disease_relation(work, relation), "nt_code"].dropna().astype(str).unique().tolist()
            )
            if not matched_nt_codes:
                continue
            matched_relation_count += 1
            per_nt_weight = float(relation_weight / max(len(matched_nt_codes), 1))
            for nt_code in matched_nt_codes:
                target_weights[nt_code] = float(target_weights.get(nt_code, 0.0) + desired_sign * per_nt_weight)
        l1 = float(sum(abs(value) for value in target_weights.values()))
        if l1 <= 0:
            return {}, {
                "enabled": False,
                "reason": "no_matched_microbes",
                "requested_disease_name": str(disease_name),
                "canonical_disease_name": canonical_disease_name,
                "relation_count": int(len(reference)),
                "matched_relation_count": int(matched_relation_count),
            }
        normalized_profile = {key: float(value / l1) for key, value in target_weights.items() if abs(value) > 1e-12}
        return normalized_profile, {
            "enabled": True,
            "requested_disease_name": str(disease_name),
            "canonical_disease_name": canonical_disease_name,
            "relation_count": int(len(reference)),
            "matched_relation_count": int(matched_relation_count),
            "target_microbe_count": int(len(normalized_profile)),
        }

    def _profile_from_frame(self, frame: pd.DataFrame) -> dict[str, object]:
        work = self._annotate_step2_mechanism(self._annotate_amr(frame))
        candidate_diseases = self._candidate_diseases_from_frame(work)
        row = work.iloc[0]
        similar_library_drugs = self._top_similar_library_drugs(row)
        similarity_disease_context = self._similarity_disease_context(row)
        biotransform_sidecar = {
            "enabled": _safe_bool(row.get("biotransform_sidecar_enabled")),
            "mode": row.get("biotransform_sidecar_mode"),
            "reason": row.get("biotransform_sidecar_reason"),
            "support_score": _safe_float(row.get("biotransform_sidecar_support_score")),
            "similarity": _safe_float(row.get("biotransform_sidecar_similarity")),
            "confidence": row.get("biotransform_sidecar_confidence"),
            "reference_count": _safe_float(row.get("biotransform_sidecar_reference_count")),
            "weighted_product_count": _safe_float(row.get("biotransform_sidecar_weighted_product_count")),
            "weighted_fraction_in_gut": _safe_float(row.get("biotransform_sidecar_weighted_fraction_in_gut")),
            "reference_drugs": row.get("biotransform_sidecar_reference_drugs"),
            "reference_product_ids": row.get("biotransform_sidecar_reference_product_ids"),
            "reference_ec_numbers": row.get("biotransform_sidecar_reference_ec_numbers"),
            "reference_reaction_centers": row.get("biotransform_sidecar_reference_reaction_centers"),
            "primary_products": row.get("biotransform_sidecar_primary_products"),
        }

        label_column = "display_step1_predicted_effect_label"
        probability_column = "display_step1_predicted_inhibit_probability"
        score_column = "display_step1_predicted_effect_score"

        if score_column in work.columns:
            work["_step1_abs_effect"] = pd.to_numeric(work[score_column], errors="coerce").abs()
        else:
            work["_step1_abs_effect"] = np.nan

        effect_sorted = work.sort_values(
            ["_step1_abs_effect", probability_column if probability_column in work.columns else "nt_code"],
            ascending=[False, False],
        )
        top_effect_microbes = effect_sorted.head(12)
        metabolism_sorted = work.sort_values(
            ["predicted_metabolized_probability", "predicted_parent_depletion_magnitude"],
            ascending=[False, False],
        )
        top_metabolism_microbes = metabolism_sorted.head(12)

        step1_counts = {}
        if label_column in work.columns:
            step1_counts = {
                str(key): int(value)
                for key, value in work[label_column].fillna("missing").value_counts().to_dict().items()
            }
        step2_counts = {
            str(key): int(value)
            for key, value in work["predicted_metabolism_label"].fillna("missing").value_counts().to_dict().items()
        }

        top_effect_records = []
        for _, effect_row in top_effect_microbes.iterrows():
            top_effect_records.append(
                {
                    "nt_code": effect_row.get("nt_code"),
                    "microbe_label": effect_row.get("microbe_label"),
                    "genus": effect_row.get("genus"),
                    "phylum": effect_row.get("phylum"),
                    "predicted_effect_label": effect_row.get(label_column),
                    "predicted_inhibit_probability": _safe_float(effect_row.get(probability_column)),
                    "predicted_effect_score": _safe_float(effect_row.get(score_column)),
                    "predicted_promote_probability_base": _safe_float(effect_row.get("predicted_promote_probability_base")),
                    "predicted_promote_probability_refined": _safe_float(
                        effect_row.get("predicted_promote_probability_refined")
                    ),
                    "predicted_promote_support_score": _safe_float(effect_row.get("predicted_promote_support_score")),
                    "predicted_promote_support_type": effect_row.get("predicted_promote_support_type"),
                    "predicted_promote_evidence_type": effect_row.get("predicted_promote_evidence_type"),
                    "predicted_cross_feeding_reference_flag": _safe_bool(
                        effect_row.get("predicted_cross_feeding_reference_flag")
                    ),
                    "predicted_cross_feeding_support_microbe": effect_row.get("predicted_cross_feeding_support_microbe"),
                    "predicted_cross_feeding_match_mode": effect_row.get("predicted_cross_feeding_match_mode"),
                    "predicted_cross_feeding_matched_term": effect_row.get("predicted_cross_feeding_matched_term"),
                    "raw_predicted_effect_label": effect_row.get("raw_step1_predicted_effect_label"),
                    "raw_predicted_inhibit_probability": _safe_float(effect_row.get("raw_step1_predicted_inhibit_probability")),
                    "raw_predicted_effect_score": _safe_float(effect_row.get("raw_step1_predicted_effect_score")),
                    "amr_conflict_flag": bool(effect_row.get("amr_conflict_flag", False)),
                    "amr_correction_applied": bool(effect_row.get("amr_correction_applied", False)),
                    "amr_expected_phenotype": effect_row.get("amr_expected_phenotype"),
                    "amr_rule_strength": effect_row.get("amr_rule_strength"),
                    "amr_primary_drug_class": effect_row.get("amr_primary_drug_class"),
                }
            )

        panel_effect_records = []
        for _, effect_row in effect_sorted.iterrows():
            panel_effect_records.append(
                {
                    "nt_code": effect_row.get("nt_code"),
                    "microbe_label": effect_row.get("microbe_label"),
                    "genus": effect_row.get("genus"),
                    "phylum": effect_row.get("phylum"),
                    "predicted_effect_label": effect_row.get(label_column),
                    "predicted_inhibit_probability": _safe_float(effect_row.get(probability_column)),
                    "predicted_effect_score": _safe_float(effect_row.get(score_column)),
                    "predicted_promote_probability_base": _safe_float(effect_row.get("predicted_promote_probability_base")),
                    "predicted_promote_probability_refined": _safe_float(
                        effect_row.get("predicted_promote_probability_refined")
                    ),
                    "predicted_promote_support_score": _safe_float(effect_row.get("predicted_promote_support_score")),
                    "predicted_promote_support_type": effect_row.get("predicted_promote_support_type"),
                    "predicted_promote_evidence_type": effect_row.get("predicted_promote_evidence_type"),
                    "predicted_cross_feeding_reference_flag": _safe_bool(
                        effect_row.get("predicted_cross_feeding_reference_flag")
                    ),
                    "predicted_cross_feeding_support_microbe": effect_row.get("predicted_cross_feeding_support_microbe"),
                    "predicted_cross_feeding_match_mode": effect_row.get("predicted_cross_feeding_match_mode"),
                    "predicted_cross_feeding_matched_term": effect_row.get("predicted_cross_feeding_matched_term"),
                    "raw_predicted_effect_label": effect_row.get("raw_step1_predicted_effect_label"),
                    "raw_predicted_inhibit_probability": _safe_float(effect_row.get("raw_step1_predicted_inhibit_probability")),
                    "raw_predicted_effect_score": _safe_float(effect_row.get("raw_step1_predicted_effect_score")),
                    "amr_conflict_flag": bool(effect_row.get("amr_conflict_flag", False)),
                    "amr_correction_applied": bool(effect_row.get("amr_correction_applied", False)),
                    "amr_expected_phenotype": effect_row.get("amr_expected_phenotype"),
                    "amr_rule_strength": effect_row.get("amr_rule_strength"),
                    "amr_primary_drug_class": effect_row.get("amr_primary_drug_class"),
                }
            )

        def _reaction_class_display(row: pd.Series) -> object:
            primary = row.get("predicted_reaction_class")
            if primary is not None and not pd.isna(primary) and str(primary).strip():
                return primary
            fallback = row.get("predicted_enzyme_reaction_classes")
            if fallback is not None and not pd.isna(fallback) and str(fallback).strip():
                return fallback
            return primary

        top_metabolism_records = []
        for _, metabolism_row in top_metabolism_microbes.iterrows():
            top_metabolism_records.append(
                {
                    "nt_code": metabolism_row.get("nt_code"),
                    "microbe_label": metabolism_row.get("microbe_label"),
                    "genus": metabolism_row.get("genus"),
                    "phylum": metabolism_row.get("phylum"),
                    "predicted_metabolism_label": metabolism_row.get("predicted_metabolism_label"),
                    "predicted_metabolized_probability": _safe_float(metabolism_row.get("predicted_metabolized_probability")),
                    "predicted_parent_depletion_fraction": _safe_float(
                        metabolism_row.get("predicted_parent_depletion_fraction")
                    ),
                    "applicability_flag": _safe_bool(metabolism_row.get("applicability_flag")),
                    "predicted_mechanism_projection_flag": _safe_bool(
                        metabolism_row.get("predicted_mechanism_projection_flag")
                    ),
                    "predicted_reaction_class": _reaction_class_display(metabolism_row),
                    "predicted_reaction_confidence": _safe_float(metabolism_row.get("predicted_reaction_confidence")),
                    "predicted_candidate_product_ids": metabolism_row.get("predicted_candidate_product_ids"),
                    "predicted_candidate_product_count": _safe_float(metabolism_row.get("predicted_candidate_product_count")),
                    "predicted_evidence_gene_ids": metabolism_row.get("predicted_evidence_gene_ids"),
                    "predicted_evidence_gene_count": _safe_float(metabolism_row.get("predicted_evidence_gene_count")),
                    "predicted_enzyme_prior_flag": _safe_bool(metabolism_row.get("predicted_enzyme_prior_flag")),
                    "predicted_enzyme_match_count": _safe_float(metabolism_row.get("predicted_enzyme_match_count")),
                    "predicted_enzyme_ids": metabolism_row.get("predicted_enzyme_ids"),
                    "predicted_enzyme_names": metabolism_row.get("predicted_enzyme_names"),
                    "predicted_enzyme_reaction_classes": metabolism_row.get("predicted_enzyme_reaction_classes"),
                    "predicted_enzyme_bond_targets": metabolism_row.get("predicted_enzyme_bond_targets"),
                    "predicted_enzyme_presence_score": _safe_float(metabolism_row.get("predicted_enzyme_presence_score")),
                    "predicted_enzyme_support_score": _safe_float(metabolism_row.get("predicted_enzyme_support_score")),
                    "predicted_enzyme_step1_promote_support_score": _safe_float(
                        metabolism_row.get("predicted_enzyme_step1_promote_support_score")
                    ),
                    "predicted_enzyme_step1_inhibit_risk_score": _safe_float(
                        metabolism_row.get("predicted_enzyme_step1_inhibit_risk_score")
                    ),
                    "biotransform_sidecar_enabled": _safe_bool(metabolism_row.get("biotransform_sidecar_enabled")),
                    "biotransform_sidecar_mode": metabolism_row.get("biotransform_sidecar_mode"),
                    "biotransform_sidecar_support_score": _safe_float(
                        metabolism_row.get("biotransform_sidecar_support_score")
                    ),
                    "biotransform_sidecar_similarity": _safe_float(metabolism_row.get("biotransform_sidecar_similarity")),
                    "biotransform_sidecar_confidence": metabolism_row.get("biotransform_sidecar_confidence"),
                    "biotransform_sidecar_weighted_product_count": _safe_float(
                        metabolism_row.get("biotransform_sidecar_weighted_product_count")
                    ),
                    "biotransform_sidecar_reference_drugs": metabolism_row.get("biotransform_sidecar_reference_drugs"),
                }
            )

        panel_metabolism_records = []
        for _, metabolism_row in metabolism_sorted.iterrows():
            panel_metabolism_records.append(
                {
                    "nt_code": metabolism_row.get("nt_code"),
                    "microbe_label": metabolism_row.get("microbe_label"),
                    "genus": metabolism_row.get("genus"),
                    "phylum": metabolism_row.get("phylum"),
                    "predicted_metabolism_label": metabolism_row.get("predicted_metabolism_label"),
                    "predicted_metabolized_probability": _safe_float(metabolism_row.get("predicted_metabolized_probability")),
                    "predicted_parent_depletion_fraction": _safe_float(
                        metabolism_row.get("predicted_parent_depletion_fraction")
                    ),
                    "applicability_flag": _safe_bool(metabolism_row.get("applicability_flag")),
                    "predicted_mechanism_projection_flag": _safe_bool(
                        metabolism_row.get("predicted_mechanism_projection_flag")
                    ),
                    "predicted_reaction_class": _reaction_class_display(metabolism_row),
                    "predicted_reaction_confidence": _safe_float(metabolism_row.get("predicted_reaction_confidence")),
                    "predicted_candidate_product_ids": metabolism_row.get("predicted_candidate_product_ids"),
                    "predicted_candidate_product_count": _safe_float(metabolism_row.get("predicted_candidate_product_count")),
                    "predicted_evidence_gene_ids": metabolism_row.get("predicted_evidence_gene_ids"),
                    "predicted_evidence_gene_count": _safe_float(metabolism_row.get("predicted_evidence_gene_count")),
                    "predicted_enzyme_prior_flag": _safe_bool(metabolism_row.get("predicted_enzyme_prior_flag")),
                    "predicted_enzyme_match_count": _safe_float(metabolism_row.get("predicted_enzyme_match_count")),
                    "predicted_enzyme_ids": metabolism_row.get("predicted_enzyme_ids"),
                    "predicted_enzyme_names": metabolism_row.get("predicted_enzyme_names"),
                    "predicted_enzyme_reaction_classes": metabolism_row.get("predicted_enzyme_reaction_classes"),
                    "predicted_enzyme_bond_targets": metabolism_row.get("predicted_enzyme_bond_targets"),
                    "predicted_enzyme_presence_score": _safe_float(metabolism_row.get("predicted_enzyme_presence_score")),
                    "predicted_enzyme_support_score": _safe_float(metabolism_row.get("predicted_enzyme_support_score")),
                    "predicted_enzyme_step1_promote_support_score": _safe_float(
                        metabolism_row.get("predicted_enzyme_step1_promote_support_score")
                    ),
                    "predicted_enzyme_step1_inhibit_risk_score": _safe_float(
                        metabolism_row.get("predicted_enzyme_step1_inhibit_risk_score")
                    ),
                    "biotransform_sidecar_enabled": _safe_bool(metabolism_row.get("biotransform_sidecar_enabled")),
                    "biotransform_sidecar_mode": metabolism_row.get("biotransform_sidecar_mode"),
                    "biotransform_sidecar_support_score": _safe_float(
                        metabolism_row.get("biotransform_sidecar_support_score")
                    ),
                    "biotransform_sidecar_similarity": _safe_float(metabolism_row.get("biotransform_sidecar_similarity")),
                    "biotransform_sidecar_confidence": metabolism_row.get("biotransform_sidecar_confidence"),
                    "biotransform_sidecar_weighted_product_count": _safe_float(
                        metabolism_row.get("biotransform_sidecar_weighted_product_count")
                    ),
                    "biotransform_sidecar_reference_drugs": metabolism_row.get("biotransform_sidecar_reference_drugs"),
                }
            )

        aggregated = {
            "step1_counts": step1_counts,
            "step2_counts": step2_counts,
            "mean_predicted_effect_score": _safe_float(
                work[score_column].mean() if score_column in work.columns else np.nan
            ),
            "mean_predicted_inhibit_probability": _safe_float(
                work[probability_column].mean() if probability_column in work.columns else np.nan
            ),
            "mean_predicted_promote_probability_refined": _safe_float(
                work["predicted_promote_probability_refined"].mean()
                if "predicted_promote_probability_refined" in work.columns
                else np.nan
            ),
            "mean_predicted_metabolized_probability": _safe_float(work["predicted_metabolized_probability"].mean()),
            "applicability_rate": _safe_float(work["applicability_flag"].fillna(False).mean()),
            "mechanism_projection_rate": _safe_float(work["predicted_mechanism_projection_flag"].fillna(False).mean()),
            "enzyme_prior_support_rate": _safe_float(
                work.get("predicted_enzyme_prior_flag", pd.Series(False, index=work.index)).fillna(False).mean()
            ),
            "mean_enzyme_support_score": _safe_float(
                work.get("predicted_enzyme_support_score", pd.Series(np.nan, index=work.index)).mean()
            ),
            "reaction_projection_pairs": int(work["predicted_reaction_class"].fillna("").astype(str).str.strip().ne("").sum()),
            "gene_projection_pairs": int(
                work["predicted_evidence_gene_ids"].fillna("").astype(str).str.strip().ne("").sum()
            ),
            "enzyme_prior_supported_pairs": int(
                work.get("predicted_enzyme_prior_flag", pd.Series(False, index=work.index)).fillna(False).sum()
            ),
            "amr_conflict_pairs": int(work["amr_conflict_flag"].fillna(False).sum()),
            "amr_corrected_pairs": int(work["amr_correction_applied"].fillna(False).sum()),
            "metabolism_supported_promote_pairs": int(
                (
                    work.get("predicted_promote_support_type", pd.Series(np.nan, index=work.index)).eq("self_metabolism_supported")
                    & work.get(label_column, pd.Series(np.nan, index=work.index)).eq("promote")
                ).sum()
            ),
            "cross_feeding_supported_promote_pairs": int(
                (
                    work.get("predicted_cross_feeding_reference_flag", pd.Series(False, index=work.index)).fillna(False).astype(bool)
                    & work.get(label_column, pd.Series(np.nan, index=work.index)).eq("promote")
                ).sum()
            ),
            "biotransform_sidecar_enabled": _safe_bool(row.get("biotransform_sidecar_enabled")),
            "biotransform_sidecar_support_score": _safe_float(row.get("biotransform_sidecar_support_score")),
            "biotransform_sidecar_similarity": _safe_float(row.get("biotransform_sidecar_similarity")),
            "biotransform_sidecar_confidence": row.get("biotransform_sidecar_confidence"),
            "biotransform_sidecar_weighted_product_count": _safe_float(
                row.get("biotransform_sidecar_weighted_product_count")
            ),
            "biotransform_sidecar_weighted_fraction_in_gut": _safe_float(
                row.get("biotransform_sidecar_weighted_fraction_in_gut")
            ),
        }
        confidence_payload = evaluate_prediction_confidence(
            effect_frame=work,
            step1_label_column=label_column,
            step1_probability_column=probability_column,
            step1_score_column=score_column,
            drug_profile=self._infer_drug_profile_from_frame(work, row),
            molecular_weight=_first_numeric_value(row, ["molecular_weight", "rdkit_exact_mol_wt"]),
            xlogp=_first_numeric_value(row, ["xlogp", "rdkit_logp"]),
            mw_bounds=self.mw_ood_bounds,
            xlogp_bounds=self.xlogp_ood_bounds,
        )
        aggregated.update(confidence_payload)

        return _clean_json(
            {
                "drug": self._drug_metadata(row),
                "aggregated": aggregated,
                "confidence_score": confidence_payload["confidence_score"],
                "confidence_tier": confidence_payload["confidence_tier"],
                "warning_flags": confidence_payload["warning_flags"],
                "confidence_breakdown": confidence_payload["confidence_breakdown"],
                "confidence_explanation": confidence_payload["confidence_explanation"],
                "candidate_diseases": candidate_diseases,
                "marketed_disease_context": self._marketed_disease_context(row.get("chemical_name")),
                "similarity_disease_context": similarity_disease_context,
                "similar_library_drugs": similar_library_drugs,
                "biotransform_sidecar": biotransform_sidecar,
                "top_effect_microbes": top_effect_records,
                "panel_effect_microbes": panel_effect_records,
                "top_metabolism_microbes": top_metabolism_records,
                "panel_metabolism_microbes": panel_metabolism_records,
            }
        )

    def _pair_payload_from_row(self, row: pd.Series) -> dict[str, object]:
        step1 = {
            "predicted_effect_label": row.get("display_step1_predicted_effect_label")
            if "display_step1_predicted_effect_label" in row.index
            else (row.get(self.step1_label_column) if self.step1_label_column else None),
            "predicted_binary_effect_label": row.get(self.step1_binary_column) if self.step1_binary_column else None,
            "predicted_inhibit_probability": _safe_float(
                row.get("display_step1_predicted_inhibit_probability")
                if "display_step1_predicted_inhibit_probability" in row.index
                else (row.get(self.step1_probability_column) if self.step1_probability_column else None)
            ),
            "predicted_effect_score": _safe_float(
                row.get("display_step1_predicted_effect_score")
                if "display_step1_predicted_effect_score" in row.index
                else (row.get(self.step1_score_column) if self.step1_score_column else None)
            ),
            "predicted_effect_magnitude": _safe_float(
                row.get(self.step1_magnitude_column) if self.step1_magnitude_column else None
            ),
            "predicted_promote_probability_base": _safe_float(row.get("predicted_promote_probability_base")),
            "predicted_promote_probability_refined": _safe_float(row.get("predicted_promote_probability_refined")),
            "predicted_promote_support_score": _safe_float(row.get("predicted_promote_support_score")),
            "predicted_promote_support_type": row.get("predicted_promote_support_type"),
            "predicted_promote_evidence_type": row.get("predicted_promote_evidence_type"),
            "predicted_cross_feeding_reference_flag": _safe_bool(row.get("predicted_cross_feeding_reference_flag")),
            "predicted_cross_feeding_support_microbe": row.get("predicted_cross_feeding_support_microbe"),
            "predicted_cross_feeding_match_mode": row.get("predicted_cross_feeding_match_mode"),
            "predicted_cross_feeding_matched_term": row.get("predicted_cross_feeding_matched_term"),
            "raw_predicted_effect_label": row.get("raw_step1_predicted_effect_label"),
            "raw_predicted_inhibit_probability": _safe_float(row.get("raw_step1_predicted_inhibit_probability")),
            "raw_predicted_effect_score": _safe_float(row.get("raw_step1_predicted_effect_score")),
            "observed_effect_label": row.get(self.step1_observed_label_column) if self.step1_observed_label_column else None,
            "observed_binary_effect_label": row.get(self.step1_observed_binary_column)
            if self.step1_observed_binary_column
            else None,
            "observed_effect_score": _safe_float(
                row.get(self.step1_observed_score_column) if self.step1_observed_score_column else None
            ),
            "amr_conflict_flag": bool(row.get("amr_conflict_flag", False)),
            "amr_correction_applied": bool(row.get("amr_correction_applied", False)),
            "amr_expected_phenotype": row.get("amr_expected_phenotype"),
            "amr_rule_id": row.get("amr_rule_id"),
            "amr_rule_strength": row.get("amr_rule_strength"),
            "amr_rule_level": row.get("amr_rule_level"),
            "amr_primary_drug_class": row.get("amr_primary_drug_class"),
            "amr_mechanism_hint": row.get("amr_mechanism_hint"),
            "amr_source_name": row.get("amr_source_name"),
            "amr_source_url": row.get("amr_source_url"),
        }
        step2 = {
            "predicted_metabolism_label": row.get("predicted_metabolism_label"),
            "predicted_metabolized_probability": _safe_float(row.get("predicted_metabolized_probability")),
            "predicted_parent_depletion_fraction": _safe_float(row.get("predicted_parent_depletion_fraction")),
            "predicted_parent_depletion_magnitude": _safe_float(row.get("predicted_parent_depletion_magnitude")),
            "drug_max_fingerprint_jaccard": _safe_float(row.get("drug_max_fingerprint_jaccard")),
            "scaffold_seen_in_training": _safe_bool(row.get("scaffold_seen_in_training")),
            "microbe_genus_seen_in_training": _safe_bool(row.get("microbe_genus_seen_in_training")),
            "microbe_phylum_seen_in_training": _safe_bool(row.get("microbe_phylum_seen_in_training")),
            "applicability_flag": _safe_bool(row.get("applicability_flag")),
            "predicted_mechanism_projection_flag": _safe_bool(row.get("predicted_mechanism_projection_flag")),
            "predicted_reaction_class": row.get("predicted_reaction_class"),
            "predicted_reaction_confidence": _safe_float(row.get("predicted_reaction_confidence")),
            "predicted_reaction_support_pairs": _safe_float(row.get("predicted_reaction_support_pairs")),
            "predicted_mechanism_support_score": _safe_float(row.get("predicted_mechanism_support_score")),
            "predicted_mechanism_support_scopes": row.get("predicted_mechanism_support_scopes"),
            "predicted_candidate_product_ids": row.get("predicted_candidate_product_ids"),
            "predicted_candidate_product_count": _safe_float(row.get("predicted_candidate_product_count")),
            "predicted_evidence_gene_ids": row.get("predicted_evidence_gene_ids"),
            "predicted_evidence_gene_count": _safe_float(row.get("predicted_evidence_gene_count")),
            "predicted_enzyme_prior_flag": _safe_bool(row.get("predicted_enzyme_prior_flag")),
            "predicted_enzyme_match_count": _safe_float(row.get("predicted_enzyme_match_count")),
            "predicted_enzyme_ids": row.get("predicted_enzyme_ids"),
            "predicted_enzyme_names": row.get("predicted_enzyme_names"),
            "predicted_enzyme_reaction_classes": row.get("predicted_enzyme_reaction_classes"),
            "predicted_enzyme_bond_targets": row.get("predicted_enzyme_bond_targets"),
            "predicted_enzyme_presence_score": _safe_float(row.get("predicted_enzyme_presence_score")),
            "predicted_enzyme_support_score": _safe_float(row.get("predicted_enzyme_support_score")),
            "predicted_enzyme_step1_promote_support_score": _safe_float(
                row.get("predicted_enzyme_step1_promote_support_score")
            ),
            "predicted_enzyme_step1_inhibit_risk_score": _safe_float(
                row.get("predicted_enzyme_step1_inhibit_risk_score")
            ),
            "biotransform_sidecar_enabled": _safe_bool(row.get("biotransform_sidecar_enabled")),
            "biotransform_sidecar_mode": row.get("biotransform_sidecar_mode"),
            "biotransform_sidecar_support_score": _safe_float(row.get("biotransform_sidecar_support_score")),
            "biotransform_sidecar_similarity": _safe_float(row.get("biotransform_sidecar_similarity")),
            "biotransform_sidecar_confidence": row.get("biotransform_sidecar_confidence"),
            "biotransform_sidecar_weighted_product_count": _safe_float(
                row.get("biotransform_sidecar_weighted_product_count")
            ),
            "biotransform_sidecar_weighted_fraction_in_gut": _safe_float(
                row.get("biotransform_sidecar_weighted_fraction_in_gut")
            ),
            "biotransform_sidecar_reference_drugs": row.get("biotransform_sidecar_reference_drugs"),
            "biotransform_sidecar_reference_product_ids": row.get("biotransform_sidecar_reference_product_ids"),
            "biotransform_sidecar_reference_ec_numbers": row.get("biotransform_sidecar_reference_ec_numbers"),
        }
        return _clean_json(
            {
                "drug": self._drug_metadata(row),
                "microbe": self._microbe_metadata(row),
                "step1": step1,
                "step2": step2,
            }
        )

    def _build_custom_drug_table(
        self,
        drug_name: str,
        smiles: str,
        drug_id: str,
        therapeutic_class: str | None = None,
        therapeutic_effect: str | None = None,
        target_species: str = "human",
        human_use: bool = True,
        veterinary: bool = False,
    ) -> pd.DataFrame:
        drug = pd.DataFrame(
            [
                {
                    "prestwick_id": drug_id,
                    "chemical_name": drug_name,
                    "cid_flat": np.nan,
                    "cid_active": np.nan,
                    "cid_main": np.nan,
                    "main_component_smiles": smiles,
                    "smiles": smiles,
                    "molecular_formula": np.nan,
                    "molecular_weight": np.nan,
                    "xlogp": np.nan,
                    "tpsa": np.nan,
                    "complexity": np.nan,
                    "volume3d": np.nan,
                    "therapeutic_class": therapeutic_class,
                    "therapeutic_effect": therapeutic_effect,
                    "atc_codes": np.nan,
                    "atc_primary_l1": np.nan,
                    "atc_primary_l3": np.nan,
                    "atc_primary_l4": np.nan,
                    "target_species": target_species,
                    "human_use": human_use,
                    "veterinary": veterinary,
                    "dose_umol": np.nan,
                    "estimated_intestine_concentration_um": np.nan,
                    "plasma_concentration_um": np.nan,
                    "fraction_excreted_in_feces": np.nan,
                    "fraction_excreted_in_urine": np.nan,
                    "estimated_colon_concentration_um": np.nan,
                    "screen_conc_20_um_as_ug_ml": np.nan,
                }
            ]
        )
        drug = _compute_smiles_descriptors(drug)
        drug = enrich_drug_table_with_rdkit(drug, smiles_columns=["main_component_smiles", "smiles"])
        drug = drug.loc[:, ~drug.columns.duplicated()].copy()
        rdkit_valid = pd.to_numeric(drug.get("rdkit_valid_smiles"), errors="coerce").fillna(0.0)
        if rdkit_valid.empty or float(rdkit_valid.iloc[0]) <= 0:
            raise ValueError(f"SMILES 格式无效，无法解析: {smiles}")
        drug = annotate_compound_semantics(drug)
        if "molecular_formula" in drug.columns:
            drug["molecular_formula"] = drug["molecular_formula"].fillna(drug.get("rdkit_formula"))
        if "molecular_weight" in drug.columns:
            drug["molecular_weight"] = drug["molecular_weight"].fillna(drug.get("rdkit_exact_mol_wt"))
        if "xlogp" in drug.columns:
            drug["xlogp"] = drug["xlogp"].fillna(drug.get("rdkit_logp"))
        if "tpsa" in drug.columns:
            drug["tpsa"] = drug["tpsa"].fillna(drug.get("rdkit_tpsa"))
        if "complexity" not in drug.columns:
            drug["complexity"] = np.nan
        if "volume3d" not in drug.columns:
            drug["volume3d"] = np.nan
        return drug

    def _build_custom_pair_table(
        self,
        drug_name: str,
        smiles: str,
        drug_id: str,
        therapeutic_class: str | None = None,
        therapeutic_effect: str | None = None,
        target_species: str = "human",
        human_use: bool = True,
        veterinary: bool = False,
    ) -> pd.DataFrame:
        drug_table = self._build_custom_drug_table(
            drug_name=drug_name,
            smiles=smiles,
            drug_id=drug_id,
            therapeutic_class=therapeutic_class,
            therapeutic_effect=therapeutic_effect,
            target_species=target_species,
            human_use=human_use,
            veterinary=veterinary,
        )
        microbes = self.microbe_feature_table.copy()
        microbes["_merge_key"] = 1
        drug_table["_merge_key"] = 1
        pairs = microbes.merge(drug_table, on="_merge_key", how="inner").drop(columns=["_merge_key"])
        pairs["pair_id"] = pairs["prestwick_id"].astype(str) + "::" + pairs["nt_code"].astype(str)
        pairs["effect_label"] = np.nan
        pairs["binary_effect_label"] = np.nan
        pairs["effect_score"] = np.nan
        pairs["source_dataset"] = "custom_input"
        pairs["label_tier"] = "inference"
        return pairs

    def _create_custom_prediction_session(
        self,
        drug_name: str,
        smiles: str,
        therapeutic_class: str | None = None,
        therapeutic_effect: str | None = None,
        target_species: str = "human",
        human_use: bool = True,
        veterinary: bool = False,
    ) -> tuple[str, pd.DataFrame]:
        session_id = uuid.uuid4().hex[:12]
        drug_id = f"Custom-{session_id}"
        session_dir = self.temp_root / f"custom_session_{session_id}"
        session_dir.mkdir(parents=True, exist_ok=True)

        pair_table = self._build_custom_pair_table(
            drug_name=drug_name,
            smiles=smiles,
            drug_id=drug_id,
            therapeutic_class=therapeutic_class,
            therapeutic_effect=therapeutic_effect,
            target_species=target_species,
            human_use=human_use,
            veterinary=veterinary,
        )
        custom_step1_input_path = session_dir / "custom_step1_pairs.csv"
        pair_table.to_csv(custom_step1_input_path, index=False)

        step1_output_dir = session_dir / "step1"
        predict_step1_hybrid(
            input_table_path=custom_step1_input_path,
            output_dir=step1_output_dir,
            classification_prepare_dir=self.step1_chemprop_prepare_dir,
            chemprop_model_path=self.step1_chemprop_model_path,
            regressor_path=self.step1_regressor_path,
            regressor_metrics_path=self.step1_regressor_metrics_path,
        )

        step2_input_dir = session_dir / "step2_inputs"
        build_step2_input_tables(
            step1_predictions_path=step1_output_dir / "predictions.csv",
            output_dir=step2_input_dir,
            step2_label_table_paths=None,
        )

        step2_output_dir = session_dir / "step2"
        predict_step2_baseline(
            input_table_path=step2_input_dir / "step2_candidate_pairs_full.csv",
            output_dir=step2_output_dir,
            classifier_path=self.step2_classifier_path,
            regressor_path=self.step2_regressor_path,
            metrics_path=self.step2_metrics_path,
            applicability_reference_path=self.step2_applicability_reference_path,
            mechanism_reference_path=self.step2_mechanism_reference_path,
            enzyme_microbe_panel_path=self.step2_enzyme_microbe_panel_path,
            enzyme_function_catalog_path=self.step2_enzyme_function_catalog_path,
        )

        integrated_frame = pd.read_csv(step2_output_dir / "predictions.csv", low_memory=False)
        # Step 3 custom simulation resolves one drug subset by prestwick_id / chemical_name.
        # Some custom Step 2 outputs may drop the original name columns, so restore them here.
        if "prestwick_id" not in integrated_frame.columns:
            integrated_frame["prestwick_id"] = drug_id
        else:
            integrated_frame["prestwick_id"] = integrated_frame["prestwick_id"].fillna(drug_id).replace("", drug_id)
        if "chemical_name" not in integrated_frame.columns:
            integrated_frame["chemical_name"] = drug_name
        else:
            integrated_frame["chemical_name"] = integrated_frame["chemical_name"].fillna(drug_name).replace("", drug_name)
        integrated_frame = refine_step1_promote_with_step2(
            integrated_frame,
            promote_classifier_path=self.step1_promote_classifier_path,
            promote_metrics_path=self.step1_promote_metrics_path,
            cross_feeding_reference_path=self.cross_feeding_reference_path,
        )
        integrated_frame = self._annotate_biotransform_sidecar(integrated_frame)
        integrated_frame.to_csv(step2_output_dir / "predictions.csv", index=False)
        self.custom_sessions[session_id] = {
            "session_id": session_id,
            "session_dir": session_dir,
            "drug_id": drug_id,
            "drug_name": drug_name,
            "smiles": smiles,
            "frame": integrated_frame,
            "integrated_predictions_path": step2_output_dir / "predictions.csv",
        }
        return session_id, integrated_frame

    def bootstrap(self) -> dict[str, object]:
        step1_counts = {}
        if self.step1_label_column is not None:
            step1_counts = {
                str(key): int(value)
                for key, value in self.frame[self.step1_label_column].fillna("missing").value_counts().to_dict().items()
            }
        step2_counts = {
            str(key): int(value)
            for key, value in self.frame["predicted_metabolism_label"].fillna("missing").value_counts().to_dict().items()
        }

        drugs = [
            _clean_json(
                {
                    "prestwick_id": row["prestwick_id"],
                    "chemical_name": row["chemical_name"],
                    "label": f"{row['chemical_name']} ({row['prestwick_id']})",
                    "therapeutic_class": row.get("therapeutic_class"),
                    "therapeutic_effect": row.get("therapeutic_effect"),
                    "atc_primary_l1": row.get("atc_primary_l1"),
                    "murcko_scaffold": row.get("murcko_scaffold"),
                }
            )
            for _, row in self.drug_table.iterrows()
        ]
        microbes = [
            _clean_json(
                {
                    "nt_code": row["nt_code"],
                    "microbe_label": row["microbe_label"],
                    "label": f"{row['microbe_label']} ({row['nt_code']})",
                    "species_label": row.get("species_label"),
                    "genus": row.get("genus"),
                    "phylum": row.get("phylum"),
                }
            )
            for _, row in self.library_microbe_table.iterrows()
        ]
        custom_microbes = [
            _clean_json(
                {
                    "nt_code": row["nt_code"],
                    "microbe_label": row["microbe_label"],
                    "label": f"{row['microbe_label']} ({row['nt_code']})",
                    "species_label": row.get("species_label"),
                    "genus": row.get("genus"),
                    "phylum": row.get("phylum"),
                }
            )
            for _, row in self.microbe_table.iterrows()
        ]
        scenarios = [
            {
                "scenario_name": name,
                "description": str(payload.get("description", "")),
            }
            for name, payload in sorted(BUILTIN_SCENARIOS.items())
        ]
        summary = {
            "n_pairs": int(len(self.frame)),
            "n_drugs": int(self.drug_table["prestwick_id"].nunique()),
            "n_microbes": int(self.library_microbe_table["nt_code"].nunique()),
            "n_custom_microbes": int(self.microbe_table["nt_code"].nunique()),
            "step1_counts": step1_counts,
            "step2_counts": step2_counts,
            "n_applicable_pairs": int(self.frame.get("applicability_flag", pd.Series(dtype=bool)).fillna(False).sum()),
        }
        return _clean_json(
            {
                "summary": summary,
                "drugs": drugs,
                "microbes": microbes,
                "custom_microbes": custom_microbes,
                "diseases": self.disease_catalog,
                "scenarios": scenarios,
                "cohort_communities": self.cohort_communities,
                "demo_candidates": self.demo_ranking,
            }
        )

    def get_drug_profile(self, drug_query: str) -> dict[str, object]:
        frame = self._drug_frame(drug_query)
        return self._profile_from_frame(frame)

    def get_pair_prediction(self, drug_query: str, microbe_query: str) -> dict[str, object]:
        row = self._pair_row(drug_query, microbe_query)
        return self._pair_payload_from_row(row)

    def predict_custom_drug(
        self,
        drug_name: str | None,
        smiles: str,
        microbe_query: str | None = None,
        therapeutic_class: str | None = None,
        therapeutic_effect: str | None = None,
        target_species: str = "human",
        human_use: bool = True,
        veterinary: bool = False,
    ) -> dict[str, object]:
        if not str(smiles).strip():
            raise ValueError("smiles 不能为空。")
        normalized_smiles = str(smiles).strip()
        normalized_name = str(drug_name or "").strip()
        if not normalized_name:
            normalized_name = f"Custom drug [{normalized_smiles[:18]}]"

        session_id, frame = self._create_custom_prediction_session(
            drug_name=normalized_name,
            smiles=normalized_smiles,
            therapeutic_class=therapeutic_class,
            therapeutic_effect=therapeutic_effect,
            target_species=target_species,
            human_use=human_use,
            veterinary=veterinary,
        )
        microbe_query = microbe_query or str(self.microbe_table.iloc[0]["nt_code"])
        pair_payload = self._pair_payload_from_row(self._pair_row_from_frame(frame, microbe_query))
        profile = self._profile_from_frame(frame)
        return _clean_json(
            {
                "session_id": session_id,
                "confidence_score": profile.get("confidence_score"),
                "confidence_tier": profile.get("confidence_tier"),
                "warning_flags": profile.get("warning_flags", []),
                "confidence_explanation": profile.get("confidence_explanation"),
                "profile": profile,
                "selected_pair": pair_payload,
            }
        )

    def get_custom_pair_prediction(self, session_id: str, microbe_query: str) -> dict[str, object]:
        frame = self._custom_frame(session_id)
        row = self._pair_row_from_frame(frame, microbe_query)
        return self._pair_payload_from_row(row)

    def simulate_step3(
        self,
        drug_query: str,
        scenario_name: str = "healthy_reference",
        community_table_path: str | None = None,
        disease_name: str | None = None,
        n_steps: int = 14,
        initial_dose: float = 1.0,
        repeat_dose: float = 1.0,
        dosing_interval: int = 1,
        drug_clearance_rate: float = 0.12,
        product_clearance_rate: float = 0.18,
        metabolism_scale: float = 0.85,
        effect_scale: float = 0.55,
        ecology_strength: float = 0.20,
        experimental_multi_product_enabled: bool = False,
        experimental_branching_scale: float = 0.35,
        experimental_secondary_metabolism_rate: float = 0.10,
    ) -> dict[str, object]:
        drug_id = self._resolve_drug_id(drug_query)
        resolved_community_path = self._resolve_community_table_path(community_table_path)
        if disease_name and str(disease_name).strip():
            disease_target_profile, disease_target_metadata = self._build_disease_target_profile_from_frame(
                self._drug_frame(drug_query),
                disease_name=disease_name,
            )
        else:
            disease_target_profile, disease_target_metadata = {}, {"enabled": False, "reason": "disease_not_specified"}

        with tempfile.TemporaryDirectory(prefix="gut_step3_web_", dir=self.temp_root) as temp_dir_name:
            if resolved_community_path is None and disease_name and str(disease_name).strip():
                resolved_community_path = self._write_disease_adjusted_community(
                    Path(temp_dir_name) / f"{_canonicalize_key(disease_name) or 'disease'}_community.csv",
                    disease_name=str(disease_name).strip(),
                )
            if resolved_community_path is None and scenario_name not in BUILTIN_SCENARIOS:
                raise ValueError(f"不支持的 scenario_name: {scenario_name}")
            active_dir = Path(temp_dir_name) / "active"
            placebo_dir = Path(temp_dir_name) / "placebo"
            active_summary = run_step3_simulation(
                integrated_predictions_path=self.integrated_predictions_path,
                output_dir=active_dir,
                drug_query=drug_id,
                scenario_name=scenario_name,
                community_table_path=resolved_community_path,
                tcg_proxy_mapping_path=self.step3_health_signature_proxy_path,
                n_steps=int(n_steps),
                initial_dose=float(initial_dose),
                repeat_dose=float(repeat_dose),
                dosing_interval=int(dosing_interval),
                drug_clearance_rate=float(drug_clearance_rate),
                product_clearance_rate=float(product_clearance_rate),
                metabolism_scale=float(metabolism_scale),
                effect_scale=float(effect_scale),
                ecology_strength=float(ecology_strength),
                disease_target_profile=disease_target_profile or None,
                experimental_multi_product_enabled=bool(experimental_multi_product_enabled),
                experimental_branching_scale=float(experimental_branching_scale),
                experimental_secondary_metabolism_rate=float(experimental_secondary_metabolism_rate),
            )
            placebo_summary = run_step3_simulation(
                integrated_predictions_path=self.integrated_predictions_path,
                output_dir=placebo_dir,
                drug_query=drug_id,
                scenario_name=scenario_name,
                community_table_path=resolved_community_path,
                tcg_proxy_mapping_path=self.step3_health_signature_proxy_path,
                n_steps=int(n_steps),
                initial_dose=0.0,
                repeat_dose=0.0,
                dosing_interval=int(dosing_interval),
                drug_clearance_rate=float(drug_clearance_rate),
                product_clearance_rate=float(product_clearance_rate),
                metabolism_scale=float(metabolism_scale),
                effect_scale=float(effect_scale),
                ecology_strength=float(ecology_strength),
                disease_target_profile=disease_target_profile or None,
                experimental_multi_product_enabled=bool(experimental_multi_product_enabled),
                experimental_branching_scale=float(experimental_branching_scale),
                experimental_secondary_metabolism_rate=float(experimental_secondary_metabolism_rate),
            )
            summary = _attach_placebo_summary_deltas(active_summary, placebo_summary)
            summary["disease_target_profile_metadata"] = disease_target_metadata
            trajectory_metrics = pd.read_csv(active_summary["trajectory_metrics_path"], low_memory=False)
            placebo_trajectory_metrics = pd.read_csv(placebo_summary["trajectory_metrics_path"], low_memory=False)
            trajectory_metrics = _inject_placebo_deltas(trajectory_metrics, placebo_trajectory_metrics)
            top_microbe_changes = pd.read_csv(active_summary["top_microbe_changes_path"], low_memory=False)

        keep_metric_columns = [
            column
            for column in [
                "timepoint",
                "health_index",
                "parent_retention_ratio",
                "aggregate_metabolite_pool",
                "experimental_aggregate_metabolite_pool",
                "development_score",
                "experimental_development_score",
                "development_score_legacy",
                "development_score_balance",
                "experimental_development_score_balance",
                "benefit_subscore",
                "risk_subscore",
                "experimental_risk_subscore",
                "dysbiosis_penalty",
                "interaction_dysbiosis_penalty",
                "uncertainty_penalty",
                "metabolite_burden_penalty",
                "experimental_metabolite_burden_penalty",
                "disease_target_alignment_score",
                "disease_target_coverage",
                "mean_applicability",
                "diversity",
                "beneficial_fraction",
                "risk_fraction",
                "health_index_legacy",
                "interaction_component",
                "interaction_balance_rho",
                "interaction_balance_shift",
                "positive_interaction_strength",
                "negative_interaction_strength",
                "tcg_health_index",
                "tcg_guild_1_fraction",
                "tcg_guild_2_fraction",
                "tcg_mapped_fraction",
                "development_score_delta_vs_placebo",
                "development_score_normalized_vs_placebo",
                "experimental_development_score_delta_vs_placebo",
                "development_score_balance_delta_vs_placebo",
                "experimental_development_score_balance_delta_vs_placebo",
                "health_index_delta_vs_placebo",
                "benefit_subscore_delta_vs_placebo",
                "risk_subscore_delta_vs_placebo",
                "experimental_risk_subscore_delta_vs_placebo",
                "parent_retention_ratio_delta_vs_placebo",
                "experimental_aggregate_metabolite_pool_delta_vs_placebo",
                "disease_target_alignment_score_delta_vs_placebo",
            ]
            if column in trajectory_metrics.columns
        ]
        keep_change_columns = [
            column
            for column in [
                "nt_code",
                "species_label",
                "initial_abundance",
                "final_abundance",
                "delta_abundance",
                "fold_change",
            ]
            if column in top_microbe_changes.columns
        ]

        return _clean_json(
            {
                "summary": {
                    key: value
                    for key, value in summary.items()
                    if not str(key).endswith("_path")
                },
                "disease_name": disease_name,
                "trajectory_metrics": trajectory_metrics.loc[:, keep_metric_columns],
                "top_microbe_changes": top_microbe_changes.loc[:, keep_change_columns],
            }
        )

    def simulate_custom_step3(
        self,
        session_id: str,
        scenario_name: str = "healthy_reference",
        community_table_path: str | None = None,
        disease_name: str | None = None,
        n_steps: int = 14,
        initial_dose: float = 1.0,
        repeat_dose: float = 1.0,
        dosing_interval: int = 1,
        drug_clearance_rate: float = 0.12,
        product_clearance_rate: float = 0.18,
        metabolism_scale: float = 0.85,
        effect_scale: float = 0.55,
        ecology_strength: float = 0.20,
        experimental_multi_product_enabled: bool = False,
        experimental_branching_scale: float = 0.35,
        experimental_secondary_metabolism_rate: float = 0.10,
    ) -> dict[str, object]:
        session = self._custom_session(session_id)
        integrated_predictions_path = Path(session["integrated_predictions_path"])  # type: ignore[arg-type]
        drug_id = str(session["drug_id"])
        resolved_community_path = self._resolve_community_table_path(community_table_path)
        if disease_name and str(disease_name).strip():
            disease_target_profile, disease_target_metadata = self._build_disease_target_profile_from_frame(
                self._custom_frame(session_id),
                disease_name=disease_name,
            )
        else:
            disease_target_profile, disease_target_metadata = {}, {"enabled": False, "reason": "disease_not_specified"}

        with tempfile.TemporaryDirectory(prefix="gut_step3_custom_", dir=self.temp_root) as temp_dir_name:
            if resolved_community_path is None and disease_name and str(disease_name).strip():
                resolved_community_path = self._write_disease_adjusted_community(
                    Path(temp_dir_name) / f"{_canonicalize_key(disease_name) or 'disease'}_community.csv",
                    disease_name=str(disease_name).strip(),
                )
            if resolved_community_path is None and scenario_name not in BUILTIN_SCENARIOS:
                raise ValueError(f"不支持的 scenario_name: {scenario_name}")
            active_dir = Path(temp_dir_name) / "active"
            placebo_dir = Path(temp_dir_name) / "placebo"
            active_summary = run_step3_simulation(
                integrated_predictions_path=integrated_predictions_path,
                output_dir=active_dir,
                drug_query=drug_id,
                scenario_name=scenario_name,
                community_table_path=resolved_community_path,
                tcg_proxy_mapping_path=self.step3_health_signature_proxy_path,
                n_steps=int(n_steps),
                initial_dose=float(initial_dose),
                repeat_dose=float(repeat_dose),
                dosing_interval=int(dosing_interval),
                drug_clearance_rate=float(drug_clearance_rate),
                product_clearance_rate=float(product_clearance_rate),
                metabolism_scale=float(metabolism_scale),
                effect_scale=float(effect_scale),
                ecology_strength=float(ecology_strength),
                disease_target_profile=disease_target_profile or None,
                experimental_multi_product_enabled=bool(experimental_multi_product_enabled),
                experimental_branching_scale=float(experimental_branching_scale),
                experimental_secondary_metabolism_rate=float(experimental_secondary_metabolism_rate),
            )
            placebo_summary = run_step3_simulation(
                integrated_predictions_path=integrated_predictions_path,
                output_dir=placebo_dir,
                drug_query=drug_id,
                scenario_name=scenario_name,
                community_table_path=resolved_community_path,
                tcg_proxy_mapping_path=self.step3_health_signature_proxy_path,
                n_steps=int(n_steps),
                initial_dose=0.0,
                repeat_dose=0.0,
                dosing_interval=int(dosing_interval),
                drug_clearance_rate=float(drug_clearance_rate),
                product_clearance_rate=float(product_clearance_rate),
                metabolism_scale=float(metabolism_scale),
                effect_scale=float(effect_scale),
                ecology_strength=float(ecology_strength),
                disease_target_profile=disease_target_profile or None,
                experimental_multi_product_enabled=bool(experimental_multi_product_enabled),
                experimental_branching_scale=float(experimental_branching_scale),
                experimental_secondary_metabolism_rate=float(experimental_secondary_metabolism_rate),
            )
            summary = _attach_placebo_summary_deltas(active_summary, placebo_summary)
            summary["disease_target_profile_metadata"] = disease_target_metadata
            trajectory_metrics = pd.read_csv(active_summary["trajectory_metrics_path"], low_memory=False)
            placebo_trajectory_metrics = pd.read_csv(placebo_summary["trajectory_metrics_path"], low_memory=False)
            trajectory_metrics = _inject_placebo_deltas(trajectory_metrics, placebo_trajectory_metrics)
            top_microbe_changes = pd.read_csv(active_summary["top_microbe_changes_path"], low_memory=False)

        keep_metric_columns = [
            column
            for column in [
                "timepoint",
                "health_index",
                "parent_retention_ratio",
                "aggregate_metabolite_pool",
                "experimental_aggregate_metabolite_pool",
                "development_score",
                "experimental_development_score",
                "development_score_legacy",
                "development_score_balance",
                "experimental_development_score_balance",
                "benefit_subscore",
                "risk_subscore",
                "experimental_risk_subscore",
                "dysbiosis_penalty",
                "interaction_dysbiosis_penalty",
                "uncertainty_penalty",
                "metabolite_burden_penalty",
                "experimental_metabolite_burden_penalty",
                "disease_target_alignment_score",
                "disease_target_coverage",
                "mean_applicability",
                "diversity",
                "beneficial_fraction",
                "risk_fraction",
                "health_index_legacy",
                "interaction_component",
                "interaction_balance_rho",
                "interaction_balance_shift",
                "positive_interaction_strength",
                "negative_interaction_strength",
                "tcg_health_index",
                "tcg_guild_1_fraction",
                "tcg_guild_2_fraction",
                "tcg_mapped_fraction",
                "development_score_delta_vs_placebo",
                "development_score_normalized_vs_placebo",
                "experimental_development_score_delta_vs_placebo",
                "development_score_balance_delta_vs_placebo",
                "experimental_development_score_balance_delta_vs_placebo",
                "health_index_delta_vs_placebo",
                "benefit_subscore_delta_vs_placebo",
                "risk_subscore_delta_vs_placebo",
                "experimental_risk_subscore_delta_vs_placebo",
                "parent_retention_ratio_delta_vs_placebo",
                "experimental_aggregate_metabolite_pool_delta_vs_placebo",
                "disease_target_alignment_score_delta_vs_placebo",
            ]
            if column in trajectory_metrics.columns
        ]
        keep_change_columns = [
            column
            for column in [
                "nt_code",
                "species_label",
                "initial_abundance",
                "final_abundance",
                "delta_abundance",
                "fold_change",
            ]
            if column in top_microbe_changes.columns
        ]
        return _clean_json(
            {
                "summary": {
                    key: value
                    for key, value in summary.items()
                    if not str(key).endswith("_path")
                },
                "disease_name": disease_name,
                "trajectory_metrics": trajectory_metrics.loc[:, keep_metric_columns],
                "top_microbe_changes": top_microbe_changes.loc[:, keep_change_columns],
                "session_id": session_id,
            }
        )

    def scenario_grid_step3(
        self,
        drug_query: str,
        community_table_path: str | None = None,
        disease_name: str | None = None,
        n_steps: int = 14,
        initial_dose: float = 1.0,
        repeat_dose: float = 1.0,
        dosing_interval: int = 1,
        drug_clearance_rate: float = 0.12,
        product_clearance_rate: float = 0.18,
        metabolism_scale: float = 0.85,
        effect_scale: float = 0.55,
        ecology_strength: float = 0.20,
        experimental_multi_product_enabled: bool = False,
        experimental_branching_scale: float = 0.35,
        experimental_secondary_metabolism_rate: float = 0.10,
    ) -> dict[str, object]:
        summaries: list[dict[str, object]] = []
        if disease_name or self._resolve_community_table_path(community_table_path) is not None:
            result = self.simulate_step3(
                drug_query=drug_query,
                community_table_path=community_table_path,
                disease_name=disease_name,
                n_steps=n_steps,
                initial_dose=initial_dose,
                repeat_dose=repeat_dose,
                dosing_interval=dosing_interval,
                drug_clearance_rate=drug_clearance_rate,
                product_clearance_rate=product_clearance_rate,
                metabolism_scale=metabolism_scale,
                effect_scale=effect_scale,
                ecology_strength=ecology_strength,
                experimental_multi_product_enabled=experimental_multi_product_enabled,
                experimental_branching_scale=experimental_branching_scale,
                experimental_secondary_metabolism_rate=experimental_secondary_metabolism_rate,
            )
            summaries.append(dict(result["summary"]))
        else:
            for scenario_name in sorted(BUILTIN_SCENARIOS):
                result = self.simulate_step3(
                    drug_query=drug_query,
                    scenario_name=scenario_name,
                    n_steps=n_steps,
                    initial_dose=initial_dose,
                    repeat_dose=repeat_dose,
                    dosing_interval=dosing_interval,
                    drug_clearance_rate=drug_clearance_rate,
                    product_clearance_rate=product_clearance_rate,
                    metabolism_scale=metabolism_scale,
                    effect_scale=effect_scale,
                    ecology_strength=ecology_strength,
                    experimental_multi_product_enabled=experimental_multi_product_enabled,
                    experimental_branching_scale=experimental_branching_scale,
                    experimental_secondary_metabolism_rate=experimental_secondary_metabolism_rate,
                )
                summary = dict(result["summary"])
                summaries.append(summary)

        return _clean_json(
            {
                "drug_query": self._resolve_drug_id(drug_query),
                "scenario_summaries": summaries,
            }
        )

    def scenario_grid_custom_step3(
        self,
        session_id: str,
        community_table_path: str | None = None,
        disease_name: str | None = None,
        n_steps: int = 14,
        initial_dose: float = 1.0,
        repeat_dose: float = 1.0,
        dosing_interval: int = 1,
        drug_clearance_rate: float = 0.12,
        product_clearance_rate: float = 0.18,
        metabolism_scale: float = 0.85,
        effect_scale: float = 0.55,
        ecology_strength: float = 0.20,
        experimental_multi_product_enabled: bool = False,
        experimental_branching_scale: float = 0.35,
        experimental_secondary_metabolism_rate: float = 0.10,
    ) -> dict[str, object]:
        summaries: list[dict[str, object]] = []
        session = self._custom_session(session_id)
        if disease_name or self._resolve_community_table_path(community_table_path) is not None:
            result = self.simulate_custom_step3(
                session_id=session_id,
                community_table_path=community_table_path,
                disease_name=disease_name,
                n_steps=n_steps,
                initial_dose=initial_dose,
                repeat_dose=repeat_dose,
                dosing_interval=dosing_interval,
                drug_clearance_rate=drug_clearance_rate,
                product_clearance_rate=product_clearance_rate,
                metabolism_scale=metabolism_scale,
                effect_scale=effect_scale,
                ecology_strength=ecology_strength,
                experimental_multi_product_enabled=experimental_multi_product_enabled,
                experimental_branching_scale=experimental_branching_scale,
                experimental_secondary_metabolism_rate=experimental_secondary_metabolism_rate,
            )
            summaries.append(dict(result["summary"]))
        else:
            for scenario_name in sorted(BUILTIN_SCENARIOS):
                result = self.simulate_custom_step3(
                    session_id=session_id,
                    scenario_name=scenario_name,
                    n_steps=n_steps,
                    initial_dose=initial_dose,
                    repeat_dose=repeat_dose,
                    dosing_interval=dosing_interval,
                    drug_clearance_rate=drug_clearance_rate,
                    product_clearance_rate=product_clearance_rate,
                    metabolism_scale=metabolism_scale,
                    effect_scale=effect_scale,
                    ecology_strength=ecology_strength,
                    experimental_multi_product_enabled=experimental_multi_product_enabled,
                    experimental_branching_scale=experimental_branching_scale,
                    experimental_secondary_metabolism_rate=experimental_secondary_metabolism_rate,
                )
                summaries.append(dict(result["summary"]))

        return _clean_json(
            {
                "session_id": session_id,
                "drug_query": str(session["drug_id"]),
                "scenario_summaries": summaries,
            }
        )
