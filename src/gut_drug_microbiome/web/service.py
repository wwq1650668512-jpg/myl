from __future__ import annotations

import math
import tempfile
import uuid
from pathlib import Path

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
DEFAULT_STEP3_COHORT_ROOT = ROOT / "data/processed/step3/cohorts"
DEFAULT_STEP3_HEALTH_SIGNATURE_PROXY_PATH = ROOT / "data/processed/health_signature/microbe_tcg_proxy_mapping.csv"
DEFAULT_DISEASE_MICROBE_REFERENCE_PATH = ROOT / "data/reference/disease_microbe_dictionary.csv"
DEFAULT_DISEASE_DRUG_REFERENCE_PATH = ROOT / "data/reference/disease_marketed_drug_catalog.csv"

REQUIRED_DISEASE_CATALOG_ENTRIES = [
    "肠易激综合征（IBS）",
    "肠易激综合征-腹泻型（IBS-D）",
    "肠易激综合征-便秘型（IBS-C）",
]

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
        step3_cohort_root: str | Path = DEFAULT_STEP3_COHORT_ROOT,
        step3_health_signature_proxy_path: str | Path = DEFAULT_STEP3_HEALTH_SIGNATURE_PROXY_PATH,
        disease_microbe_reference_path: str | Path | None = DEFAULT_DISEASE_MICROBE_REFERENCE_PATH,
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
        self.step3_cohort_root = Path(step3_cohort_root)
        self.step3_health_signature_proxy_path = Path(step3_health_signature_proxy_path)
        self.disease_microbe_reference_path = (
            None if disease_microbe_reference_path is None else Path(disease_microbe_reference_path)
        )
        self.disease_drug_reference_path = None if disease_drug_reference_path is None else Path(disease_drug_reference_path)

        usecols = _existing_usecols(self.integrated_predictions_path, DESIRED_COLUMNS)
        self.frame = pd.read_csv(self.integrated_predictions_path, usecols=usecols, low_memory=False)
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
        self.drug_table["search_key"] = self.drug_table["prestwick_id"].map(_canonicalize_key)
        self.drug_table["name_key"] = self.drug_table["chemical_name"].map(_canonicalize_key)

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
        self.disease_microbe_reference = self._load_optional_reference(self.disease_microbe_reference_path)
        self.disease_drug_reference = self._load_optional_reference(self.disease_drug_reference_path)
        self.disease_catalog = self._build_disease_catalog()
        self.custom_sessions: dict[str, dict[str, object]] = {}

    def _load_optional_reference(self, path: Path | None) -> pd.DataFrame:
        """Load an optional CSV reference table, returning an empty frame when unavailable."""
        if path is None or not path.exists():
            return pd.DataFrame()
        return pd.read_csv(path, low_memory=False)

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

        catalog: list[dict[str, object]] = []
        for disease_name in sorted(disease_names):
            disease_key = _canonicalize_key(disease_name)
            microbe_count = 0
            if not self.disease_microbe_reference.empty:
                microbe_count = int(
                    self.disease_microbe_reference["disease_name"].map(_canonicalize_key).eq(disease_key).sum()
                )
            marketed_count = 0
            if not self.disease_drug_reference.empty:
                marketed_count = int(
                    self.disease_drug_reference["disease_name"].map(_canonicalize_key).eq(disease_key).sum()
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
            return work["species_label"].map(_canonicalize_key).eq(microbe_key) | work["microbe_label"].map(_canonicalize_key).eq(microbe_key)
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

    def _candidate_diseases_from_frame(self, work: pd.DataFrame) -> list[dict[str, object]]:
        """Score diseases whose curated microbe patterns are directionally consistent with Step 1 outputs."""
        if self.disease_microbe_reference.empty:
            return []

        disease_records: list[dict[str, object]] = []
        level_weights = {"species": 1.0, "genus": 0.75, "family": 0.45, "phylum": 0.30, "class": 0.25, "order": 0.25}
        source_weights = {"microbe_to_disease": 1.0, "disease_to_microbe": 0.7}
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
            relation_scores: list[float] = []
            evidence_rows: list[dict[str, object]] = []
            mechanism_rows: list[pd.DataFrame] = []
            for _, relation in group.iterrows():
                desired_effect = str(relation.get("desired_step1_effect", "unknown"))
                if desired_effect not in {"promote", "inhibit"}:
                    continue
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
                weight = level_weights.get(str(relation.get("taxon_level", "unknown")), 0.2) * source_weights.get(
                    str(relation.get("source_sheet", "")),
                    0.5,
                )
                weighted_score = relation_score * weight
                relation_scores.append(weighted_score)
                matched = matched.copy()
                matched["relation_weight"] = float(weight)
                matched["desired_step1_effect"] = desired_effect
                mechanism_rows.append(matched)
                strongest = matched.iloc[0]
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
                        "taxon_level": relation.get("taxon_level"),
                        "source_sheet": relation.get("source_sheet"),
                    }
                )

            if not relation_scores:
                continue

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
            default_fusion_mode = "weighted_0.65_0.35"
            disease_score_mechanism = fuse_disease_scores(
                raw_score=raw_microbe_score,
                mechanism_score=mechanism_balance,
                fusion_mode=default_fusion_mode,
            )
            marketed_examples = []
            if not self.disease_drug_reference.empty:
                disease_key = _canonicalize_key(disease_name)
                marketed_examples = (
                    self.disease_drug_reference.loc[
                        self.disease_drug_reference["disease_name"].map(_canonicalize_key).eq(disease_key),
                        "marketed_drug_name_raw",
                    ]
                    .dropna()
                    .astype(str)
                    .head(5)
                    .tolist()
                )
            disease_records.append(
                {
                    "disease_name": disease_name,
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
                    "marketed_drug_examples": marketed_examples,
                    "evidence_examples": evidence_rows[:5],
                }
            )

        disease_records.sort(
            key=lambda item: (
                item["disease_score_mechanism"],
                item["disease_score_raw_only"],
                item["matched_relation_count"],
            ),
            reverse=True,
        )
        return disease_records[:8]

    def _write_disease_adjusted_community(self, output_path: Path, disease_name: str) -> Path:
        """Create a temporary community table reflecting one curated disease profile."""
        if self.disease_microbe_reference.empty:
            raise ValueError("当前没有可用的疾病-微生物参考表。")
        community = build_disease_adjusted_community(
            microbe_metadata=self.microbe_feature_table,
            disease_name=disease_name,
            disease_microbe_reference=self.disease_microbe_reference,
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        community.to_csv(output_path, index=False)
        return output_path

    def _profile_from_frame(self, frame: pd.DataFrame) -> dict[str, object]:
        work = self._annotate_step2_mechanism(self._annotate_amr(frame))
        candidate_diseases = self._candidate_diseases_from_frame(work)
        row = work.iloc[0]

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
                    "predicted_reaction_class": metabolism_row.get("predicted_reaction_class"),
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
                    "predicted_reaction_class": metabolism_row.get("predicted_reaction_class"),
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
        }

        return _clean_json(
            {
                "drug": self._drug_metadata(row),
                "aggregated": aggregated,
                "candidate_diseases": candidate_diseases,
                "marketed_disease_context": self._marketed_disease_context(row.get("chemical_name")),
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
        )

        integrated_frame = pd.read_csv(step2_output_dir / "predictions.csv", low_memory=False)
        integrated_frame = refine_step1_promote_with_step2(
            integrated_frame,
            promote_classifier_path=self.step1_promote_classifier_path,
            promote_metrics_path=self.step1_promote_metrics_path,
            cross_feeding_reference_path=self.cross_feeding_reference_path,
        )
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
    ) -> dict[str, object]:
        drug_id = self._resolve_drug_id(drug_query)
        resolved_community_path = self._resolve_community_table_path(community_table_path)

        with tempfile.TemporaryDirectory(prefix="gut_step3_web_", dir=self.temp_root) as temp_dir_name:
            if resolved_community_path is None and disease_name and str(disease_name).strip():
                resolved_community_path = self._write_disease_adjusted_community(
                    Path(temp_dir_name) / f"{_canonicalize_key(disease_name) or 'disease'}_community.csv",
                    disease_name=str(disease_name).strip(),
                )
            if resolved_community_path is None and scenario_name not in BUILTIN_SCENARIOS:
                raise ValueError(f"不支持的 scenario_name: {scenario_name}")
            summary = run_step3_simulation(
                integrated_predictions_path=self.integrated_predictions_path,
                output_dir=Path(temp_dir_name),
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
            )
            trajectory_metrics = pd.read_csv(summary["trajectory_metrics_path"], low_memory=False)
            top_microbe_changes = pd.read_csv(summary["top_microbe_changes_path"], low_memory=False)

        keep_metric_columns = [
            column
            for column in [
                "timepoint",
                "health_index",
                "parent_retention_ratio",
                "aggregate_metabolite_pool",
                "development_score",
                "development_score_legacy",
                "development_score_balance",
                "benefit_subscore",
                "risk_subscore",
                "dysbiosis_penalty",
                "interaction_dysbiosis_penalty",
                "uncertainty_penalty",
                "metabolite_burden_penalty",
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
    ) -> dict[str, object]:
        session = self._custom_session(session_id)
        integrated_predictions_path = Path(session["integrated_predictions_path"])  # type: ignore[arg-type]
        drug_id = str(session["drug_id"])
        resolved_community_path = self._resolve_community_table_path(community_table_path)

        with tempfile.TemporaryDirectory(prefix="gut_step3_custom_", dir=self.temp_root) as temp_dir_name:
            if resolved_community_path is None and disease_name and str(disease_name).strip():
                resolved_community_path = self._write_disease_adjusted_community(
                    Path(temp_dir_name) / f"{_canonicalize_key(disease_name) or 'disease'}_community.csv",
                    disease_name=str(disease_name).strip(),
                )
            if resolved_community_path is None and scenario_name not in BUILTIN_SCENARIOS:
                raise ValueError(f"不支持的 scenario_name: {scenario_name}")
            summary = run_step3_simulation(
                integrated_predictions_path=integrated_predictions_path,
                output_dir=Path(temp_dir_name),
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
            )
            trajectory_metrics = pd.read_csv(summary["trajectory_metrics_path"], low_memory=False)
            top_microbe_changes = pd.read_csv(summary["top_microbe_changes_path"], low_memory=False)

        keep_metric_columns = [
            column
            for column in [
                "timepoint",
                "health_index",
                "parent_retention_ratio",
                "aggregate_metabolite_pool",
                "development_score",
                "development_score_legacy",
                "development_score_balance",
                "benefit_subscore",
                "risk_subscore",
                "dysbiosis_penalty",
                "interaction_dysbiosis_penalty",
                "uncertainty_penalty",
                "metabolite_burden_penalty",
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
                )
                summaries.append(dict(result["summary"]))

        return _clean_json(
            {
                "session_id": session_id,
                "drug_query": str(session["drug_id"]),
                "scenario_summaries": summaries,
            }
        )
