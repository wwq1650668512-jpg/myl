from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

from gut_drug_microbiome.utils.text import canonicalize_key as _canonicalize_key
from gut_drug_microbiome.utils.text import normalize_whitespace as _canonicalize_text

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_TCG_PROXY_MAPPING_PATH = ROOT / "data/processed/health_signature/microbe_tcg_proxy_mapping.csv"
DEFAULT_CROSS_FEEDING_EDGES_PATH = ROOT / "data/reference/cross_feeding_edges.csv"
DEFAULT_ENZYME_FUNCTION_CATALOG_PATH = ROOT / "data/reference/step2_enzyme_function_catalog.csv"


BENEFICIAL_GENERA = {
    "Akkermansia",
    "Anaerostipes",
    "Bifidobacterium",
    "Blautia",
    "Coprococcus",
    "Eubacterium",
    "Faecalibacterium",
    "Parabacteroides",
    "Roseburia",
    "Ruminococcus",
}

RISK_GENERA = {
    "Clostridium",
    "Enterococcus",
    "Escherichia",
    "Fusobacterium",
    "Klebsiella",
    "Proteus",
    "Streptococcus",
}


BUILTIN_SCENARIOS: dict[str, dict[str, object]] = {
    "healthy_reference": {
        "description": "Balanced reference community mapped onto the available Step 1 strain panel.",
        "default_weight": 1.0,
        "concentration_power": 1.0,
        "genus_multipliers": {
            "Bifidobacterium": 2.6,
            "Blautia": 2.1,
            "Faecalibacterium": 2.8,
            "Roseburia": 2.7,
            "Ruminococcus": 1.8,
            "Akkermansia": 2.0,
            "Escherichia": 0.4,
            "Proteus": 0.4,
            "Klebsiella": 0.4,
        },
        "phylum_multipliers": {
            "Firmicutes": 1.3,
            "Actinobacteria": 1.2,
            "Proteobacteria": 0.5,
        },
    },
    "high_fiber": {
        "description": "Fiber-rich community proxy enriched for SCFA-associated commensals.",
        "default_weight": 1.0,
        "concentration_power": 1.05,
        "genus_multipliers": {
            "Bifidobacterium": 3.0,
            "Blautia": 2.4,
            "Faecalibacterium": 3.2,
            "Roseburia": 3.0,
            "Ruminococcus": 2.1,
            "Parabacteroides": 1.6,
            "Escherichia": 0.35,
            "Proteus": 0.35,
        },
        "phylum_multipliers": {
            "Firmicutes": 1.4,
            "Actinobacteria": 1.4,
            "Proteobacteria": 0.4,
        },
    },
    "high_fat": {
        "description": "High-fat westernized community proxy with lower beneficial commensal support.",
        "default_weight": 1.0,
        "concentration_power": 1.15,
        "genus_multipliers": {
            "Bacteroides": 1.8,
            "Alistipes": 1.6,
            "Bilophila": 1.8,
            "Escherichia": 1.4,
            "Bifidobacterium": 0.6,
            "Faecalibacterium": 0.65,
            "Roseburia": 0.65,
            "Akkermansia": 0.9,
        },
        "phylum_multipliers": {
            "Bacteroidetes": 1.3,
            "Proteobacteria": 1.2,
            "Actinobacteria": 0.7,
        },
    },
    "antibiotic_perturbed": {
        "description": "Low-diversity dysbiotic proxy after broad antibiotic exposure.",
        "default_weight": 1.0,
        "concentration_power": 1.35,
        "genus_multipliers": {
            "Escherichia": 2.8,
            "Proteus": 2.6,
            "Klebsiella": 2.6,
            "Enterococcus": 2.4,
            "Clostridium": 1.6,
            "Bifidobacterium": 0.4,
            "Faecalibacterium": 0.35,
            "Roseburia": 0.35,
            "Blautia": 0.45,
        },
        "phylum_multipliers": {
            "Proteobacteria": 2.0,
            "Firmicutes": 0.9,
            "Actinobacteria": 0.5,
        },
    },
}

def _normalize_abundances(abundances: dict[str, float], minimum: float = 0.0) -> dict[str, float]:
    """Clip, renormalize, and return a valid relative-abundance distribution."""
    clipped = {key: max(float(value), minimum) for key, value in abundances.items() if float(value) > 0 or minimum > 0}
    total = sum(clipped.values())
    if total <= 0:
        return {}
    return {key: value / total for key, value in clipped.items()}


def _parse_multi_value(value: object) -> list[str]:
    """Split a semicolon-delimited field into a list of non-empty items."""
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    return [item.strip() for item in text.split(";") if item.strip()]


def _parse_flexible_multi_value(value: object) -> list[str]:
    """Split a flexible delimited field into a list of non-empty items."""
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    return [item.strip() for item in re.split(r"[;|]", text) if item.strip()]


def _coerce_float(value: object, default: float = 0.0) -> float:
    """Convert a scalar-like value to float, falling back to a default."""
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(numeric):
        return float(default)
    return float(numeric)


def _coerce_bool_text(value: object) -> bool:
    """Interpret common truthy strings and booleans as a Python bool."""
    if pd.isna(value):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return _canonicalize_key(value) in {"1", "true", "yes", "y", "ready"}


def _clip01(value: float) -> float:
    """Clip a scalar into the closed interval [0, 1]."""
    return float(np.clip(float(value), 0.0, 1.0))


def _optional_column_map(frame: pd.DataFrame, candidates: list[str]) -> dict[str, object]:
    """Return the first available nt_code-indexed value map from a list of candidate columns."""
    for column in candidates:
        if column in frame.columns:
            return frame.set_index("nt_code")[column].to_dict()
    return {}


def _evidence_level_weight(value: object) -> float:
    """Map qualitative evidence tiers onto a small numeric scale."""
    key = _canonicalize_key(value)
    if key == "high":
        return 1.0
    if key == "medium":
        return 0.7
    if key == "low":
        return 0.45
    return 0.3


def _jaccard_similarity(left: set[str], right: set[str]) -> float:
    """Compute a Jaccard similarity score between two finite sets."""
    if not left or not right:
        return 0.0
    union = left | right
    if not union:
        return 0.0
    return float(len(left & right) / len(union))


def _first_nonempty_text(frame: pd.DataFrame, column: str, default: str) -> str:
    """Return the first non-empty text value from a column, else a default string."""
    if column not in frame.columns:
        return _canonicalize_text(default)
    values = frame[column].dropna().astype(str).map(str.strip)
    values = values[values.ne("")]
    if values.empty:
        return _canonicalize_text(default)
    return _canonicalize_text(values.iloc[0])


def _scenario_names() -> list[str]:
    """Return the sorted names of all builtin Step 3 simulation scenarios."""
    return sorted(BUILTIN_SCENARIOS)


def _build_drug_context_keys(drug_frame: pd.DataFrame) -> set[str]:
    """Collect flexible compound-name and semantic-family keys for reference matching."""
    keys: set[str] = set()
    for column in [
        "chemical_name",
        "compound_name_normalized",
        "compound_semantic_family",
        "compound_semantic_aliases",
        "compound_semantic_keywords",
    ]:
        if column not in drug_frame.columns:
            continue
        for value in drug_frame[column].dropna().astype(str):
            key = _canonicalize_key(value)
            if key:
                keys.add(key)
            for item in _parse_flexible_multi_value(value):
                item_key = _canonicalize_key(item)
                if item_key:
                    keys.add(item_key)
    return keys


def _build_microbe_metadata(frame: pd.DataFrame) -> pd.DataFrame:
    """Extract one-row-per-microbe metadata used throughout the simulation pipeline."""
    columns = [
        column
        for column in [
            "nt_code",
            "microbe_label",
            "species_label",
            "species_name",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "gram_stain",
        ]
        if column in frame.columns
    ]
    metadata = frame.loc[:, columns].drop_duplicates(subset=["nt_code"]).reset_index(drop=True)
    return metadata


def _build_builtin_community(
    microbe_metadata: pd.DataFrame,
    scenario_name: str,
) -> tuple[dict[str, float], dict[str, object]]:
    """Create a normalized starting community from one of the builtin scenario templates."""
    if scenario_name not in BUILTIN_SCENARIOS:
        raise ValueError(f"Unsupported scenario_name: {scenario_name}. Available: {_scenario_names()}")

    profile = BUILTIN_SCENARIOS[scenario_name]
    default_weight = float(profile.get("default_weight", 1.0))
    concentration_power = float(profile.get("concentration_power", 1.0))
    genus_multipliers = {str(key): float(value) for key, value in profile.get("genus_multipliers", {}).items()}
    phylum_multipliers = {str(key): float(value) for key, value in profile.get("phylum_multipliers", {}).items()}

    weights: dict[str, float] = {}
    for _, row in microbe_metadata.iterrows():
        nt_code = str(row["nt_code"])
        genus = _canonicalize_text(row.get("genus"))
        phylum = _canonicalize_text(row.get("phylum"))
        weight = default_weight
        if genus:
            weight *= genus_multipliers.get(genus, 1.0)
        if phylum:
            weight *= phylum_multipliers.get(phylum, 1.0)
        weights[nt_code] = max(weight, 1e-6)

    normalized = _normalize_abundances(weights)
    if concentration_power != 1.0:
        normalized = _normalize_abundances(
            {key: value**concentration_power for key, value in normalized.items()},
            minimum=1e-8,
        )

    metadata = {
        "scenario_name": scenario_name,
        "scenario_description": str(profile.get("description", "")),
        "source": "builtin_panel_proxy",
        "n_microbes": int(len(normalized)),
    }
    return normalized, metadata


def _load_custom_community(
    community_table_path: str | Path,
    microbe_metadata: pd.DataFrame,
) -> tuple[dict[str, float], dict[str, object]]:
    """Load a user-provided community table and map it onto the available microbe panel."""
    community_table_path = Path(community_table_path)
    raw = pd.read_csv(community_table_path, low_memory=False)
    alias_candidates = {
        "microbe_key": ["nt_code", "microbe_id", "microbe_name", "species_label", "microbe_label"],
        "abundance": ["abundance", "relative_abundance", "biomass", "weight"],
    }
    matched: dict[str, str] = {}
    for standard_name, aliases in alias_candidates.items():
        normalized_map = {_canonicalize_key(column): column for column in raw.columns}
        for alias in aliases:
            key = _canonicalize_key(alias)
            if key in normalized_map:
                matched[standard_name] = normalized_map[key]
                break
    if {"microbe_key", "abundance"} - set(matched):
        raise ValueError(
            f"Community table must contain columns compatible with {alias_candidates}; got {raw.columns.tolist()}"
        )

    lookup: dict[str, str] = {}
    for _, row in microbe_metadata.iterrows():
        nt_code = str(row["nt_code"])
        for value in [
            nt_code,
            row.get("microbe_label"),
            row.get("species_label"),
            row.get("species_name"),
        ]:
            key = _canonicalize_key(value)
            if key:
                lookup[key] = nt_code

    abundances: dict[str, float] = {}
    for _, row in raw.iterrows():
        key = _canonicalize_key(row[matched["microbe_key"]])
        if key not in lookup:
            continue
        abundance = pd.to_numeric(pd.Series([row[matched["abundance"]]]), errors="coerce").iloc[0]
        if pd.isna(abundance) or float(abundance) <= 0:
            continue
        abundances[lookup[key]] = abundances.get(lookup[key], 0.0) + float(abundance)

    if not abundances:
        raise RuntimeError("No custom community rows could be mapped onto the available Step 1 microbe panel.")

    return _normalize_abundances(abundances), {
        "scenario_name": community_table_path.stem,
        "scenario_description": "custom community table",
        "source": str(community_table_path),
        "n_microbes": int(len(abundances)),
    }


def _resolve_drug_subset(frame: pd.DataFrame, drug_query: str) -> pd.DataFrame:
    """Resolve a user drug query to the matching subset of prediction rows."""
    available_columns = [column for column in ["prestwick_id", "chemical_name"] if column in frame.columns]
    if not available_columns:
        raise ValueError("Prediction frame must contain at least one of: prestwick_id, chemical_name")

    candidates = frame.loc[:, available_columns].drop_duplicates().reset_index(drop=True)
    query_key = _canonicalize_key(drug_query)
    exact_mask = pd.Series(False, index=candidates.index)
    if "prestwick_id" in candidates.columns:
        exact_mask = exact_mask | candidates["prestwick_id"].map(_canonicalize_key).eq(query_key)
    if "chemical_name" in candidates.columns:
        exact_mask = exact_mask | candidates["chemical_name"].map(_canonicalize_key).eq(query_key)
    if exact_mask.any():
        selected = candidates.loc[exact_mask].iloc[0]
        if "prestwick_id" in frame.columns and "prestwick_id" in selected.index:
            return frame[frame["prestwick_id"] == selected["prestwick_id"]].copy()
        return frame[frame["chemical_name"] == selected["chemical_name"]].copy()

    if "chemical_name" in candidates.columns:
        contains_mask = candidates["chemical_name"].astype(str).str.contains(drug_query, case=False, na=False)
        if contains_mask.sum() == 1:
            selected = candidates.loc[contains_mask].iloc[0]
            if "prestwick_id" in frame.columns and "prestwick_id" in selected.index:
                return frame[frame["prestwick_id"] == selected["prestwick_id"]].copy()
            return frame[frame["chemical_name"] == selected["chemical_name"]].copy()
        if contains_mask.sum() > 1:
            option_columns = [column for column in ["prestwick_id", "chemical_name"] if column in candidates.columns]
            options = candidates.loc[contains_mask, option_columns].head(10).to_dict(orient="records")
            raise ValueError(f"drug_query is ambiguous; examples: {options}")

    raise ValueError(f"Could not resolve drug_query={drug_query!r}")


def _shannon_diversity(abundances: dict[str, float]) -> float:
    """Compute Shannon diversity for a relative-abundance distribution."""
    values = np.asarray([value for value in abundances.values() if value > 0], dtype=float)
    if values.size == 0:
        return 0.0
    return float(-(values * np.log(values)).sum())


def _normalized_shannon(abundances: dict[str, float]) -> float:
    """Compute Shannon diversity normalized by the maximum value for the observed richness."""
    values = np.asarray([value for value in abundances.values() if value > 0], dtype=float)
    if values.size <= 1:
        return 0.0
    shannon = _shannon_diversity(abundances)
    return float(shannon / math.log(values.size))


def _weighted_fraction(
    abundances: dict[str, float],
    microbe_metadata: pd.DataFrame,
    target_genera: set[str],
) -> float:
    """Sum the abundance fraction belonging to a target set of genera."""
    genus_map = microbe_metadata.set_index("nt_code")["genus"].to_dict()
    total = 0.0
    for nt_code, abundance in abundances.items():
        genus = _canonicalize_text(genus_map.get(nt_code))
        if genus in target_genera:
            total += float(abundance)
    return float(total)


def _load_tcg_proxy_mapping(
    tcg_proxy_mapping_path: str | Path | None,
    microbe_metadata: pd.DataFrame,
) -> tuple[dict[str, str], dict[str, object]]:
    """Load optional guild mappings used to compute TCG-style health summaries."""
    empty_metadata = {
        "tcg_proxy_mapping_path": None,
        "tcg_source_name": None,
        "tcg_source_url": None,
        "tcg_ready_microbe_count": 0,
        "tcg_mapped_microbe_count": 0,
        "tcg_panel_microbe_count": int(len(microbe_metadata)),
    }
    if tcg_proxy_mapping_path is None:
        return {}, empty_metadata

    path = Path(tcg_proxy_mapping_path)
    if not path.exists():
        return {}, empty_metadata

    raw = pd.read_csv(path, low_memory=False)
    if raw.empty or "nt_code" not in raw.columns:
        metadata = empty_metadata.copy()
        metadata["tcg_proxy_mapping_path"] = str(path)
        return {}, metadata

    panel_nt_codes = {str(value) for value in microbe_metadata["nt_code"].astype(str)}
    work = raw.copy()
    work["nt_code"] = work["nt_code"].astype(str).str.strip()
    work = work[work["nt_code"].isin(panel_nt_codes)].copy()

    membership_map: dict[str, str] = {}
    for _, row in work.iterrows():
        nt_code = str(row["nt_code"])
        membership_key = _canonicalize_key(row.get("tcg_membership"))
        if membership_key in {"guild1", "guild_1"}:
            canonical_membership = "guild_1"
        elif membership_key in {"guild2", "guild_2"}:
            canonical_membership = "guild_2"
        else:
            continue
        if not _coerce_bool_text(row.get("ready_for_step3")):
            continue
        membership_map[nt_code] = canonical_membership

    source_name = None
    if "tcg_source_name" in work.columns:
        values = work["tcg_source_name"].dropna().astype(str).str.strip()
        values = values[values.ne("")]
        if not values.empty:
            source_name = values.iloc[0]

    source_url = None
    if "tcg_source_url" in work.columns:
        values = work["tcg_source_url"].dropna().astype(str).str.strip()
        values = values[values.ne("")]
        if not values.empty:
            source_url = values.iloc[0]

    metadata = {
        "tcg_proxy_mapping_path": str(path),
        "tcg_source_name": source_name,
        "tcg_source_url": source_url,
        "tcg_ready_microbe_count": int(len(membership_map)),
        "tcg_mapped_microbe_count": int(
            work["tcg_membership"]
            .fillna("")
            .astype(str)
            .str.strip()
            .str.lower()
            .isin(["guild_1", "guild_2", "guild1", "guild2"])
            .sum()
        ),
        "tcg_panel_microbe_count": int(len(microbe_metadata)),
    }
    return membership_map, metadata


def _stability_score(
    abundances: dict[str, float],
    baseline_abundances: dict[str, float],
) -> float:
    """Measure community stability as one minus the L1 distance from baseline."""
    keys = sorted(set(abundances) | set(baseline_abundances))
    l1_distance = sum(abs(abundances.get(key, 0.0) - baseline_abundances.get(key, 0.0)) for key in keys)
    return float(max(0.0, 1.0 - 0.5 * l1_distance))


def _health_index(
    abundances: dict[str, float],
    baseline_abundances: dict[str, float],
    microbe_metadata: pd.DataFrame,
    tcg_membership_map: dict[str, str] | None = None,
    interaction_metrics: dict[str, object] | None = None,
    baseline_interaction_rho: float = 0.0,
) -> dict[str, float]:
    """Compute composite health metrics from diversity, balance, and baseline stability."""
    diversity = _normalized_shannon(abundances)
    beneficial_fraction = _weighted_fraction(abundances, microbe_metadata, BENEFICIAL_GENERA)
    risk_fraction = _weighted_fraction(abundances, microbe_metadata, RISK_GENERA)
    stability = _stability_score(abundances, baseline_abundances)
    balance = max(0.0, min(1.0, 0.5 + beneficial_fraction - risk_fraction))
    health_index_legacy = 100.0 * (0.40 * diversity + 0.35 * balance + 0.25 * stability)
    interaction_metrics = interaction_metrics or {}
    interaction_balance_rho = float(np.clip(_coerce_float(interaction_metrics.get("interaction_balance_rho"), default=0.0), -1.0, 1.0))
    interaction_balance_shift = float(interaction_balance_rho - float(baseline_interaction_rho))
    interaction_coverage = _clip01(_coerce_float(interaction_metrics.get("interaction_coverage"), default=0.0))
    interaction_competitiveness = float(max(0.0, min(1.0, 0.5 - 0.5 * interaction_balance_rho)))
    interaction_shift_score = float(max(0.0, min(1.0, 1.0 - max(0.0, interaction_balance_shift))))
    interaction_component = float(
        interaction_coverage * (0.65 * interaction_competitiveness + 0.35 * interaction_shift_score)
        + (1.0 - interaction_coverage) * 0.5
    )
    health_index = 100.0 * (
        0.30 * diversity
        + 0.20 * balance
        + 0.20 * stability
        + 0.30 * interaction_component
    )
    tcg_membership_map = tcg_membership_map or {}
    tcg_guild_1_fraction = float(
        sum(float(abundances.get(nt_code, 0.0)) for nt_code, membership in tcg_membership_map.items() if membership == "guild_1")
    )
    tcg_guild_2_fraction = float(
        sum(float(abundances.get(nt_code, 0.0)) for nt_code, membership in tcg_membership_map.items() if membership == "guild_2")
    )
    tcg_mapped_fraction = float(max(0.0, tcg_guild_1_fraction + tcg_guild_2_fraction))
    if tcg_mapped_fraction > 0:
        tcg_guild_1_share = float(tcg_guild_1_fraction / tcg_mapped_fraction)
        tcg_guild_2_share = float(tcg_guild_2_fraction / tcg_mapped_fraction)
        tcg_balance = float(max(0.0, min(1.0, 0.5 + tcg_guild_1_share - tcg_guild_2_share)))
        tcg_balance_coverage_adjusted = float(tcg_mapped_fraction * tcg_balance + (1.0 - tcg_mapped_fraction) * 0.5)
        tcg_health_index = float(
            100.0
            * (
                0.28 * diversity
                + 0.22 * tcg_balance_coverage_adjusted
                + 0.20 * stability
                + 0.30 * interaction_component
            )
        )
    else:
        tcg_guild_1_share = float("nan")
        tcg_guild_2_share = float("nan")
        tcg_balance = float("nan")
        tcg_balance_coverage_adjusted = float("nan")
        tcg_health_index = float("nan")
    return {
        "health_index": float(max(0.0, min(100.0, health_index))),
        "health_index_legacy": float(max(0.0, min(100.0, health_index_legacy))),
        "diversity": float(diversity),
        "beneficial_fraction": float(beneficial_fraction),
        "risk_fraction": float(risk_fraction),
        "stability": float(stability),
        "interaction_component": float(interaction_component),
        "interaction_coverage": float(interaction_coverage),
        "interaction_balance_rho": float(interaction_balance_rho),
        "interaction_balance_shift": float(interaction_balance_shift),
        "interaction_competitiveness": float(interaction_competitiveness),
        "interaction_shift_score": float(interaction_shift_score),
        "positive_interaction_strength": float(
            _coerce_float(interaction_metrics.get("positive_interaction_strength"), default=0.0)
        ),
        "negative_interaction_strength": float(
            _coerce_float(interaction_metrics.get("negative_interaction_strength"), default=0.0)
        ),
        "interaction_positive_share": float(
            _coerce_float(interaction_metrics.get("interaction_positive_share"), default=0.0)
        ),
        "interaction_negative_share": float(
            _coerce_float(interaction_metrics.get("interaction_negative_share"), default=0.0)
        ),
        "tcg_guild_1_fraction": float(tcg_guild_1_fraction),
        "tcg_guild_2_fraction": float(tcg_guild_2_fraction),
        "tcg_mapped_fraction": float(tcg_mapped_fraction),
        "tcg_unmapped_fraction": float(max(0.0, 1.0 - tcg_mapped_fraction)),
        "tcg_guild_1_share_within_mapped": float(tcg_guild_1_share),
        "tcg_guild_2_share_within_mapped": float(tcg_guild_2_share),
        "tcg_balance": float(tcg_balance),
        "tcg_balance_coverage_adjusted": float(tcg_balance_coverage_adjusted),
        "tcg_health_index": float(tcg_health_index),
    }


def _sigmoid_score(value: float, center: float = 45.0, scale: float = 8.0) -> float:
    """Map an unbounded score balance into a 0-100 development score."""
    safe_scale = max(float(scale), 1e-6)
    return float(100.0 / (1.0 + math.exp(-(float(value) - float(center)) / safe_scale)))


def _top_microbe_changes(
    baseline_abundances: dict[str, float],
    final_abundances: dict[str, float],
    microbe_metadata: pd.DataFrame,
    top_n: int = 10,
) -> pd.DataFrame:
    """Return the microbes with the largest absolute abundance changes."""
    label_map = microbe_metadata.set_index("nt_code")["species_label"].to_dict()
    rows = []
    for nt_code in sorted(set(baseline_abundances) | set(final_abundances)):
        start_value = float(baseline_abundances.get(nt_code, 0.0))
        end_value = float(final_abundances.get(nt_code, 0.0))
        rows.append(
            {
                "nt_code": nt_code,
                "species_label": label_map.get(nt_code, nt_code),
                "initial_abundance": start_value,
                "final_abundance": end_value,
                "delta_abundance": end_value - start_value,
                "fold_change": np.nan if start_value <= 0 else end_value / start_value,
            }
        )
    frame = pd.DataFrame(rows)
    return frame.reindex(frame["delta_abundance"].abs().sort_values(ascending=False).index).head(top_n).reset_index(drop=True)


def _mean_applicability(
    abundances: dict[str, float],
    applicability_map: dict[str, float],
) -> float:
    """Compute the abundance-weighted mean applicability score across microbes."""
    score = 0.0
    for nt_code, abundance in abundances.items():
        score += float(abundance) * float(applicability_map.get(nt_code, 0.0))
    return float(score)


def _normalize_disease_target_profile(disease_target_profile: dict[str, float] | None) -> dict[str, float]:
    """Normalize a disease-target profile into signed nt_code weights with unit L1 norm."""
    if not disease_target_profile:
        return {}
    cleaned: dict[str, float] = {}
    for key, value in disease_target_profile.items():
        nt_code = str(key).strip()
        if not nt_code:
            continue
        weight = _coerce_float(value, default=0.0)
        if abs(weight) <= 1e-12:
            continue
        cleaned[nt_code] = float(weight)
    l1 = float(sum(abs(value) for value in cleaned.values()))
    if l1 <= 0:
        return {}
    return {key: float(value / l1) for key, value in cleaned.items()}


def _disease_target_alignment(
    abundances: dict[str, float],
    baseline_abundances: dict[str, float],
    disease_target_profile: dict[str, float],
    alignment_scale: float = 0.20,
) -> dict[str, float]:
    """Measure how much current abundance shifts align with disease-target directions."""
    if not disease_target_profile:
        return {
            "disease_target_alignment_raw": 0.0,
            "disease_target_alignment_score": 50.0,
            "disease_target_coverage": 0.0,
        }
    total_weight = float(sum(abs(weight) for weight in disease_target_profile.values()))
    if total_weight <= 0:
        return {
            "disease_target_alignment_raw": 0.0,
            "disease_target_alignment_score": 50.0,
            "disease_target_coverage": 0.0,
        }
    matched_weight = 0.0
    alignment_raw = 0.0
    for nt_code, weight in disease_target_profile.items():
        baseline = float(baseline_abundances.get(nt_code, 0.0))
        current = float(abundances.get(nt_code, 0.0))
        if baseline <= 0 and current <= 0:
            continue
        matched_weight += abs(weight)
        alignment_raw += float(weight) * (current - baseline)
    coverage = float(np.clip(matched_weight / total_weight, 0.0, 1.0))
    alignment_01 = _clip01(0.5 + alignment_raw / max(float(alignment_scale), 1e-6))
    effective_alignment = float(coverage * alignment_01 + (1.0 - coverage) * 0.5)
    return {
        "disease_target_alignment_raw": float(alignment_raw),
        "disease_target_alignment_score": float(100.0 * effective_alignment),
        "disease_target_coverage": float(coverage),
    }


def _build_applicability_map(drug_frame: pd.DataFrame) -> dict[str, float]:
    """Build per-microbe applicability scores from Step 2 inference support signals."""
    applicability: dict[str, float] = {}
    for _, row in drug_frame.iterrows():
        nt_code = str(row["nt_code"])
        hard_flag = 1.0 if bool(row.get("applicability_flag", False)) else 0.0
        drug_similarity = float(np.clip(_coerce_float(row.get("drug_max_fingerprint_jaccard"), default=0.0), 0.0, 1.0))
        scaffold_seen = 1.0 if bool(row.get("scaffold_seen_in_training", False)) else 0.0
        genus_seen = 1.0 if bool(row.get("microbe_genus_seen_in_training", False)) else 0.0
        phylum_seen = 1.0 if bool(row.get("microbe_phylum_seen_in_training", False)) else 0.0
        microbe_support = max(genus_seen, 0.5 * phylum_seen)
        soft_score = 0.45 * drug_similarity + 0.25 * scaffold_seen + 0.30 * microbe_support
        applicability[nt_code] = float(max(hard_flag, min(1.0, soft_score)))
    return applicability


def _load_cross_feeding_reference(
    cross_feeding_reference_path: str | Path | None,
    compound_context_keys: set[str],
) -> list[dict[str, object]]:
    """Load compound-aware curated cross-feeding edges relevant to the current drug context."""
    if cross_feeding_reference_path is None:
        return []
    path = Path(cross_feeding_reference_path)
    if not path.exists():
        return []

    raw = pd.read_csv(path, low_memory=False)
    if raw.empty:
        return []

    matched_rows: list[dict[str, object]] = []
    for _, row in raw.iterrows():
        row_keys: set[str] = set()
        for column in [
            "compound_name_raw",
            "compound_name_normalized",
            "compound_aliases",
            "compound_family",
            "match_keywords",
        ]:
            key = _canonicalize_key(row.get(column))
            if key:
                row_keys.add(key)
            for item in _parse_flexible_multi_value(row.get(column)):
                item_key = _canonicalize_key(item)
                if item_key:
                    row_keys.add(item_key)
        if compound_context_keys and not (row_keys & compound_context_keys):
            continue
        matched_rows.append(
            {
                "producer_key": _canonicalize_key(row.get("producer_microbe_label")),
                "consumer_key": _canonicalize_key(row.get("consumer_microbe_label")),
                "evidence_level": row.get("evidence_level"),
                "evidence_type": row.get("evidence_type"),
            }
        )
    return matched_rows


def _load_cross_feeding_enzyme_weights(
    enzyme_function_catalog_path: str | Path | None,
) -> dict[str, float]:
    """Load enzyme-weight priors for functions likely to create cross-feeding opportunities."""
    if enzyme_function_catalog_path is None:
        return {}
    path = Path(enzyme_function_catalog_path)
    if not path.exists():
        return {}

    raw = pd.read_csv(path, low_memory=False)
    if raw.empty:
        return {}

    support_roles = {
        _canonicalize_key("cross_feeding_support"),
        _canonicalize_key("metabolism_supported_promote"),
        _canonicalize_key("nutrient_release_support"),
        _canonicalize_key("community_context_modulation"),
    }
    weights: dict[str, float] = {}
    for _, row in raw.iterrows():
        enzyme_id = _canonicalize_text(row.get("enzyme_id"))
        if not enzyme_id:
            continue
        role_key = _canonicalize_key(row.get("step1_feedback_role"))
        if role_key not in support_roles:
            continue
        weights[enzyme_id] = _clip01(
            max(
                _coerce_float(row.get("step1_promote_weight"), default=0.0),
                0.7 * _coerce_float(row.get("metabolism_weight"), default=0.0),
            )
        )
    return weights


def _build_interaction_model(
    drug_frame: pd.DataFrame,
    microbe_metadata: pd.DataFrame,
    cross_feeding_reference_path: str | Path | None,
    enzyme_function_catalog_path: str | Path | None,
) -> dict[str, object]:
    """Build a lightweight positive-vs-negative interaction model for one drug scenario."""
    compound_context_keys = _build_drug_context_keys(drug_frame)
    curated_edges = _load_cross_feeding_reference(
        cross_feeding_reference_path=cross_feeding_reference_path,
        compound_context_keys=compound_context_keys,
    )
    cross_feeding_enzyme_weights = _load_cross_feeding_enzyme_weights(
        enzyme_function_catalog_path=enzyme_function_catalog_path
    )

    species_lookup: dict[str, str] = {}
    for _, row in microbe_metadata.iterrows():
        nt_code = str(row["nt_code"])
        for candidate in [
            row.get("microbe_label"),
            row.get("species_label"),
            row.get("species_name"),
        ]:
            key = _canonicalize_key(candidate)
            if key:
                species_lookup[key] = nt_code

    curated_positive_edges: dict[tuple[str, str], float] = {}
    for item in curated_edges:
        producer_nt = species_lookup.get(str(item["producer_key"]))
        consumer_nt = species_lookup.get(str(item["consumer_key"]))
        if not producer_nt or not consumer_nt or producer_nt == consumer_nt:
            continue
        curated_weight = _evidence_level_weight(item.get("evidence_level"))
        if _canonicalize_key(item.get("evidence_type")) == _canonicalize_key("cross_feeding_supported_promote"):
            curated_weight += 0.15
        edge_key = (producer_nt, consumer_nt)
        curated_positive_edges[edge_key] = max(curated_positive_edges.get(edge_key, 0.0), curated_weight)

    metadata_map = microbe_metadata.set_index("nt_code").to_dict(orient="index")
    microbe_profiles: dict[str, dict[str, object]] = {}
    for _, row in drug_frame.iterrows():
        nt_code = str(row["nt_code"])
        effect_score = _coerce_float(
            row.get("step1_predicted_effect_score", row.get("predicted_effect_score")),
            default=0.0,
        )
        positive_effect_score = _clip01(max(effect_score, 0.0))
        metabolized_probability = _clip01(_coerce_float(row.get("predicted_metabolized_probability"), default=0.0))
        depletion_strength = _clip01(
            max(0.0, -_coerce_float(row.get("predicted_parent_depletion_fraction"), default=0.0))
        )
        mechanism_support = _clip01(_coerce_float(row.get("predicted_mechanism_support_score"), default=0.0))
        promote_support = _clip01(_coerce_float(row.get("predicted_enzyme_step1_promote_support_score"), default=0.0))
        enzyme_presence = _clip01(_coerce_float(row.get("predicted_enzyme_presence_score"), default=0.0))
        candidate_product_score = _clip01(_coerce_float(row.get("predicted_candidate_product_count"), default=0.0) / 3.0)
        enzyme_ids = {_canonicalize_text(item) for item in _parse_multi_value(row.get("predicted_enzyme_ids"))}
        enzyme_ids = {item for item in enzyme_ids if item}
        reaction_classes = {_canonicalize_text(item) for item in _parse_multi_value(row.get("predicted_reaction_class"))}
        reaction_classes = {item for item in reaction_classes if item}
        cross_feeding_enzyme_score = _clip01(
            sum(cross_feeding_enzyme_weights.get(enzyme_id, 0.0) for enzyme_id in enzyme_ids)
        )
        producer_score = _clip01(
            0.25 * promote_support
            + 0.15 * enzyme_presence
            + 0.15 * candidate_product_score
            + 0.20 * cross_feeding_enzyme_score
            + 0.10 * mechanism_support
            + 0.15 * metabolized_probability
            + 0.10 * positive_effect_score
        )
        consumer_score = _clip01(
            0.40 * promote_support
            + 0.25 * positive_effect_score
            + 0.15 * cross_feeding_enzyme_score
            + 0.10 * candidate_product_score
            + 0.10 * metabolized_probability
        )
        competition_score = _clip01(
            0.50 * metabolized_probability
            + 0.35 * depletion_strength
            + 0.15 * max(mechanism_support, enzyme_presence)
        )
        microbe_profiles[nt_code] = {
            "producer_score": producer_score,
            "consumer_score": consumer_score,
            "competition_score": competition_score,
            "enzyme_ids": enzyme_ids,
            "reaction_classes": reaction_classes,
            "genus": _canonicalize_text(metadata_map.get(nt_code, {}).get("genus")),
            "family": _canonicalize_text(metadata_map.get(nt_code, {}).get("family")),
            "phylum": _canonicalize_text(metadata_map.get(nt_code, {}).get("phylum")),
        }

    positive_edges: dict[tuple[str, str], float] = {}
    negative_pairs: dict[tuple[str, str], float] = {}
    nt_codes = sorted(microbe_profiles)
    for source_nt in nt_codes:
        source_profile = microbe_profiles[source_nt]
        for target_nt in nt_codes:
            if source_nt == target_nt:
                continue
            target_profile = microbe_profiles[target_nt]
            generic_positive = math.sqrt(
                float(source_profile["producer_score"]) * float(target_profile["consumer_score"])
            )
            if source_profile["genus"] and source_profile["genus"] == target_profile["genus"]:
                generic_positive *= 0.8
            elif source_profile["family"] and source_profile["family"] == target_profile["family"]:
                generic_positive *= 0.9
            curated_bonus = curated_positive_edges.get((source_nt, target_nt), 0.0)
            positive_weight = min(1.5, 0.55 * generic_positive + curated_bonus)
            if positive_weight >= 0.05:
                positive_edges[(source_nt, target_nt)] = float(positive_weight)

        for target_nt in nt_codes:
            if source_nt >= target_nt:
                continue
            target_profile = microbe_profiles[target_nt]
            taxonomic_overlap = 0.0
            if source_profile["genus"] and source_profile["genus"] == target_profile["genus"]:
                taxonomic_overlap = 0.25
            elif source_profile["family"] and source_profile["family"] == target_profile["family"]:
                taxonomic_overlap = 0.12
            elif source_profile["phylum"] and source_profile["phylum"] == target_profile["phylum"]:
                taxonomic_overlap = 0.05
            enzyme_overlap = _jaccard_similarity(
                set(source_profile["enzyme_ids"]),
                set(target_profile["enzyme_ids"]),
            )
            reaction_overlap = _jaccard_similarity(
                set(source_profile["reaction_classes"]),
                set(target_profile["reaction_classes"]),
            )
            competition_weight = min(
                1.5,
                0.55 * math.sqrt(float(source_profile["competition_score"]) * float(target_profile["competition_score"]))
                + 0.20 * enzyme_overlap
                + 0.10 * reaction_overlap
                + taxonomic_overlap,
            )
            if competition_weight >= 0.08:
                negative_pairs[(source_nt, target_nt)] = float(competition_weight)

    involved_nt_codes = {
        nt_code
        for pair in list(positive_edges) + list(negative_pairs)
        for nt_code in pair
    }
    return {
        "compound_context_keys": sorted(compound_context_keys),
        "curated_reference_edge_count": int(len(curated_positive_edges)),
        "positive_edge_count": int(len(positive_edges)),
        "negative_edge_count": int(len(negative_pairs)),
        "involved_nt_codes": involved_nt_codes,
        "positive_edges": positive_edges,
        "negative_pairs": negative_pairs,
    }


def _summarize_interaction_state(
    abundances: dict[str, float],
    interaction_model: dict[str, object],
) -> dict[str, object]:
    """Summarize the current interaction regime as positive-vs-negative network balance."""
    positive_edges = interaction_model.get("positive_edges", {})
    negative_pairs = interaction_model.get("negative_pairs", {})
    positive_incoming_pressure = {nt_code: 0.0 for nt_code in abundances}
    negative_incoming_pressure = {nt_code: 0.0 for nt_code in abundances}

    positive_strength = 0.0
    for (source_nt, target_nt), weight in positive_edges.items():
        source_abundance = float(abundances.get(source_nt, 0.0))
        target_abundance = float(abundances.get(target_nt, 0.0))
        if source_abundance <= 0 or target_abundance <= 0:
            continue
        positive_strength += source_abundance * target_abundance * float(weight)
        positive_incoming_pressure[target_nt] = positive_incoming_pressure.get(target_nt, 0.0) + source_abundance * float(weight)

    negative_strength = 0.0
    for (left_nt, right_nt), weight in negative_pairs.items():
        left_abundance = float(abundances.get(left_nt, 0.0))
        right_abundance = float(abundances.get(right_nt, 0.0))
        if left_abundance <= 0 or right_abundance <= 0:
            continue
        negative_strength += 2.0 * left_abundance * right_abundance * float(weight)
        negative_incoming_pressure[left_nt] = negative_incoming_pressure.get(left_nt, 0.0) + right_abundance * float(weight)
        negative_incoming_pressure[right_nt] = negative_incoming_pressure.get(right_nt, 0.0) + left_abundance * float(weight)

    total_strength = positive_strength + negative_strength
    interaction_balance_rho = 0.0 if total_strength <= 0 else (positive_strength - negative_strength) / total_strength
    involved_nt_codes = set(interaction_model.get("involved_nt_codes", set()))
    interaction_coverage = float(sum(float(abundances.get(nt_code, 0.0)) for nt_code in involved_nt_codes))
    return {
        "positive_interaction_strength": float(positive_strength),
        "negative_interaction_strength": float(negative_strength),
        "total_interaction_strength": float(total_strength),
        "interaction_balance_rho": float(np.clip(interaction_balance_rho, -1.0, 1.0)),
        "interaction_positive_share": float(0.0 if total_strength <= 0 else positive_strength / total_strength),
        "interaction_negative_share": float(0.0 if total_strength <= 0 else negative_strength / total_strength),
        "interaction_coverage": float(np.clip(interaction_coverage, 0.0, 1.0)),
        "positive_incoming_pressure": positive_incoming_pressure,
        "negative_incoming_pressure": negative_incoming_pressure,
    }


def run_step3_simulation(
    integrated_predictions_path: str | Path,
    output_dir: str | Path,
    drug_query: str,
    scenario_name: str = "healthy_reference",
    community_table_path: str | Path | None = None,
    tcg_proxy_mapping_path: str | Path | None = DEFAULT_TCG_PROXY_MAPPING_PATH,
    cross_feeding_reference_path: str | Path | None = DEFAULT_CROSS_FEEDING_EDGES_PATH,
    enzyme_function_catalog_path: str | Path | None = DEFAULT_ENZYME_FUNCTION_CATALOG_PATH,
    n_steps: int = 14,
    initial_dose: float = 1.0,
    repeat_dose: float = 1.0,
    dosing_interval: int = 1,
    drug_clearance_rate: float = 0.12,
    product_clearance_rate: float = 0.18,
    metabolism_scale: float = 0.85,
    effect_scale: float = 0.55,
    ecology_strength: float = 0.20,
    interaction_scale: float = 0.18,
    cooperative_metabolism_boost: float = 0.35,
    abundance_floor: float = 1e-6,
    disease_target_profile: dict[str, float] | None = None,
    experimental_multi_product_enabled: bool = False,
    experimental_branching_scale: float = 0.35,
    experimental_secondary_metabolism_rate: float = 0.10,
) -> dict[str, object]:
    """Run the Step 3 community simulation and export trajectories plus summaries.

    Args:
        integrated_predictions_path: Integrated Step 1 and Step 2 prediction table.
        output_dir: Directory where trajectory outputs and summary JSON are written.
        drug_query: Prestwick ID or drug name used to select one drug for simulation.
        scenario_name: Builtin scenario used when no custom community table is supplied.
        community_table_path: Optional custom community abundance table.
        tcg_proxy_mapping_path: Optional guild mapping file used for extra health metrics.
        cross_feeding_reference_path: Optional curated positive-interaction reference table.
        enzyme_function_catalog_path: Optional enzyme catalog used to infer interaction roles.
        n_steps: Number of discrete simulation time steps.
        initial_dose: Initial parent-drug dose at time zero.
        repeat_dose: Dose amount added at each dosing interval after the first step.
        dosing_interval: Number of time steps between repeated doses.
        drug_clearance_rate: Fraction of parent drug cleared per step.
        product_clearance_rate: Fraction of metabolite pool cleared per step.
        metabolism_scale: Global scale factor on microbial metabolism pressure.
        effect_scale: Global scale factor on Step 1 drug-pressure effects.
        ecology_strength: Strength of pull back toward the baseline community.
        interaction_scale: Strength of the interaction-aware feedback term.
        cooperative_metabolism_boost: Extra metabolism gain in cross-feeding-dominant states.
        abundance_floor: Minimum abundance floor applied during renormalization.
        disease_target_profile: Optional signed nt_code profile (+ promote target, - inhibit target)
            used to compute disease-target alignment rewards.
        experimental_multi_product_enabled: When True, emit a parallel experimental score that
            inflates metabolite burden for drugs with multiple predicted products.
        experimental_branching_scale: Extra direct product burden scale under the experimental mode.
        experimental_secondary_metabolism_rate: Fraction of the experimental product pool allowed
            to recursively generate downstream burden per step.

    Returns:
        A summary dictionary describing the simulation outputs and final metrics.
    """
    integrated_predictions_path = Path(integrated_predictions_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    frame = pd.read_csv(integrated_predictions_path, low_memory=False)
    drug_frame = _resolve_drug_subset(frame, drug_query=drug_query)
    if drug_frame.empty:
        raise RuntimeError(f"No rows found for drug_query={drug_query!r}")

    microbe_metadata = _build_microbe_metadata(drug_frame)
    tcg_membership_map, tcg_mapping_metadata = _load_tcg_proxy_mapping(
        tcg_proxy_mapping_path=tcg_proxy_mapping_path,
        microbe_metadata=microbe_metadata,
    )
    if community_table_path is None:
        template_abundances, scenario_metadata = _build_builtin_community(microbe_metadata, scenario_name=scenario_name)
    else:
        template_abundances, scenario_metadata = _load_custom_community(community_table_path, microbe_metadata=microbe_metadata)

    effect_map = (
        drug_frame.set_index("nt_code")["step1_predicted_effect_score"]
        if "step1_predicted_effect_score" in drug_frame.columns
        else drug_frame.set_index("nt_code")["predicted_effect_score"]
    ).to_dict()
    inhibit_prob_map = (
        drug_frame.set_index("nt_code")["step1_predicted_inhibit_probability"]
        if "step1_predicted_inhibit_probability" in drug_frame.columns
        else drug_frame.set_index("nt_code")["predicted_inhibit_probability"]
    ).to_dict()
    metabolism_prob_map = drug_frame.set_index("nt_code")["predicted_metabolized_probability"].to_dict()
    depletion_map = drug_frame.set_index("nt_code")["predicted_parent_depletion_fraction"].to_dict()
    experimental_product_count_map = _optional_column_map(
        drug_frame,
        ["experimental_biotransform_product_count", "predicted_candidate_product_count"],
    )
    experimental_fraction_in_gut_map = _optional_column_map(
        drug_frame,
        ["experimental_biotransform_fraction_in_gut"],
    )
    applicability_map = _build_applicability_map(drug_frame)
    interaction_model = _build_interaction_model(
        drug_frame=drug_frame,
        microbe_metadata=microbe_metadata,
        cross_feeding_reference_path=cross_feeding_reference_path,
        enzyme_function_catalog_path=enzyme_function_catalog_path,
    )
    normalized_disease_target_profile = _normalize_disease_target_profile(disease_target_profile)
    disease_target_enabled = bool(normalized_disease_target_profile)
    prestwick_id = _first_nonempty_text(drug_frame, "prestwick_id", default=drug_query)
    chemical_name = _first_nonempty_text(drug_frame, "chemical_name", default=prestwick_id or drug_query)

    current_abundances = template_abundances.copy()
    baseline_abundances = template_abundances.copy()
    cumulative_dose = float(initial_dose)
    current_parent_drug = float(initial_dose)
    current_metabolite_pool = 0.0
    current_experimental_metabolite_pool = 0.0

    history_rows: list[dict[str, object]] = []
    abundance_rows: list[dict[str, object]] = []

    baseline_interaction_metrics = _summarize_interaction_state(
        abundances=baseline_abundances,
        interaction_model=interaction_model,
    )
    baseline_interaction_rho = float(baseline_interaction_metrics["interaction_balance_rho"])
    initial_health = _health_index(
        current_abundances,
        baseline_abundances,
        microbe_metadata=microbe_metadata,
        tcg_membership_map=tcg_membership_map,
        interaction_metrics=baseline_interaction_metrics,
        baseline_interaction_rho=baseline_interaction_rho,
    )

    def _record_timepoint(timepoint: int) -> None:
        nonlocal history_rows, abundance_rows
        interaction_metrics = _summarize_interaction_state(
            abundances=current_abundances,
            interaction_model=interaction_model,
        )
        health = _health_index(
            current_abundances,
            baseline_abundances,
            microbe_metadata=microbe_metadata,
            tcg_membership_map=tcg_membership_map,
            interaction_metrics=interaction_metrics,
            baseline_interaction_rho=baseline_interaction_rho,
        )
        mean_applicability = _mean_applicability(current_abundances, applicability_map)
        parent_retention = 0.0 if cumulative_dose <= 0 else current_parent_drug / cumulative_dose
        dysbiosis_penalty = max(0.0, initial_health["health_index"] - health["health_index"])
        interaction_dysbiosis_penalty = 100.0 * max(0.0, health["interaction_balance_shift"]) * max(
            0.25, health["interaction_coverage"]
        )
        uncertainty_penalty = 100.0 * max(0.0, 1.0 - mean_applicability)
        efficacy_proxy = 100.0 * parent_retention
        stability_score = 100.0 * health["stability"]
        community_preservation_score = 0.65 * health["health_index"] + 0.35 * stability_score
        metabolite_burden_ratio = 0.0 if cumulative_dose <= 0 else current_metabolite_pool / cumulative_dose
        metabolite_burden_penalty = 100.0 * max(0.0, min(1.0, metabolite_burden_ratio))
        active_experimental_pool = (
            current_experimental_metabolite_pool if experimental_multi_product_enabled else current_metabolite_pool
        )
        experimental_metabolite_burden_ratio = 0.0 if cumulative_dose <= 0 else active_experimental_pool / cumulative_dose
        experimental_metabolite_burden_penalty = 100.0 * max(0.0, min(1.0, experimental_metabolite_burden_ratio))
        disease_target_metrics = _disease_target_alignment(
            abundances=current_abundances,
            baseline_abundances=baseline_abundances,
            disease_target_profile=normalized_disease_target_profile,
        )
        disease_target_alignment_score = float(disease_target_metrics["disease_target_alignment_score"])
        if disease_target_enabled:
            benefit_subscore = (
                0.45 * efficacy_proxy
                + 0.35 * community_preservation_score
                + 0.20 * disease_target_alignment_score
            )
        else:
            benefit_subscore = 0.55 * efficacy_proxy + 0.45 * community_preservation_score
        risk_subscore = (
            0.38 * dysbiosis_penalty
            + 0.22 * interaction_dysbiosis_penalty
            + 0.22 * uncertainty_penalty
            + 0.18 * metabolite_burden_penalty
        )
        experimental_risk_subscore = (
            0.38 * dysbiosis_penalty
            + 0.22 * interaction_dysbiosis_penalty
            + 0.22 * uncertainty_penalty
            + 0.18 * experimental_metabolite_burden_penalty
        )
        development_score_balance = benefit_subscore - risk_subscore
        experimental_development_score_balance = benefit_subscore - experimental_risk_subscore
        development_score_raw = development_score_balance
        development_score = _sigmoid_score(development_score_balance, center=45.0, scale=8.0)
        experimental_development_score_raw = experimental_development_score_balance
        experimental_development_score = _sigmoid_score(
            experimental_development_score_balance,
            center=45.0,
            scale=8.0,
        )
        development_score_legacy_raw = efficacy_proxy - 0.7 * dysbiosis_penalty - 0.3 * uncertainty_penalty
        development_score_legacy = max(0.0, min(100.0, development_score_legacy_raw))
        history_rows.append(
            {
                "timepoint": timepoint,
                "parent_drug_concentration": float(current_parent_drug),
                "aggregate_metabolite_pool": float(current_metabolite_pool),
                "experimental_aggregate_metabolite_pool": float(active_experimental_pool),
                "cumulative_dose": float(cumulative_dose),
                "parent_retention_ratio": float(parent_retention),
                "health_index": health["health_index"],
                "health_index_legacy": health["health_index_legacy"],
                "diversity": health["diversity"],
                "beneficial_fraction": health["beneficial_fraction"],
                "risk_fraction": health["risk_fraction"],
                "stability": health["stability"],
                "interaction_component": health["interaction_component"],
                "interaction_coverage": health["interaction_coverage"],
                "positive_interaction_strength": health["positive_interaction_strength"],
                "negative_interaction_strength": health["negative_interaction_strength"],
                "interaction_positive_share": health["interaction_positive_share"],
                "interaction_negative_share": health["interaction_negative_share"],
                "interaction_balance_rho": health["interaction_balance_rho"],
                "interaction_balance_shift": health["interaction_balance_shift"],
                "interaction_competitiveness": health["interaction_competitiveness"],
                "interaction_shift_score": health["interaction_shift_score"],
                "tcg_health_index": health["tcg_health_index"],
                "tcg_guild_1_fraction": health["tcg_guild_1_fraction"],
                "tcg_guild_2_fraction": health["tcg_guild_2_fraction"],
                "tcg_mapped_fraction": health["tcg_mapped_fraction"],
                "tcg_unmapped_fraction": health["tcg_unmapped_fraction"],
                "tcg_guild_1_share_within_mapped": health["tcg_guild_1_share_within_mapped"],
                "tcg_guild_2_share_within_mapped": health["tcg_guild_2_share_within_mapped"],
                "tcg_balance": health["tcg_balance"],
                "tcg_balance_coverage_adjusted": health["tcg_balance_coverage_adjusted"],
                "stability_score": float(stability_score),
                "mean_applicability": mean_applicability,
                "dysbiosis_penalty": float(dysbiosis_penalty),
                "interaction_dysbiosis_penalty": float(interaction_dysbiosis_penalty),
                "uncertainty_penalty": float(uncertainty_penalty),
                "metabolite_burden_ratio": float(metabolite_burden_ratio),
                "metabolite_burden_penalty": float(metabolite_burden_penalty),
                "experimental_metabolite_burden_ratio": float(experimental_metabolite_burden_ratio),
                "experimental_metabolite_burden_penalty": float(experimental_metabolite_burden_penalty),
                "efficacy_proxy": float(efficacy_proxy),
                "community_preservation_score": float(community_preservation_score),
                "disease_target_alignment_raw": float(disease_target_metrics["disease_target_alignment_raw"]),
                "disease_target_alignment_score": float(disease_target_alignment_score),
                "disease_target_coverage": float(disease_target_metrics["disease_target_coverage"]),
                "disease_target_reward_enabled": bool(disease_target_enabled),
                "benefit_subscore": float(benefit_subscore),
                "risk_subscore": float(risk_subscore),
                "experimental_risk_subscore": float(experimental_risk_subscore),
                "development_score_balance": float(development_score_balance),
                "experimental_development_score_balance": float(experimental_development_score_balance),
                "development_score_legacy_raw": float(development_score_legacy_raw),
                "development_score_legacy": float(development_score_legacy),
                "development_score_raw": float(development_score_raw),
                "development_score": float(development_score),
                "experimental_development_score_raw": float(experimental_development_score_raw),
                "experimental_development_score": float(experimental_development_score),
                "experimental_multi_product_enabled": bool(experimental_multi_product_enabled),
            }
        )
        for nt_code, abundance in sorted(current_abundances.items()):
            abundance_rows.append(
                {
                    "timepoint": timepoint,
                    "nt_code": nt_code,
                    "abundance": float(abundance),
                }
            )

    _record_timepoint(timepoint=0)

    for timepoint in range(1, n_steps + 1):
        dose_input = 0.0
        if dosing_interval > 0 and repeat_dose > 0 and timepoint > 1 and (timepoint - 1) % dosing_interval == 0:
            dose_input = float(repeat_dose)

        drug_exposure = current_parent_drug / max(initial_dose, 1e-8)
        interaction_state = _summarize_interaction_state(
            abundances=current_abundances,
            interaction_model=interaction_model,
        )
        raw_abundances: dict[str, float] = {}
        weighted_metabolism = 0.0
        weighted_branching_signal = 0.0
        for nt_code, current_abundance in current_abundances.items():
            effect_score = _coerce_float(effect_map.get(nt_code), default=0.0)
            ecology_pull = ecology_strength * (baseline_abundances.get(nt_code, 0.0) - current_abundance)
            positive_pressure = _coerce_float(
                interaction_state["positive_incoming_pressure"].get(nt_code),
                default=0.0,
            )
            negative_pressure = _coerce_float(
                interaction_state["negative_incoming_pressure"].get(nt_code),
                default=0.0,
            )
            interaction_delta = float(
                np.clip(
                    interaction_scale * (0.20 + 0.80 * drug_exposure) * (positive_pressure - negative_pressure),
                    -0.75,
                    0.75,
                )
            )
            # Clip the exponent to keep custom-SMILES simulations numerically stable
            # even when a regressor produces outlier effect scores.
            growth_log_delta = float(
                np.clip(
                    effect_scale * effect_score * drug_exposure + ecology_pull + interaction_delta,
                    -12.0,
                    12.0,
                )
            )
            growth_multiplier = math.exp(growth_log_delta)
            next_value = max(abundance_floor, current_abundance * growth_multiplier)
            raw_abundances[nt_code] = next_value

            metabolized_probability = _coerce_float(metabolism_prob_map.get(nt_code), default=0.0)
            depletion_strength = float(
                np.clip(
                    -_coerce_float(depletion_map.get(nt_code), default=0.0),
                    0.0,
                    1.0,
                )
            )
            metabolism_pressure = (
                current_abundance * metabolized_probability * max(depletion_strength, 0.05 * metabolized_probability)
            )
            weighted_metabolism += metabolism_pressure

            product_count = max(_coerce_float(experimental_product_count_map.get(nt_code), default=0.0), 0.0)
            extra_product_score = _clip01(max(product_count - 1.0, 0.0) / 4.0)
            fraction_in_gut = _clip01(_coerce_float(experimental_fraction_in_gut_map.get(nt_code), default=0.0))
            branching_signal = extra_product_score * (0.55 + 0.45 * fraction_in_gut)
            weighted_branching_signal += metabolism_pressure * branching_signal

        current_abundances = _normalize_abundances(raw_abundances, minimum=abundance_floor)
        cooperative_metabolism_multiplier = 1.0 + cooperative_metabolism_boost * max(
            0.0, _coerce_float(interaction_state.get("interaction_balance_rho"), default=0.0)
        )
        parent_consumed = min(
            current_parent_drug,
            current_parent_drug * metabolism_scale * cooperative_metabolism_multiplier * weighted_metabolism,
        )
        parent_cleared = current_parent_drug * drug_clearance_rate
        current_parent_drug = max(0.0, current_parent_drug - parent_consumed - parent_cleared + dose_input)
        current_metabolite_pool = max(
            0.0,
            current_metabolite_pool * (1.0 - product_clearance_rate) + parent_consumed,
        )
        branching_ratio = 0.0 if weighted_metabolism <= 1e-8 else _clip01(weighted_branching_signal / weighted_metabolism)
        if experimental_multi_product_enabled:
            experimental_direct_generation = parent_consumed * (1.0 + experimental_branching_scale * branching_ratio)
            experimental_secondary_generation = (
                current_experimental_metabolite_pool
                * max(float(experimental_secondary_metabolism_rate), 0.0)
                * branching_ratio
            )
            current_experimental_metabolite_pool = max(
                0.0,
                current_experimental_metabolite_pool * (1.0 - product_clearance_rate)
                + experimental_direct_generation
                + experimental_secondary_generation,
            )
        else:
            current_experimental_metabolite_pool = current_metabolite_pool
        cumulative_dose += dose_input
        _record_timepoint(timepoint=timepoint)

    trajectory_metrics = pd.DataFrame(history_rows)
    trajectory_abundances = pd.DataFrame(abundance_rows)
    trajectory_wide = trajectory_abundances.pivot(index="timepoint", columns="nt_code", values="abundance").reset_index()
    top_changes = _top_microbe_changes(
        baseline_abundances=baseline_abundances,
        final_abundances=current_abundances,
        microbe_metadata=microbe_metadata,
    )
    final_row = trajectory_metrics.iloc[-1]
    summary = {
        "integrated_predictions_path": str(integrated_predictions_path),
        "output_dir": str(output_dir),
        "drug_query": drug_query,
        "prestwick_id": prestwick_id,
        "chemical_name": chemical_name,
        "scenario_name": scenario_metadata["scenario_name"],
        "scenario_description": scenario_metadata["scenario_description"],
        "community_source": scenario_metadata["source"],
        "n_steps": int(n_steps),
        "initial_dose": float(initial_dose),
        "repeat_dose": float(repeat_dose),
        "dosing_interval": int(dosing_interval),
        "drug_clearance_rate": float(drug_clearance_rate),
        "product_clearance_rate": float(product_clearance_rate),
        "metabolism_scale": float(metabolism_scale),
        "effect_scale": float(effect_scale),
        "ecology_strength": float(ecology_strength),
        "interaction_scale": float(interaction_scale),
        "cooperative_metabolism_boost": float(cooperative_metabolism_boost),
        "experimental_multi_product_enabled": bool(experimental_multi_product_enabled),
        "experimental_branching_scale": float(experimental_branching_scale),
        "experimental_secondary_metabolism_rate": float(experimental_secondary_metabolism_rate),
        "experimental_product_annotation_pairs": int(
            pd.Series(list(experimental_product_count_map.values())).map(lambda value: _coerce_float(value, default=0.0)).gt(0).sum()
        ),
        "disease_target_profile_size": int(len(normalized_disease_target_profile)),
        "disease_target_reward_enabled": bool(disease_target_enabled),
        "cross_feeding_reference_path": None if cross_feeding_reference_path is None else str(cross_feeding_reference_path),
        "enzyme_function_catalog_path": None if enzyme_function_catalog_path is None else str(enzyme_function_catalog_path),
        "interaction_reference_edge_count": int(interaction_model["curated_reference_edge_count"]),
        "interaction_positive_edge_count": int(interaction_model["positive_edge_count"]),
        "interaction_negative_edge_count": int(interaction_model["negative_edge_count"]),
        "health_signature_mode": "tcg_proxy_secondary",
        "health_signature_source_path": tcg_mapping_metadata["tcg_proxy_mapping_path"],
        "health_signature_source_name": tcg_mapping_metadata["tcg_source_name"],
        "health_signature_source_url": tcg_mapping_metadata["tcg_source_url"],
        "health_signature_ready_microbe_count": int(tcg_mapping_metadata["tcg_ready_microbe_count"]),
        "health_signature_mapped_microbe_count": int(tcg_mapping_metadata["tcg_mapped_microbe_count"]),
        "health_signature_panel_microbe_count": int(tcg_mapping_metadata["tcg_panel_microbe_count"]),
        "initial_health_index": float(trajectory_metrics.iloc[0]["health_index"]),
        "initial_health_index_legacy": float(trajectory_metrics.iloc[0]["health_index_legacy"]),
        "final_health_index": float(final_row["health_index"]),
        "final_health_index_legacy": float(final_row["health_index_legacy"]),
        "initial_tcg_health_index": float(trajectory_metrics.iloc[0]["tcg_health_index"]),
        "final_tcg_health_index": float(final_row["tcg_health_index"]),
        "initial_interaction_balance_rho": float(trajectory_metrics.iloc[0]["interaction_balance_rho"]),
        "final_interaction_balance_rho": float(final_row["interaction_balance_rho"]),
        "final_interaction_balance_shift": float(final_row["interaction_balance_shift"]),
        "final_interaction_component": float(final_row["interaction_component"]),
        "final_positive_interaction_strength": float(final_row["positive_interaction_strength"]),
        "final_negative_interaction_strength": float(final_row["negative_interaction_strength"]),
        "interaction_dysbiosis_penalty_final": float(final_row["interaction_dysbiosis_penalty"]),
        "final_tcg_guild_1_fraction": float(final_row["tcg_guild_1_fraction"]),
        "final_tcg_guild_2_fraction": float(final_row["tcg_guild_2_fraction"]),
        "final_tcg_mapped_fraction": float(final_row["tcg_mapped_fraction"]),
        "final_tcg_balance": float(final_row["tcg_balance"]),
        "initial_parent_drug_concentration": float(trajectory_metrics.iloc[0]["parent_drug_concentration"]),
        "final_parent_drug_concentration": float(final_row["parent_drug_concentration"]),
        "final_parent_retention_ratio": float(final_row["parent_retention_ratio"]),
        "final_aggregate_metabolite_pool": float(final_row["aggregate_metabolite_pool"]),
        "final_experimental_aggregate_metabolite_pool": float(final_row["experimental_aggregate_metabolite_pool"]),
        "mean_applicability_final": float(final_row["mean_applicability"]),
        "final_stability": float(final_row["stability"]),
        "efficacy_proxy_final": float(final_row["efficacy_proxy"]),
        "community_preservation_final": float(final_row["community_preservation_score"]),
        "disease_target_alignment_raw_final": float(final_row["disease_target_alignment_raw"]),
        "disease_target_alignment_score_final": float(final_row["disease_target_alignment_score"]),
        "disease_target_coverage_final": float(final_row["disease_target_coverage"]),
        "benefit_subscore_final": float(final_row["benefit_subscore"]),
        "risk_subscore_final": float(final_row["risk_subscore"]),
        "experimental_risk_subscore_final": float(final_row["experimental_risk_subscore"]),
        "dysbiosis_penalty_final": float(final_row["dysbiosis_penalty"]),
        "uncertainty_penalty_final": float(final_row["uncertainty_penalty"]),
        "metabolite_burden_penalty_final": float(final_row["metabolite_burden_penalty"]),
        "experimental_metabolite_burden_penalty_final": float(final_row["experimental_metabolite_burden_penalty"]),
        "development_score_balance": float(final_row["development_score_balance"]),
        "experimental_development_score_balance": float(final_row["experimental_development_score_balance"]),
        "development_score_legacy": float(final_row["development_score_legacy"]),
        "development_score_legacy_raw": float(final_row["development_score_legacy_raw"]),
        "development_score": float(final_row["development_score"]),
        "experimental_development_score": float(final_row["experimental_development_score"]),
        "experimental_development_score_raw": float(final_row["experimental_development_score_raw"]),
        "development_score_raw": float(final_row["development_score_raw"]),
        "trajectory_metrics_path": str(output_dir / "trajectory_metrics.csv"),
        "trajectory_abundances_path": str(output_dir / "trajectory_abundances.csv"),
        "trajectory_abundances_wide_path": str(output_dir / "trajectory_abundances_wide.csv"),
        "top_microbe_changes_path": str(output_dir / "top_microbe_changes.csv"),
        "top_microbe_changes": top_changes.to_dict(orient="records"),
    }

    trajectory_metrics.to_csv(output_dir / "trajectory_metrics.csv", index=False)
    trajectory_abundances.to_csv(output_dir / "trajectory_abundances.csv", index=False)
    trajectory_wide.to_csv(output_dir / "trajectory_abundances_wide.csv", index=False)
    top_changes.to_csv(output_dir / "top_microbe_changes.csv", index=False)
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    return summary
