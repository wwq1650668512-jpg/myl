from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_TCG_PROXY_MAPPING_PATH = ROOT / "data/processed/health_signature/microbe_tcg_proxy_mapping.csv"


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


def _canonicalize_text(value: object) -> str:
    """Normalize free text by stripping and collapsing repeated whitespace."""
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value).strip())


def _canonicalize_key(value: object) -> str:
    """Build a lowercase alphanumeric key used for fuzzy matching."""
    return re.sub(r"[^a-z0-9]+", "", _canonicalize_text(value).lower())


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
    candidates = frame.loc[:, ["prestwick_id", "chemical_name"]].drop_duplicates().reset_index(drop=True)
    query_key = _canonicalize_key(drug_query)
    exact_mask = (
        candidates["prestwick_id"].map(_canonicalize_key).eq(query_key)
        | candidates["chemical_name"].map(_canonicalize_key).eq(query_key)
    )
    if exact_mask.any():
        selected = candidates.loc[exact_mask].iloc[0]
        return frame[frame["prestwick_id"] == selected["prestwick_id"]].copy()

    contains_mask = candidates["chemical_name"].astype(str).str.contains(drug_query, case=False, na=False)
    if contains_mask.sum() == 1:
        selected = candidates.loc[contains_mask].iloc[0]
        return frame[frame["prestwick_id"] == selected["prestwick_id"]].copy()
    if contains_mask.sum() > 1:
        options = candidates.loc[contains_mask, ["prestwick_id", "chemical_name"]].head(10).to_dict(orient="records")
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
) -> dict[str, float]:
    """Compute composite health metrics from diversity, balance, and baseline stability."""
    diversity = _normalized_shannon(abundances)
    beneficial_fraction = _weighted_fraction(abundances, microbe_metadata, BENEFICIAL_GENERA)
    risk_fraction = _weighted_fraction(abundances, microbe_metadata, RISK_GENERA)
    stability = _stability_score(abundances, baseline_abundances)
    balance = max(0.0, min(1.0, 0.5 + beneficial_fraction - risk_fraction))
    health_index = 100.0 * (0.40 * diversity + 0.35 * balance + 0.25 * stability)
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
        tcg_health_index = float(100.0 * (0.40 * diversity + 0.35 * tcg_balance_coverage_adjusted + 0.25 * stability))
    else:
        tcg_guild_1_share = float("nan")
        tcg_guild_2_share = float("nan")
        tcg_balance = float("nan")
        tcg_balance_coverage_adjusted = float("nan")
        tcg_health_index = float("nan")
    return {
        "health_index": float(max(0.0, min(100.0, health_index))),
        "diversity": float(diversity),
        "beneficial_fraction": float(beneficial_fraction),
        "risk_fraction": float(risk_fraction),
        "stability": float(stability),
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


def run_step3_simulation(
    integrated_predictions_path: str | Path,
    output_dir: str | Path,
    drug_query: str,
    scenario_name: str = "healthy_reference",
    community_table_path: str | Path | None = None,
    tcg_proxy_mapping_path: str | Path | None = DEFAULT_TCG_PROXY_MAPPING_PATH,
    n_steps: int = 14,
    initial_dose: float = 1.0,
    repeat_dose: float = 1.0,
    dosing_interval: int = 1,
    drug_clearance_rate: float = 0.12,
    product_clearance_rate: float = 0.18,
    metabolism_scale: float = 0.85,
    effect_scale: float = 0.55,
    ecology_strength: float = 0.20,
    abundance_floor: float = 1e-6,
) -> dict[str, object]:
    """Run the Step 3 community simulation and export trajectories plus summaries.

    Args:
        integrated_predictions_path: Integrated Step 1 and Step 2 prediction table.
        output_dir: Directory where trajectory outputs and summary JSON are written.
        drug_query: Prestwick ID or drug name used to select one drug for simulation.
        scenario_name: Builtin scenario used when no custom community table is supplied.
        community_table_path: Optional custom community abundance table.
        tcg_proxy_mapping_path: Optional guild mapping file used for extra health metrics.
        n_steps: Number of discrete simulation time steps.
        initial_dose: Initial parent-drug dose at time zero.
        repeat_dose: Dose amount added at each dosing interval after the first step.
        dosing_interval: Number of time steps between repeated doses.
        drug_clearance_rate: Fraction of parent drug cleared per step.
        product_clearance_rate: Fraction of metabolite pool cleared per step.
        metabolism_scale: Global scale factor on microbial metabolism pressure.
        effect_scale: Global scale factor on Step 1 drug-pressure effects.
        ecology_strength: Strength of pull back toward the baseline community.
        abundance_floor: Minimum abundance floor applied during renormalization.

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
    applicability_map = _build_applicability_map(drug_frame)
    prestwick_id = _first_nonempty_text(drug_frame, "prestwick_id", default=drug_query)
    chemical_name = _first_nonempty_text(drug_frame, "chemical_name", default=prestwick_id or drug_query)

    current_abundances = template_abundances.copy()
    baseline_abundances = template_abundances.copy()
    cumulative_dose = float(initial_dose)
    current_parent_drug = float(initial_dose)
    current_metabolite_pool = 0.0

    history_rows: list[dict[str, object]] = []
    abundance_rows: list[dict[str, object]] = []

    initial_health = _health_index(
        current_abundances,
        baseline_abundances,
        microbe_metadata=microbe_metadata,
        tcg_membership_map=tcg_membership_map,
    )

    def _record_timepoint(timepoint: int) -> None:
        nonlocal history_rows, abundance_rows
        health = _health_index(
            current_abundances,
            baseline_abundances,
            microbe_metadata=microbe_metadata,
            tcg_membership_map=tcg_membership_map,
        )
        mean_applicability = _mean_applicability(current_abundances, applicability_map)
        parent_retention = 0.0 if cumulative_dose <= 0 else current_parent_drug / cumulative_dose
        dysbiosis_penalty = max(0.0, initial_health["health_index"] - health["health_index"])
        uncertainty_penalty = 100.0 * max(0.0, 1.0 - mean_applicability)
        efficacy_proxy = 100.0 * parent_retention
        stability_score = 100.0 * health["stability"]
        community_preservation_score = 0.65 * health["health_index"] + 0.35 * stability_score
        metabolite_burden_ratio = 0.0 if cumulative_dose <= 0 else current_metabolite_pool / cumulative_dose
        metabolite_burden_penalty = 100.0 * max(0.0, min(1.0, metabolite_burden_ratio))
        benefit_subscore = 0.55 * efficacy_proxy + 0.45 * community_preservation_score
        risk_subscore = (
            0.50 * dysbiosis_penalty
            + 0.30 * uncertainty_penalty
            + 0.20 * metabolite_burden_penalty
        )
        development_score_balance = benefit_subscore - risk_subscore
        development_score_raw = development_score_balance
        development_score = _sigmoid_score(development_score_balance, center=45.0, scale=8.0)
        development_score_legacy_raw = efficacy_proxy - 0.7 * dysbiosis_penalty - 0.3 * uncertainty_penalty
        development_score_legacy = max(0.0, min(100.0, development_score_legacy_raw))
        history_rows.append(
            {
                "timepoint": timepoint,
                "parent_drug_concentration": float(current_parent_drug),
                "aggregate_metabolite_pool": float(current_metabolite_pool),
                "cumulative_dose": float(cumulative_dose),
                "parent_retention_ratio": float(parent_retention),
                "health_index": health["health_index"],
                "diversity": health["diversity"],
                "beneficial_fraction": health["beneficial_fraction"],
                "risk_fraction": health["risk_fraction"],
                "stability": health["stability"],
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
                "uncertainty_penalty": float(uncertainty_penalty),
                "metabolite_burden_ratio": float(metabolite_burden_ratio),
                "metabolite_burden_penalty": float(metabolite_burden_penalty),
                "efficacy_proxy": float(efficacy_proxy),
                "community_preservation_score": float(community_preservation_score),
                "benefit_subscore": float(benefit_subscore),
                "risk_subscore": float(risk_subscore),
                "development_score_balance": float(development_score_balance),
                "development_score_legacy_raw": float(development_score_legacy_raw),
                "development_score_legacy": float(development_score_legacy),
                "development_score_raw": float(development_score_raw),
                "development_score": float(development_score),
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
        raw_abundances: dict[str, float] = {}
        weighted_metabolism = 0.0
        for nt_code, current_abundance in current_abundances.items():
            effect_score = _coerce_float(effect_map.get(nt_code), default=0.0)
            ecology_pull = ecology_strength * (baseline_abundances.get(nt_code, 0.0) - current_abundance)
            # Clip the exponent to keep custom-SMILES simulations numerically stable
            # even when a regressor produces outlier effect scores.
            growth_log_delta = float(np.clip(effect_scale * effect_score * drug_exposure + ecology_pull, -12.0, 12.0))
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
            weighted_metabolism += current_abundance * metabolized_probability * max(depletion_strength, 0.05 * metabolized_probability)

        current_abundances = _normalize_abundances(raw_abundances, minimum=abundance_floor)
        parent_consumed = min(current_parent_drug, current_parent_drug * metabolism_scale * weighted_metabolism)
        parent_cleared = current_parent_drug * drug_clearance_rate
        current_parent_drug = max(0.0, current_parent_drug - parent_consumed - parent_cleared + dose_input)
        current_metabolite_pool = max(
            0.0,
            current_metabolite_pool * (1.0 - product_clearance_rate) + parent_consumed,
        )
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
        "health_signature_mode": "tcg_proxy_secondary",
        "health_signature_source_path": tcg_mapping_metadata["tcg_proxy_mapping_path"],
        "health_signature_source_name": tcg_mapping_metadata["tcg_source_name"],
        "health_signature_source_url": tcg_mapping_metadata["tcg_source_url"],
        "health_signature_ready_microbe_count": int(tcg_mapping_metadata["tcg_ready_microbe_count"]),
        "health_signature_mapped_microbe_count": int(tcg_mapping_metadata["tcg_mapped_microbe_count"]),
        "health_signature_panel_microbe_count": int(tcg_mapping_metadata["tcg_panel_microbe_count"]),
        "initial_health_index": float(trajectory_metrics.iloc[0]["health_index"]),
        "final_health_index": float(final_row["health_index"]),
        "initial_tcg_health_index": float(trajectory_metrics.iloc[0]["tcg_health_index"]),
        "final_tcg_health_index": float(final_row["tcg_health_index"]),
        "final_tcg_guild_1_fraction": float(final_row["tcg_guild_1_fraction"]),
        "final_tcg_guild_2_fraction": float(final_row["tcg_guild_2_fraction"]),
        "final_tcg_mapped_fraction": float(final_row["tcg_mapped_fraction"]),
        "final_tcg_balance": float(final_row["tcg_balance"]),
        "initial_parent_drug_concentration": float(trajectory_metrics.iloc[0]["parent_drug_concentration"]),
        "final_parent_drug_concentration": float(final_row["parent_drug_concentration"]),
        "final_parent_retention_ratio": float(final_row["parent_retention_ratio"]),
        "final_aggregate_metabolite_pool": float(final_row["aggregate_metabolite_pool"]),
        "mean_applicability_final": float(final_row["mean_applicability"]),
        "final_stability": float(final_row["stability"]),
        "efficacy_proxy_final": float(final_row["efficacy_proxy"]),
        "community_preservation_final": float(final_row["community_preservation_score"]),
        "benefit_subscore_final": float(final_row["benefit_subscore"]),
        "risk_subscore_final": float(final_row["risk_subscore"]),
        "dysbiosis_penalty_final": float(final_row["dysbiosis_penalty"]),
        "uncertainty_penalty_final": float(final_row["uncertainty_penalty"]),
        "metabolite_burden_penalty_final": float(final_row["metabolite_burden_penalty"]),
        "development_score_balance": float(final_row["development_score_balance"]),
        "development_score_legacy": float(final_row["development_score_legacy"]),
        "development_score_legacy_raw": float(final_row["development_score_legacy_raw"]),
        "development_score": float(final_row["development_score"]),
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
