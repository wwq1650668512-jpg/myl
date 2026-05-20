from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd

from gut_drug_microbiome.utils.text import canonicalize_key as _canonicalize_key


MECHANISM_SCORE_KEYS = [
    "anti_inflammatory_score",
    "pro_inflammatory_score",
    "butyrate_support_score",
    "barrier_protection_score",
    "toxin_risk_score",
    "mucus_degradation_score",
    "pathobiont_load",
    "competition_vs_crossfeeding_proxy",
    "mechanism_benefit_score",
    "mechanism_risk_score",
    "mechanism_balance_score",
]

FUSION_MODE_WEIGHTS: dict[str, tuple[float, float]] = {
    "raw_only": (1.0, 0.0),
    "mechanism_only": (0.0, 1.0),
    "weighted_0.3_0.7": (0.3, 0.7),
    # Backward-compatible default used by existing endpoints.
    "weighted_0.65_0.35": (0.65, 0.35),
}


# Genus-level priors distilled from high-confidence literature patterns.
GENUS_TRAIT_PRIORS: dict[str, dict[str, float]] = {
    "faecalibacterium": {
        "anti_inflammatory_weight": 1.0,
        "butyrate_weight": 1.0,
        "barrier_weight": 0.9,
    },
    "roseburia": {
        "anti_inflammatory_weight": 0.9,
        "butyrate_weight": 1.0,
        "barrier_weight": 0.8,
    },
    "bifidobacterium": {
        "anti_inflammatory_weight": 0.85,
        "butyrate_weight": 0.45,
        "barrier_weight": 0.8,
    },
    "anaerostipes": {
        "anti_inflammatory_weight": 0.8,
        "butyrate_weight": 0.9,
        "barrier_weight": 0.7,
    },
    "eubacterium": {
        "anti_inflammatory_weight": 0.7,
        "butyrate_weight": 0.8,
        "barrier_weight": 0.6,
    },
    "coprococcus": {
        "anti_inflammatory_weight": 0.7,
        "butyrate_weight": 0.75,
    },
    "butyrivibrio": {
        "anti_inflammatory_weight": 0.6,
        "butyrate_weight": 0.75,
    },
    "akkermansia": {
        "anti_inflammatory_weight": 0.7,
        "barrier_weight": 0.9,
        "mucus_degradation_weight": 0.35,
    },
    "parabacteroides": {
        "anti_inflammatory_weight": 0.45,
        "barrier_weight": 0.35,
    },
    "blautia": {
        "anti_inflammatory_weight": 0.5,
        "butyrate_weight": 0.4,
        "barrier_weight": 0.35,
    },
    "escherichia": {
        "pro_inflammatory_weight": 1.0,
        "toxin_weight": 0.75,
        "pathobiont_weight": 1.0,
    },
    "klebsiella": {
        "pro_inflammatory_weight": 0.95,
        "pathobiont_weight": 0.95,
        "toxin_weight": 0.45,
    },
    "enterococcus": {
        "pro_inflammatory_weight": 0.85,
        "pathobiont_weight": 0.9,
    },
    "fusobacterium": {
        "pro_inflammatory_weight": 1.0,
        "toxin_weight": 1.0,
        "pathobiont_weight": 1.0,
    },
    "proteus": {
        "pro_inflammatory_weight": 0.95,
        "pathobiont_weight": 0.95,
        "toxin_weight": 0.35,
    },
    "bilophila": {
        "pro_inflammatory_weight": 0.9,
        "pathobiont_weight": 0.9,
        "mucus_degradation_weight": 0.35,
    },
    "desulfovibrio": {
        "pro_inflammatory_weight": 0.8,
        "pathobiont_weight": 0.75,
    },
    "clostridium": {
        "pro_inflammatory_weight": 0.6,
        "pathobiont_weight": 0.6,
        "toxin_weight": 0.55,
    },
    "bacteroides": {
        "mucus_degradation_weight": 0.25,
        "toxin_weight": 0.35,
        "pathobiont_weight": 0.35,
    },
    "ruminococcus": {
        "mucus_degradation_weight": 0.45,
        "pro_inflammatory_weight": 0.25,
    },
}


# Species-level overrides where evidence is much stronger than genus averages.
SPECIES_TRAIT_PRIORS: dict[str, dict[str, float]] = {
    "faecalibacteriumprausnitzii": {
        "anti_inflammatory_weight": 1.0,
        "butyrate_weight": 1.0,
        "barrier_weight": 1.0,
    },
    "roseburiahominis": {
        "anti_inflammatory_weight": 1.0,
        "butyrate_weight": 1.0,
        "barrier_weight": 0.9,
    },
    "akkermansiamuciniphila": {
        "anti_inflammatory_weight": 0.85,
        "barrier_weight": 0.95,
        "mucus_degradation_weight": 0.35,
    },
    "ruminococcusgnavus": {
        "pro_inflammatory_weight": 0.95,
        "mucus_degradation_weight": 1.0,
        "pathobiont_weight": 0.85,
    },
    "ruminococcustorques": {
        "pro_inflammatory_weight": 0.85,
        "mucus_degradation_weight": 0.95,
        "pathobiont_weight": 0.7,
    },
    "escherichiacoli": {
        "pro_inflammatory_weight": 1.0,
        "toxin_weight": 0.85,
        "pathobiont_weight": 1.0,
    },
    "fusobacteriumnucleatum": {
        "pro_inflammatory_weight": 1.0,
        "toxin_weight": 1.0,
        "pathobiont_weight": 1.0,
    },
    "bacteroidesfragilis": {
        "pro_inflammatory_weight": 0.55,
        "toxin_weight": 0.8,
        "pathobiont_weight": 0.7,
    },
    "clostridiumdifficile": {
        "pro_inflammatory_weight": 1.0,
        "toxin_weight": 1.0,
        "pathobiont_weight": 1.0,
    },
    "bilophilawadsworthia": {
        "pro_inflammatory_weight": 0.95,
        "pathobiont_weight": 0.95,
        "mucus_degradation_weight": 0.35,
    },
}


def _coerce_float(value: object, default: float = 0.0) -> float:
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(numeric):
        return float(default)
    return float(numeric)


def _clip01(value: float) -> float:
    return float(np.clip(float(value), 0.0, 1.0))


def _empty_trait_weights() -> dict[str, float]:
    return {
        "anti_inflammatory_weight": 0.0,
        "pro_inflammatory_weight": 0.0,
        "butyrate_weight": 0.0,
        "barrier_weight": 0.0,
        "toxin_weight": 0.0,
        "mucus_degradation_weight": 0.0,
        "pathobiont_weight": 0.0,
    }


def _merge_weights(base: dict[str, float], overlay: dict[str, float]) -> dict[str, float]:
    merged = dict(base)
    for key, value in overlay.items():
        merged[key] = max(_coerce_float(merged.get(key), default=0.0), _coerce_float(value, default=0.0))
    return merged


def infer_microbe_trait_priors(frame: pd.DataFrame) -> pd.DataFrame:
    """Attach literature-derived mechanism priors per microbe row.

    Scientific rationale:
    - We use conservative genus priors as defaults.
    - Species-level overrides are applied only where evidence is stronger and more specific.
    """
    result = frame.copy()
    trait_rows = []
    for _, row in result.iterrows():
        species_key = _canonicalize_key(row.get("species_label")) or _canonicalize_key(row.get("microbe_label"))
        genus_key = _canonicalize_key(row.get("genus"))

        weights = _empty_trait_weights()
        if genus_key in GENUS_TRAIT_PRIORS:
            weights = _merge_weights(weights, GENUS_TRAIT_PRIORS[genus_key])
        if species_key in SPECIES_TRAIT_PRIORS:
            weights = _merge_weights(weights, SPECIES_TRAIT_PRIORS[species_key])
        trait_rows.append(weights)

    traits = pd.DataFrame(trait_rows, index=result.index)
    return pd.concat([result, traits], axis=1)


def attach_action_signals(
    frame: pd.DataFrame,
    score_column: str,
    inhibit_probability_column: str,
    promote_probability_column: str = "predicted_promote_probability_refined",
) -> pd.DataFrame:
    """Convert model outputs into signed promote/inhibit action signals.

    Scientific rationale:
    - Step1 score captures direction/magnitude of growth shift.
    - Probabilities calibrate uncertainty and soften hard-threshold artifacts.
    """
    result = frame.copy()
    score = pd.to_numeric(result.get(score_column, pd.Series(np.nan, index=result.index)), errors="coerce").fillna(0.0)
    inhibit_prob = pd.to_numeric(
        result.get(inhibit_probability_column, pd.Series(np.nan, index=result.index)),
        errors="coerce",
    ).fillna(0.0)
    promote_prob = pd.to_numeric(
        result.get(promote_probability_column, pd.Series(np.nan, index=result.index)),
        errors="coerce",
    ).fillna(0.0)

    result["mechanism_promote_signal"] = score.clip(lower=0.0) + 0.75 * promote_prob.clip(lower=0.0, upper=1.0)
    result["mechanism_inhibit_signal"] = (-score).clip(lower=0.0) + 0.75 * inhibit_prob.clip(lower=0.0, upper=1.0)
    result["mechanism_net_signal"] = result["mechanism_promote_signal"] - result["mechanism_inhibit_signal"]
    return result


def _empty_mechanism_scores() -> dict[str, float]:
    return {
        "anti_inflammatory_score": 0.0,
        "pro_inflammatory_score": 0.0,
        "butyrate_support_score": 0.0,
        "barrier_protection_score": 0.0,
        "toxin_risk_score": 0.0,
        "mucus_degradation_score": 0.0,
        "pathobiont_load": 0.0,
        "competition_vs_crossfeeding_proxy": 0.5,
        "mechanism_benefit_score": 0.0,
        "mechanism_risk_score": 0.0,
        "mechanism_balance_score": 0.0,
    }


def combine_disease_scores(
    raw_microbe_score: float,
    mechanism_balance_score: float,
    raw_weight: float = 0.65,
) -> float:
    """Blend legacy raw-microbe score and mechanism-layer balance for ablation-ready scoring."""
    raw_weight = _clip01(raw_weight)
    mechanism_weight = 1.0 - raw_weight
    return float(raw_weight * raw_microbe_score + mechanism_weight * mechanism_balance_score)


def fusion_weights(fusion_mode: str) -> tuple[float, float]:
    """Return (raw_weight, mechanism_weight) for one fusion mode."""
    mode = str(fusion_mode or "").strip()
    if mode not in FUSION_MODE_WEIGHTS:
        available = ", ".join(sorted(FUSION_MODE_WEIGHTS))
        raise ValueError(f"Unsupported fusion_mode={fusion_mode!r}. Available: {available}")
    return FUSION_MODE_WEIGHTS[mode]


def fuse_disease_scores(
    raw_score: float,
    mechanism_score: float,
    fusion_mode: str = "weighted_0.65_0.35",
) -> float:
    """Fuse raw and mechanism scores under a named mode."""
    raw_weight, mechanism_weight = fusion_weights(fusion_mode)
    return float(raw_weight * float(raw_score) + mechanism_weight * float(mechanism_score))


def compute_mechanism_layer(
    matched_rows: pd.DataFrame,
    relation_weight_column: str = "relation_weight",
    top_n_contributors: int = 5,
) -> dict[str, object]:
    """Aggregate drug->microbe effects into interpretable mechanism-layer scores.

    Scientific rationale:
    - Benefit channels reward promotion of protective/butyrate/barrier taxa and inhibition of inflammatory taxa.
    - Risk channels quantify promotion pressure on pathobiont/toxin/mucus-degradation axes.
    - Competition-vs-crossfeeding proxy summarizes interaction regime from existing Step1/Step2 signals.
    """
    if matched_rows.empty:
        return {"scores": _empty_mechanism_scores(), "top_contributors": {key: [] for key in MECHANISM_SCORE_KEYS}}

    sums = defaultdict(float)
    contributors: dict[str, list[dict[str, object]]] = defaultdict(list)
    weight_sum = 0.0
    competition_pressure = 0.0
    crossfeeding_pressure = 0.0

    for _, row in matched_rows.iterrows():
        row_weight = max(_coerce_float(row.get(relation_weight_column), default=1.0), 1e-6)
        weight_sum += row_weight

        promote_signal = max(0.0, _coerce_float(row.get("mechanism_promote_signal"), default=0.0))
        inhibit_signal = max(0.0, _coerce_float(row.get("mechanism_inhibit_signal"), default=0.0))
        net_signal = _coerce_float(row.get("mechanism_net_signal"), default=0.0)
        positive_shift = max(0.0, net_signal)
        negative_shift = max(0.0, -net_signal)

        anti_w = _clip01(_coerce_float(row.get("anti_inflammatory_weight"), default=0.0))
        pro_w = _clip01(_coerce_float(row.get("pro_inflammatory_weight"), default=0.0))
        butyrate_w = _clip01(_coerce_float(row.get("butyrate_weight"), default=0.0))
        barrier_w = _clip01(_coerce_float(row.get("barrier_weight"), default=0.0))
        toxin_w = _clip01(_coerce_float(row.get("toxin_weight"), default=0.0))
        mucus_w = _clip01(_coerce_float(row.get("mucus_degradation_weight"), default=0.0))
        pathobiont_w = _clip01(_coerce_float(row.get("pathobiont_weight"), default=0.0))

        anti_contrib = row_weight * (anti_w * positive_shift + 0.65 * pro_w * negative_shift)
        pro_contrib = row_weight * (pro_w * positive_shift + 0.35 * anti_w * negative_shift)
        butyrate_contrib = row_weight * butyrate_w * positive_shift
        barrier_contrib = row_weight * barrier_w * positive_shift
        toxin_contrib = row_weight * toxin_w * positive_shift
        mucus_contrib = row_weight * mucus_w * positive_shift
        pathobiont_contrib = row_weight * pathobiont_w * positive_shift

        sums["anti_inflammatory_score"] += anti_contrib
        sums["pro_inflammatory_score"] += pro_contrib
        sums["butyrate_support_score"] += butyrate_contrib
        sums["barrier_protection_score"] += barrier_contrib
        sums["toxin_risk_score"] += toxin_contrib
        sums["mucus_degradation_score"] += mucus_contrib
        sums["pathobiont_load"] += pathobiont_contrib

        cross_ref = 1.0 if bool(row.get("predicted_cross_feeding_reference_flag", False)) else 0.0
        enzyme_promote_support = _clip01(_coerce_float(row.get("predicted_enzyme_step1_promote_support_score"), default=0.0))
        enzyme_inhibit_risk = _clip01(_coerce_float(row.get("predicted_enzyme_step1_inhibit_risk_score"), default=0.0))
        mechanism_support = _clip01(_coerce_float(row.get("predicted_mechanism_support_score"), default=0.0))
        depletion = _clip01(max(0.0, -_coerce_float(row.get("predicted_parent_depletion_fraction"), default=0.0)))
        applicability = 1.0 if bool(row.get("applicability_flag", False)) else 0.0

        crossfeeding_pressure += row_weight * (
            0.45 * cross_ref + 0.30 * _clip01(promote_signal) + 0.15 * enzyme_promote_support + 0.10 * mechanism_support
        )
        competition_pressure += row_weight * (
            0.50 * _clip01(inhibit_signal) + 0.20 * enzyme_inhibit_risk + 0.20 * depletion + 0.10 * applicability
        )

        row_identity = {
            "nt_code": row.get("nt_code"),
            "microbe_label": row.get("microbe_label"),
            "genus": row.get("genus"),
            "species_label": row.get("species_label"),
            "desired_step1_effect": row.get("desired_step1_effect"),
            "relation_weight": float(row_weight),
            "net_signal": float(net_signal),
        }
        for key, contribution in [
            ("anti_inflammatory_score", anti_contrib),
            ("pro_inflammatory_score", pro_contrib),
            ("butyrate_support_score", butyrate_contrib),
            ("barrier_protection_score", barrier_contrib),
            ("toxin_risk_score", toxin_contrib),
            ("mucus_degradation_score", mucus_contrib),
            ("pathobiont_load", pathobiont_contrib),
        ]:
            if contribution <= 0:
                continue
            payload = dict(row_identity)
            payload["contribution"] = float(contribution)
            contributors[key].append(payload)

    normalization = max(weight_sum, 1e-8)
    scores = _empty_mechanism_scores()
    for key in [
        "anti_inflammatory_score",
        "pro_inflammatory_score",
        "butyrate_support_score",
        "barrier_protection_score",
        "toxin_risk_score",
        "mucus_degradation_score",
        "pathobiont_load",
    ]:
        scores[key] = _clip01(sums[key] / normalization)

    scores["competition_vs_crossfeeding_proxy"] = _clip01(
        competition_pressure / max(competition_pressure + crossfeeding_pressure, 1e-8)
    )
    scores["mechanism_benefit_score"] = _clip01(
        0.30 * scores["anti_inflammatory_score"]
        + 0.25 * scores["butyrate_support_score"]
        + 0.20 * scores["barrier_protection_score"]
        + 0.25 * scores["competition_vs_crossfeeding_proxy"]
    )
    scores["mechanism_risk_score"] = _clip01(
        0.35 * scores["pro_inflammatory_score"]
        + 0.30 * scores["pathobiont_load"]
        + 0.20 * scores["toxin_risk_score"]
        + 0.15 * scores["mucus_degradation_score"]
    )
    scores["mechanism_balance_score"] = float(scores["mechanism_benefit_score"] - scores["mechanism_risk_score"])

    top_contributors: dict[str, list[dict[str, object]]] = {}
    for key in MECHANISM_SCORE_KEYS:
        rows = contributors.get(key, [])
        rows = sorted(rows, key=lambda item: item["contribution"], reverse=True)
        top_contributors[key] = rows[:top_n_contributors]

    return {
        "scores": scores,
        "top_contributors": top_contributors,
        "weight_sum": float(weight_sum),
        "competition_pressure": float(competition_pressure),
        "crossfeeding_pressure": float(crossfeeding_pressure),
        "n_rows": int(len(matched_rows)),
    }
