from __future__ import annotations

import pandas as pd

from gut_drug_microbiome.step1.hybrid import refine_step1_promote_with_step2


def test_refine_step1_promote_with_step2_uses_enzyme_support() -> None:
    base = pd.DataFrame(
        [
            {
                "step1_predicted_inhibit_probability": 0.10,
                "step1_predicted_effect_score": 0.25,
                "step1_predicted_effect_label_hybrid": "promote",
                "predicted_metabolized_probability": 0.62,
                "predicted_parent_depletion_fraction": 0.35,
                "predicted_mechanism_projection_flag": False,
                "predicted_reaction_confidence": 0.0,
                "predicted_mechanism_support_score": 0.0,
                "drug_max_fingerprint_jaccard": 0.30,
                "applicability_flag": True,
                "scaffold_seen_in_training": True,
                "microbe_genus_seen_in_training": True,
                "microbe_phylum_seen_in_training": True,
            }
        ]
    )
    with_enzyme = base.copy()
    with_enzyme["predicted_enzyme_prior_flag"] = True
    with_enzyme["predicted_enzyme_support_score"] = 0.75
    with_enzyme["predicted_enzyme_step1_promote_support_score"] = 0.60
    with_enzyme["predicted_enzyme_step1_inhibit_risk_score"] = 0.0

    without_enzyme = base.copy()
    without_enzyme["predicted_enzyme_prior_flag"] = False
    without_enzyme["predicted_enzyme_support_score"] = 0.0
    without_enzyme["predicted_enzyme_step1_promote_support_score"] = 0.0
    without_enzyme["predicted_enzyme_step1_inhibit_risk_score"] = 0.0

    enriched = refine_step1_promote_with_step2(with_enzyme)
    plain = refine_step1_promote_with_step2(without_enzyme)

    assert float(enriched.loc[0, "predicted_promote_probability_refined"]) > float(
        plain.loc[0, "predicted_promote_probability_refined"]
    )
