from __future__ import annotations

import pandas as pd

from gut_drug_microbiome.web.service import _attach_placebo_summary_deltas
from gut_drug_microbiome.web.service import _inject_placebo_deltas


def test_inject_placebo_deltas_adds_delta_columns() -> None:
    active = pd.DataFrame(
        [
            {"timepoint": 0, "development_score": 20.0, "experimental_development_score": 19.0, "health_index": 80.0},
            {"timepoint": 1, "development_score": 30.0, "experimental_development_score": 27.0, "health_index": 75.0},
        ]
    )
    placebo = pd.DataFrame(
        [
            {"timepoint": 0, "development_score": 20.0, "experimental_development_score": 18.0, "health_index": 80.0},
            {"timepoint": 1, "development_score": 10.0, "experimental_development_score": 9.5, "health_index": 78.0},
        ]
    )
    merged = _inject_placebo_deltas(active, placebo)

    assert "development_score_delta_vs_placebo" in merged.columns
    assert "experimental_development_score_delta_vs_placebo" in merged.columns
    assert "health_index_delta_vs_placebo" in merged.columns
    assert "development_score_normalized_vs_placebo" in merged.columns
    assert float(merged.loc[merged["timepoint"].eq(1), "development_score_delta_vs_placebo"].iloc[0]) == 20.0
    assert float(merged.loc[merged["timepoint"].eq(1), "experimental_development_score_delta_vs_placebo"].iloc[0]) == 17.5
    assert float(merged.loc[merged["timepoint"].eq(1), "development_score_normalized_vs_placebo"].iloc[0]) > 50.0


def test_attach_placebo_summary_deltas_reports_placebo_baseline() -> None:
    active_summary = {
        "development_score": 42.0,
        "experimental_development_score": 39.0,
        "development_score_balance": 12.0,
        "experimental_development_score_balance": 9.0,
        "final_health_index": 70.0,
        "final_parent_retention_ratio": 0.4,
        "final_experimental_aggregate_metabolite_pool": 0.6,
        "benefit_subscore_final": 55.0,
        "risk_subscore_final": 43.0,
        "experimental_risk_subscore_final": 46.0,
        "disease_target_alignment_score_final": 62.0,
    }
    placebo_summary = {
        "development_score": 30.0,
        "experimental_development_score": 29.0,
        "development_score_balance": 5.0,
        "experimental_development_score_balance": 4.0,
        "final_health_index": 72.0,
        "final_parent_retention_ratio": 0.0,
        "final_experimental_aggregate_metabolite_pool": 0.1,
        "benefit_subscore_final": 48.0,
        "risk_subscore_final": 43.0,
        "experimental_risk_subscore_final": 43.5,
        "disease_target_alignment_score_final": 50.0,
    }
    merged = _attach_placebo_summary_deltas(active_summary, placebo_summary)

    assert merged["development_score_delta_vs_placebo"] == 12.0
    assert merged["experimental_development_score_delta_vs_placebo"] == 10.0
    assert merged["development_score_balance_delta_vs_placebo"] == 7.0
    assert merged["experimental_development_score_balance_delta_vs_placebo"] == 5.0
    assert merged["health_index_delta_vs_placebo"] == -2.0
    assert merged["experimental_aggregate_metabolite_pool_delta_vs_placebo"] == 0.5
    assert merged["experimental_risk_subscore_delta_vs_placebo"] == 2.5
    assert merged["development_score_normalized_vs_placebo"] > 50.0
    assert "placebo_baseline" in merged
    assert merged["placebo_baseline"]["final_development_score"] == 30.0
    assert merged["placebo_baseline"]["final_experimental_development_score"] == 29.0
