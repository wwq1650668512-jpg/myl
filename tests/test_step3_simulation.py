from __future__ import annotations

import pandas as pd

from gut_drug_microbiome.step3 import run_step3_simulation


def test_run_step3_simulation_writes_expected_outputs(tmp_path) -> None:
    prediction_table = pd.DataFrame(
        [
            {
                "prestwick_id": "P1",
                "chemical_name": "Drug One",
                "nt_code": "NT1",
                "microbe_label": "Microbe A",
                "species_label": "Bifidobacterium longum",
                "genus": "Bifidobacterium",
                "phylum": "Actinobacteria",
                "step1_predicted_effect_score": 0.1,
                "step1_predicted_inhibit_probability": 0.2,
                "predicted_metabolized_probability": 0.4,
                "predicted_parent_depletion_fraction": -0.3,
                "drug_max_fingerprint_jaccard": 0.7,
                "scaffold_seen_in_training": True,
                "microbe_genus_seen_in_training": True,
                "microbe_phylum_seen_in_training": True,
                "applicability_flag": True,
            },
            {
                "prestwick_id": "P1",
                "chemical_name": "Drug One",
                "nt_code": "NT2",
                "microbe_label": "Microbe B",
                "species_label": "Escherichia coli",
                "genus": "Escherichia",
                "phylum": "Proteobacteria",
                "step1_predicted_effect_score": -0.2,
                "step1_predicted_inhibit_probability": 0.8,
                "predicted_metabolized_probability": 0.6,
                "predicted_parent_depletion_fraction": -0.5,
                "drug_max_fingerprint_jaccard": 0.5,
                "scaffold_seen_in_training": False,
                "microbe_genus_seen_in_training": True,
                "microbe_phylum_seen_in_training": True,
                "applicability_flag": True,
            },
        ]
    )
    predictions_path = tmp_path / "predictions.csv"
    prediction_table.to_csv(predictions_path, index=False)

    summary = run_step3_simulation(
        integrated_predictions_path=predictions_path,
        output_dir=tmp_path / "simulation",
        drug_query="P1",
        scenario_name="healthy_reference",
        tcg_proxy_mapping_path=None,
        n_steps=4,
    )

    metrics = pd.read_csv(tmp_path / "simulation" / "trajectory_metrics.csv")
    abundances = pd.read_csv(tmp_path / "simulation" / "trajectory_abundances.csv")

    assert summary["prestwick_id"] == "P1"
    assert summary["n_steps"] == 4
    assert summary["interaction_positive_edge_count"] >= 0
    assert "final_interaction_balance_rho" in summary
    assert len(metrics) == 5
    assert set(abundances["nt_code"]) == {"NT1", "NT2"}
    per_timepoint_sum = abundances.groupby("timepoint")["abundance"].sum().round(6)
    assert per_timepoint_sum.eq(1.0).all()
    assert {"interaction_balance_rho", "interaction_balance_shift", "health_index_legacy"} <= set(metrics.columns)


def test_run_step3_simulation_uses_curated_cross_feeding_reference(tmp_path) -> None:
    prediction_table = pd.DataFrame(
        [
            {
                "prestwick_id": "P2",
                "chemical_name": "QRR",
                "nt_code": "NTA",
                "microbe_label": "Producer A",
                "species_label": "Bacteroides uniformis",
                "species_name": "Bacteroides uniformis",
                "genus": "Bacteroides",
                "family": "Bacteroidaceae",
                "phylum": "Bacteroidetes",
                "step1_predicted_effect_score": 0.0,
                "step1_predicted_inhibit_probability": 0.2,
                "predicted_metabolized_probability": 0.7,
                "predicted_parent_depletion_fraction": -0.6,
                "drug_max_fingerprint_jaccard": 0.8,
                "scaffold_seen_in_training": True,
                "microbe_genus_seen_in_training": True,
                "microbe_phylum_seen_in_training": True,
                "applicability_flag": True,
                "compound_name_normalized": "qrr",
                "compound_semantic_family": "flavonoid_glycoside",
                "predicted_enzyme_step1_promote_support_score": 0.7,
                "predicted_candidate_product_count": 2,
                "predicted_enzyme_ids": "ENZ008",
            },
            {
                "prestwick_id": "P2",
                "chemical_name": "QRR",
                "nt_code": "NTB",
                "microbe_label": "Consumer B",
                "species_label": "Akkermansia muciniphila",
                "species_name": "Akkermansia muciniphila",
                "genus": "Akkermansia",
                "family": "Akkermansiaceae",
                "phylum": "Verrucomicrobia",
                "step1_predicted_effect_score": 0.5,
                "step1_predicted_inhibit_probability": 0.2,
                "predicted_metabolized_probability": 0.3,
                "predicted_parent_depletion_fraction": -0.2,
                "drug_max_fingerprint_jaccard": 0.8,
                "scaffold_seen_in_training": True,
                "microbe_genus_seen_in_training": True,
                "microbe_phylum_seen_in_training": True,
                "applicability_flag": True,
                "compound_name_normalized": "qrr",
                "compound_semantic_family": "flavonoid_glycoside",
                "predicted_enzyme_step1_promote_support_score": 0.8,
                "predicted_candidate_product_count": 1,
                "predicted_enzyme_ids": "ENZ015",
            },
        ]
    )
    predictions_path = tmp_path / "predictions.csv"
    prediction_table.to_csv(predictions_path, index=False)

    summary = run_step3_simulation(
        integrated_predictions_path=predictions_path,
        output_dir=tmp_path / "simulation",
        drug_query="P2",
        scenario_name="healthy_reference",
        tcg_proxy_mapping_path=None,
        n_steps=3,
    )

    metrics = pd.read_csv(tmp_path / "simulation" / "trajectory_metrics.csv")

    assert summary["interaction_reference_edge_count"] >= 1
    assert summary["interaction_positive_edge_count"] >= 1
    assert metrics["positive_interaction_strength"].max() > 0


def test_run_step3_simulation_supports_disease_target_reward(tmp_path) -> None:
    prediction_table = pd.DataFrame(
        [
            {
                "prestwick_id": "P3",
                "chemical_name": "Drug Target",
                "nt_code": "NTA",
                "microbe_label": "Bifido A",
                "species_label": "Bifidobacterium longum",
                "genus": "Bifidobacterium",
                "phylum": "Actinobacteria",
                "step1_predicted_effect_score": 1.2,
                "step1_predicted_inhibit_probability": 0.1,
                "predicted_metabolized_probability": 0.3,
                "predicted_parent_depletion_fraction": -0.2,
                "drug_max_fingerprint_jaccard": 0.8,
                "scaffold_seen_in_training": True,
                "microbe_genus_seen_in_training": True,
                "microbe_phylum_seen_in_training": True,
                "applicability_flag": True,
            },
            {
                "prestwick_id": "P3",
                "chemical_name": "Drug Target",
                "nt_code": "NTB",
                "microbe_label": "Escherichia B",
                "species_label": "Escherichia coli",
                "genus": "Escherichia",
                "phylum": "Proteobacteria",
                "step1_predicted_effect_score": -0.2,
                "step1_predicted_inhibit_probability": 0.7,
                "predicted_metabolized_probability": 0.4,
                "predicted_parent_depletion_fraction": -0.4,
                "drug_max_fingerprint_jaccard": 0.8,
                "scaffold_seen_in_training": True,
                "microbe_genus_seen_in_training": True,
                "microbe_phylum_seen_in_training": True,
                "applicability_flag": True,
            },
        ]
    )
    predictions_path = tmp_path / "predictions.csv"
    prediction_table.to_csv(predictions_path, index=False)

    promote_summary = run_step3_simulation(
        integrated_predictions_path=predictions_path,
        output_dir=tmp_path / "sim_promote",
        drug_query="P3",
        scenario_name="healthy_reference",
        tcg_proxy_mapping_path=None,
        n_steps=4,
        disease_target_profile={"NTA": 1.0},
    )
    inhibit_summary = run_step3_simulation(
        integrated_predictions_path=predictions_path,
        output_dir=tmp_path / "sim_inhibit",
        drug_query="P3",
        scenario_name="healthy_reference",
        tcg_proxy_mapping_path=None,
        n_steps=4,
        disease_target_profile={"NTA": -1.0},
    )
    promote_metrics = pd.read_csv(tmp_path / "sim_promote" / "trajectory_metrics.csv")

    assert promote_summary["disease_target_reward_enabled"] is True
    assert promote_summary["disease_target_profile_size"] == 1
    assert "disease_target_alignment_score" in promote_metrics.columns
    assert promote_summary["disease_target_alignment_score_final"] > 50.0
    assert promote_summary["benefit_subscore_final"] > inhibit_summary["benefit_subscore_final"]


def test_run_step3_simulation_exposes_additive_multi_product_experiment(tmp_path) -> None:
    prediction_table = pd.DataFrame(
        [
            {
                "prestwick_id": "P4",
                "chemical_name": "Branch Drug",
                "nt_code": "NTA",
                "microbe_label": "Metabolizer A",
                "species_label": "Bifidobacterium longum",
                "genus": "Bifidobacterium",
                "phylum": "Actinobacteria",
                "step1_predicted_effect_score": 0.0,
                "step1_predicted_inhibit_probability": 0.2,
                "predicted_metabolized_probability": 0.8,
                "predicted_parent_depletion_fraction": -0.6,
                "drug_max_fingerprint_jaccard": 0.9,
                "scaffold_seen_in_training": True,
                "microbe_genus_seen_in_training": True,
                "microbe_phylum_seen_in_training": True,
                "applicability_flag": True,
                "experimental_biotransform_product_count": 4,
                "experimental_biotransform_fraction_in_gut": 0.6,
            },
            {
                "prestwick_id": "P4",
                "chemical_name": "Branch Drug",
                "nt_code": "NTB",
                "microbe_label": "Metabolizer B",
                "species_label": "Escherichia coli",
                "genus": "Escherichia",
                "phylum": "Proteobacteria",
                "step1_predicted_effect_score": 0.0,
                "step1_predicted_inhibit_probability": 0.2,
                "predicted_metabolized_probability": 0.7,
                "predicted_parent_depletion_fraction": -0.5,
                "drug_max_fingerprint_jaccard": 0.8,
                "scaffold_seen_in_training": True,
                "microbe_genus_seen_in_training": True,
                "microbe_phylum_seen_in_training": True,
                "applicability_flag": True,
                "experimental_biotransform_product_count": 4,
                "experimental_biotransform_fraction_in_gut": 0.6,
            },
        ]
    )
    predictions_path = tmp_path / "predictions.csv"
    prediction_table.to_csv(predictions_path, index=False)

    summary = run_step3_simulation(
        integrated_predictions_path=predictions_path,
        output_dir=tmp_path / "simulation",
        drug_query="P4",
        scenario_name="healthy_reference",
        tcg_proxy_mapping_path=None,
        n_steps=4,
        experimental_multi_product_enabled=True,
        experimental_branching_scale=0.5,
        experimental_secondary_metabolism_rate=0.2,
    )
    metrics = pd.read_csv(tmp_path / "simulation" / "trajectory_metrics.csv")

    assert summary["experimental_multi_product_enabled"] is True
    assert summary["experimental_product_annotation_pairs"] == 2
    assert summary["final_experimental_aggregate_metabolite_pool"] >= summary["final_aggregate_metabolite_pool"]
    assert summary["experimental_metabolite_burden_penalty_final"] >= summary["metabolite_burden_penalty_final"]
    assert summary["experimental_development_score"] <= summary["development_score"]
    assert "experimental_aggregate_metabolite_pool" in metrics.columns
    assert "experimental_development_score" in metrics.columns
