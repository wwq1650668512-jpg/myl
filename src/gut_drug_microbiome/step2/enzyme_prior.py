from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

from gut_drug_microbiome.step1.compound_semantics import annotate_compound_semantics
from gut_drug_microbiome.utils.text import canonicalize_key
from gut_drug_microbiome.utils.text import normalize_whitespace


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ENZYME_CATALOG_PATH = ROOT / "data/reference/step2_enzyme_function_catalog.csv"
DEFAULT_MICROBE_ENZYME_PANEL_PATH = ROOT / "data/reference/step2_microbe_enzyme_prior_long.csv"
DEFAULT_MICROBE_ENZYME_MATRIX_PATH = ROOT / "data/reference/step2_microbe_enzyme_prior_matrix.csv"
DEFAULT_MICROBE_ENZYME_EVIDENCE_LEDGER_PATH = ROOT / "data/reference/step2_microbe_enzyme_evidence_ledger.csv"
DEFAULT_MICROBE_ENZYME_CURATION_TEMPLATE_PATH = ROOT / "data/reference/step2_microbe_enzyme_curation_template.csv"
DEFAULT_ENZYME_SUMMARY_PATH = ROOT / "data/reference/step2_enzyme_prior_summary.json"

PRESENCE_WEIGHT_MAP = {
    "curated_present": 1.00,
    "likely_present": 0.80,
    "genus_prior": 0.60,
    "weak_prior": 0.35,
    "absent": 0.00,
    "unknown": math.nan,
}

EVIDENCE_SCOPE_PRIORITY = {
    "strainliterature": 6.0,
    "straingenome": 5.5,
    "speciesliterature": 5.0,
    "speciesgenome": 4.5,
    "genusliterature": 3.5,
    "manualcuration": 3.0,
    "genuspriorseed": 1.0,
}

CURATION_STATUS_PRIORITY = {
    "reviewedaccepted": 0.8,
    "curatedliterature": 0.65,
    "manualcurated": 0.55,
    "starterseed": 0.0,
    "needscuration": -0.2,
}

STRAIN_MATCH_PRIORITY = {
    "samestrain": 0.7,
    "typestrain": 0.6,
    "samespecies": 0.4,
    "samegenus": 0.1,
}

ENZYME_CATALOG_SEED: tuple[dict[str, object], ...] = (
    {
        "enzyme_id": "ENZ001",
        "enzyme_name": "beta_glucuronidase",
        "enzyme_family": "glycosidase",
        "ec_number": "EC 3.2.1.31",
        "reaction_class": "deconjugation",
        "bond_target": "glucuronide glycosidic C-O bond",
        "substrate_scope": "host_phase_ii_conjugates_and_xenobiotics",
        "compound_semantic_families": "",
        "substrate_keywords": "glucuronide;glucuronidated;beta-d-glucuronide",
        "likely_products_or_outcomes": "aglycone release;re-activation of conjugated metabolites",
        "step2_mechanistic_role": "direct_drug_or_metabolite_deconjugation",
        "step1_feedback_role": "mixed_detox_and_reactivation",
        "metabolism_weight": 0.95,
        "step1_promote_weight": 0.20,
        "step1_inhibit_weight": 0.35,
        "notes": "Important for deconjugation of host/drug glucuronides and enterohepatic recycling.",
    },
    {
        "enzyme_id": "ENZ002",
        "enzyme_name": "sulfatase",
        "enzyme_family": "hydrolase",
        "ec_number": "EC 3.1.6.-",
        "reaction_class": "deconjugation",
        "bond_target": "sulfate ester bond",
        "substrate_scope": "sulfated host molecules and sulfated xenobiotics",
        "compound_semantic_families": "",
        "substrate_keywords": "sulfate;sulfated;sulfate ester",
        "likely_products_or_outcomes": "desulfation;aglycone release",
        "step2_mechanistic_role": "direct_deconjugation",
        "step1_feedback_role": "mixed_detox_and_reactivation",
        "metabolism_weight": 0.80,
        "step1_promote_weight": 0.18,
        "step1_inhibit_weight": 0.22,
        "notes": "Gut sulfatases can release phenolic or host-derived sulfated compounds.",
    },
    {
        "enzyme_id": "ENZ003",
        "enzyme_name": "azoreductase",
        "enzyme_family": "reductase",
        "ec_number": "EC 1.7.-.-",
        "reaction_class": "reduction",
        "bond_target": "azo N=N bond",
        "substrate_scope": "azo dyes and azo-linked prodrugs",
        "compound_semantic_families": "",
        "substrate_keywords": "azo;diazo;sulfasalazine;prontosil",
        "likely_products_or_outcomes": "azo bond cleavage;amine release",
        "step2_mechanistic_role": "direct_drug_activation_or_cleavage",
        "step1_feedback_role": "bioactivation_risk",
        "metabolism_weight": 0.95,
        "step1_promote_weight": 0.05,
        "step1_inhibit_weight": 0.35,
        "notes": "A classic gut-drug metabolism mechanism for azo compounds.",
    },
    {
        "enzyme_id": "ENZ004",
        "enzyme_name": "nitroreductase",
        "enzyme_family": "reductase",
        "ec_number": "EC 1.7.1.-",
        "reaction_class": "reduction",
        "bond_target": "nitro group",
        "substrate_scope": "nitroaromatics and nitroheterocycles",
        "compound_semantic_families": "",
        "substrate_keywords": "nitro;nitrofuran;nitroimidazole;nitroaromatic",
        "likely_products_or_outcomes": "nitro reduction to hydroxylamine or amine",
        "step2_mechanistic_role": "direct_drug_reduction",
        "step1_feedback_role": "bioactivation_risk",
        "metabolism_weight": 0.90,
        "step1_promote_weight": 0.05,
        "step1_inhibit_weight": 0.30,
        "notes": "Can contribute to activation or detoxification depending on the compound.",
    },
    {
        "enzyme_id": "ENZ005",
        "enzyme_name": "carboxylesterase",
        "enzyme_family": "hydrolase",
        "ec_number": "EC 3.1.1.-",
        "reaction_class": "hydrolysis",
        "bond_target": "carboxylic ester bond",
        "substrate_scope": "ester-containing drugs, prodrugs, and lipophilic esters",
        "compound_semantic_families": "",
        "substrate_keywords": "ester;acetate ester;prodrug ester",
        "likely_products_or_outcomes": "ester hydrolysis;acid and alcohol release",
        "step2_mechanistic_role": "direct_drug_hydrolysis",
        "step1_feedback_role": "detox_or_nutrient_release",
        "metabolism_weight": 0.85,
        "step1_promote_weight": 0.22,
        "step1_inhibit_weight": 0.12,
        "notes": "A broad class useful for ester cleavage priors.",
    },
    {
        "enzyme_id": "ENZ006",
        "enzyme_name": "amidase",
        "enzyme_family": "hydrolase",
        "ec_number": "EC 3.5.1.-",
        "reaction_class": "hydrolysis",
        "bond_target": "amide bond",
        "substrate_scope": "amide-containing drugs and simple amides",
        "compound_semantic_families": "",
        "substrate_keywords": "amide;acetamide;lactam",
        "likely_products_or_outcomes": "amide hydrolysis;acid and amine release",
        "step2_mechanistic_role": "direct_drug_hydrolysis",
        "step1_feedback_role": "detox_or_nutrient_release",
        "metabolism_weight": 0.78,
        "step1_promote_weight": 0.18,
        "step1_inhibit_weight": 0.10,
        "notes": "Broad hydrolase prior for drug amides and related bonds.",
    },
    {
        "enzyme_id": "ENZ007",
        "enzyme_name": "beta_glucosidase",
        "enzyme_family": "glycosidase",
        "ec_number": "EC 3.2.1.21",
        "reaction_class": "deconjugation",
        "bond_target": "beta-glucosidic bond",
        "substrate_scope": "plant glycosides, flavonoid glycosides, and oligosaccharides",
        "compound_semantic_families": "flavonoid_glycoside;plant_polysaccharide",
        "substrate_keywords": "glucoside;glycoside;cellobiose;arbutin;salicin",
        "likely_products_or_outcomes": "aglycone release;simple sugar release",
        "step2_mechanistic_role": "dietary_xenobiotic_and_glycoside_processing",
        "step1_feedback_role": "cross_feeding_support",
        "metabolism_weight": 0.60,
        "step1_promote_weight": 0.75,
        "step1_inhibit_weight": 0.05,
        "notes": "Important for releasing aglycones and fermentable sugars from glycosides.",
    },
    {
        "enzyme_id": "ENZ008",
        "enzyme_name": "alpha_rhamnosidase",
        "enzyme_family": "glycosidase",
        "ec_number": "EC 3.2.1.40",
        "reaction_class": "deconjugation",
        "bond_target": "rhamnoside bond",
        "substrate_scope": "rhamnose-containing flavonoid glycosides",
        "compound_semantic_families": "flavonoid_glycoside",
        "substrate_keywords": "rhamnoside;rutinoside;rhamnose;rutinose;rutin",
        "likely_products_or_outcomes": "deglycosylation;aglycone release",
        "step2_mechanistic_role": "glycoside_processing",
        "step1_feedback_role": "cross_feeding_support",
        "metabolism_weight": 0.68,
        "step1_promote_weight": 0.80,
        "step1_inhibit_weight": 0.04,
        "notes": "Useful for flavonoid rutinosides and similar plant-derived glycosides.",
    },
    {
        "enzyme_id": "ENZ009",
        "enzyme_name": "beta_galactosidase",
        "enzyme_family": "glycosidase",
        "ec_number": "EC 3.2.1.23",
        "reaction_class": "deconjugation",
        "bond_target": "beta-galactosidic bond",
        "substrate_scope": "galactosides and galactan-containing substrates",
        "compound_semantic_families": "plant_polysaccharide",
        "substrate_keywords": "galactoside;galactan;lactose;arabinogalactan",
        "likely_products_or_outcomes": "simple sugar release",
        "step2_mechanistic_role": "carbohydrate_processing",
        "step1_feedback_role": "cross_feeding_support",
        "metabolism_weight": 0.52,
        "step1_promote_weight": 0.72,
        "step1_inhibit_weight": 0.03,
        "notes": "Supports utilization of galactose-containing glycans and dietary substrates.",
    },
    {
        "enzyme_id": "ENZ010",
        "enzyme_name": "bile_salt_hydrolase",
        "enzyme_family": "hydrolase",
        "ec_number": "EC 3.5.1.24",
        "reaction_class": "deconjugation",
        "bond_target": "amide bond in conjugated bile salts",
        "substrate_scope": "host conjugated bile acids",
        "compound_semantic_families": "",
        "substrate_keywords": "taurocholate;glycocholate;bile salt;conjugated bile acid",
        "likely_products_or_outcomes": "bile acid deconjugation",
        "step2_mechanistic_role": "host_metabolite_processing",
        "step1_feedback_role": "community_context_modulation",
        "metabolism_weight": 0.35,
        "step1_promote_weight": 0.42,
        "step1_inhibit_weight": 0.06,
        "notes": "Primarily host-bile related, but relevant to community-level ecology.",
    },
    {
        "enzyme_id": "ENZ011",
        "enzyme_name": "o_demethylase",
        "enzyme_family": "transferase_complex",
        "ec_number": "NA",
        "reaction_class": "demethylation",
        "bond_target": "aryl O-CH3 bond",
        "substrate_scope": "methoxylated aromatics and polyphenols",
        "compound_semantic_families": "polyphenol;catechin_gallate",
        "substrate_keywords": "methoxy;guaiacol;anisole;vanillate;ferulate",
        "likely_products_or_outcomes": "O-demethylation;catechol-like intermediates",
        "step2_mechanistic_role": "aromatic_xenobiotic_transformation",
        "step1_feedback_role": "metabolism_supported_promote",
        "metabolism_weight": 0.72,
        "step1_promote_weight": 0.30,
        "step1_inhibit_weight": 0.12,
        "notes": "A useful prior for anaerobic polyphenol biotransformation.",
    },
    {
        "enzyme_id": "ENZ012",
        "enzyme_name": "dehydroxylase",
        "enzyme_family": "lyase_or_reductive_complex",
        "ec_number": "NA",
        "reaction_class": "dehydroxylation",
        "bond_target": "aryl C-OH functionality",
        "substrate_scope": "polyphenols, catechols, and steroid-like substrates",
        "compound_semantic_families": "polyphenol;catechin_gallate",
        "substrate_keywords": "catechol;polyphenol;hydroxylated aromatic",
        "likely_products_or_outcomes": "dehydroxylated aromatic products",
        "step2_mechanistic_role": "anaerobic_aromatic_transformation",
        "step1_feedback_role": "metabolism_supported_promote",
        "metabolism_weight": 0.68,
        "step1_promote_weight": 0.28,
        "step1_inhibit_weight": 0.10,
        "notes": "Represents anaerobic reduction/removal of hydroxyl substituents.",
    },
    {
        "enzyme_id": "ENZ013",
        "enzyme_name": "deacetylase",
        "enzyme_family": "hydrolase",
        "ec_number": "EC 3.1.1.-",
        "reaction_class": "deacetylation",
        "bond_target": "acetyl ester or acetamide group",
        "substrate_scope": "acetylated carbohydrates and acetylated xenobiotics",
        "compound_semantic_families": "plant_polysaccharide",
        "substrate_keywords": "acetyl;acetate;acetylated",
        "likely_products_or_outcomes": "deacetylation;acetate release",
        "step2_mechanistic_role": "deprotection_and_carbohydrate_access",
        "step1_feedback_role": "cross_feeding_support",
        "metabolism_weight": 0.62,
        "step1_promote_weight": 0.46,
        "step1_inhibit_weight": 0.08,
        "notes": "Useful for loosening acetylated glycans or acetylated small molecules.",
    },
    {
        "enzyme_id": "ENZ014",
        "enzyme_name": "phosphatase",
        "enzyme_family": "hydrolase",
        "ec_number": "EC 3.1.3.-",
        "reaction_class": "hydrolysis",
        "bond_target": "phosphate ester bond",
        "substrate_scope": "phosphorylated metabolites and phosphate esters",
        "compound_semantic_families": "",
        "substrate_keywords": "phosphate;phospho;phosphorylated",
        "likely_products_or_outcomes": "dephosphorylation",
        "step2_mechanistic_role": "cofactor_or_metabolite_processing",
        "step1_feedback_role": "nutrient_release_support",
        "metabolism_weight": 0.45,
        "step1_promote_weight": 0.24,
        "step1_inhibit_weight": 0.04,
        "notes": "A broad hydrolysis prior for phosphate-containing substrates.",
    },
    {
        "enzyme_id": "ENZ015",
        "enzyme_name": "mucin_glycosidase",
        "enzyme_family": "glycosidase_complex",
        "ec_number": "NA",
        "reaction_class": "deconjugation",
        "bond_target": "host glycan glycosidic bonds",
        "substrate_scope": "mucin and host-derived glycans",
        "compound_semantic_families": "",
        "substrate_keywords": "mucin;host glycan;fucose;sialic acid",
        "likely_products_or_outcomes": "host glycan liberation;community cross-feeding",
        "step2_mechanistic_role": "host_glycan_processing",
        "step1_feedback_role": "cross_feeding_support",
        "metabolism_weight": 0.30,
        "step1_promote_weight": 0.58,
        "step1_inhibit_weight": 0.06,
        "notes": "Especially relevant to mucin-degrading taxa such as Akkermansia.",
    },
    {
        "enzyme_id": "ENZ016",
        "enzyme_name": "polysaccharide_lyase",
        "enzyme_family": "lyase",
        "ec_number": "EC 4.2.2.-",
        "reaction_class": "ring_cleavage",
        "bond_target": "complex polysaccharide backbone",
        "substrate_scope": "dietary polysaccharides and glycan polymers",
        "compound_semantic_families": "plant_polysaccharide;mannan_oligosaccharide",
        "substrate_keywords": "galactomannan;mannan;oligosaccharide;pectin;alginate",
        "likely_products_or_outcomes": "oligosaccharide release;cross-feeding substrates",
        "step2_mechanistic_role": "substrate_liberation_for_community_metabolism",
        "step1_feedback_role": "cross_feeding_support",
        "metabolism_weight": 0.40,
        "step1_promote_weight": 0.82,
        "step1_inhibit_weight": 0.02,
        "notes": "A proxy for PUL-like carbohydrate access functions in glycan degraders.",
    },
)

GENUS_ENZYME_PRIORS: dict[str, list[tuple[str, str]]] = {
    "akkermansia": [("ENZ002", "likely_present"), ("ENZ009", "likely_present"), ("ENZ015", "curated_present")],
    "alistipes": [("ENZ001", "genus_prior"), ("ENZ002", "genus_prior"), ("ENZ007", "weak_prior")],
    "bacteroides": [
        ("ENZ001", "likely_present"),
        ("ENZ002", "likely_present"),
        ("ENZ005", "genus_prior"),
        ("ENZ006", "genus_prior"),
        ("ENZ007", "likely_present"),
        ("ENZ008", "genus_prior"),
        ("ENZ009", "likely_present"),
        ("ENZ014", "genus_prior"),
        ("ENZ016", "likely_present"),
    ],
    "bifidobacterium": [("ENZ007", "likely_present"), ("ENZ008", "likely_present"), ("ENZ009", "likely_present"), ("ENZ010", "genus_prior")],
    "bilophila": [("ENZ002", "weak_prior"), ("ENZ010", "weak_prior")],
    "blautia": [("ENZ001", "weak_prior"), ("ENZ011", "genus_prior"), ("ENZ012", "weak_prior"), ("ENZ013", "weak_prior")],
    "butyrivibrio": [("ENZ007", "genus_prior"), ("ENZ011", "weak_prior"), ("ENZ013", "weak_prior"), ("ENZ016", "genus_prior")],
    "citrobacter": [("ENZ001", "likely_present"), ("ENZ003", "genus_prior"), ("ENZ004", "likely_present"), ("ENZ006", "genus_prior")],
    "clostridium": [("ENZ001", "genus_prior"), ("ENZ006", "weak_prior"), ("ENZ011", "likely_present"), ("ENZ012", "likely_present"), ("ENZ013", "genus_prior")],
    "coprococcus": [("ENZ011", "weak_prior"), ("ENZ012", "weak_prior"), ("ENZ013", "weak_prior")],
    "dialister": [("ENZ004", "weak_prior")],
    "dorea": [("ENZ011", "weak_prior"), ("ENZ012", "weak_prior"), ("ENZ013", "weak_prior")],
    "enterococcus": [("ENZ007", "genus_prior"), ("ENZ009", "likely_present"), ("ENZ010", "likely_present")],
    "escherichia": [("ENZ001", "curated_present"), ("ENZ003", "likely_present"), ("ENZ004", "curated_present"), ("ENZ005", "genus_prior"), ("ENZ006", "genus_prior")],
    "eubacterium": [("ENZ001", "genus_prior"), ("ENZ011", "likely_present"), ("ENZ012", "likely_present"), ("ENZ013", "genus_prior")],
    "faecalibacterium": [("ENZ011", "genus_prior"), ("ENZ012", "weak_prior"), ("ENZ013", "weak_prior")],
    "haemophilus": [("ENZ004", "weak_prior")],
    "klebsiella": [("ENZ001", "genus_prior"), ("ENZ003", "genus_prior"), ("ENZ004", "likely_present"), ("ENZ006", "genus_prior")],
    "lactobacillus": [("ENZ005", "weak_prior"), ("ENZ007", "likely_present"), ("ENZ009", "likely_present"), ("ENZ010", "likely_present")],
    "lactiplantibacillus": [("ENZ005", "weak_prior"), ("ENZ007", "likely_present"), ("ENZ009", "likely_present"), ("ENZ010", "likely_present")],
    "limosilactobacillus": [("ENZ005", "weak_prior"), ("ENZ007", "likely_present"), ("ENZ009", "likely_present"), ("ENZ010", "likely_present")],
    "parabacteroides": [("ENZ001", "genus_prior"), ("ENZ002", "genus_prior"), ("ENZ007", "weak_prior"), ("ENZ009", "genus_prior"), ("ENZ016", "genus_prior")],
    "peptoclostridium": [("ENZ001", "weak_prior"), ("ENZ011", "weak_prior"), ("ENZ012", "weak_prior")],
    "prevotella": [("ENZ002", "genus_prior"), ("ENZ007", "genus_prior"), ("ENZ008", "weak_prior"), ("ENZ009", "genus_prior"), ("ENZ014", "weak_prior"), ("ENZ016", "genus_prior")],
    "proteus": [("ENZ003", "genus_prior"), ("ENZ004", "likely_present"), ("ENZ006", "weak_prior")],
    "pseudoflavonifractor": [("ENZ011", "likely_present"), ("ENZ012", "likely_present"), ("ENZ013", "weak_prior")],
    "roseburia": [("ENZ011", "genus_prior"), ("ENZ012", "weak_prior"), ("ENZ013", "weak_prior")],
    "ruminococcus": [("ENZ007", "genus_prior"), ("ENZ009", "genus_prior"), ("ENZ013", "weak_prior"), ("ENZ016", "likely_present")],
    "streptococcus": [("ENZ007", "weak_prior"), ("ENZ009", "likely_present"), ("ENZ010", "genus_prior")],
    "veillonella": [("ENZ004", "weak_prior")],
}


def _starter_enzyme_catalog() -> pd.DataFrame:
    return pd.DataFrame(ENZYME_CATALOG_SEED)


def _split_multi_value(value: object) -> list[str]:
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    parts = re.split(r"[;|,]", text)
    return [part.strip() for part in parts if part.strip()]


def _infer_genus(row: pd.Series) -> str:
    for column in ["genus", "microbe_label", "species_label", "species_name"]:
        value = normalize_whitespace(row.get(column))
        if not value:
            continue
        if column == "genus":
            return value
        return value.split()[0]
    return ""


def _normalize_presence_call(value: object) -> str:
    key = canonicalize_key(value)
    if key in {"curatedpresent", "present", "highconfidence"}:
        return "curated_present"
    if key in {"likelypresent", "likely"}:
        return "likely_present"
    if key in {"genusprior", "prior"}:
        return "genus_prior"
    if key in {"weakprior", "weak"}:
        return "weak_prior"
    if key in {"absent", "notpresent"}:
        return "absent"
    return "unknown"


def _ensure_presence_weight(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    if "presence_call" not in result.columns:
        result["presence_call"] = "unknown"
    result["presence_call"] = result["presence_call"].map(_normalize_presence_call)
    if "presence_weight" not in result.columns:
        result["presence_weight"] = result["presence_call"].map(PRESENCE_WEIGHT_MAP)
    else:
        result["presence_weight"] = pd.to_numeric(result["presence_weight"], errors="coerce")
        missing_mask = result["presence_weight"].isna()
        result.loc[missing_mask, "presence_weight"] = result.loc[missing_mask, "presence_call"].map(PRESENCE_WEIGHT_MAP)
    return result


def _starter_long_table(microbes: pd.DataFrame, catalog: pd.DataFrame) -> pd.DataFrame:
    long_rows: list[dict[str, object]] = []
    for _, row in microbes.iterrows():
        genus = _infer_genus(row)
        genus_key = canonicalize_key(genus)
        priors = GENUS_ENZYME_PRIORS.get(genus_key, [])
        for enzyme_id, presence_call in priors:
            enzyme_record = catalog.loc[catalog["enzyme_id"].eq(enzyme_id)].iloc[0]
            long_rows.append(
                {
                    "nt_code": row.get("nt_code"),
                    "microbe_label": row.get("microbe_label"),
                    "species_label": row.get("species_label"),
                    "species_name": row.get("species_name"),
                    "species": row.get("species"),
                    "strain": row.get("strain"),
                    "genus": genus or row.get("genus"),
                    "family": row.get("family"),
                    "phylum": row.get("phylum"),
                    "enzyme_id": enzyme_id,
                    "enzyme_name": enzyme_record["enzyme_name"],
                    "presence_call": presence_call,
                    "presence_weight": PRESENCE_WEIGHT_MAP[presence_call],
                    "evidence_scope": "genus_prior_seed",
                    "evidence_source": "project_starter_prior",
                    "source_database": "",
                    "literature_citation": "",
                    "pmid": "",
                    "doi": "",
                    "genome_accession": "",
                    "strain_match_level": "same_genus",
                    "evidence_note": (
                        f"Starter prior seeded from genus-level gut microbiome enzyme knowledge "
                        f"for {genus or row.get('species_label')}."
                    ),
                    "curation_status": "starter_seed",
                }
            )
    long_table = pd.DataFrame(long_rows)
    if long_table.empty:
        return _ensure_presence_weight(long_table)
    return _ensure_presence_weight(long_table).sort_values(["nt_code", "enzyme_id"]).reset_index(drop=True)


def _build_curation_template_frame(
    microbes: pd.DataFrame,
    catalog: pd.DataFrame,
    starter_long: pd.DataFrame,
    resolved_evidence: pd.DataFrame | None = None,
) -> pd.DataFrame:
    starter_lookup: dict[tuple[str, str], dict[str, object]] = {}
    if not starter_long.empty:
        starter_lookup = {
            (str(row["nt_code"]), str(row["enzyme_id"])): row.to_dict()
            for _, row in starter_long.iterrows()
        }
    resolved_lookup: dict[tuple[str, str], dict[str, object]] = {}
    if resolved_evidence is not None and not resolved_evidence.empty:
        resolved_lookup = {
            (str(row["nt_code"]), str(row["enzyme_id"])): row.to_dict()
            for _, row in resolved_evidence.iterrows()
        }

    rows: list[dict[str, object]] = []
    for _, microbe in microbes.iterrows():
        for _, enzyme in catalog.iterrows():
            key = (str(microbe.get("nt_code") or ""), str(enzyme.get("enzyme_id") or ""))
            starter_row = starter_lookup.get(key, {})
            resolved_row = resolved_lookup.get(key, {})
            rows.append(
                {
                    "nt_code": microbe.get("nt_code"),
                    "microbe_label": microbe.get("microbe_label"),
                    "species_label": microbe.get("species_label"),
                    "species_name": microbe.get("species_name"),
                    "species": microbe.get("species"),
                    "strain": microbe.get("strain"),
                    "genus": microbe.get("genus"),
                    "family": microbe.get("family"),
                    "phylum": microbe.get("phylum"),
                    "enzyme_id": enzyme.get("enzyme_id"),
                    "enzyme_name": enzyme.get("enzyme_name"),
                    "enzyme_family": enzyme.get("enzyme_family"),
                    "reaction_class": enzyme.get("reaction_class"),
                    "bond_target": enzyme.get("bond_target"),
                    "starter_presence_call": starter_row.get("presence_call", ""),
                    "starter_presence_weight": starter_row.get("presence_weight", ""),
                    "starter_evidence_scope": starter_row.get("evidence_scope", ""),
                    "starter_evidence_source": starter_row.get("evidence_source", ""),
                    "starter_evidence_note": starter_row.get("evidence_note", ""),
                    "curated_presence_call": resolved_row.get("presence_call", ""),
                    "curated_presence_weight": resolved_row.get("presence_weight", ""),
                    "curated_evidence_scope": resolved_row.get("evidence_scope", ""),
                    "curated_evidence_source": resolved_row.get("evidence_source", ""),
                    "curated_source_database": resolved_row.get("source_database", ""),
                    "curated_literature_citation": resolved_row.get("literature_citation", ""),
                    "curated_pmid": resolved_row.get("pmid", ""),
                    "curated_doi": resolved_row.get("doi", ""),
                    "curated_genome_accession": resolved_row.get("genome_accession", ""),
                    "curated_strain_match_level": resolved_row.get("strain_match_level", ""),
                    "curated_evidence_note": resolved_row.get("evidence_note", ""),
                    "curated_curation_status": resolved_row.get("curation_status", "needs_curation"),
                    "review_priority": (
                        "review"
                        if resolved_row
                        else ("high" if starter_row else "medium")
                    ),
                }
            )
    return pd.DataFrame(rows).sort_values(["nt_code", "enzyme_id"]).reset_index(drop=True)


def build_step2_enzyme_curation_template(
    microbe_table_path: str | Path | pd.DataFrame,
    output_path: str | Path = DEFAULT_MICROBE_ENZYME_CURATION_TEMPLATE_PATH,
) -> dict[str, object]:
    """Build an editable species/strain-level curation worksheet for enzyme evidence."""
    if isinstance(microbe_table_path, pd.DataFrame):
        microbes = microbe_table_path.drop_duplicates(subset=["nt_code"]).reset_index(drop=True).copy()
        microbe_table_repr = "<dataframe>"
    else:
        microbes = pd.read_csv(microbe_table_path, low_memory=False).drop_duplicates(subset=["nt_code"]).reset_index(drop=True)
        microbe_table_repr = str(microbe_table_path)

    catalog = _starter_enzyme_catalog()
    starter_long = _starter_long_table(microbes, catalog)
    template = _build_curation_template_frame(microbes, catalog, starter_long)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    template.to_csv(output_path, index=False)
    return {
        "microbe_table_path": microbe_table_repr,
        "output_path": str(output_path),
        "n_microbes": int(microbes["nt_code"].nunique()),
        "n_enzymes": int(len(catalog)),
        "n_template_rows": int(len(template)),
    }


def _normalize_literature_evidence_table(
    evidence: pd.DataFrame,
    microbes: pd.DataFrame,
    catalog: pd.DataFrame,
) -> pd.DataFrame:
    if evidence.empty:
        return pd.DataFrame()

    result = evidence.copy()
    curated_columns = {
        "presence_call": "curated_presence_call",
        "presence_weight": "curated_presence_weight",
        "evidence_scope": "curated_evidence_scope",
        "evidence_source": "curated_evidence_source",
        "source_database": "curated_source_database",
        "literature_citation": "curated_literature_citation",
        "pmid": "curated_pmid",
        "doi": "curated_doi",
        "genome_accession": "curated_genome_accession",
        "strain_match_level": "curated_strain_match_level",
        "evidence_note": "curated_evidence_note",
        "curation_status": "curated_curation_status",
    }
    if any(column in result.columns for column in curated_columns.values()):
        for target, source in curated_columns.items():
            if source in result.columns:
                result[target] = result[source]

    required_defaults = {
        "nt_code": "",
        "microbe_label": "",
        "species_label": "",
        "species_name": "",
        "species": "",
        "strain": "",
        "genus": "",
        "family": "",
        "phylum": "",
        "enzyme_id": "",
        "enzyme_name": "",
        "presence_call": "",
        "presence_weight": np.nan,
        "evidence_scope": "",
        "evidence_source": "",
        "source_database": "",
        "literature_citation": "",
        "pmid": "",
        "doi": "",
        "genome_accession": "",
        "strain_match_level": "",
        "evidence_note": "",
        "curation_status": "",
    }
    for column, default in required_defaults.items():
        if column not in result.columns:
            result[column] = default

    microbes_lookup = microbes.drop_duplicates(subset=["nt_code"]).copy()
    for column in ["microbe_label", "species_label", "species_name", "species", "strain", "genus", "family", "phylum"]:
        if column not in microbes_lookup.columns:
            microbes_lookup[column] = ""
    result = result.merge(
        microbes_lookup[
            ["nt_code", "microbe_label", "species_label", "species_name", "species", "strain", "genus", "family", "phylum"]
        ],
        on="nt_code",
        how="left",
        suffixes=("", "_microbe"),
    )
    for column in ["microbe_label", "species_label", "species_name", "species", "strain", "genus", "family", "phylum"]:
        microbe_column = f"{column}_microbe"
        if microbe_column in result.columns:
            missing_mask = result[column].isna() | result[column].astype(str).str.strip().eq("")
            result.loc[missing_mask, column] = result.loc[missing_mask, microbe_column]
            result = result.drop(columns=[microbe_column])

    catalog_lookup = catalog.copy()
    result = result.merge(
        catalog_lookup[["enzyme_id", "enzyme_name"]],
        on="enzyme_id",
        how="left",
        suffixes=("", "_catalog"),
    )
    if "enzyme_name_catalog" in result.columns:
        missing_mask = result["enzyme_name"].isna() | result["enzyme_name"].astype(str).str.strip().eq("")
        result.loc[missing_mask, "enzyme_name"] = result.loc[missing_mask, "enzyme_name_catalog"]
        result = result.drop(columns=["enzyme_name_catalog"])

    result["nt_code"] = result["nt_code"].map(normalize_whitespace)
    result["enzyme_id"] = result["enzyme_id"].map(normalize_whitespace)
    result["presence_call"] = result["presence_call"].fillna("").astype(str)
    has_content = (
        result["presence_call"].str.strip().ne("")
        | result["evidence_scope"].fillna("").astype(str).str.strip().ne("")
        | result["evidence_source"].fillna("").astype(str).str.strip().ne("")
        | result["source_database"].fillna("").astype(str).str.strip().ne("")
        | result["literature_citation"].fillna("").astype(str).str.strip().ne("")
        | result["pmid"].fillna("").astype(str).str.strip().ne("")
        | result["doi"].fillna("").astype(str).str.strip().ne("")
        | result["genome_accession"].fillna("").astype(str).str.strip().ne("")
        | result["evidence_note"].fillna("").astype(str).str.strip().ne("")
    )
    result = result.loc[result["nt_code"].ne("") & result["enzyme_id"].ne("") & has_content].copy()
    if result.empty:
        return result

    result["evidence_scope"] = result["evidence_scope"].replace("", np.nan)
    missing_scope_mask = result["evidence_scope"].isna()
    literature_mask = (
        result["pmid"].fillna("").astype(str).str.strip().ne("")
        | result["doi"].fillna("").astype(str).str.strip().ne("")
        | result["literature_citation"].fillna("").astype(str).str.strip().ne("")
    )
    result.loc[missing_scope_mask & literature_mask, "evidence_scope"] = "species_literature"
    result.loc[missing_scope_mask & ~literature_mask, "evidence_scope"] = "manual_curation"
    result["evidence_source"] = result["evidence_source"].replace("", np.nan)
    result["evidence_source"] = result["evidence_source"].fillna(result["source_database"].replace("", np.nan))
    result["evidence_source"] = result["evidence_source"].fillna("manual_curation")
    result["curation_status"] = result["curation_status"].replace("", np.nan).fillna("curated_literature")
    result["strain_match_level"] = result["strain_match_level"].replace("", np.nan).fillna("same_species")
    return _ensure_presence_weight(result)


def _rank_evidence_row(row: pd.Series) -> float:
    scope_score = EVIDENCE_SCOPE_PRIORITY.get(canonicalize_key(row.get("evidence_scope")), 0.5)
    status_score = CURATION_STATUS_PRIORITY.get(canonicalize_key(row.get("curation_status")), 0.0)
    strain_score = STRAIN_MATCH_PRIORITY.get(canonicalize_key(row.get("strain_match_level")), 0.0)
    presence_weight = pd.to_numeric(pd.Series([row.get("presence_weight")]), errors="coerce").iloc[0]
    presence_score = 0.0 if pd.isna(presence_weight) else float(presence_weight) * 0.1
    citation_score = 0.15 if normalize_whitespace(row.get("pmid") or row.get("doi") or row.get("literature_citation")) else 0.0
    return float(scope_score + status_score + strain_score + presence_score + citation_score)


def _resolve_microbe_enzyme_evidence(evidence_ledger: pd.DataFrame) -> pd.DataFrame:
    if evidence_ledger.empty:
        return evidence_ledger.copy()
    ranked = _ensure_presence_weight(evidence_ledger)
    ranked["evidence_rank_score"] = ranked.apply(_rank_evidence_row, axis=1)
    ranked = ranked.sort_values(
        ["nt_code", "enzyme_id", "evidence_rank_score", "presence_weight"],
        ascending=[True, True, False, False],
        na_position="last",
    )
    return ranked.drop_duplicates(subset=["nt_code", "enzyme_id"], keep="first").reset_index(drop=True)


def build_step2_enzyme_reference_tables(
    microbe_table_path: str | Path | pd.DataFrame,
    enzyme_catalog_path: str | Path = DEFAULT_ENZYME_CATALOG_PATH,
    microbe_enzyme_long_path: str | Path = DEFAULT_MICROBE_ENZYME_PANEL_PATH,
    microbe_enzyme_matrix_path: str | Path = DEFAULT_MICROBE_ENZYME_MATRIX_PATH,
    evidence_ledger_path: str | Path | None = DEFAULT_MICROBE_ENZYME_EVIDENCE_LEDGER_PATH,
    curation_template_path: str | Path | None = DEFAULT_MICROBE_ENZYME_CURATION_TEMPLATE_PATH,
    literature_evidence_path: str | Path | pd.DataFrame | None = None,
    summary_path: str | Path | None = DEFAULT_ENZYME_SUMMARY_PATH,
) -> dict[str, object]:
    """Build enzyme-reference tables for the 83-microbe panel with optional species-level overrides."""
    if isinstance(microbe_table_path, pd.DataFrame):
        microbe_table = microbe_table_path.copy()
        microbe_table_repr = "<dataframe>"
    else:
        microbe_table = pd.read_csv(microbe_table_path, low_memory=False)
        microbe_table_repr = str(microbe_table_path)

    catalog = _starter_enzyme_catalog()
    microbes = microbe_table.drop_duplicates(subset=["nt_code"]).reset_index(drop=True).copy()
    starter_long = _starter_long_table(microbes, catalog)

    curated_evidence = _normalize_literature_evidence_table(
        _load_catalog(literature_evidence_path),
        microbes=microbes,
        catalog=catalog,
    )
    evidence_ledger = pd.concat([starter_long, curated_evidence], ignore_index=True, sort=False)
    evidence_ledger = _ensure_presence_weight(evidence_ledger)
    long_table = _resolve_microbe_enzyme_evidence(evidence_ledger)

    matrix_input = long_table.copy()
    for column in ["nt_code", "microbe_label", "species_label", "species_name", "genus", "family", "phylum"]:
        if column in matrix_input.columns:
            matrix_input[column] = matrix_input[column].fillna("")
    matrix_source = matrix_input.pivot_table(
        index=["nt_code", "microbe_label", "species_label", "species_name", "genus", "family", "phylum"],
        columns="enzyme_id",
        values="presence_weight",
        aggfunc="max",
    ).reset_index()
    if not matrix_source.empty:
        matrix_source.columns = [
            column if isinstance(column, str) else column[0]
            for column in matrix_source.columns
        ]

    enzyme_catalog_path = Path(enzyme_catalog_path)
    microbe_enzyme_long_path = Path(microbe_enzyme_long_path)
    microbe_enzyme_matrix_path = Path(microbe_enzyme_matrix_path)
    output_paths = [enzyme_catalog_path, microbe_enzyme_long_path, microbe_enzyme_matrix_path]
    if evidence_ledger_path is not None:
        output_paths.append(Path(evidence_ledger_path))
    if curation_template_path is not None:
        output_paths.append(Path(curation_template_path))
    for path in output_paths:
        path.parent.mkdir(parents=True, exist_ok=True)

    catalog.to_csv(enzyme_catalog_path, index=False)
    long_table.to_csv(microbe_enzyme_long_path, index=False)
    matrix_source.to_csv(microbe_enzyme_matrix_path, index=False)
    if evidence_ledger_path is not None:
        Path(evidence_ledger_path).parent.mkdir(parents=True, exist_ok=True)
        evidence_ledger.sort_values(["nt_code", "enzyme_id"]).to_csv(evidence_ledger_path, index=False)
    if curation_template_path is not None:
        template = _build_curation_template_frame(microbes, catalog, starter_long, resolved_evidence=long_table)
        Path(curation_template_path).parent.mkdir(parents=True, exist_ok=True)
        template.to_csv(curation_template_path, index=False)

    summary = {
        "microbe_table_path": microbe_table_repr,
        "enzyme_catalog_path": str(enzyme_catalog_path),
        "microbe_enzyme_long_path": str(microbe_enzyme_long_path),
        "microbe_enzyme_matrix_path": str(microbe_enzyme_matrix_path),
        "evidence_ledger_path": str(evidence_ledger_path) if evidence_ledger_path is not None else None,
        "curation_template_path": str(curation_template_path) if curation_template_path is not None else None,
        "literature_evidence_path": "<dataframe>"
        if isinstance(literature_evidence_path, pd.DataFrame)
        else (str(literature_evidence_path) if literature_evidence_path is not None else None),
        "n_microbes": int(microbes["nt_code"].nunique()),
        "n_seeded_microbe_enzyme_rows": int(len(starter_long)),
        "n_curated_evidence_rows": int(len(curated_evidence)),
        "n_evidence_ledger_rows": int(len(evidence_ledger)),
        "n_resolved_microbe_enzyme_rows": int(len(long_table)),
        "n_enzyme_families": int(len(catalog)),
        "seeded_genera": sorted(set(starter_long["genus"].dropna().astype(str))) if not starter_long.empty else [],
        "resolved_present_rows": int(
            long_table.get("presence_weight", pd.Series(dtype=float)).fillna(0).gt(0).sum()
        ),
    }
    if summary_path is not None:
        summary_path = Path(summary_path)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
        summary["summary_path"] = str(summary_path)
    return summary


def _load_catalog(path: str | Path | pd.DataFrame | None) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame()
    if isinstance(path, pd.DataFrame):
        return path.copy()
    candidate = Path(path)
    if not candidate.exists():
        return pd.DataFrame()
    return pd.read_csv(candidate, low_memory=False)


def _build_panel_lookup(panel: pd.DataFrame, catalog: pd.DataFrame) -> dict[str, list[dict[str, object]]]:
    if panel.empty or catalog.empty:
        return {}
    panel = _ensure_presence_weight(panel)
    merged = panel.merge(catalog, on=["enzyme_id", "enzyme_name"], how="left")
    lookup: dict[str, list[dict[str, object]]] = {}
    for _, row in merged.iterrows():
        nt_code = normalize_whitespace(row.get("nt_code"))
        if not nt_code:
            continue
        lookup.setdefault(nt_code, []).append(row.to_dict())
    return lookup


def _normalized_text_blob(row: pd.Series) -> str:
    fields = [
        row.get("chemical_name"),
        row.get("therapeutic_class"),
        row.get("therapeutic_effect"),
        row.get("compound_semantic_aliases"),
        row.get("compound_semantic_keywords"),
        row.get("compound_name_normalized"),
    ]
    return " ".join(normalize_whitespace(value).lower() for value in fields if normalize_whitespace(value))


def _reaction_key_from_row(row: pd.Series) -> str:
    for column in ["predicted_reaction_class", "step2_reaction_class", "reaction_class"]:
        key = canonicalize_key(row.get(column))
        if key:
            return key
    return ""


def annotate_step2_with_enzyme_priors(
    frame: pd.DataFrame,
    microbe_enzyme_panel_path: str | Path | pd.DataFrame | None = DEFAULT_MICROBE_ENZYME_PANEL_PATH,
    enzyme_catalog_path: str | Path | pd.DataFrame | None = DEFAULT_ENZYME_CATALOG_PATH,
) -> pd.DataFrame:
    """Annotate Step 2 candidate or prediction rows with enzyme-prior support features."""
    result = frame.copy()
    defaults = {
        "predicted_enzyme_prior_flag": False,
        "predicted_enzyme_match_count": 0,
        "predicted_enzyme_ids": "",
        "predicted_enzyme_names": "",
        "predicted_enzyme_reaction_classes": "",
        "predicted_enzyme_bond_targets": "",
        "predicted_enzyme_presence_score": np.nan,
        "predicted_enzyme_support_score": np.nan,
        "predicted_enzyme_step1_promote_support_score": np.nan,
        "predicted_enzyme_step1_inhibit_risk_score": np.nan,
    }
    for column, default in defaults.items():
        result[column] = default

    panel = _load_catalog(microbe_enzyme_panel_path)
    catalog = _load_catalog(enzyme_catalog_path)
    if result.empty or panel.empty or catalog.empty:
        return result

    if "compound_semantic_family" not in result.columns:
        result = annotate_compound_semantics(result)

    panel_lookup = _build_panel_lookup(panel, catalog)
    if not panel_lookup:
        return result

    for row_index, row in result.iterrows():
        nt_code = normalize_whitespace(row.get("nt_code"))
        enzyme_records = panel_lookup.get(nt_code, [])
        if not enzyme_records:
            continue

        family_key = canonicalize_key(row.get("compound_semantic_family"))
        text_blob = _normalized_text_blob(row)
        reaction_key = _reaction_key_from_row(row)

        matched_entries: list[dict[str, object]] = []
        for enzyme in enzyme_records:
            presence_weight = pd.to_numeric(pd.Series([enzyme.get("presence_weight")]), errors="coerce").iloc[0]
            if pd.isna(presence_weight) or float(presence_weight) <= 0:
                continue

            family_tokens = {canonicalize_key(value) for value in _split_multi_value(enzyme.get("compound_semantic_families"))}
            keyword_tokens = {normalize_whitespace(value).lower() for value in _split_multi_value(enzyme.get("substrate_keywords"))}
            enzyme_reaction_key = canonicalize_key(enzyme.get("reaction_class"))

            family_match = bool(family_key and family_key in family_tokens)
            keyword_match = bool(keyword_tokens and any(token and token in text_blob for token in keyword_tokens))
            reaction_match = bool(reaction_key and enzyme_reaction_key and reaction_key == enzyme_reaction_key)

            match_strength = min(
                1.0,
                (0.55 if family_match else 0.0)
                + (0.35 if keyword_match else 0.0)
                + (0.20 if reaction_match else 0.0),
            )
            if match_strength <= 0:
                continue

            metabolism_weight = float(pd.to_numeric(pd.Series([enzyme.get("metabolism_weight")]), errors="coerce").iloc[0] or 0.0)
            promote_weight = float(pd.to_numeric(pd.Series([enzyme.get("step1_promote_weight")]), errors="coerce").iloc[0] or 0.0)
            inhibit_weight = float(pd.to_numeric(pd.Series([enzyme.get("step1_inhibit_weight")]), errors="coerce").iloc[0] or 0.0)
            matched_entries.append(
                {
                    "enzyme_id": str(enzyme.get("enzyme_id") or ""),
                    "enzyme_name": str(enzyme.get("enzyme_name") or ""),
                    "reaction_class": str(enzyme.get("reaction_class") or ""),
                    "bond_target": str(enzyme.get("bond_target") or ""),
                    "presence_weight": float(presence_weight),
                    "support": float(presence_weight) * metabolism_weight * match_strength,
                    "promote": float(presence_weight) * promote_weight * match_strength,
                    "inhibit": float(presence_weight) * inhibit_weight * match_strength,
                }
            )

        if not matched_entries:
            continue

        matched_entries = sorted(matched_entries, key=lambda item: item["support"], reverse=True)
        result.at[row_index, "predicted_enzyme_prior_flag"] = True
        result.at[row_index, "predicted_enzyme_match_count"] = int(len(matched_entries))
        result.at[row_index, "predicted_enzyme_ids"] = ";".join(dict.fromkeys(item["enzyme_id"] for item in matched_entries if item["enzyme_id"]))
        result.at[row_index, "predicted_enzyme_names"] = ";".join(dict.fromkeys(item["enzyme_name"] for item in matched_entries if item["enzyme_name"]))
        result.at[row_index, "predicted_enzyme_reaction_classes"] = ";".join(
            dict.fromkeys(item["reaction_class"] for item in matched_entries if item["reaction_class"])
        )
        result.at[row_index, "predicted_enzyme_bond_targets"] = ";".join(
            dict.fromkeys(item["bond_target"] for item in matched_entries if item["bond_target"])
        )
        result.at[row_index, "predicted_enzyme_presence_score"] = float(
            max(item["presence_weight"] for item in matched_entries)
        )
        result.at[row_index, "predicted_enzyme_support_score"] = float(
            min(1.0, sum(item["support"] for item in matched_entries[:3]))
        )
        result.at[row_index, "predicted_enzyme_step1_promote_support_score"] = float(
            min(1.0, sum(item["promote"] for item in matched_entries[:3]))
        )
        result.at[row_index, "predicted_enzyme_step1_inhibit_risk_score"] = float(
            min(1.0, sum(item["inhibit"] for item in matched_entries[:3]))
        )
    return result
