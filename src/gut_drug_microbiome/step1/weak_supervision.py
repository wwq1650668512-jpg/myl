from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .chem_features import enrich_drug_table_with_rdkit


MDIPID_DOWNLOADS = {
    "drug_info.csv": "https://mdipid.idrblab.net/sites/default/files/download/2.General%20Information%20of%20Drug%20and%20Substance.csv",
    "microbe_info.csv": "https://mdipid.idrblab.net/sites/default/files/download/1.General%20Information%20of%20Microbiota.csv",
    "deim.csv": "https://mdipid.idrblab.net/sites/default/files/download/6.Data%20of%20Drug%20or%20other%20exogenous%20substances%20impact%20on%20microbiota%20%28DEIM%29.csv",
}

MASI_INFERRED_DOWNLOADS = {
    "MASI_v2.0_download_microbeInfo.xlsx": "https://www.aiddlab.com/MASI2025/download/MASI_v2.0_download_microbeInfo.xlsx",
    "MASI_v2.0_download_substanceInfo.xlsx": "https://www.aiddlab.com/MASI2025/download/MASI_v2.0_download_substanceInfo.xlsx",
    "MASI_v2.0_download_microbeSubstanceInteractionRecords.xlsx": "https://www.aiddlab.com/MASI2025/download/MASI_v2.0_download_microbeSubstanceInteractionRecords.xlsx",
    "MASI_v2.0_download_DrugReactionRecords.xlsx": "https://www.aiddlab.com/MASI2025/download/MASI_v2.0_download_DrugReactionRecords.xlsx",
    "MASI_v2.0_download_DietReactionRecords.xlsx": "https://www.aiddlab.com/MASI2025/download/MASI_v2.0_download_DietReactionRecords.xlsx",
}

DRUG_LIKE_TOKENS = (
    "drug",
    "drug-like",
    "pharmaceutical",
    "medicine",
    "medication",
    "therapeutic substance",
    "antibiotic",
    "antidepressant",
    "antipsychotic",
    "statin",
    "nsaid",
    "anti-inflammatory",
    "analgesic",
    "antiviral",
    "antifungal",
    "antidiabetic",
    "proton pump inhibitor",
)


def _snake_case(text: str) -> str:
    """Normalize external weak-supervision column names into snake_case."""
    value = str(text).strip().lower()
    value = value.replace("µ", "u").replace("μ", "u").replace("/", "_").replace("-", "_").replace(".", "_")
    value = re.sub(r"[^a-z0-9]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    return value or "unnamed"


def _download_file(url: str, destination: Path) -> None:
    """Download one weak-supervision source file with curl."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["curl", "-L", "--max-time", "120", "-o", str(destination), url],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def _write_manifest(path: Path, payload: dict) -> None:
    """Write a JSON manifest or summary payload to disk."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def _download_catalog(
    catalog: dict[str, str],
    raw_dir: Path,
    overwrite: bool = False,
    allow_partial: bool = False,
) -> dict[str, dict[str, str]]:
    """Download a catalog of URLs and return per-file status metadata."""
    manifest: dict[str, dict[str, str]] = {}
    failures: list[str] = []
    for filename, url in catalog.items():
        destination = raw_dir / filename
        status = "cached"
        error = ""
        if overwrite or not destination.exists():
            try:
                _download_file(url, destination)
                status = "downloaded"
            except subprocess.CalledProcessError as exc:
                status = "failed"
                error = str(exc)
                failures.append(filename)
        manifest[filename] = {
            "path": str(destination),
            "url": url,
            "status": status,
            "error": error,
        }
    if failures and not allow_partial:
        failed = ", ".join(sorted(failures))
        raise RuntimeError(f"Failed to download required files: {failed}")
    return manifest


def download_mdipid_data(raw_dir: str | Path, overwrite: bool = False) -> dict[str, str]:
    """Download the MDIPID files used for Step 1 silver-label construction."""
    raw_dir = Path(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)

    manifest: dict[str, str] = {}
    for filename, url in MDIPID_DOWNLOADS.items():
        destination = raw_dir / filename
        if overwrite or not destination.exists():
            _download_file(url, destination)
        manifest[filename] = str(destination)

    manifest["masi_inferred_downloads_note"] = (
        "MASI v2.0 download URLs are inferred from the public download page filenames "
        "and may require manual retrieval if direct shell access to aiddlab.com is unavailable."
    )
    _write_manifest(
        raw_dir / "download_manifest.json",
        {"mdipid": manifest, "masi_inferred_downloads": MASI_INFERRED_DOWNLOADS},
    )
    return manifest


def download_masi_data(
    raw_dir: str | Path,
    overwrite: bool = False,
    allow_partial: bool = True,
) -> dict[str, object]:
    """Download inferred MASI files and record which files succeeded or failed."""
    raw_dir = Path(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)

    manifest = _download_catalog(MASI_INFERRED_DOWNLOADS, raw_dir, overwrite=overwrite, allow_partial=allow_partial)
    summary = {
        "dataset": "masi_v2_inferred",
        "note": (
            "MASI URLs are inferred from the public MASI v2.0 download page. "
            "If downloads fail, manually place the original files in this directory and rerun normalization."
        ),
        "files": manifest,
    }
    _write_manifest(raw_dir / "download_manifest.json", summary)
    return summary


def _extract_pubchem_cid(value: str | float) -> str | float:
    """Extract the first numeric PubChem CID token from a raw text field."""
    if not isinstance(value, str):
        return np.nan
    match = re.search(r"(\d+)", value)
    if not match:
        return np.nan
    return match.group(1)


def _infer_effect_label(text: str | float) -> str | float:
    """Infer a coarse Step 1 label from free-text literature descriptions."""
    if not isinstance(text, str):
        return np.nan
    normalized = text.lower()
    neutral_patterns = [
        "no significant change",
        "not significantly changed",
        "unchanged",
        "no change",
    ]
    promote_patterns = [
        "increase",
        "increase the relative abundance",
        "increased the relative abundance",
        "increase abundance",
        "increased abundance",
        "elevate",
        "enrich",
        "promote",
        "higher abundance",
        "positive association",
    ]
    inhibit_patterns = [
        "decrease",
        "decrease the relative abundance",
        "decreased the relative abundance",
        "decrease abundance",
        "decreased abundance",
        "reduce the relative abundance",
        "reduced the relative abundance",
        "reduce abundance",
        "reduced abundance",
        "lower abundance",
        "deplete",
        "suppress",
        "negative association",
        "inhibit",
    ]
    if any(pattern in normalized for pattern in neutral_patterns):
        return "no_effect"
    if any(pattern in normalized for pattern in promote_patterns):
        return "promote"
    if any(pattern in normalized for pattern in inhibit_patterns):
        return "inhibit"
    return np.nan


def _coalesce_series(frame: pd.DataFrame, candidates: Iterable[str]) -> pd.Series:
    """Combine multiple optional columns into one series by first non-null value."""
    result = pd.Series(np.nan, index=frame.index, dtype=object)
    for column in candidates:
        if column not in frame.columns:
            continue
        values = frame[column]
        mask = result.isna() & values.notna()
        result.loc[mask] = values.loc[mask]
    return result


def _read_table_auto(path: Path) -> pd.DataFrame:
    """Read a CSV, TSV, or Excel weak-supervision table based on file suffix."""
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep=None, engine="python")
    if suffix == ".csv":
        return pd.read_csv(path)
    raise ValueError(f"Unsupported table format: {path}")


def _find_best_file(raw_dir: Path, name_fragment: str) -> Path | None:
    """Find the best matching local file for a required MASI artifact name fragment."""
    candidates = sorted(
        [
            path
            for path in raw_dir.iterdir()
            if path.is_file() and name_fragment.lower() in path.name.lower()
        ],
        key=lambda path: (path.suffix.lower() not in {".xlsx", ".xls"}, path.name.lower()),
    )
    return candidates[0] if candidates else None


def _normalize_text_series(series: pd.Series) -> pd.Series:
    """Convert a series to stripped strings with empty-text fallbacks."""
    return series.fillna("").astype(str).str.strip()


def _infer_drug_like_scope(frame: pd.DataFrame) -> pd.Series:
    """Heuristically classify MASI substances as drug-like, non-drug-like, or unknown."""
    category = _normalize_text_series(
        _coalesce_series(frame, ["substance_category", "category", "active_substance_category"])
    ).str.lower()
    subcategory = _normalize_text_series(
        _coalesce_series(frame, ["substance_subcategory", "subcategory", "active_substance_subcategory"])
    ).str.lower()
    details = _normalize_text_series(
        _coalesce_series(frame, ["substance_details", "substance_detail", "details", "drug_detail"])
    ).str.lower()
    substance_name = _normalize_text_series(
        _coalesce_series(frame, ["substance_name", "compound_name", "drug_name", "substance"])
    ).str.lower()

    joined = category + " | " + subcategory + " | " + details + " | " + substance_name
    is_drug_like = joined.apply(lambda value: any(token in value for token in DRUG_LIKE_TOKENS))
    has_scope_signal = joined.str.len() > 6
    scope = pd.Series("unknown", index=frame.index, dtype=object)
    scope.loc[has_scope_signal & ~is_drug_like] = "non_drug_like"
    scope.loc[is_drug_like] = "drug_like"
    return scope


def build_mdipid_silver_table(
    raw_dir: str | Path,
    output_dir: str | Path,
    gut_only: bool = True,
) -> dict:
    """Normalize MDIPID records into a Step 1 silver-label modeling table.

    Args:
        raw_dir: Directory containing downloaded MDIPID source files.
        output_dir: Directory where the normalized silver table and summary are written.
        gut_only: Whether to keep only gut microbiota records.

    Returns:
        A summary dictionary with record counts and output paths.
    """
    raw_dir = Path(raw_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    deim = pd.read_csv(raw_dir / "deim.csv")
    deim = deim.rename(columns={column: _snake_case(column) for column in deim.columns})
    drug_info = pd.read_csv(raw_dir / "drug_info.csv")
    drug_info = drug_info.rename(columns={column: _snake_case(column) for column in drug_info.columns})
    microbe_info = pd.read_csv(raw_dir / "microbe_info.csv")
    microbe_info = microbe_info.rename(columns={column: _snake_case(column) for column in microbe_info.columns})

    silver = deim.merge(drug_info, left_on="drugid", right_on="drugid", how="left", suffixes=("", "_drug"))
    silver = silver.merge(microbe_info, on="mic_id", how="left", suffixes=("", "_microbe"))

    if gut_only and "microbiota_site" in silver.columns:
        silver = silver[silver["microbiota_site"].fillna("").str.lower() == "gut"].copy()

    silver["effect_label"] = silver["results_detail"].apply(_infer_effect_label)
    silver = silver.dropna(subset=["effect_label"]).copy()
    silver["binary_effect_label"] = silver["effect_label"]
    silver["effect_score"] = np.nan
    silver["chemical_name"] = silver["drug_name"]
    silver["pubchem_cid"] = silver["pubchem_cid"].apply(_extract_pubchem_cid)
    silver["canonical_smiles"] = silver["canonical_smiles"].replace(".", np.nan)
    silver["iso_smiles"] = silver["iso_smiles"].replace(".", np.nan)
    silver["main_component_smiles"] = silver["canonical_smiles"].fillna(silver["iso_smiles"])
    silver["molecular_formula"] = silver.get("formula")
    silver["molecular_weight"] = pd.to_numeric(silver.get("mw"), errors="coerce")
    silver["tpsa"] = pd.to_numeric(silver.get("topological_polar_surface_area"), errors="coerce")
    silver["complexity"] = pd.to_numeric(silver.get("complexity"), errors="coerce")
    silver["xlogp"] = pd.to_numeric(silver.get("xlogp"), errors="coerce")
    silver["hbond_donor_count"] = pd.to_numeric(silver.get("hbonddonor"), errors="coerce")
    silver["hbond_acceptor_count"] = pd.to_numeric(silver.get("hbondacc"), errors="coerce")
    silver["rotatable_bond_count"] = pd.to_numeric(silver.get("rotbonds"), errors="coerce")
    silver["species_label"] = silver["speciesname"]
    silver["phylum"] = silver["phylum_name"]
    silver["class"] = silver["class_name"]
    silver["order"] = silver["order_name"]
    silver["family"] = silver["family_name"]
    silver["genus"] = silver["genus_name"]
    silver["target_species"] = "human"
    silver["human_use"] = True
    silver["veterinary"] = False
    silver["pair_id"] = silver["drugid"].astype(str) + "::" + silver["mic_id"].astype(str)
    silver["source_dataset"] = "mdipid_deim"
    silver["label_tier"] = "silver"
    silver["source_record_id"] = silver["pair_id"]
    silver["supporting_text"] = silver["results_detail"]

    silver = enrich_drug_table_with_rdkit(silver, smiles_columns=["main_component_smiles", "canonical_smiles", "iso_smiles"])

    columns_to_keep = [
        "pair_id",
        "source_dataset",
        "label_tier",
        "source_record_id",
        "effect_label",
        "binary_effect_label",
        "effect_score",
        "drugid",
        "chemical_name",
        "pubchem_cid",
        "main_component_smiles",
        "canonical_smiles",
        "iso_smiles",
        "molecular_formula",
        "molecular_weight",
        "xlogp",
        "tpsa",
        "complexity",
        "hbond_donor_count",
        "hbond_acceptor_count",
        "rotatable_bond_count",
        "species_label",
        "speciesname",
        "speciesid",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "microbiota_site",
        "drug_detail",
        "experimental_materials",
        "sample_source",
        "experiment_methods",
        "supporting_text",
    ]
    rdkit_columns = [column for column in silver.columns if column.startswith("rdkit_") or column.startswith("morgan_fp_")]
    existing_columns = [column for column in columns_to_keep if column in silver.columns]
    silver = silver[existing_columns + sorted(rdkit_columns)].copy()

    output_path = output_dir / "step1_silver_mdipid.csv"
    silver.to_csv(output_path, index=False)

    summary = {
        "raw_dir": str(raw_dir),
        "output_path": str(output_path),
        "gut_only": gut_only,
        "n_records": int(len(silver)),
        "n_drugs": int(silver["drugid"].nunique()),
        "n_microbes": int(silver["species_label"].nunique()),
        "label_counts": {key: int(value) for key, value in silver["effect_label"].value_counts().to_dict().items()},
    }
    _write_manifest(output_dir / "step1_silver_mdipid_summary.json", summary)
    return summary


def build_masi_silver_table(
    raw_dir: str | Path,
    output_dir: str | Path,
    drug_like_only: bool = False,
) -> dict:
    """Normalize MASI interactions into a Step 1 silver-label table.

    Args:
        raw_dir: Directory containing the MASI source workbooks.
        output_dir: Directory where the normalized silver table and summary are written.
        drug_like_only: Whether to keep only substance records inferred as drug-like.

    Returns:
        A summary dictionary with scope counts, label counts, and output paths.
    """
    raw_dir = Path(raw_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not raw_dir.exists():
        raise FileNotFoundError(
            f"MASI raw directory not found: {raw_dir}. "
            "Run scripts/download_step1_weak_supervision.py --dataset masi or place the original MASI files there first."
        )

    interaction_path = _find_best_file(raw_dir, "microbeSubstanceInteractionRecords")
    if interaction_path is None:
        raise FileNotFoundError(
            "No MASI interaction file found. Expected a file containing "
            "'microbeSubstanceInteractionRecords' in its filename."
        )

    microbe_info_path = _find_best_file(raw_dir, "microbeInfo")
    substance_info_path = _find_best_file(raw_dir, "substanceInfo")

    interactions = _read_table_auto(interaction_path)
    interactions = interactions.rename(columns={column: _snake_case(column) for column in interactions.columns})

    microbe_info = None
    if microbe_info_path is not None:
        microbe_info = _read_table_auto(microbe_info_path)
        microbe_info = microbe_info.rename(columns={column: _snake_case(column) for column in microbe_info.columns})

    substance_info = None
    if substance_info_path is not None:
        substance_info = _read_table_auto(substance_info_path)
        substance_info = substance_info.rename(columns={column: _snake_case(column) for column in substance_info.columns})

    silver = interactions.copy()

    if microbe_info is not None:
        microbe_key = next((column for column in ["microbe_id", "microbeid", "mid", "pmdbm_id"] if column in silver.columns), None)
        microbe_info_key = next(
            (column for column in ["microbe_id", "microbeid", "mid", "pmdbm_id"] if column in microbe_info.columns),
            None,
        )
        if microbe_key and microbe_info_key:
            silver = silver.merge(microbe_info, left_on=microbe_key, right_on=microbe_info_key, how="left", suffixes=("", "_microbe"))

    if substance_info is not None:
        substance_key = next(
            (column for column in ["substance_id", "substanceid", "drug_id", "compound_id", "suid"] if column in silver.columns),
            None,
        )
        substance_info_key = next(
            (column for column in ["substance_id", "substanceid", "drug_id", "compound_id", "suid"] if column in substance_info.columns),
            None,
        )
        if substance_key and substance_info_key:
            silver = silver.merge(
                substance_info,
                left_on=substance_key,
                right_on=substance_info_key,
                how="left",
                suffixes=("", "_substance"),
            )

    interaction_category = _normalize_text_series(_coalesce_series(silver, ["interaction_category"])).str.lower()
    microbe_change = _normalize_text_series(_coalesce_series(silver, ["microbe_change"])).str.lower()
    keep_mask = interaction_category.eq("substances alter microbe abundance")
    keep_mask = keep_mask | (
        interaction_category.eq("n_a")
        & microbe_change.isin(["increase", "decrease", "no significant change"])
    )
    silver = silver[keep_mask].copy()

    silver["microbe_id_std"] = _coalesce_series(silver, ["microbe_id", "microbeid", "mid", "pmdbm_id"])
    silver["substance_id_std"] = _coalesce_series(silver, ["substance_id", "substanceid", "drug_id", "compound_id", "suid"])
    silver["microbe_id_std"] = silver["microbe_id_std"].fillna(_coalesce_series(silver, ["masi_microbe_id"]))
    silver["substance_id_std"] = silver["substance_id_std"].fillna(_coalesce_series(silver, ["masi_substance_chemicald"]))
    silver["species_label"] = _coalesce_series(
        silver,
        ["species_name", "species_label", "microbe_name", "microbe", "microbiota_name"],
    )
    silver["chemical_name"] = _coalesce_series(
        silver,
        ["substance_name", "compound_name", "drug_name", "substance", "substance_name_substance"],
    )
    silver["substance_category"] = _coalesce_series(
        silver,
        ["substance_category", "category", "active_substance_category"],
    )
    silver["substance_subcategory"] = _coalesce_series(
        silver,
        ["substance_subcategory", "subcategory", "active_substance_subcategory"],
    )
    silver["substance_details"] = _coalesce_series(
        silver,
        ["substance_details", "substance_detail", "substance_exposure_details", "details", "drug_detail"],
    )
    silver["effect_on_microbe"] = _coalesce_series(
        silver,
        ["effect_on_microbe", "microbe_change", "effect", "alteration_type", "influence_on_microbe"],
    )
    silver["effect_strength"] = _coalesce_series(silver, ["effect_strength", "strength", "effect_size"])
    silver["experimental_system"] = _coalesce_series(
        silver,
        ["experimental_system", "system", "experiment_system"],
    )
    silver["experimental_organism"] = _coalesce_series(
        silver,
        ["experimental_organism", "organism", "host"],
    )
    silver["experimental_disease_condition"] = _coalesce_series(
        silver,
        ["experimental_disease_condition", "disease_condition", "condition"],
    )
    silver["reference_pubmed_id"] = _coalesce_series(
        silver,
        ["reference_pubmed_id", "pubmed_id", "pmid", "reference"],
    )
    silver["reference_pubmed_id"] = silver["reference_pubmed_id"].fillna(_coalesce_series(silver, ["reference_id"]))
    silver["canonical_smiles"] = _coalesce_series(
        silver,
        ["canonical_smiles", "canonical_smiles_substance", "smiles", "smiles_substance", "iso_smiles"],
    )
    silver["pubchem_cid"] = _coalesce_series(
        silver,
        ["pubchem_cid", "pubchem_id", "cid", "id_pubchem", "pubchem_cid_substance"],
    )
    silver["molecular_formula"] = _coalesce_series(
        silver,
        ["molecular_formula", "formula", "chemical_formula"],
    )
    silver["therapeutic_class"] = _coalesce_series(silver, ["therapeutic_class"])
    silver["therapeutic_effect"] = _coalesce_series(silver, ["substance_subcategory"])
    silver["microbiota_site"] = _coalesce_series(silver, ["microbiota_site"])
    silver["phylum"] = _coalesce_series(silver, ["phylum", "phylum_name"])
    silver["class"] = _coalesce_series(silver, ["class", "class_name"])
    silver["order"] = _coalesce_series(silver, ["order", "order_name"])
    silver["family"] = _coalesce_series(silver, ["family", "family_name"])
    silver["genus"] = _coalesce_series(silver, ["genus", "genus_name"])

    silver["effect_label"] = (
        _coalesce_series(silver, ["effect_on_microbe", "microbe_change_statistics", "effect_strength", "substance_details"]).apply(
            _infer_effect_label
        )
    )
    silver = silver.dropna(subset=["effect_label", "species_label", "chemical_name"]).copy()

    silver["substance_scope"] = _infer_drug_like_scope(silver)
    if drug_like_only:
        silver = silver[silver["substance_scope"] == "drug_like"].copy()

    silver["binary_effect_label"] = silver["effect_label"]
    silver["effect_score"] = np.nan
    silver["main_component_smiles"] = silver["canonical_smiles"]
    silver["molecular_weight"] = pd.to_numeric(_coalesce_series(silver, ["molecular_weight", "mw"]), errors="coerce")
    silver["xlogp"] = pd.to_numeric(_coalesce_series(silver, ["xlogp", "logp"]), errors="coerce")
    silver["tpsa"] = pd.to_numeric(_coalesce_series(silver, ["tpsa", "topological_polar_surface_area"]), errors="coerce")
    silver["complexity"] = pd.to_numeric(_coalesce_series(silver, ["complexity"]), errors="coerce")
    silver["target_species"] = _coalesce_series(silver, ["experimental_organism", "experiment_model_species"])
    silver["human_use"] = silver["substance_scope"].eq("drug_like")
    silver["veterinary"] = False
    silver["pair_id"] = silver["substance_id_std"].astype(str) + "::" + silver["microbe_id_std"].astype(str)
    silver["source_dataset"] = "masi_v2"
    silver["label_tier"] = "silver"
    silver["source_record_id"] = _coalesce_series(silver, ["interation_record_id", "interaction_record_id"]).fillna(silver["pair_id"])
    silver["supporting_text"] = _coalesce_series(
        silver,
        ["microbe_change_statistics", "effect_on_microbe", "substance_details", "effect_strength", "outcome"],
    )

    silver = enrich_drug_table_with_rdkit(silver, smiles_columns=["main_component_smiles", "canonical_smiles"])

    columns_to_keep = [
        "pair_id",
        "source_dataset",
        "label_tier",
        "source_record_id",
        "effect_label",
        "binary_effect_label",
        "effect_score",
        "substance_id_std",
        "chemical_name",
        "pubchem_cid",
        "main_component_smiles",
        "canonical_smiles",
        "molecular_formula",
        "molecular_weight",
        "xlogp",
        "tpsa",
        "complexity",
        "therapeutic_class",
        "therapeutic_effect",
        "target_species",
        "human_use",
        "veterinary",
        "species_label",
        "microbe_id_std",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "microbiota_site",
        "substance_category",
        "substance_subcategory",
        "substance_scope",
        "microbe_change_statistics",
        "effect_on_microbe",
        "effect_strength",
        "experimental_system",
        "experimental_organism",
        "experimental_disease_condition",
        "reference_pubmed_id",
        "supporting_text",
    ]
    rdkit_columns = [column for column in silver.columns if column.startswith("rdkit_") or column.startswith("morgan_fp_")]
    existing_columns = [column for column in columns_to_keep if column in silver.columns]
    silver = silver[existing_columns + sorted(rdkit_columns)].copy()

    output_path = output_dir / "step1_silver_masi.csv"
    silver.to_csv(output_path, index=False)

    summary = {
        "raw_dir": str(raw_dir),
        "interaction_file": str(interaction_path),
        "output_path": str(output_path),
        "drug_like_only": drug_like_only,
        "n_records": int(len(silver)),
        "n_drug_like_records": int(silver["substance_scope"].eq("drug_like").sum()) if "substance_scope" in silver.columns else 0,
        "n_substances": int(silver["chemical_name"].nunique()),
        "n_microbes": int(silver["species_label"].nunique()),
        "label_counts": {key: int(value) for key, value in silver["effect_label"].value_counts().to_dict().items()},
        "substance_scope_counts": {key: int(value) for key, value in silver["substance_scope"].value_counts().to_dict().items()},
    }
    _write_manifest(output_dir / "step1_silver_masi_summary.json", summary)
    return summary
