from __future__ import annotations

import json
import math
import re
from difflib import SequenceMatcher
from pathlib import Path

import numpy as np
import pandas as pd

from gut_drug_microbiome.step1.chem_features import enrich_drug_table_with_rdkit
from gut_drug_microbiome.utils.text import canonicalize_key as _canonicalize_key
from gut_drug_microbiome.utils.text import normalize_whitespace as _canonicalize_whitespace


SOURCE_DATASET = "zimmermann_2019"
CONTROL_HEADERS = {"Control pH 7", "Control pH 6", "Control pH 5", "Control pH 4"}

def _snake_case(value: object) -> str:
    """Convert a raw spreadsheet header into a normalized snake_case field name."""
    text = _canonicalize_whitespace(value)
    text = text.replace("%", " pct ")
    text = text.replace("(", " ").replace(")", " ")
    text = text.replace("/", " ")
    text = re.sub(r"[^a-zA-Z0-9]+", "_", text).strip("_").lower()
    return text


def _make_unique_columns(columns: list[str]) -> list[str]:
    """Ensure parsed spreadsheet column names are unique after normalization."""
    counts: dict[str, int] = {}
    result: list[str] = []
    for column in columns:
        key = column or "unnamed"
        count = counts.get(key, 0)
        if count == 0:
            result.append(key)
        else:
            result.append(f"{key}_{count}")
        counts[key] = count + 1
    return result


def _read_raw(workbook_path: str | Path, sheet_name: str) -> pd.DataFrame:
    """Read an Excel sheet without headers so custom header parsing can be applied."""
    return pd.read_excel(workbook_path, sheet_name=sheet_name, header=None)


def _read_with_header(workbook_path: str | Path, sheet_name: str, header_row: int) -> pd.DataFrame:
    """Read an Excel sheet using the specified header row."""
    return pd.read_excel(workbook_path, sheet_name=sheet_name, header=header_row)


def _extract_group_columns(raw: pd.DataFrame, top_row_idx: int, metric_row_idx: int, start_col: int) -> list[dict[str, object]]:
    """Extract grouped metric-column metadata from multi-row spreadsheet headers."""
    top_row = raw.iloc[top_row_idx].tolist()
    metric_row = raw.iloc[metric_row_idx].tolist()
    current_group = None
    columns: list[dict[str, object]] = []
    for index in range(start_col, len(top_row)):
        top_value = top_row[index]
        if pd.notna(top_value):
            current_group = _canonicalize_whitespace(top_value)
        metric_value = metric_row[index]
        columns.append(
            {
                "column_index": index,
                "group_name": current_group,
                "metric_name": _canonicalize_whitespace(metric_value),
            }
        )
    return columns


def _group_metric_blocks(columns: list[dict[str, object]]) -> dict[str, list[dict[str, object]]]:
    """Group parsed column descriptors by their top-level header name."""
    grouped: dict[str, list[dict[str, object]]] = {}
    for record in columns:
        group_name = str(record["group_name"])
        grouped.setdefault(group_name, []).append(record)
    return grouped


def _table3_metric_key(metric_name: str) -> str:
    """Map Supplementary Table 3 metric headers to normalized field names."""
    metric_name = _canonicalize_whitespace(metric_name)
    if metric_name.startswith("% consumed STD "):
        return "percent_consumed_std"
    if metric_name.startswith("% consumed "):
        return "percent_consumed"
    if metric_name.startswith("FC STD "):
        return "fold_change_std"
    if metric_name.startswith("FC "):
        return "fold_change"
    if "p(FDR)" in metric_name:
        return "p_fdr"
    raise ValueError(f"Unsupported Supplementary Table 3 metric: {metric_name}")


def _table6_metric_key(metric_name: str) -> str:
    """Map Supplementary Table 6 metric headers to normalized field names."""
    metric_name = _canonicalize_whitespace(metric_name)
    mapping = {
        "FC (log2) drug vs no-drug pool t=0h": "fc_log2_drug_vs_no_drug_t0",
        "FDR drug vs no-drug pool t=0h": "fdr_drug_vs_no_drug_t0",
        "FC (log2) drug vs no-drug pool t=12h": "fc_log2_drug_vs_no_drug_t12",
        "FDR drug vs no-drug pool t=12h": "fdr_drug_vs_no_drug_t12",
        "Intensity mean t=0h": "intensity_mean_t0",
        "Intensity STD t=0h": "intensity_std_t0",
        "Intensity mean t=12h": "intensity_mean_t12",
        "Intensity STD t=12h": "intensity_std_t12",
        "FC t=12h vs t=0h (log2)": "fc_log2_t12_vs_t0",
        "FC STD t=12h vs t=0h": "fc_std_t12_vs_t0",
        "p(FDR) t=12h vs t=0h": "p_fdr_t12_vs_t0",
    }
    for prefix, key in mapping.items():
        if metric_name.startswith(prefix):
            return key
    raise ValueError(f"Unsupported Supplementary Table 6 metric: {metric_name}")


def _pair_identifier(drug_name: str, microbe_name: str) -> str:
    """Build the project-standard pair identifier for one drug-microbe pair."""
    return f"{drug_name}::{microbe_name}"


def _normalize_parent_depletion_fraction(percent_consumed: object) -> float | None:
    """Convert Zimmermann percent-consumed values into signed depletion fractions."""
    if pd.isna(percent_consumed):
        return None
    try:
        value = float(percent_consumed)
    except (TypeError, ValueError):
        return None
    return -value / 100.0


def _determine_metabolism_label(percent_consumed: object, adaptive_threshold_pct: object, p_fdr: object) -> str:
    """Derive the metabolized/not_metabolized label from Zimmermann assay statistics."""
    if pd.isna(percent_consumed) or pd.isna(adaptive_threshold_pct) or pd.isna(p_fdr):
        return "uncertain"
    try:
        consumed = float(percent_consumed)
        threshold = float(adaptive_threshold_pct)
        significance = float(p_fdr)
    except (TypeError, ValueError):
        return "uncertain"
    if consumed >= threshold and significance < 0.05:
        return "metabolized"
    return "not_metabolized"


def _make_candidate_metabolite_id(index_value: object) -> str:
    """Create a stable candidate metabolite identifier from a worksheet index value."""
    if pd.isna(index_value):
        return "zimmermann_candidate_metabolite_unknown"
    try:
        return f"zimmermann_candidate_metabolite_{int(index_value):06d}"
    except (TypeError, ValueError):
        return f"zimmermann_candidate_metabolite_{_snake_case(index_value)}"


def _prefer_non_empty_strings(primary: pd.Series, fallback: pd.Series) -> pd.Series:
    """Choose non-empty values from the primary series, otherwise fall back."""
    primary = primary.fillna("").astype(str)
    fallback = fallback.fillna("").astype(str)
    return primary.where(primary.str.strip().ne(""), fallback)


def parse_zimmermann_drug_table(workbook_path: str | Path, keep_drug_names: set[str] | None = None) -> pd.DataFrame:
    """Parse Zimmermann Supplementary Table 2 into a normalized drug feature table."""
    table = _read_with_header(workbook_path, "Supplementary Table 2", header_row=3)
    table.columns = _make_unique_columns([_snake_case(column) for column in table.columns])
    table = table.dropna(subset=["molename"]).reset_index(drop=True)
    table = table.copy()
    table["drug_name"] = table["molename"].map(_canonicalize_whitespace)
    table["drug_name_key"] = table["drug_name"].map(_canonicalize_key)
    if keep_drug_names is not None:
        table = table[table["drug_name"].isin(keep_drug_names)].reset_index(drop=True)
    table["drug_id"] = table["drug_name"]
    if "therapeuticindication" in table.columns:
        table["therapeutic_indication"] = table["therapeuticindication"]
    if "estimatedcolonconcentrationmaier2018um" in table.columns:
        table["estimated_colon_concentration_um"] = pd.to_numeric(
            table["estimatedcolonconcentrationmaier2018um"],
            errors="coerce",
        )
    if "molecular_formula" not in table.columns and "molformula" in table.columns:
        table["molecular_formula"] = table["molformula"]
    if "molecular_weight" not in table.columns and "molweight" in table.columns:
        table["molecular_weight"] = pd.to_numeric(table["molweight"], errors="coerce")
    if "xlogp" not in table.columns and "mollogp" in table.columns:
        table["xlogp"] = pd.to_numeric(table["mollogp"], errors="coerce")
    if "smiles" in table.columns:
        table["smiles"] = table["smiles"].map(_canonicalize_whitespace)

    table = table.copy()
    table = enrich_drug_table_with_rdkit(table, smiles_columns=["smiles"])
    table["source_dataset"] = SOURCE_DATASET
    return table


def _match_table1_metadata_to_strain(table1: pd.DataFrame, strain_name: str) -> dict[str, object]:
    """Match a strain name to the closest metadata row in Zimmermann Table 1."""
    strain_key = _canonicalize_key(strain_name)
    best_match: dict[str, object] | None = None
    best_score = -1

    for _, row in table1.iterrows():
        name = _canonicalize_whitespace(row.get("name"))
        description = _canonicalize_whitespace(row.get("phylum_genotype_or_description"))
        reference = _canonicalize_whitespace(row.get("reference"))
        if not name or not reference:
            continue
        candidates = [
            name,
            f"{name} {reference}",
            f"{name} {reference.replace(' ', '')}",
        ]
        for candidate in candidates:
            candidate_key = _canonicalize_key(candidate)
            if not candidate_key:
                continue
            if strain_key.startswith(candidate_key) or candidate_key.startswith(strain_key):
                if len(candidate_key) > best_score:
                    best_match = {
                        "species_name": name,
                        "phylum_or_description": description,
                        "reference": reference,
                    }
                    best_score = len(candidate_key)
                continue
            reference_key = _canonicalize_key(reference)
            reference_matches = bool(reference_key and reference_key in strain_key)
            if reference_matches:
                similarity = SequenceMatcher(None, strain_key, candidate_key).ratio()
                fuzzy_score = similarity + 1.0
                if similarity >= 0.9 and fuzzy_score > best_score:
                    best_match = {
                        "species_name": name,
                        "phylum_or_description": description,
                        "reference": reference,
                    }
                    best_score = fuzzy_score
    return {} if best_match is None else best_match


def parse_zimmermann_microbe_table(workbook_path: str | Path, strain_names: list[str]) -> pd.DataFrame:
    """Parse Zimmermann strain metadata for the requested list of strain names."""
    table1 = _read_with_header(workbook_path, "Supplementary Table 1", header_row=3)
    table1 = table1.rename(columns={column: _snake_case(column) for column in table1.columns})
    table1 = table1.dropna(subset=["name"]).reset_index(drop=True)

    records = []
    for strain_name in strain_names:
        metadata = _match_table1_metadata_to_strain(table1, strain_name)
        record = {
            "microbe_id": strain_name,
            "microbe_name": strain_name,
            "species_name": metadata.get("species_name"),
            "species_label": metadata.get("species_name"),
            "phylum_or_description": metadata.get("phylum_or_description"),
            "phylum": metadata.get("phylum_or_description"),
            "reference": metadata.get("reference"),
            "source_dataset": SOURCE_DATASET,
        }
        if metadata.get("species_name"):
            species_text = str(metadata["species_name"])
            parts = species_text.split()
            record["genus"] = parts[0] if parts else np.nan
            record["species_epithet"] = parts[1] if len(parts) > 1 else np.nan
        else:
            record["genus"] = np.nan
            record["species_epithet"] = np.nan
        records.append(record)
    return pd.DataFrame(records)


def parse_zimmermann_parent_screen(workbook_path: str | Path) -> pd.DataFrame:
    """Parse Zimmermann parent-drug depletion screens into a normalized label table."""
    raw = _read_raw(workbook_path, "Supplementary Table 3")
    columns = _extract_group_columns(raw, top_row_idx=3, metric_row_idx=4, start_col=2)
    grouped = _group_metric_blocks(columns)

    data = raw.iloc[5:].reset_index(drop=True).copy()
    data = data[data.iloc[:, 0].notna()].reset_index(drop=True)
    data = data.rename(columns={0: "drug_name", 1: "drug_adaptive_fc_threshold_pct"})
    data["drug_name"] = data["drug_name"].map(_canonicalize_whitespace)
    data["drug_name_key"] = data["drug_name"].map(_canonicalize_key)
    data["drug_id"] = data["drug_name"]

    records = []
    for _, row in data.iterrows():
        drug_name = row["drug_name"]
        adaptive_threshold = pd.to_numeric(pd.Series([row["drug_adaptive_fc_threshold_pct"]]), errors="coerce").iloc[0]
        for group_name, block in grouped.items():
            if group_name in CONTROL_HEADERS:
                continue
            metric_values: dict[str, object] = {}
            for item in block:
                metric_key = _table3_metric_key(str(item["metric_name"]))
                metric_values[metric_key] = row.iloc[int(item["column_index"])]

            percent_consumed = pd.to_numeric(pd.Series([metric_values.get("percent_consumed")]), errors="coerce").iloc[0]
            p_fdr = pd.to_numeric(pd.Series([metric_values.get("p_fdr")]), errors="coerce").iloc[0]
            parent_depletion_fraction = _normalize_parent_depletion_fraction(percent_consumed)
            metabolism_label = _determine_metabolism_label(percent_consumed, adaptive_threshold, p_fdr)
            microbe_name = group_name

            records.append(
                {
                    "pair_id": _pair_identifier(drug_name, microbe_name),
                    "prestwick_id": drug_name,
                    "nt_code": microbe_name,
                    "drug_id": drug_name,
                    "microbe_id": microbe_name,
                    "drug_name": drug_name,
                    "microbe_name": microbe_name,
                    "drug_name_key": _canonicalize_key(drug_name),
                    "drug_adaptive_fc_threshold_pct": adaptive_threshold,
                    "percent_consumed": percent_consumed,
                    "percent_consumed_std": pd.to_numeric(pd.Series([metric_values.get("percent_consumed_std")]), errors="coerce").iloc[0],
                    "fold_change": pd.to_numeric(pd.Series([metric_values.get("fold_change")]), errors="coerce").iloc[0],
                    "fold_change_std": pd.to_numeric(pd.Series([metric_values.get("fold_change_std")]), errors="coerce").iloc[0],
                    "p_fdr": p_fdr,
                    "parent_depletion_fraction": parent_depletion_fraction,
                    "metabolism_label": metabolism_label,
                    "reaction_class": np.nan,
                    "product_ids": "",
                    "evidence_gene_ids": "",
                    "source_dataset": SOURCE_DATASET,
                    "label_tier": "gold",
                    "source_scope": "isolate",
                    "source_record_id": _pair_identifier(drug_name, microbe_name),
                    "raw_metabolism_label": (
                        f"percent_consumed={percent_consumed};adaptive_threshold_pct={adaptive_threshold};p_fdr={p_fdr}"
                    ),
                    "raw_reaction_class": np.nan,
                }
            )
    return pd.DataFrame(records)


def parse_zimmermann_metabolite_candidates(workbook_path: str | Path) -> pd.DataFrame:
    """Parse Zimmermann candidate metabolite metadata from Supplementary Table 5."""
    table = _read_with_header(workbook_path, "Supplementary Table 5", header_row=3)
    table = table.rename(columns={column: _snake_case(column) for column in table.columns})
    table = table.dropna(subset=["parentdrug", "index"]).reset_index(drop=True)
    table["drug_name"] = table["parentdrug"].map(_canonicalize_whitespace)
    table["drug_name_key"] = table["drug_name"].map(_canonicalize_key)
    table["candidate_metabolite_id"] = table["index"].map(_make_candidate_metabolite_id)
    table["mz"] = pd.to_numeric(table["mz"], errors="coerce")
    table["rt"] = pd.to_numeric(table["rt"], errors="coerce")
    table["mz_delta"] = pd.to_numeric(table["drugmassdelta"], errors="coerce")
    table["good_filter"] = pd.to_numeric(table["goodfilter"], errors="coerce").fillna(0).astype(int)
    table["drug_consumed_flag"] = pd.to_numeric(table["drugconsumedflag"], errors="coerce").fillna(0).astype(int)
    table["source_dataset"] = SOURCE_DATASET
    return table


def parse_zimmermann_metabolite_long(
    workbook_path: str | Path,
    metabolite_candidates: pd.DataFrame,
) -> pd.DataFrame:
    """Parse long-form Zimmermann metabolite evidence from Supplementary Table 6."""
    raw = _read_raw(workbook_path, "Supplementary Table 6")
    columns = _extract_group_columns(raw, top_row_idx=3, metric_row_idx=4, start_col=4)
    grouped = _group_metric_blocks(columns)

    data = raw.iloc[5:].reset_index(drop=True).copy()
    data = data[data.iloc[:, 0].notna()].reset_index(drop=True)
    data = data.rename(columns={0: "drug_name", 1: "mz", 2: "rt", 3: "mz_delta"})
    data["drug_name"] = data["drug_name"].map(_canonicalize_whitespace)
    data["drug_name_key"] = data["drug_name"].map(_canonicalize_key)
    data["mz"] = pd.to_numeric(data["mz"], errors="coerce")
    data["rt"] = pd.to_numeric(data["rt"], errors="coerce")
    data["mz_delta"] = pd.to_numeric(data["mz_delta"], errors="coerce")

    candidate_lookup = metabolite_candidates.loc[:, ["candidate_metabolite_id", "drug_name", "mz", "rt", "good_filter"]].copy()
    candidate_lookup["join_key"] = (
        candidate_lookup["drug_name"].astype(str)
        + "||"
        + candidate_lookup["mz"].round(3).astype(str)
        + "||"
        + candidate_lookup["rt"].round(3).astype(str)
    )

    records = []
    for _, row in data.iterrows():
        join_key = f"{row['drug_name']}||{round(float(row['mz']), 3)}||{round(float(row['rt']), 3)}"
        matches = candidate_lookup[candidate_lookup["join_key"] == join_key]
        candidate_id = matches["candidate_metabolite_id"].iloc[0] if not matches.empty else _make_candidate_metabolite_id(join_key)
        good_filter = int(matches["good_filter"].iloc[0]) if not matches.empty else 0

        for group_name, block in grouped.items():
            if group_name in CONTROL_HEADERS:
                continue
            metric_values: dict[str, object] = {}
            for item in block:
                metric_key = _table6_metric_key(str(item["metric_name"]))
                metric_values[metric_key] = row.iloc[int(item["column_index"])]
            fc_t12_vs_t0 = pd.to_numeric(pd.Series([metric_values.get("fc_log2_t12_vs_t0")]), errors="coerce").iloc[0]
            p_fdr_t12_vs_t0 = pd.to_numeric(pd.Series([metric_values.get("p_fdr_t12_vs_t0")]), errors="coerce").iloc[0]
            product_detected = bool(
                good_filter == 1
                and pd.notna(fc_t12_vs_t0)
                and pd.notna(p_fdr_t12_vs_t0)
                and float(fc_t12_vs_t0) > 1.0
                and float(p_fdr_t12_vs_t0) < 0.05
            )
            records.append(
                {
                    "pair_id": _pair_identifier(row["drug_name"], group_name),
                    "drug_name": row["drug_name"],
                    "drug_name_key": row["drug_name_key"],
                    "microbe_name": group_name,
                    "candidate_metabolite_id": candidate_id,
                    "mz": row["mz"],
                    "rt": row["rt"],
                    "mz_delta": row["mz_delta"],
                    "good_filter": good_filter,
                    "fc_log2_drug_vs_no_drug_t0": pd.to_numeric(pd.Series([metric_values.get("fc_log2_drug_vs_no_drug_t0")]), errors="coerce").iloc[0],
                    "fdr_drug_vs_no_drug_t0": pd.to_numeric(pd.Series([metric_values.get("fdr_drug_vs_no_drug_t0")]), errors="coerce").iloc[0],
                    "fc_log2_drug_vs_no_drug_t12": pd.to_numeric(pd.Series([metric_values.get("fc_log2_drug_vs_no_drug_t12")]), errors="coerce").iloc[0],
                    "fdr_drug_vs_no_drug_t12": pd.to_numeric(pd.Series([metric_values.get("fdr_drug_vs_no_drug_t12")]), errors="coerce").iloc[0],
                    "intensity_mean_t0": pd.to_numeric(pd.Series([metric_values.get("intensity_mean_t0")]), errors="coerce").iloc[0],
                    "intensity_std_t0": pd.to_numeric(pd.Series([metric_values.get("intensity_std_t0")]), errors="coerce").iloc[0],
                    "intensity_mean_t12": pd.to_numeric(pd.Series([metric_values.get("intensity_mean_t12")]), errors="coerce").iloc[0],
                    "intensity_std_t12": pd.to_numeric(pd.Series([metric_values.get("intensity_std_t12")]), errors="coerce").iloc[0],
                    "fc_log2_t12_vs_t0": fc_t12_vs_t0,
                    "fc_std_t12_vs_t0": pd.to_numeric(pd.Series([metric_values.get("fc_std_t12_vs_t0")]), errors="coerce").iloc[0],
                    "p_fdr_t12_vs_t0": p_fdr_t12_vs_t0,
                    "product_detected": product_detected,
                    "source_dataset": SOURCE_DATASET,
                }
            )
    return pd.DataFrame(records)


def parse_zimmermann_gene_links(workbook_path: str | Path) -> pd.DataFrame:
    """Parse Zimmermann gene-drug association evidence from Supplementary Table 13."""
    raw = _read_raw(workbook_path, "Supplementary Table 13")
    top_row = raw.iloc[3].tolist()
    sub_row = raw.iloc[4].tolist()
    data = raw.iloc[5:].reset_index(drop=True).copy()
    data = data[data.iloc[:, 0].notna()].reset_index(drop=True)
    data = data.rename(
        columns={
            0: "refseq_locus_tag",
            1: "patric_id",
            2: "product",
            3: "protein_id",
        }
    )

    records = []
    for column_index in range(4, len(top_row)):
        drug_name = _canonicalize_whitespace(top_row[column_index])
        sub_header = _canonicalize_whitespace(sub_row[column_index])
        if not drug_name or sub_header != "Parent drug":
            continue
        indicator = pd.to_numeric(data.iloc[:, column_index], errors="coerce").fillna(0)
        hits = data[indicator == 1].copy()
        for _, row in hits.iterrows():
            records.append(
                {
                    "drug_name": drug_name,
                    "drug_name_key": _canonicalize_key(drug_name),
                    "refseq_locus_tag": _canonicalize_whitespace(row["refseq_locus_tag"]),
                    "patric_id": _canonicalize_whitespace(row["patric_id"]),
                    "product": _canonicalize_whitespace(row["product"]),
                    "protein_id": _canonicalize_whitespace(row["protein_id"]),
                    "source_dataset": SOURCE_DATASET,
                }
            )
    return pd.DataFrame(records)


def normalize_zimmermann_2019(
    input_path: str | Path,
    output_dir: str | Path,
) -> dict[str, object]:
    """Normalize the Zimmermann 2019 workbook into Step 2 tables and summaries.

    Args:
        input_path: Path to the Zimmermann supplementary workbook.
        output_dir: Directory where normalized tables and summary JSON are written.

    Returns:
        A summary dictionary with counts and generated file paths.
    """
    input_path = Path(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    parent_screen = parse_zimmermann_parent_screen(input_path)
    drug_names = set(parent_screen["drug_name"].dropna().unique().tolist())
    strain_names = sorted(parent_screen["microbe_name"].dropna().unique().tolist())

    drug_table = parse_zimmermann_drug_table(input_path, keep_drug_names=drug_names)
    microbe_table = parse_zimmermann_microbe_table(input_path, strain_names=strain_names)
    metabolite_candidates = parse_zimmermann_metabolite_candidates(input_path)
    metabolite_long = parse_zimmermann_metabolite_long(input_path, metabolite_candidates=metabolite_candidates)
    gene_links = parse_zimmermann_gene_links(input_path)

    gene_summary = (
        gene_links.groupby("drug_name_key")["refseq_locus_tag"]
        .apply(lambda values: ";".join(sorted({value for value in values if value})))
        .reset_index()
        .rename(columns={"refseq_locus_tag": "evidence_gene_ids"})
    )
    product_summary = (
        metabolite_long[metabolite_long["product_detected"]]
        .groupby("pair_id")["candidate_metabolite_id"]
        .apply(lambda values: ";".join(sorted({value for value in values if value})))
        .reset_index()
        .rename(columns={"candidate_metabolite_id": "product_ids"})
    )

    label_table = parent_screen.merge(gene_summary, on="drug_name_key", how="left", suffixes=("", "_new"))
    if "evidence_gene_ids_new" in label_table.columns:
        label_table["evidence_gene_ids"] = _prefer_non_empty_strings(
            label_table["evidence_gene_ids_new"],
            label_table["evidence_gene_ids"],
        )
        label_table = label_table.drop(columns=["evidence_gene_ids_new"])
    label_table = label_table.merge(product_summary, on="pair_id", how="left", suffixes=("", "_new"))
    if "product_ids_new" in label_table.columns:
        label_table["product_ids"] = _prefer_non_empty_strings(
            label_table["product_ids_new"],
            label_table["product_ids"],
        )
        label_table = label_table.drop(columns=["product_ids_new"])

    label_table["reaction_class"] = label_table["reaction_class"].astype(object)
    unresolved_mask = label_table["metabolism_label"].eq("metabolized") & label_table["product_ids"].fillna("").eq("")
    label_table.loc[unresolved_mask, "reaction_class"] = "bioaccumulation_or_unresolved_depletion"

    modeling_table = label_table.merge(
        drug_table,
        on="drug_name",
        how="left",
        suffixes=("", "_drug"),
    )
    modeling_table = modeling_table.merge(
        microbe_table,
        on="microbe_name",
        how="left",
        suffixes=("", "_microbe"),
    )
    modeling_table["n_product_candidates"] = modeling_table["product_ids"].fillna("").map(
        lambda value: 0 if not str(value).strip() else len([item for item in str(value).split(";") if item.strip()])
    )
    modeling_table["n_evidence_genes"] = modeling_table["evidence_gene_ids"].fillna("").map(
        lambda value: 0 if not str(value).strip() else len([item for item in str(value).split(";") if item.strip()])
    )
    modeling_table["step2_label_available"] = True
    modeling_table["step2_has_resolved_products"] = modeling_table["product_ids"].fillna("").astype(str).str.strip().ne("")

    label_table_path = output_dir / "zimmermann_2019_label_table.csv"
    drug_table_path = output_dir / "zimmermann_2019_drug_table.csv"
    microbe_table_path = output_dir / "zimmermann_2019_microbe_table.csv"
    metabolite_candidates_path = output_dir / "zimmermann_2019_metabolite_candidates.csv"
    metabolite_long_path = output_dir / "zimmermann_2019_metabolite_long.csv"
    gene_links_path = output_dir / "zimmermann_2019_gene_links.csv"
    modeling_table_path = output_dir / "zimmermann_2019_modeling_table.csv"

    label_table.to_csv(label_table_path, index=False)
    drug_table.to_csv(drug_table_path, index=False)
    microbe_table.to_csv(microbe_table_path, index=False)
    metabolite_candidates.to_csv(metabolite_candidates_path, index=False)
    metabolite_long.to_csv(metabolite_long_path, index=False)
    gene_links.to_csv(gene_links_path, index=False)
    modeling_table.to_csv(modeling_table_path, index=False)

    summary = {
        "input_path": str(input_path),
        "output_dir": str(output_dir),
        "n_label_rows": int(len(label_table)),
        "n_drugs": int(label_table["drug_name"].nunique()),
        "n_microbes": int(label_table["microbe_name"].nunique()),
        "n_metabolized": int(label_table["metabolism_label"].eq("metabolized").sum()),
        "n_not_metabolized": int(label_table["metabolism_label"].eq("not_metabolized").sum()),
        "n_uncertain": int(label_table["metabolism_label"].eq("uncertain").sum()),
        "n_pairs_with_candidate_products": int(label_table["product_ids"].fillna("").astype(str).str.strip().ne("").sum()),
        "n_unique_candidate_metabolites": int(metabolite_candidates["candidate_metabolite_id"].nunique()),
        "n_gene_links": int(len(gene_links)),
        "n_pairs_with_gene_evidence": int(label_table["evidence_gene_ids"].fillna("").astype(str).str.strip().ne("").sum()),
        "label_table_path": str(label_table_path),
        "drug_table_path": str(drug_table_path),
        "microbe_table_path": str(microbe_table_path),
        "metabolite_candidates_path": str(metabolite_candidates_path),
        "metabolite_long_path": str(metabolite_long_path),
        "gene_links_path": str(gene_links_path),
        "modeling_table_path": str(modeling_table_path),
        "metabolism_label_counts": {
            str(key): int(value)
            for key, value in label_table["metabolism_label"].value_counts(dropna=False).to_dict().items()
        },
        "reaction_class_counts": {
            str(key): int(value)
            for key, value in label_table["reaction_class"].fillna("missing").value_counts(dropna=False).to_dict().items()
        },
    }
    (output_dir / "zimmermann_2019_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return summary
