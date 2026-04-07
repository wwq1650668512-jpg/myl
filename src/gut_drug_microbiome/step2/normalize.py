from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


NORMALIZED_STEP2_COLUMNS = [
    "pair_id",
    "prestwick_id",
    "nt_code",
    "drug_name",
    "microbe_name",
    "metabolism_label",
    "reaction_class",
    "parent_depletion_fraction",
    "product_ids",
    "evidence_gene_ids",
    "source_dataset",
    "label_tier",
    "source_scope",
    "source_record_id",
    "raw_metabolism_label",
    "raw_reaction_class",
]

COLUMN_ALIASES = {
    "pair_id": ["pair_id", "drug_microbe_pair_id"],
    "prestwick_id": ["prestwick_id", "drug_id", "compound_id", "drug_identifier"],
    "nt_code": ["nt_code", "microbe_id", "bacterial_isolate", "isolate_id", "strain_id", "taxonomy_id"],
    "drug_name": ["drug_name", "compound_name", "drug", "compound", "parent_drug"],
    "microbe_name": ["microbe_name", "bacteria", "bacterial_isolate_name", "strain_name", "microbe"],
    "metabolism_label": [
        "metabolism_label",
        "depletion_or_metabolism_label",
        "metabolized",
        "label",
        "metabolism",
        "depletion_label",
    ],
    "reaction_class": [
        "reaction_class",
        "reaction_type",
        "biotransformation_type",
        "reaction_annotation",
        "reaction",
    ],
    "parent_depletion_fraction": [
        "parent_depletion_fraction",
        "parent_fraction_change",
        "parent_drug_fraction_change",
        "drug_depletion_fraction",
        "depletion_fraction",
        "percent_parent_remaining_change",
        "parent_change",
    ],
    "product_ids": [
        "product_ids",
        "metabolite_identity",
        "metabolite_id",
        "metabolite_ids",
        "product_id",
        "product_name",
        "products",
    ],
    "evidence_gene_ids": [
        "evidence_gene_ids",
        "gene_association",
        "gene_ids",
        "genes",
        "associated_genes",
    ],
    "source_record_id": ["source_record_id", "record_id", "sample_id", "assay_id", "reference_id"],
}


def _canonicalize_column_name(value: str) -> str:
    """Normalize a raw column label into a lowercase matching key."""
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")


def read_table_auto(input_path: str | Path, sheet_name: str | int | None = None) -> pd.DataFrame:
    """Read a CSV, TSV, or Excel table based on its suffix."""
    input_path = Path(input_path)
    suffix = input_path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(input_path, low_memory=False)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(input_path, sep="\t", low_memory=False)
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(input_path, sheet_name=sheet_name if sheet_name is not None else 0)
    raise ValueError(f"Unsupported table format: {input_path}")


def _find_matching_column(frame: pd.DataFrame, aliases: list[str]) -> str | None:
    """Find the first column in a frame that matches one of the allowed aliases."""
    normalized_map = {_canonicalize_column_name(column): column for column in frame.columns}
    for alias in aliases:
        key = _canonicalize_column_name(alias)
        if key in normalized_map:
            return normalized_map[key]
    return None


def _normalize_fraction(value: object) -> float | None:
    """Normalize a numeric fraction field, accepting percentages or 0-1 values."""
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text:
        return None
    text = text.replace("%", "")
    try:
        numeric = float(text)
    except ValueError:
        return None
    if abs(numeric) > 1.0 and abs(numeric) <= 100.0:
        return numeric / 100.0
    return numeric


def _normalize_multi_value(value: object) -> str:
    """Split and deduplicate a multi-value field into a semicolon-delimited string."""
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    parts = re.split(r"[;,|]", text)
    normalized = []
    seen = set()
    for part in parts:
        item = part.strip()
        if not item:
            continue
        lowered = item.lower()
        if lowered in seen:
            continue
        seen.add(lowered)
        normalized.append(item)
    return ";".join(normalized)


def _normalize_metabolism_label(value: object) -> str:
    """Map heterogeneous label text or booleans into the Step 2 label vocabulary."""
    if pd.isna(value):
        return "uncertain"
    if isinstance(value, (bool, np.bool_)):
        return "metabolized" if bool(value) else "not_metabolized"

    text = str(value).strip().lower()
    if not text:
        return "uncertain"
    if text in {"1", "yes", "y", "true", "metabolized", "metabolism", "depleted", "transformed", "converted"}:
        return "metabolized"
    if text in {"0", "no", "n", "false", "not_metabolized", "not metabolized", "stable", "unchanged"}:
        return "not_metabolized"
    if any(token in text for token in ["uncertain", "ambiguous", "mixed", "partial", "possible"]):
        return "uncertain"
    if any(token in text for token in ["deplet", "metabol", "product", "transform", "convert"]):
        return "metabolized"
    if any(token in text for token in ["no change", "not detect", "none"]):
        return "not_metabolized"
    return "uncertain"


def _normalize_reaction_class(value: object) -> str | None:
    """Map free-text reaction descriptions into coarse Step 2 reaction classes."""
    if pd.isna(value):
        return None
    text = str(value).strip().lower()
    if not text:
        return None
    if "deacetyl" in text:
        return "deacetylation"
    if "dehydroxyl" in text:
        return "dehydroxylation"
    if "demethyl" in text:
        return "demethylation"
    if "deconjug" in text or "deglucuron" in text or "deglycos" in text or "desulfat" in text:
        return "deconjugation"
    if "ring cleav" in text or "ring-open" in text or "ring open" in text:
        return "ring_cleavage"
    if "hydrolys" in text:
        return "hydrolysis"
    if "reduct" in text:
        return "reduction"
    if "bioaccum" in text or "unresolved" in text or "depletion" in text:
        return "bioaccumulation_or_unresolved_depletion"
    if text in {"other", "misc", "unknown"}:
        return "other"
    return "other"


def normalize_step2_label_table(
    input_path: str | Path,
    output_path: str | Path,
    source_dataset: str,
    label_tier: str = "gold",
    source_scope: str = "isolate",
    sheet_name: str | int | None = None,
) -> dict[str, object]:
    """Normalize an external Step 2 label table into the project-wide schema.

    Args:
        input_path: Raw input table path.
        output_path: Destination CSV for the normalized table.
        source_dataset: Dataset identifier written into the normalized output.
        label_tier: Label quality tier such as gold or silver.
        source_scope: Experimental scope such as isolate or community.
        sheet_name: Optional Excel sheet selector.

    Returns:
        A summary dictionary with selected columns, counts, and output metadata.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    raw = read_table_auto(input_path, sheet_name=sheet_name)
    frame = raw.copy()
    selected_columns: dict[str, str] = {}
    for standard_name, aliases in COLUMN_ALIASES.items():
        matched = _find_matching_column(frame, aliases)
        if matched is not None:
            selected_columns[standard_name] = matched

    normalized = pd.DataFrame(index=frame.index)
    if "pair_id" in selected_columns:
        normalized["pair_id"] = frame[selected_columns["pair_id"]]
    if "prestwick_id" in selected_columns:
        normalized["prestwick_id"] = frame[selected_columns["prestwick_id"]]
    if "nt_code" in selected_columns:
        normalized["nt_code"] = frame[selected_columns["nt_code"]]
    if "drug_name" in selected_columns:
        normalized["drug_name"] = frame[selected_columns["drug_name"]]
    if "microbe_name" in selected_columns:
        normalized["microbe_name"] = frame[selected_columns["microbe_name"]]

    if "pair_id" not in normalized.columns:
        if {"prestwick_id", "nt_code"}.issubset(normalized.columns):
            normalized["pair_id"] = (
                normalized["prestwick_id"].fillna("unknown_drug").astype(str)
                + "::"
                + normalized["nt_code"].fillna("unknown_microbe").astype(str)
            )
        else:
            normalized["pair_id"] = [f"{source_dataset}_row_{index}" for index in range(len(normalized))]

    raw_metabolism_series = frame[selected_columns["metabolism_label"]] if "metabolism_label" in selected_columns else pd.Series(np.nan, index=frame.index)
    raw_reaction_series = frame[selected_columns["reaction_class"]] if "reaction_class" in selected_columns else pd.Series(np.nan, index=frame.index)

    normalized["raw_metabolism_label"] = raw_metabolism_series
    normalized["raw_reaction_class"] = raw_reaction_series
    normalized["metabolism_label"] = raw_metabolism_series.map(_normalize_metabolism_label)
    normalized["reaction_class"] = raw_reaction_series.map(_normalize_reaction_class)

    if "parent_depletion_fraction" in selected_columns:
        normalized["parent_depletion_fraction"] = frame[selected_columns["parent_depletion_fraction"]].map(_normalize_fraction)
    else:
        normalized["parent_depletion_fraction"] = np.nan

    if "product_ids" in selected_columns:
        normalized["product_ids"] = frame[selected_columns["product_ids"]].map(_normalize_multi_value)
    else:
        normalized["product_ids"] = ""

    if "evidence_gene_ids" in selected_columns:
        normalized["evidence_gene_ids"] = frame[selected_columns["evidence_gene_ids"]].map(_normalize_multi_value)
    else:
        normalized["evidence_gene_ids"] = ""

    if "source_record_id" in selected_columns:
        normalized["source_record_id"] = frame[selected_columns["source_record_id"]]
    else:
        normalized["source_record_id"] = [f"{source_dataset}_{index}" for index in range(len(normalized))]

    normalized["source_dataset"] = source_dataset
    normalized["label_tier"] = label_tier
    normalized["source_scope"] = source_scope

    unresolved_mask = (
        normalized["metabolism_label"].eq("metabolized")
        & normalized["reaction_class"].isna()
        & normalized["product_ids"].eq("")
    )
    normalized.loc[unresolved_mask, "reaction_class"] = "bioaccumulation_or_unresolved_depletion"

    for column in NORMALIZED_STEP2_COLUMNS:
        if column not in normalized.columns:
            normalized[column] = np.nan

    normalized = normalized.loc[:, NORMALIZED_STEP2_COLUMNS].copy()
    normalized.to_csv(output_path, index=False)

    summary = {
        "input_path": str(input_path),
        "output_path": str(output_path),
        "source_dataset": source_dataset,
        "label_tier": label_tier,
        "source_scope": source_scope,
        "n_rows": int(len(normalized)),
        "n_unique_pairs": int(normalized["pair_id"].nunique(dropna=True)),
        "metabolism_label_counts": {
            str(key): int(value)
            for key, value in normalized["metabolism_label"].value_counts(dropna=False).to_dict().items()
        },
        "reaction_class_counts": {
            str(key): int(value)
            for key, value in normalized["reaction_class"].fillna("missing").value_counts(dropna=False).to_dict().items()
        },
        "selected_columns": selected_columns,
    }
    summary_path = output_path.with_name(output_path.stem + "_summary.json")
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    return summary
