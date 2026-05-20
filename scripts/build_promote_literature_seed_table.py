from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.step1.chem_features import enrich_drug_table_with_rdkit


REQUIRED_COLUMNS = [
    "record_id",
    "compound_name_raw",
    "compound_name_normalized",
    "smiles",
    "inchikey",
    "microbe_name_raw",
    "microbe_label_normalized",
    "strain_label",
    "effect_direction",
    "effect_type",
    "effect_score_proxy",
    "effect_score_proxy_type",
    "dose_value",
    "dose_unit",
    "culture_context",
    "supporting_microbe",
    "source_pmid",
    "source_title",
    "evidence_level",
    "notes",
]

VALID_EFFECT_DIRECTIONS = {"promote", "inhibit", "no_effect"}
VALID_EFFECT_TYPES = {"direct_promote", "metabolism_supported_promote", "functional_promote", "direct_inhibit", "other"}
VALID_EVIDENCE_LEVELS = {"high", "medium", "low"}
VALID_PROXY_TYPES = {"od_ratio", "cfu_change", "relative_abundance_change", "qualitative_increase", "qualitative_decrease"}


def _ensure_required_columns(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    for column in REQUIRED_COLUMNS:
        if column not in result.columns:
            result[column] = np.nan
    return result.loc[:, REQUIRED_COLUMNS].copy()


def _clean_text(value: object) -> object:
    if pd.isna(value):
        return np.nan
    text = str(value).strip()
    return text if text else np.nan


def _normalize_frame(frame: pd.DataFrame) -> pd.DataFrame:
    result = _ensure_required_columns(frame)
    for column in result.columns:
        if column in {"effect_score_proxy", "dose_value", "source_pmid"}:
            continue
        result[column] = result[column].map(_clean_text)

    result["effect_score_proxy"] = pd.to_numeric(result["effect_score_proxy"], errors="coerce")
    result["dose_value"] = pd.to_numeric(result["dose_value"], errors="coerce")
    result["source_pmid"] = pd.to_numeric(result["source_pmid"], errors="coerce").astype("Int64")

    result["effect_direction"] = result["effect_direction"].str.lower()
    result["effect_type"] = result["effect_type"].str.lower()
    result["evidence_level"] = result["evidence_level"].str.lower()
    result["effect_score_proxy_type"] = result["effect_score_proxy_type"].str.lower()
    return result


def _validate_frame(frame: pd.DataFrame) -> list[str]:
    issues: list[str] = []
    if frame["record_id"].isna().any():
        issues.append("record_id contains missing values")
    if frame["record_id"].duplicated().any():
        duplicated = frame.loc[frame["record_id"].duplicated(), "record_id"].astype(str).tolist()
        issues.append(f"record_id duplicated: {duplicated[:5]}")
    invalid_direction = sorted(set(frame["effect_direction"].dropna()) - VALID_EFFECT_DIRECTIONS)
    if invalid_direction:
        issues.append(f"invalid effect_direction values: {invalid_direction}")
    invalid_type = sorted(set(frame["effect_type"].dropna()) - VALID_EFFECT_TYPES)
    if invalid_type:
        issues.append(f"invalid effect_type values: {invalid_type}")
    invalid_level = sorted(set(frame["evidence_level"].dropna()) - VALID_EVIDENCE_LEVELS)
    if invalid_level:
        issues.append(f"invalid evidence_level values: {invalid_level}")
    invalid_proxy_type = sorted(set(frame["effect_score_proxy_type"].dropna()) - VALID_PROXY_TYPES)
    if invalid_proxy_type:
        issues.append(f"invalid effect_score_proxy_type values: {invalid_proxy_type}")
    if frame["compound_name_raw"].isna().any():
        issues.append("compound_name_raw contains missing values")
    if frame["microbe_name_raw"].isna().any():
        issues.append("microbe_name_raw contains missing values")
    return issues


def build_promote_literature_seed_table(
    input_path: str | Path,
    output_path: str | Path,
    summary_path: str | Path | None = None,
) -> dict[str, object]:
    input_path = Path(input_path)
    output_path = Path(output_path)
    summary_path = output_path.with_suffix(".summary.json") if summary_path is None else Path(summary_path)

    raw = pd.read_csv(input_path, low_memory=False)
    normalized = _normalize_frame(raw)
    issues = _validate_frame(normalized)
    if issues:
        raise ValueError("Promote literature seed table validation failed: " + "; ".join(issues))

    chem_input = pd.DataFrame(
        {
            "record_id": normalized["record_id"],
            "main_component_smiles": normalized["smiles"],
            "smiles": normalized["smiles"],
        }
    )
    chem_enriched = enrich_drug_table_with_rdkit(chem_input, smiles_columns=["main_component_smiles", "smiles"])
    keep_chem_columns = [
        column
        for column in [
            "record_id",
            "canonical_smiles_rdkit",
            "inchikey",
            "murcko_scaffold",
            "rdkit_formula",
            "rdkit_valid_smiles",
            "rdkit_exact_mol_wt",
            "rdkit_logp",
            "rdkit_tpsa",
        ]
        if column in chem_enriched.columns
    ]
    output = normalized.merge(chem_enriched.loc[:, keep_chem_columns], on="record_id", how="left", suffixes=("", "_rdkit"))

    if "inchikey_rdkit" in output.columns:
        if "inchikey" in output.columns:
            output["inchikey"] = output["inchikey"].combine_first(output["inchikey_rdkit"])
        else:
            output["inchikey"] = output["inchikey_rdkit"]
        output = output.drop(columns=["inchikey_rdkit"])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, index=False)

    summary = {
        "input_path": str(input_path),
        "output_path": str(output_path),
        "n_rows": int(len(output)),
        "effect_direction_counts": {
            str(key): int(value) for key, value in output["effect_direction"].fillna("missing").value_counts().to_dict().items()
        },
        "effect_type_counts": {
            str(key): int(value) for key, value in output["effect_type"].fillna("missing").value_counts().to_dict().items()
        },
        "evidence_level_counts": {
            str(key): int(value) for key, value in output["evidence_level"].fillna("missing").value_counts().to_dict().items()
        },
        "n_rows_with_smiles": int(output["smiles"].notna().sum()),
        "n_rows_with_valid_rdkit_smiles": int(pd.to_numeric(output.get("rdkit_valid_smiles"), errors="coerce").fillna(0).astype(bool).sum()),
    }
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Normalize a manually curated promote literature seed table.")
    parser.add_argument(
        "--input",
        default=ROOT / "data/reference/promote_literature_seed_template.csv",
        type=Path,
        help="Manual literature curation template CSV.",
    )
    parser.add_argument(
        "--output",
        default=ROOT / "data/reference/promote_literature_seed_table.csv",
        type=Path,
        help="Normalized literature seed table output path.",
    )
    parser.add_argument(
        "--summary",
        default=None,
        type=Path,
        help="Optional JSON summary path.",
    )
    args = parser.parse_args()
    summary = build_promote_literature_seed_table(args.input, args.output, args.summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
