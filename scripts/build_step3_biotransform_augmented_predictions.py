from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_PATH = ROOT / "predictions/step2/baseline_scaffold_v1_83/predictions.csv"
DEFAULT_OUTPUT_PATH = ROOT / "predictions/step2/baseline_scaffold_v1_83/predictions_biotransform_experimental.csv"
DEFAULT_SUMMARY_PATH = ROOT / "predictions/step2/baseline_scaffold_v1_83/predictions_biotransform_experimental.summary.json"
DEFAULT_PROXIMAL_OUTPUT_ROOT = (
    ROOT
    / "external/papers/Computational_analysis_gut_microbiota_drug_metabolism_code/PROXIMAL2_Supplementary/output/set1_subset2_products"
)
DEFAULT_PROXIMAL_FRACTION_PATH = (
    ROOT
    / "external/papers/Computational_analysis_gut_microbiota_drug_metabolism_code/PROXIMAL2_Supplementary/sup_data/supplementary_data_1e.csv"
)


def _canonicalize_key(value: object) -> str:
    text = str(value or "").strip().lower()
    return re.sub(r"[^a-z0-9]+", "", text)


def _collect_proximal_product_annotations(output_root: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for summary_path in sorted(output_root.glob("*/*_output_summary.csv")):
        frame = pd.read_csv(summary_path, low_memory=False)
        if frame.empty:
            continue
        drug_name = str(frame.iloc[0].get("DrugName") or "").strip()
        match_key = _canonicalize_key(drug_name)
        if not match_key:
            continue
        product_ids = []
        for value in frame.get("ProductInChiKey", pd.Series(dtype=object)).fillna("").astype(str):
            item = value.strip()
            if item and item not in product_ids:
                product_ids.append(item)
        ec_numbers = []
        for value in frame.get("EC", pd.Series(dtype=object)).fillna("").astype(str):
            item = value.strip()
            if item and item not in ec_numbers:
                ec_numbers.append(item)
        reaction_centers = []
        for value in frame.get("ReactionCenter", pd.Series(dtype=object)).fillna("").astype(str):
            item = value.strip()
            if item and item not in reaction_centers:
                reaction_centers.append(item)
        rows.append(
            {
                "biotransform_match_key": match_key,
                "experimental_biotransform_drug_name": drug_name,
                "experimental_biotransform_product_count": int(len(product_ids)),
                "experimental_biotransform_product_ids": ";".join(product_ids[:16]),
                "experimental_biotransform_ec_numbers": ";".join(ec_numbers[:16]),
                "experimental_biotransform_reaction_centers": ";".join(reaction_centers[:16]),
                "experimental_biotransform_support_rows": int(len(frame)),
                "experimental_biotransform_summary_path": str(summary_path.relative_to(ROOT)),
            }
        )
    return pd.DataFrame(rows)


def _collect_proximal_fraction_annotations(fraction_path: Path) -> pd.DataFrame:
    frame = pd.read_csv(fraction_path, low_memory=False)
    if frame.empty:
        return pd.DataFrame(columns=["biotransform_match_key", "experimental_biotransform_fraction_in_gut"])
    frame["biotransform_match_key"] = frame["DrugName"].map(_canonicalize_key)
    grouped = (
        frame.groupby("biotransform_match_key", as_index=False)
        .agg(
            experimental_biotransform_fraction_in_gut=("FractionInGut", "max"),
            experimental_biotransform_primary_product_id=("ProductID", "first"),
            experimental_biotransform_primary_product_name=("ProductName", "first"),
            experimental_biotransform_drugbank_id=("DrugID", "first"),
        )
    )
    return grouped


def build_augmented_predictions(
    *,
    input_path: Path,
    output_path: Path,
    summary_path: Path,
    proximal_output_root: Path,
    proximal_fraction_path: Path,
) -> dict[str, object]:
    predictions = pd.read_csv(input_path, low_memory=False).copy()
    predictions["biotransform_match_key"] = predictions["chemical_name"].map(_canonicalize_key)

    product_annotations = _collect_proximal_product_annotations(proximal_output_root)
    fraction_annotations = _collect_proximal_fraction_annotations(proximal_fraction_path)

    augmented = predictions.merge(product_annotations, on="biotransform_match_key", how="left")
    augmented = augmented.merge(fraction_annotations, on="biotransform_match_key", how="left")
    augmented.drop(columns=["biotransform_match_key"], inplace=True)

    if "experimental_biotransform_product_count" in augmented.columns:
        augmented["experimental_biotransform_product_count"] = (
            pd.to_numeric(augmented["experimental_biotransform_product_count"], errors="coerce").fillna(0).astype(int)
        )
    for column in [
        "experimental_biotransform_product_ids",
        "experimental_biotransform_ec_numbers",
        "experimental_biotransform_reaction_centers",
        "experimental_biotransform_summary_path",
        "experimental_biotransform_primary_product_id",
        "experimental_biotransform_primary_product_name",
        "experimental_biotransform_drugbank_id",
    ]:
        if column in augmented.columns:
            augmented[column] = augmented[column].fillna("").astype(str)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    augmented.to_csv(output_path, index=False)

    matched = augmented[augmented["experimental_biotransform_product_count"].gt(0)].copy()
    matched_drugs = (
        matched.loc[:, ["prestwick_id", "chemical_name", "experimental_biotransform_product_count"]]
        .drop_duplicates()
        .sort_values(["experimental_biotransform_product_count", "prestwick_id"], ascending=[False, True])
    )
    summary = {
        "input_path": str(input_path),
        "output_path": str(output_path),
        "n_rows": int(len(augmented)),
        "n_rows_with_product_annotations": int(matched.shape[0]),
        "n_drugs_with_product_annotations": int(matched_drugs.shape[0]),
        "matched_drugs": matched_drugs.head(25).to_dict(orient="records"),
        "proximal_output_root": str(proximal_output_root),
        "proximal_fraction_path": str(proximal_fraction_path),
    }
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Build an additive Step 3 BioTransformer/PROXIMAL-style experimental predictions copy.")
    parser.add_argument("--input-path", default=DEFAULT_INPUT_PATH, type=Path)
    parser.add_argument("--output-path", default=DEFAULT_OUTPUT_PATH, type=Path)
    parser.add_argument("--summary-path", default=DEFAULT_SUMMARY_PATH, type=Path)
    parser.add_argument("--proximal-output-root", default=DEFAULT_PROXIMAL_OUTPUT_ROOT, type=Path)
    parser.add_argument("--proximal-fraction-path", default=DEFAULT_PROXIMAL_FRACTION_PATH, type=Path)
    args = parser.parse_args()
    summary = build_augmented_predictions(
        input_path=args.input_path,
        output_path=args.output_path,
        summary_path=args.summary_path,
        proximal_output_root=args.proximal_output_root,
        proximal_fraction_path=args.proximal_fraction_path,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
