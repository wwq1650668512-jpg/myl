from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MICROBE_TABLE = ROOT / "data/processed/step1/step1_microbe_table.csv"

MICROBE_COLUMN_ALIASES = [
    "nt_code",
    "microbe_id",
    "microbe_name",
    "microbe_label",
    "species_label",
    "species_name",
    "species",
    "taxon",
    "taxon_name",
    "genus",
]
ABUNDANCE_COLUMN_ALIASES = [
    "abundance",
    "relative_abundance",
    "relative_frequency",
    "weight",
    "biomass",
]
SAMPLE_COLUMN_ALIASES = [
    "sample_id",
    "sample",
    "sample_name",
    "subject_id",
]


def _canonicalize_key(value: object) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def _pick_column(frame: pd.DataFrame, aliases: list[str], explicit: str | None) -> str:
    if explicit:
        if explicit not in frame.columns:
            raise ValueError(f"Column {explicit!r} not found. Available: {frame.columns.tolist()}")
        return explicit

    normalized = {_canonicalize_key(column): column for column in frame.columns}
    for alias in aliases:
        key = _canonicalize_key(alias)
        if key in normalized:
            return normalized[key]
    raise ValueError(f"Could not match any alias from {aliases}; available columns: {frame.columns.tolist()}")


def _build_lookup(microbe_table: pd.DataFrame) -> tuple[dict[str, list[dict[str, object]]], dict[str, list[dict[str, object]]]]:
    exact_lookup: dict[str, dict[str, dict[str, object]]] = {}
    genus_lookup: dict[str, dict[str, dict[str, object]]] = {}
    for _, row in microbe_table.iterrows():
        record = {
            "nt_code": str(row["nt_code"]),
            "microbe_label": str(row.get("microbe_label", "") or ""),
            "species_label": str(row.get("species_label", "") or ""),
            "species_name": str(row.get("species_name", "") or ""),
            "genus": str(row.get("genus", "") or ""),
        }
        for value in [
            record["nt_code"],
            record["microbe_label"],
            record["species_label"],
            record["species_name"],
        ]:
            key = _canonicalize_key(value)
            if key:
                exact_lookup.setdefault(key, {})[record["nt_code"]] = record
        genus_key = _canonicalize_key(record["genus"])
        if genus_key:
            genus_lookup.setdefault(genus_key, {})[record["nt_code"]] = record
    exact_result = {key: list(value.values()) for key, value in exact_lookup.items()}
    genus_result = {key: list(value.values()) for key, value in genus_lookup.items()}
    return exact_result, genus_result


def main() -> None:
    parser = argparse.ArgumentParser(description="Map a cohort abundance table onto the Step 3 community panel format.")
    parser.add_argument("--input-table", required=True, type=Path, help="Raw cohort abundance table.")
    parser.add_argument("--output-path", required=True, type=Path, help="Output community table path.")
    parser.add_argument(
        "--microbe-table",
        default=DEFAULT_MICROBE_TABLE,
        type=Path,
        help="Step 1/Step 3 microbe panel table used as the mapping target.",
    )
    parser.add_argument("--sample-id", default=None, help="Optional sample id to select from a multi-sample table.")
    parser.add_argument("--sample-column", default=None, help="Explicit sample column name.")
    parser.add_argument("--microbe-column", default=None, help="Explicit taxon column name.")
    parser.add_argument("--abundance-column", default=None, help="Explicit abundance column name.")
    args = parser.parse_args()

    raw = pd.read_csv(args.input_table, low_memory=False)
    microbe_table = pd.read_csv(args.microbe_table, low_memory=False).drop_duplicates(subset=["nt_code"]).reset_index(drop=True)

    microbe_column = _pick_column(raw, MICROBE_COLUMN_ALIASES, args.microbe_column)
    abundance_column = _pick_column(raw, ABUNDANCE_COLUMN_ALIASES, args.abundance_column)

    sample_column = None
    sample_candidates = {_canonicalize_key(column): column for column in raw.columns}
    if args.sample_column is not None:
        sample_column = _pick_column(raw, SAMPLE_COLUMN_ALIASES, args.sample_column)
    else:
        for alias in SAMPLE_COLUMN_ALIASES:
            key = _canonicalize_key(alias)
            if key in sample_candidates:
                sample_column = sample_candidates[key]
                break

    selected_sample_id = args.sample_id
    work = raw.copy()
    if sample_column is not None:
        sample_values = work[sample_column].dropna().astype(str).map(str.strip)
        sample_values = sample_values[sample_values.ne("")]
        if selected_sample_id is None:
            unique_values = sorted(sample_values.unique().tolist())
            if len(unique_values) > 1:
                raise ValueError(
                    f"Detected multiple samples in column {sample_column!r}; provide --sample-id. Examples: {unique_values[:10]}"
                )
            selected_sample_id = unique_values[0] if unique_values else None
        if selected_sample_id is not None:
            work = work[work[sample_column].astype(str).str.strip() == str(selected_sample_id).strip()].copy()

    exact_lookup, genus_lookup = _build_lookup(microbe_table)

    mapped_rows: list[dict[str, object]] = []
    unmapped_rows: list[dict[str, object]] = []
    for _, row in work.iterrows():
        source_taxon = str(row[microbe_column]).strip()
        abundance = pd.to_numeric(pd.Series([row[abundance_column]]), errors="coerce").iloc[0]
        if not source_taxon or pd.isna(abundance) or float(abundance) <= 0:
            continue

        key = _canonicalize_key(source_taxon)
        matches = exact_lookup.get(key, [])
        mapping_mode = "exact"
        if not matches:
            matches = genus_lookup.get(key, [])
            mapping_mode = "genus_split"

        if not matches:
            unmapped_rows.append(
                {
                    "source_taxon": source_taxon,
                    "abundance": float(abundance),
                    "sample_id": selected_sample_id,
                }
            )
            continue

        distributed_abundance = float(abundance) / len(matches)
        for match in matches:
            mapped_rows.append(
                {
                    "sample_id": selected_sample_id,
                    "source_taxon": source_taxon,
                    "nt_code": match["nt_code"],
                    "microbe_label": match["microbe_label"],
                    "species_label": match["species_label"],
                    "genus": match["genus"],
                    "mapping_mode": mapping_mode,
                    "matched_panel_count": int(len(matches)),
                    "input_abundance": float(abundance),
                    "distributed_abundance": float(distributed_abundance),
                }
            )

    if not mapped_rows:
        raise RuntimeError("No abundance rows could be mapped onto the Step 3 panel.")

    mapping_report = pd.DataFrame(mapped_rows)
    community = (
        mapping_report.groupby(["nt_code", "microbe_label", "species_label"], as_index=False)["distributed_abundance"]
        .sum()
        .rename(columns={"distributed_abundance": "abundance"})
    )
    community["abundance"] = community["abundance"] / community["abundance"].sum()

    output_path = args.output_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_path = output_path.with_name(f"{output_path.stem}_mapping_report.csv")
    unmapped_path = output_path.with_name(f"{output_path.stem}_unmapped.csv")
    summary_path = output_path.with_name(f"{output_path.stem}_summary.json")

    community.to_csv(output_path, index=False)
    mapping_report.to_csv(report_path, index=False)
    pd.DataFrame(unmapped_rows).to_csv(unmapped_path, index=False)

    summary = {
        "input_table": str(args.input_table),
        "output_path": str(output_path),
        "mapping_report_path": str(report_path),
        "unmapped_path": str(unmapped_path),
        "sample_id": selected_sample_id,
        "n_input_rows": int(len(work)),
        "n_mapped_rows": int(len(mapping_report)),
        "n_unmapped_rows": int(len(unmapped_rows)),
        "n_output_microbes": int(len(community)),
        "mapping_modes": {
            str(key): int(value)
            for key, value in mapping_report["mapping_mode"].value_counts().to_dict().items()
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
