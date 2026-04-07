from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PANEL_PATH = ROOT / "data/processed/step1/step1_microbe_table.csv"
DEFAULT_OUTPUT_DIR = ROOT / "data/processed/food_tcm"
DEFAULT_MICROBETCM_CATALOG = Path("/tmp/microbetcm_microbe.json")
DEFAULT_MICROBETCM_ASSOC = Path("/tmp/microbetcm_all_data.json")

HIGH_FREQUENCY_MICROBES = [
    "Akkermansia muciniphila",
    "Bacteroides vulgatus",
    "Bacteroides uniformis",
    "Bacteroides thetaiotaomicron",
    "Bifidobacterium adolescentis",
    "Bilophila wadsworthia",
    "Blautia obeum",
    "Collinsella aerofaciens",
    "Eggerthella lenta",
    "Faecalibacterium prausnitzii",
    "Fusobacterium nucleatum",
    "Lactobacillus acidophilus",
    "Lactobacillus gasseri",
    "Parabacteroides distasonis",
    "Prevotella copri",
    "Roseburia hominis",
    "Roseburia intestinalis",
    "Ruminococcus bromii",
    "Ruminococcus gnavus",
]


def _normalize_taxon(text: object) -> str:
    if text is None or pd.isna(text):
        return ""
    normalized = str(text).strip().lower().replace("_", " ")
    normalized = re.sub(r"\([^)]*\)", "", normalized)
    replacements = {
        "para bacteroides": "parabacteroides",
        "parabacteroides": "parabacteroides",
        "blautiaobeum": "blautia obeum",
        "rosenburia": "roseburia",
    }
    for source, target in replacements.items():
        normalized = normalized.replace(source, target)
    normalized = re.sub(r"[^a-z0-9 ]+", " ", normalized)
    normalized = re.sub(r"\s+", " ", normalized).strip()
    return normalized


def _split_taxa_field(text: object) -> list[str]:
    if text is None or pd.isna(text):
        return []
    tokens: list[str] = []
    for part in re.split(r"[;/]", str(text)):
        token = _normalize_taxon(part)
        if token:
            tokens.append(token)
    return tokens


def _clean_entity(text: object) -> str:
    if text is None or pd.isna(text):
        return ""
    value = str(text).strip()
    if value.upper() == "NA":
        return ""
    return value


def _base_species_label(species_label: str) -> str:
    return re.sub(r"\s+\([^)]*\)$", "", species_label).strip()


def _build_panel(panel_path: Path) -> pd.DataFrame:
    panel = pd.read_csv(panel_path, low_memory=False)
    panel = panel[
        ["nt_code", "species_label", "species_name", "genus", "phylum"]
    ].drop_duplicates(subset=["nt_code"])
    panel["base_species_label"] = panel["species_label"].map(_base_species_label)
    panel["species_norm"] = panel["species_label"].map(_normalize_taxon)
    panel["base_species_norm"] = panel["base_species_label"].map(_normalize_taxon)
    panel["genus_norm"] = panel["genus"].map(_normalize_taxon)
    panel["phylum_norm"] = panel["phylum"].map(_normalize_taxon)
    return panel.sort_values(["species_label", "nt_code"], na_position="last").reset_index(drop=True)


def _load_json_frame(path: Path) -> pd.DataFrame:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return pd.DataFrame(payload)


def _best_catalog_match(panel_row: pd.Series, catalog: pd.DataFrame) -> dict[str, str]:
    species_hits = catalog[catalog["name_norm"] == panel_row["species_norm"]]
    if not species_hits.empty:
        match = species_hits.iloc[0]
        return {
            "catalog_match_type": "species_exact",
            "catalog_match_name": match["English_Name"],
            "catalog_match_rank": match["Rank"],
            "catalog_match_id": match["MicrobeTCM_ID"],
            "catalog_match_confidence": "high",
        }
    base_species_hits = catalog[catalog["name_norm"] == panel_row["base_species_norm"]]
    if not base_species_hits.empty:
        match = base_species_hits.iloc[0]
        return {
            "catalog_match_type": "base_species_exact",
            "catalog_match_name": match["English_Name"],
            "catalog_match_rank": match["Rank"],
            "catalog_match_id": match["MicrobeTCM_ID"],
            "catalog_match_confidence": "high",
        }
    genus_hits = catalog[catalog["name_norm"] == panel_row["genus_norm"]]
    if not genus_hits.empty:
        match = genus_hits.iloc[0]
        return {
            "catalog_match_type": "genus_proxy",
            "catalog_match_name": match["English_Name"],
            "catalog_match_rank": match["Rank"],
            "catalog_match_id": match["MicrobeTCM_ID"],
            "catalog_match_confidence": "medium",
        }
    return {
        "catalog_match_type": "no_match",
        "catalog_match_name": "",
        "catalog_match_rank": "",
        "catalog_match_id": "",
        "catalog_match_confidence": "none",
    }


def _choose_tcm_entity(row: pd.Series) -> tuple[str, str]:
    ordered = [
        ("formula", _clean_entity(row.get("Formula", ""))),
        ("single_herb", _clean_entity(row.get("Herb", ""))),
        ("component", _clean_entity(row.get("Ingredient", ""))),
        ("technology", _clean_entity(row.get("Technology", ""))),
        ("acupoint", _clean_entity(row.get("Acupoint", ""))),
    ]
    for entity_type, value in ordered:
        if value:
            return value, entity_type
    return "", ""


def _join_examples(values: list[str], limit: int = 3) -> str:
    counter = Counter(value for value in values if value)
    if not counter:
        return ""
    return "; ".join(value for value, _ in counter.most_common(limit))


def _build_microbetcm_tables(
    panel: pd.DataFrame,
    catalog_path: Path,
    assoc_path: Path,
    output_dir: Path,
) -> dict[str, object]:
    catalog = _load_json_frame(catalog_path)
    assoc = _load_json_frame(assoc_path)
    catalog["name_norm"] = catalog["English_Name"].map(_normalize_taxon)

    relation_rows: list[dict[str, object]] = []
    for assoc_index, assoc_row in assoc.iterrows():
        microbe_tokens = set(_split_taxa_field(assoc_row.get("Microbe", "")))
        genus_tokens = set(_split_taxa_field(assoc_row.get("Genus", "")))
        all_tokens = microbe_tokens | genus_tokens
        if not all_tokens:
            continue

        for _, panel_row in panel.iterrows():
            match_type = ""
            matched_taxon = ""
            matched_field = ""
            if panel_row["species_norm"] and panel_row["species_norm"] in microbe_tokens:
                match_type = "species_exact"
                matched_taxon = panel_row["species_label"]
                matched_field = "Microbe"
            elif panel_row["base_species_norm"] and panel_row["base_species_norm"] in microbe_tokens:
                match_type = "base_species_exact"
                matched_taxon = panel_row["base_species_label"]
                matched_field = "Microbe"
            elif panel_row["genus_norm"] and panel_row["genus_norm"] in genus_tokens:
                match_type = "genus_proxy"
                matched_taxon = panel_row["genus"]
                matched_field = "Genus"
            elif panel_row["genus_norm"] and panel_row["genus_norm"] in microbe_tokens:
                match_type = "genus_proxy"
                matched_taxon = panel_row["genus"]
                matched_field = "Microbe"

            if not match_type:
                continue

            tcm_name, tcm_type = _choose_tcm_entity(assoc_row)
            relation_rows.append(
                {
                    "nt_code": panel_row["nt_code"],
                    "species_label": panel_row["species_label"],
                    "species_name": panel_row["species_name"],
                    "genus": panel_row["genus"],
                    "phylum": panel_row["phylum"],
                    "match_type": match_type,
                    "matched_taxon": matched_taxon,
                    "matched_taxon_field": matched_field,
                    "microbe_field_raw": _clean_entity(assoc_row.get("Microbe", "")),
                    "genus_field_raw": _clean_entity(assoc_row.get("Genus", "")),
                    "disease_context": _clean_entity(assoc_row.get("Disease", "")),
                    "sample_context": _clean_entity(assoc_row.get("Sample", "")),
                    "method": _clean_entity(assoc_row.get("Method", "")),
                    "tcm_name": tcm_name,
                    "tcm_type": tcm_type,
                    "herb": _clean_entity(assoc_row.get("Herb", "")),
                    "formula": _clean_entity(assoc_row.get("Formula", "")),
                    "ingredient": _clean_entity(assoc_row.get("Ingredient", "")),
                    "technology": _clean_entity(assoc_row.get("Technology", "")),
                    "acupoint": _clean_entity(assoc_row.get("Acupoint", "")),
                    "relation_direction": _clean_entity(assoc_row.get("Correlation_with_microbe", "")),
                    "reference_pmid": str(assoc_row.get("Reference_PMID", "")).strip(),
                    "published_year": str(assoc_row.get("Published_Year", "")).strip(),
                    "source_record_id": f"MicrobeTCM_ASSOC_{assoc_index + 1}",
                    "source_name": "MicrobeTCM",
                    "source_url": "https://www.microbetcm.com/microbetcm/static/AllDataForFigure.json",
                    "curation_status": "first_pass_database_match",
                    "notes": "",
                }
            )

    relation_frame = pd.DataFrame(relation_rows)
    relation_frame = relation_frame.sort_values(
        ["species_label", "match_type", "disease_context", "reference_pmid", "source_record_id"],
        na_position="last",
    ).reset_index(drop=True)

    hit_rows: list[dict[str, object]] = []
    for _, panel_row in panel.iterrows():
        catalog_match = _best_catalog_match(panel_row, catalog)
        matched_relations = relation_frame[relation_frame["nt_code"] == panel_row["nt_code"]]
        hit_rows.append(
            {
                "nt_code": panel_row["nt_code"],
                "species_label": panel_row["species_label"],
                "base_species_label": panel_row["base_species_label"],
                "genus": panel_row["genus"],
                "phylum": panel_row["phylum"],
                **catalog_match,
                "association_row_count": int(len(matched_relations)),
                "association_species_exact_count": int(
                    matched_relations["match_type"].isin(["species_exact", "base_species_exact"]).sum()
                ),
                "association_genus_proxy_count": int((matched_relations["match_type"] == "genus_proxy").sum()),
                "example_direction": _join_examples(matched_relations["relation_direction"].tolist(), limit=2),
                "example_tcm_entities": _join_examples(matched_relations["tcm_name"].tolist(), limit=3),
                "example_pmids": _join_examples(matched_relations["reference_pmid"].tolist(), limit=3),
                "curation_status": (
                    "matched_associations"
                    if len(matched_relations) > 0
                    else ("catalog_only_match" if catalog_match["catalog_match_type"] != "no_match" else "no_match")
                ),
                "notes": "",
            }
        )

    hit_frame = pd.DataFrame(hit_rows).sort_values(
        ["association_row_count", "species_label"], ascending=[False, True]
    )

    hits_path = output_dir / "microbetcm_first_pass_hits.csv"
    relations_path = output_dir / "microbetcm_relation_first_pass.csv"
    summary_path = output_dir / "microbetcm_first_pass_summary.json"

    hit_frame.to_csv(hits_path, index=False)
    relation_frame.to_csv(relations_path, index=False)
    summary = {
        "catalog_path": str(catalog_path),
        "association_path": str(assoc_path),
        "hits_path": str(hits_path),
        "relations_path": str(relations_path),
        "n_panel_microbes": int(len(panel)),
        "n_catalog_matches": int((hit_frame["catalog_match_type"] != "no_match").sum()),
        "n_microbes_with_associations": int((hit_frame["association_row_count"] > 0).sum()),
        "n_relation_rows": int(len(relation_frame)),
    }
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def _build_fgmdi_table(panel: pd.DataFrame, output_dir: Path) -> dict[str, object]:
    grouped_rows: list[dict[str, object]] = []
    for species_label in HIGH_FREQUENCY_MICROBES:
        subset = panel[panel["species_label"] == species_label]
        if subset.empty:
            continue
        first = subset.iloc[0]
        genus = next((value for value in subset["genus"].tolist() if isinstance(value, str) and value.strip()), "")
        phylum = next((value for value in subset["phylum"].tolist() if isinstance(value, str) and value.strip()), "")
        grouped_rows.append(
            {
                "priority_rank": HIGH_FREQUENCY_MICROBES.index(species_label) + 1,
                "species_label": species_label,
                "representative_nt_code": first["nt_code"],
                "genus": genus,
                "phylum": phylum,
            }
        )
    high_frequency_panel = pd.DataFrame(grouped_rows)
    rows: list[dict[str, object]] = []
    for _, row in high_frequency_panel.sort_values("priority_rank").iterrows():
        if row["species_label"] == "Akkermansia muciniphila":
            rows.append(
                {
                    "priority_rank": int(row["priority_rank"]),
                    "representative_nt_code": row["representative_nt_code"],
                    "species_label": row["species_label"],
                    "genus": row["genus"],
                    "phylum": row["phylum"],
                    "food_name": "Berberine",
                    "food_category": "plant_compound / herbal_medicine",
                    "food_component": "berberine",
                    "relation_direction": "increase",
                    "relation_scope": "article_level_example",
                    "evidence_type": "FGMDI_article_intro_example",
                    "host_context": "mouse acute alcohol liver injury model",
                    "disease_context": "alcohol-related liver injury",
                    "source_name": "FGMDI article + cited source",
                    "source_url": "https://doi.org/10.1016/j.fbio.2024.104091",
                    "source_record_id": "FGMDI_intro_example_berberine_akkermansia",
                    "pmid_or_doi": "PMID:32790968",
                    "curation_status": "seeded_from_accessible_official_example",
                    "notes": (
                        "FGMDI row-level public export was not reachable from accessible official endpoints; "
                        "this row is seeded from the official FGMDI article introduction and its cited source."
                    ),
                }
            )
        else:
            rows.append(
                {
                    "priority_rank": int(row["priority_rank"]),
                    "representative_nt_code": row["representative_nt_code"],
                    "species_label": row["species_label"],
                    "genus": row["genus"],
                    "phylum": row["phylum"],
                    "food_name": "",
                    "food_category": "",
                    "food_component": "",
                    "relation_direction": "",
                    "relation_scope": "",
                    "evidence_type": "",
                    "host_context": "",
                    "disease_context": "",
                    "source_name": "FGMDI",
                    "source_url": "https://doi.org/10.1016/j.fbio.2024.104091",
                    "source_record_id": "",
                    "pmid_or_doi": "10.1016/j.fbio.2024.104091",
                    "curation_status": "await_public_row_level_export",
                    "notes": (
                        "The accessible FGMDI official pages provide article-level descriptions but not a public "
                        "row-level downloadable table for this species."
                    ),
                }
            )

    fgmdi_frame = pd.DataFrame(rows)
    fgmdi_path = output_dir / "fgmdi_high_frequency_food_first_pass.csv"
    fgmdi_frame.to_csv(fgmdi_path, index=False)
    return {
        "fgmdi_path": str(fgmdi_path),
        "n_high_frequency_microbes": int(len(fgmdi_frame)),
        "n_seeded_rows": int((fgmdi_frame["curation_status"] == "seeded_from_accessible_official_example").sum()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build first-pass MicrobeTCM and FGMDI lookup tables for the current 83-microbe panel."
    )
    parser.add_argument("--panel-path", type=Path, default=DEFAULT_PANEL_PATH)
    parser.add_argument("--microbetcm-catalog", type=Path, default=DEFAULT_MICROBETCM_CATALOG)
    parser.add_argument("--microbetcm-associations", type=Path, default=DEFAULT_MICROBETCM_ASSOC)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    panel = _build_panel(args.panel_path)

    microbetcm_summary = _build_microbetcm_tables(
        panel=panel,
        catalog_path=args.microbetcm_catalog,
        assoc_path=args.microbetcm_associations,
        output_dir=args.output_dir,
    )
    fgmdi_summary = _build_fgmdi_table(panel=panel, output_dir=args.output_dir)

    combined_summary = {
        "panel_path": str(args.panel_path),
        "n_panel_microbes": int(len(panel)),
        "microbetcm": microbetcm_summary,
        "fgmdi": fgmdi_summary,
    }
    summary_path = args.output_dir / "food_tcm_first_pass_build_summary.json"
    summary_path.write_text(json.dumps(combined_summary, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
