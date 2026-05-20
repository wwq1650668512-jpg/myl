from __future__ import annotations

import json
import re
import ssl
import time
import urllib.parse
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd

from gut_drug_microbiome.step2.enzyme_prior import DEFAULT_ENZYME_CATALOG_PATH
from gut_drug_microbiome.utils.text import canonicalize_key
from gut_drug_microbiome.utils.text import normalize_whitespace


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_UNIPROT_SPARQL_ENDPOINT = "https://sparql.uniprot.org/sparql"
DEFAULT_UNIPROT_EVIDENCE_OUTPUT_PATH = ROOT / "data/reference/step2_uniprot_enzyme_candidate_evidence.csv"
DEFAULT_UNIPROT_RAW_OUTPUT_PATH = ROOT / "data/reference/step2_uniprot_protein_enzyme_hits.csv"
DEFAULT_UNIPROT_UNRESOLVED_OUTPUT_PATH = ROOT / "data/reference/step2_uniprot_unresolved_taxa.csv"
DEFAULT_UNIPROT_SUMMARY_OUTPUT_PATH = ROOT / "data/reference/step2_uniprot_enzyme_fetch_summary.json"

STARTER_ENZYME_MATCH_RULES: dict[str, dict[str, tuple[str, ...]]] = {
    "ENZ001": {
        "ec_exact": ("3.2.1.31",),
        "name_keywords": ("beta-glucuronidase", "beta glucuronidase", "beta-d-glucuronidase", "beta-d glucuronidase"),
    },
    "ENZ002": {
        "ec_prefix": ("3.1.6.",),
        "name_keywords": ("sulfatase", "sulphatase", "arylsulfatase"),
    },
    "ENZ003": {
        "name_keywords": ("azoreductase", "azo reductase"),
    },
    "ENZ004": {
        "name_keywords": ("nitroreductase", "nitro reductase"),
    },
    "ENZ005": {
        "ec_prefix": ("3.1.1.",),
        "name_keywords": ("carboxylesterase", "esterase", "lipase"),
    },
    "ENZ006": {
        "ec_prefix": ("3.5.1.",),
        "name_keywords": ("amidase",),
    },
    "ENZ007": {
        "ec_exact": ("3.2.1.21",),
        "name_keywords": ("beta-glucosidase", "beta glucosidase", "beta-glucoside glucosidase"),
    },
    "ENZ008": {
        "ec_exact": ("3.2.1.40",),
        "name_keywords": ("alpha-l-rhamnosidase", "alpha rhamnosidase", "alpha-l-rhamnoside", "rhamnosidase"),
    },
    "ENZ009": {
        "ec_exact": ("3.2.1.23",),
        "name_keywords": ("beta-galactosidase", "beta galactosidase"),
    },
    "ENZ010": {
        "ec_exact": ("3.5.1.24",),
        "name_keywords": ("bile salt hydrolase", "choloylglycine hydrolase"),
    },
    "ENZ011": {
        "name_keywords": ("o-demethylase", "o demethylase", "demethylase"),
    },
    "ENZ012": {
        "name_keywords": ("dehydroxylase", "dehydroxylation"),
    },
    "ENZ013": {
        "name_keywords": ("deacetylase", "deacetylation", "acetylesterase"),
    },
    "ENZ014": {
        "ec_prefix": ("3.1.3.",),
        "name_keywords": ("phosphatase",),
    },
    "ENZ015": {
        "name_keywords": (
            "mucin glycosidase",
            "sialidase",
            "neuraminidase",
            "fucosidase",
            "hexosaminidase",
            "galactosidase",
        ),
    },
    "ENZ016": {
        "ec_prefix": ("4.2.2.",),
        "name_keywords": ("polysaccharide lyase", "pectate lyase", "alginate lyase", "hyaluronate lyase"),
    },
}


def _clean_species_text(value: object) -> str:
    text = normalize_whitespace(value)
    if not text:
        return ""
    text = re.sub(r"\([^)]*\)", " ", text)
    text = re.sub(r"\bDSM\b.*$", " ", text, flags=re.IGNORECASE)
    text = re.sub(r"\bATCC\b.*$", " ", text, flags=re.IGNORECASE)
    text = re.sub(r"\bJCM\b.*$", " ", text, flags=re.IGNORECASE)
    text = re.sub(r"\bstrain\b.*$", " ", text, flags=re.IGNORECASE)
    text = normalize_whitespace(text)
    return text


def candidate_taxon_names_for_microbe(row: pd.Series) -> list[str]:
    values = [
        row.get("species_name"),
        row.get("species_label"),
        row.get("species"),
        row.get("microbe_label"),
    ]
    candidates: list[str] = []
    for value in values:
        cleaned = _clean_species_text(value)
        if cleaned:
            candidates.append(cleaned)
            tokens = cleaned.split()
            if len(tokens) >= 2:
                candidates.append(" ".join(tokens[:2]))
    genus = normalize_whitespace(row.get("genus"))
    species = _clean_species_text(row.get("species"))
    if genus and species:
        species_tokens = species.split()
        if len(species_tokens) >= 2:
            candidates.append(f"{genus} {species_tokens[1]}")
    deduped: list[str] = []
    seen: set[str] = set()
    for candidate in candidates:
        key = canonicalize_key(candidate)
        if key and key not in seen:
            seen.add(key)
            deduped.append(candidate)
    return deduped


def build_uniprot_taxonomy_query(scientific_names: list[str]) -> str:
    values = " ".join(f'"{name}"' for name in scientific_names)
    return f"""
PREFIX up: <http://purl.uniprot.org/core/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

SELECT DISTINCT ?query ?taxon ?name ?matchField ?obsolete ?replacement ?replacementName
FROM <http://sparql.uniprot.org/taxonomy>
WHERE {{
  VALUES ?query {{ {values} }}
  ?taxon a up:Taxon .
  {{
    ?taxon up:scientificName ?name .
    BIND("scientificName" AS ?matchField)
  }}
  UNION
  {{
    ?taxon up:synonym ?name .
    BIND("synonym" AS ?matchField)
  }}
  UNION
  {{
    ?taxon up:otherName ?name .
    BIND("otherName" AS ?matchField)
  }}
  OPTIONAL {{ ?taxon up:obsolete ?obsolete . }}
  OPTIONAL {{
    ?taxon up:replacedBy ?replacement .
    OPTIONAL {{ ?replacement up:scientificName ?replacementName . }}
  }}
  FILTER(LCASE(STR(?name)) = LCASE(STR(?query)))
}}
""".strip()


def build_uniprot_enzyme_query(taxon_uri: str, reviewed_only: bool = False) -> str:
    reviewed_filter = "  ?protein up:reviewed true .\n" if reviewed_only else ""
    return f"""
PREFIX up: <http://purl.uniprot.org/core/>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

SELECT DISTINCT ?protein ?accession ?organism ?organismName ?reviewed ?proteinName ?enzyme ?ecNumber
WHERE {{
  ?protein a up:Protein ;
           up:organism ?organism .
{reviewed_filter}  ?organism rdfs:subClassOf <{taxon_uri}> ;
            up:scientificName ?organismName .
  OPTIONAL {{ ?protein up:reviewed ?reviewed . }}
  OPTIONAL {{
    ?protein up:recommendedName ?recommendedName .
    ?recommendedName up:fullName ?proteinName .
  }}
  ?protein (up:enzyme | up:domain/up:enzyme | up:component/up:enzyme) ?enzyme .
  BIND(REPLACE(STR(?protein), "^.*/", "") AS ?accession)
  BIND(REPLACE(STR(?enzyme), "^.*/", "") AS ?ecNumber)
}}
ORDER BY ?accession ?ecNumber
""".strip()


def _run_sparql_query(
    query: str,
    endpoint: str = DEFAULT_UNIPROT_SPARQL_ENDPOINT,
    timeout_seconds: int = 90,
    max_retries: int = 5,
    retry_sleep_seconds: float = 2.0,
) -> list[dict[str, str]]:
    params = urllib.parse.urlencode({"query": query, "format": "json"})
    payload: dict[str, object] | None = None
    last_error: Exception | None = None
    for attempt in range(max_retries):
        request = urllib.request.Request(
            f"{endpoint}?{params}",
            headers={
                "Accept": "application/sparql-results+json",
                "User-Agent": "gut-drug-microbiome/0.1.0",
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
                payload = json.loads(response.read().decode("utf-8"))
            break
        except (urllib.error.URLError, TimeoutError, ssl.SSLError) as exc:
            last_error = exc
            if attempt >= max_retries - 1:
                raise
            time.sleep(retry_sleep_seconds * (attempt + 1))
    if payload is None:
        if last_error is not None:
            raise last_error
        return []
    rows: list[dict[str, str]] = []
    for binding in payload.get("results", {}).get("bindings", []):
        row = {key: value.get("value", "") for key, value in binding.items()}
        rows.append(row)
    return rows


def _resolve_uniprot_taxon(row: pd.Series, endpoint: str, timeout_seconds: int) -> dict[str, Any]:
    candidates = candidate_taxon_names_for_microbe(row)
    if not candidates:
        return {"status": "no_query_name", "query_names": "", "matched_taxon_uri": "", "matched_taxon_name": ""}
    query = build_uniprot_taxonomy_query(candidates)
    hits = _run_sparql_query(query, endpoint=endpoint, timeout_seconds=timeout_seconds)
    if not hits:
        return {
            "status": "not_found",
            "query_names": ";".join(candidates),
            "matched_taxon_uri": "",
            "matched_taxon_name": "",
        }
    query_rank = {canonicalize_key(name): index for index, name in enumerate(candidates)}
    match_field_rank = {
        "scientificname": 0,
        "synonym": 1,
        "othername": 2,
    }
    hits = sorted(
        hits,
        key=lambda item: query_rank.get(canonicalize_key(item.get("query")), 999),
    )
    hits = sorted(
        hits,
        key=lambda item: (
            query_rank.get(canonicalize_key(item.get("query")), 999),
            match_field_rank.get(canonicalize_key(item.get("matchField")), 99),
            1 if str(item.get("obsolete", "")).strip().lower() == "true" else 0,
        ),
    )
    best = hits[0]
    replacement_uri = normalize_whitespace(best.get("replacement"))
    replacement_name = normalize_whitespace(best.get("replacementName"))
    matched_taxon_uri = replacement_uri or normalize_whitespace(best.get("taxon"))
    matched_name = replacement_name or normalize_whitespace(best.get("name"))
    strain_match_level = "same_species"
    if len(matched_name.split()) >= 3:
        strain_match_level = "same_strain"
    return {
        "status": "resolved",
        "query_names": ";".join(candidates),
        "matched_taxon_uri": matched_taxon_uri,
        "matched_taxon_name": matched_name,
        "matched_query_name": normalize_whitespace(best.get("query")),
        "matched_name_field": normalize_whitespace(best.get("matchField")),
        "obsolete_taxon_replaced": bool(replacement_uri),
        "strain_match_level": strain_match_level,
    }


def _match_uniprot_record_to_project_enzymes(ec_number: str, protein_name: str) -> list[tuple[str, str]]:
    matches: list[tuple[str, str]] = []
    ec_number = normalize_whitespace(ec_number)
    protein_key = normalize_whitespace(protein_name).lower()
    for enzyme_id, rule in STARTER_ENZYME_MATCH_RULES.items():
        for ec_exact in rule.get("ec_exact", ()):
            if ec_number == ec_exact:
                matches.append((enzyme_id, "exact_ec"))
                break
        else:
            for ec_prefix in rule.get("ec_prefix", ()):
                if ec_number.startswith(ec_prefix):
                    matches.append((enzyme_id, "ec_prefix"))
                    break
            else:
                for keyword in rule.get("name_keywords", ()):
                    if keyword in protein_key:
                        matches.append((enzyme_id, "name_keyword"))
                        break
    deduped: list[tuple[str, str]] = []
    seen: set[str] = set()
    for enzyme_id, match_mode in matches:
        if enzyme_id not in seen:
            seen.add(enzyme_id)
            deduped.append((enzyme_id, match_mode))
    return deduped


def _presence_call_from_match_modes(match_modes: list[str], reviewed_hits: int) -> str:
    if "exact_ec" in match_modes and reviewed_hits > 0:
        return "curated_present"
    if "exact_ec" in match_modes or "ec_prefix" in match_modes:
        return "likely_present"
    return "weak_prior"


def _write_uniprot_outputs(
    aggregated_rows: list[dict[str, Any]],
    raw_rows: list[dict[str, Any]],
    unresolved_rows: list[dict[str, Any]],
    evidence_output_path: str | Path,
    raw_output_path: str | Path,
    unresolved_output_path: str | Path,
    summary_output_path: str | Path | None,
    summary: dict[str, Any],
) -> None:
    evidence_output_path = Path(evidence_output_path)
    raw_output_path = Path(raw_output_path)
    unresolved_output_path = Path(unresolved_output_path)
    for path in [evidence_output_path, raw_output_path, unresolved_output_path]:
        path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(aggregated_rows).to_csv(evidence_output_path, index=False)
    pd.DataFrame(raw_rows).to_csv(raw_output_path, index=False)
    pd.DataFrame(unresolved_rows).to_csv(unresolved_output_path, index=False)
    if summary_output_path is not None:
        summary_output_path = Path(summary_output_path)
        summary_output_path.parent.mkdir(parents=True, exist_ok=True)
        summary_output_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")


def fetch_uniprot_enzyme_candidates(
    microbe_table_path: str | Path | pd.DataFrame,
    enzyme_catalog_path: str | Path | pd.DataFrame = DEFAULT_ENZYME_CATALOG_PATH,
    evidence_output_path: str | Path = DEFAULT_UNIPROT_EVIDENCE_OUTPUT_PATH,
    raw_output_path: str | Path = DEFAULT_UNIPROT_RAW_OUTPUT_PATH,
    unresolved_output_path: str | Path = DEFAULT_UNIPROT_UNRESOLVED_OUTPUT_PATH,
    summary_output_path: str | Path | None = DEFAULT_UNIPROT_SUMMARY_OUTPUT_PATH,
    endpoint: str = DEFAULT_UNIPROT_SPARQL_ENDPOINT,
    reviewed_only: bool = False,
    limit_microbes: int | None = None,
    sleep_seconds: float = 0.2,
    timeout_seconds: int = 90,
    checkpoint_every: int = 10,
) -> dict[str, Any]:
    if isinstance(microbe_table_path, pd.DataFrame):
        microbes = microbe_table_path.copy()
        microbe_table_repr = "<dataframe>"
    else:
        microbes = pd.read_csv(microbe_table_path, low_memory=False)
        microbe_table_repr = str(microbe_table_path)
    microbes = microbes.drop_duplicates(subset=["nt_code"]).reset_index(drop=True)
    if limit_microbes is not None:
        microbes = microbes.head(limit_microbes).copy()

    if isinstance(enzyme_catalog_path, pd.DataFrame):
        catalog = enzyme_catalog_path.copy()
    else:
        catalog = pd.read_csv(enzyme_catalog_path, low_memory=False)
    enzyme_name_lookup = {
        str(row["enzyme_id"]): str(row["enzyme_name"])
        for _, row in catalog.iterrows()
        if normalize_whitespace(row.get("enzyme_id"))
    }

    raw_rows: list[dict[str, Any]] = []
    unresolved_rows: list[dict[str, Any]] = []
    aggregated_rows: list[dict[str, Any]] = []

    for index, (_, microbe) in enumerate(microbes.iterrows(), start=1):
        resolved = _resolve_uniprot_taxon(microbe, endpoint=endpoint, timeout_seconds=timeout_seconds)
        if resolved.get("status") != "resolved":
            unresolved_rows.append(
                {
                    "nt_code": microbe.get("nt_code"),
                    "microbe_label": microbe.get("microbe_label"),
                    "species_label": microbe.get("species_label"),
                    "species_name": microbe.get("species_name"),
                    "query_names": resolved.get("query_names", ""),
                    "status": resolved.get("status", "not_found"),
                }
            )
        else:
            query = build_uniprot_enzyme_query(str(resolved["matched_taxon_uri"]), reviewed_only=reviewed_only)
            hits = _run_sparql_query(query, endpoint=endpoint, timeout_seconds=timeout_seconds)
            if sleep_seconds > 0:
                time.sleep(sleep_seconds)

            if not hits:
                unresolved_rows.append(
                    {
                        "nt_code": microbe.get("nt_code"),
                        "microbe_label": microbe.get("microbe_label"),
                        "species_label": microbe.get("species_label"),
                        "species_name": microbe.get("species_name"),
                        "query_names": resolved.get("query_names", ""),
                        "status": "no_protein_hits",
                    }
                )
            else:
                microbe_match_rows: list[dict[str, Any]] = []
                for hit in hits:
                    protein_name = normalize_whitespace(hit.get("proteinName"))
                    ec_number = normalize_whitespace(hit.get("ecNumber"))
                    match_pairs = _match_uniprot_record_to_project_enzymes(ec_number=ec_number, protein_name=protein_name)
                    if not match_pairs:
                        continue
                    reviewed_value = str(hit.get("reviewed", "")).strip().lower() == "true"
                    for enzyme_id, match_mode in match_pairs:
                        row = {
                            "nt_code": microbe.get("nt_code"),
                            "microbe_label": microbe.get("microbe_label"),
                            "species_label": microbe.get("species_label"),
                            "species_name": microbe.get("species_name"),
                            "species": microbe.get("species"),
                            "strain": microbe.get("strain"),
                            "genus": microbe.get("genus"),
                            "family": microbe.get("family"),
                            "phylum": microbe.get("phylum"),
                            "matched_taxon_uri": resolved.get("matched_taxon_uri"),
                            "matched_taxon_name": resolved.get("matched_taxon_name"),
                            "matched_query_name": resolved.get("matched_query_name"),
                            "strain_match_level": resolved.get("strain_match_level", "same_species"),
                            "uniprot_accession": hit.get("accession"),
                            "uniprot_protein_uri": hit.get("protein"),
                            "uniprot_organism_name": hit.get("organismName"),
                            "uniprot_protein_name": protein_name,
                            "ec_number": ec_number,
                            "reviewed": reviewed_value,
                            "enzyme_id": enzyme_id,
                            "enzyme_name": enzyme_name_lookup.get(enzyme_id, enzyme_id),
                            "annotation_match_mode": match_mode,
                        }
                        raw_rows.append(row)
                        microbe_match_rows.append(row)

                if not microbe_match_rows:
                    unresolved_rows.append(
                        {
                            "nt_code": microbe.get("nt_code"),
                            "microbe_label": microbe.get("microbe_label"),
                            "species_label": microbe.get("species_label"),
                            "species_name": microbe.get("species_name"),
                            "query_names": resolved.get("query_names", ""),
                            "status": "no_project_enzyme_match",
                        }
                    )
                else:
                    match_frame = pd.DataFrame(microbe_match_rows)
                    for enzyme_id, enzyme_hits in match_frame.groupby("enzyme_id"):
                        reviewed_hits = int(enzyme_hits["reviewed"].fillna(False).astype(bool).sum())
                        match_modes = enzyme_hits["annotation_match_mode"].astype(str).tolist()
                        presence_call = _presence_call_from_match_modes(match_modes, reviewed_hits)
                        aggregated_rows.append(
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
                                "enzyme_id": enzyme_id,
                                "enzyme_name": enzyme_name_lookup.get(enzyme_id, enzyme_id),
                                "presence_call": presence_call,
                                "evidence_scope": "species_genome",
                                "evidence_source": "uniprot_sparql_automated",
                                "source_database": "UniProt",
                                "literature_citation": "",
                                "pmid": "",
                                "doi": "",
                                "genome_accession": "",
                                "strain_match_level": resolved.get("strain_match_level", "same_species"),
                                "evidence_note": (
                                    f"Automated UniProt SPARQL match from taxon {resolved.get('matched_taxon_name', '')}; "
                                    f"{len(enzyme_hits)} protein hit(s), {reviewed_hits} reviewed."
                                ).strip(),
                                "curation_status": "needs_curation",
                                "matched_taxon_uri": resolved.get("matched_taxon_uri"),
                                "matched_taxon_name": resolved.get("matched_taxon_name"),
                                "matched_query_name": resolved.get("matched_query_name"),
                                "annotation_match_modes": ";".join(sorted(set(match_modes))),
                                "uniprot_accessions": ";".join(sorted(set(enzyme_hits["uniprot_accession"].astype(str)))),
                                "uniprot_protein_names": ";".join(
                                    sorted(set(name for name in enzyme_hits["uniprot_protein_name"].astype(str) if name))
                                ),
                                "ec_numbers": ";".join(sorted(set(value for value in enzyme_hits["ec_number"].astype(str) if value))),
                                "n_matching_proteins": int(len(enzyme_hits)),
                                "n_reviewed_proteins": reviewed_hits,
                            }
                        )

        if checkpoint_every > 0 and index % checkpoint_every == 0:
            checkpoint_summary = {
                "microbe_table_path": microbe_table_repr,
                "endpoint": endpoint,
                "reviewed_only": bool(reviewed_only),
                "n_input_microbes": int(len(microbes)),
                "n_processed_microbes": int(index),
                "n_resolved_taxa": int(len({row["nt_code"] for row in aggregated_rows})),
                "n_candidate_evidence_rows": int(len(aggregated_rows)),
                "n_raw_protein_hits": int(len(raw_rows)),
                "n_unresolved_microbes": int(len({row["nt_code"] for row in unresolved_rows})),
                "evidence_output_path": str(evidence_output_path),
                "raw_output_path": str(raw_output_path),
                "unresolved_output_path": str(unresolved_output_path),
            }
            _write_uniprot_outputs(
                aggregated_rows=aggregated_rows,
                raw_rows=raw_rows,
                unresolved_rows=unresolved_rows,
                evidence_output_path=evidence_output_path,
                raw_output_path=raw_output_path,
                unresolved_output_path=unresolved_output_path,
                summary_output_path=summary_output_path,
                summary=checkpoint_summary,
            )

    raw_frame = pd.DataFrame(raw_rows)
    aggregated = pd.DataFrame(aggregated_rows)
    unresolved = pd.DataFrame(unresolved_rows)

    summary = {
        "microbe_table_path": microbe_table_repr,
        "endpoint": endpoint,
        "reviewed_only": bool(reviewed_only),
        "n_input_microbes": int(len(microbes)),
        "n_processed_microbes": int(len(microbes)),
        "n_resolved_taxa": int(aggregated["nt_code"].nunique()) if not aggregated.empty else 0,
        "n_candidate_evidence_rows": int(len(aggregated)),
        "n_raw_protein_hits": int(len(raw_frame)),
        "n_unresolved_microbes": int(unresolved["nt_code"].nunique()) if not unresolved.empty else 0,
        "evidence_output_path": str(evidence_output_path),
        "raw_output_path": str(raw_output_path),
        "unresolved_output_path": str(unresolved_output_path),
    }
    _write_uniprot_outputs(
        aggregated_rows=aggregated_rows,
        raw_rows=raw_rows,
        unresolved_rows=unresolved_rows,
        evidence_output_path=evidence_output_path,
        raw_output_path=raw_output_path,
        unresolved_output_path=unresolved_output_path,
        summary_output_path=summary_output_path,
        summary=summary,
    )
    if summary_output_path is not None:
        summary["summary_output_path"] = str(summary_output_path)
    return summary
