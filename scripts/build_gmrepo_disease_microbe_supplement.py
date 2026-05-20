from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.utils.text import canonicalize_key as _canonicalize_key
from gut_drug_microbiome.utils.text import normalize_whitespace as _clean_text


HEALTH_MESH_ID = "D006262"
SOURCE_URL = "https://gmrepo.humangut.info"

# Map common GMrepo disease terms to project canonical names to reduce duplicate disease labels.
DISEASE_TERM_ALIAS_MAP = {
    _canonicalize_key("Colorectal Neoplasms"): "结直肠癌（CRC）",
    _canonicalize_key("Crohn Disease"): "克罗恩病（CD）",
    _canonicalize_key("Ulcerative Colitis"): "溃疡性结肠炎（UC）",
    _canonicalize_key("Inflammatory Bowel Diseases"): "炎症性肠病（IBD）",
    _canonicalize_key("Irritable Bowel Syndrome"): "肠易激综合征（IBS）",
    _canonicalize_key("Diarrhea"): "腹泻（Diarrhea）",
    _canonicalize_key("Constipation"): "便秘（Constipation）",
}

CONFIDENCE_RANK = {"low": 1, "medium": 2, "high": 3}


def _post_json(
    base_url: str,
    endpoint: str,
    payload: dict[str, object],
    *,
    timeout_seconds: float,
    retries: int,
    sleep_seconds: float,
    verbose: bool,
) -> dict[str, object]:
    url = f"{base_url.rstrip('/')}/{endpoint.strip('/')}/"
    body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url=url, data=body, headers=headers, method="POST")
            with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
                data = response.read().decode("utf-8")
            if sleep_seconds > 0:
                time.sleep(sleep_seconds)
            return json.loads(data)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last_error = exc
            if verbose:
                print(f"[warn] {endpoint} attempt {attempt}/{retries} failed: {exc}", file=sys.stderr)
            if attempt < retries:
                time.sleep(min(1.5 * attempt, 5.0))
    raise RuntimeError(f"Failed to fetch GMrepo endpoint: {endpoint}") from last_error


def _canonical_disease_name(mesh_id: str, term: str) -> str:
    mapped = DISEASE_TERM_ALIAS_MAP.get(_canonicalize_key(term))
    if mapped:
        return mapped
    cleaned_term = _clean_text(term)
    if cleaned_term:
        return f"{cleaned_term}（MeSH:{mesh_id}）"
    return f"MeSH:{mesh_id}"


def _infer_taxon_level(raw_value: object) -> str:
    text = _clean_text(raw_value).lower()
    if text in {"species", "genus", "family", "order", "class", "phylum"}:
        return text
    # GMrepo markers are usually species/genus; fallback to species if binomial name appears.
    return "species" if " " in text else "genus"


def _genus_hint(scientific_name: str) -> str:
    tokens = scientific_name.split()
    return tokens[0] if tokens else ""


def _confidence_from_marker(abs_lda: float, nrproj: int, conflict: int) -> str:
    if conflict:
        return "low"
    if nrproj >= 5 and abs_lda >= 3.0:
        return "high"
    if nrproj >= 2 and abs_lda >= 2.0:
        return "medium"
    return "low"


def _map_detail_row(
    row: dict[str, object],
    *,
    disease_mesh_id: str,
    disease_name: str,
    mesh1: str,
    mesh2: str,
) -> dict[str, object] | None:
    scientific_name = _clean_text(row.get("scientific_name"))
    if not scientific_name:
        return None
    try:
        lda = float(row.get("LDA", 0.0))
    except (TypeError, ValueError):
        return None
    try:
        nrproj = int(float(row.get("nrproj", 1)))
    except (TypeError, ValueError):
        nrproj = 1
    try:
        conflict = int(float(row.get("conflict", 0)))
    except (TypeError, ValueError):
        conflict = 0

    # From GMrepo comparison detail page semantics:
    # LDA < 0 => phenotype1 enriched; LDA > 0 => phenotype2 enriched.
    disease_is_mesh1 = mesh1 == disease_mesh_id
    disease_enriched = (lda < 0 and disease_is_mesh1) or (lda > 0 and not disease_is_mesh1)
    disease_effect = "increase" if disease_enriched else "decrease"
    desired_effect = "inhibit" if disease_enriched else "promote"
    microbe_role = "risk" if disease_enriched else "protective"
    confidence = _confidence_from_marker(abs(lda), nrproj, conflict)
    taxon_level = _infer_taxon_level(row.get("taxon_rank_level", ""))
    mechanism_note = (
        "GMrepo phenotype comparison marker; "
        f"mesh_pair={mesh1}|{mesh2}; "
        f"project_id={_clean_text(row.get('project_id')) or 'unknown'}; "
        f"LDA={lda:.3f}; nrproj={nrproj}; conflict={conflict}"
    )

    return {
        "source_sheet": "gmrepo_health_vs_disease",
        "disease_name": disease_name,
        "microbe_name_raw": scientific_name,
        "microbe_name_cn": "",
        "microbe_key": _canonicalize_key(scientific_name),
        "genus_hint": _genus_hint(scientific_name),
        "taxon_level": taxon_level,
        "abundance_change_raw": f"LDA={lda:.3f}",
        "abundance_shift": disease_effect,
        "microbe_role_raw": "disease_enriched" if disease_enriched else "health_enriched",
        "microbe_role_in_disease": microbe_role,
        "disease_effect_on_microbe": disease_effect,
        "desired_step1_effect": desired_effect,
        "relation_confidence": confidence,
        "mechanism_note": mechanism_note,
        "source_database": "GMrepo v3",
        "source_url": SOURCE_URL,
        "mesh_id_disease": disease_mesh_id,
        "mesh_id_health": HEALTH_MESH_ID,
        "mesh_id_comparison_1": mesh1,
        "mesh_id_comparison_2": mesh2,
        "lda_score": float(lda),
        "marker_nr_projects": int(nrproj),
        "marker_conflict_flag": int(conflict),
    }


def _load_or_fetch_details(
    *,
    cache_dir: Path | None,
    refresh_cache: bool,
    base_url: str,
    timeout_seconds: float,
    retries: int,
    sleep_seconds: float,
    verbose: bool,
    mesh1: str,
    mesh2: str,
) -> dict[str, object]:
    cache_path = None
    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path = cache_dir / f"{mesh1}__{mesh2}.json"
        if cache_path.exists() and not refresh_cache:
            with cache_path.open("r", encoding="utf-8") as handle:
                return json.load(handle)

    payload = {"mesh_id1": mesh1, "mesh_id2": mesh2}
    data = _post_json(
        base_url=base_url,
        endpoint="getPhenotypeComparisonsDetails",
        payload=payload,
        timeout_seconds=timeout_seconds,
        retries=retries,
        sleep_seconds=sleep_seconds,
        verbose=verbose,
    )
    if cache_path is not None:
        with cache_path.open("w", encoding="utf-8") as handle:
            json.dump(data, handle, ensure_ascii=False, indent=2)
    return data


def build_gmrepo_supplement(
    *,
    output_path: Path,
    summary_path: Path,
    base_url: str,
    min_abs_lda: float,
    min_nrproj: int,
    include_conflict: bool,
    max_comparisons: int | None,
    timeout_seconds: float,
    retries: int,
    sleep_seconds: float,
    cache_dir: Path | None,
    refresh_cache: bool,
    verbose: bool,
) -> dict[str, object]:
    phenotypes_payload = _post_json(
        base_url=base_url,
        endpoint="get_all_phenotypes",
        payload={},
        timeout_seconds=timeout_seconds,
        retries=retries,
        sleep_seconds=sleep_seconds,
        verbose=verbose,
    )
    term_by_mesh = {
        _clean_text(item.get("disease")): _clean_text(item.get("term"))
        for item in phenotypes_payload.get("phenotypes", [])
        if _clean_text(item.get("disease"))
    }

    comparisons_payload = _post_json(
        base_url=base_url,
        endpoint="get_all_phenotype_comparisons",
        payload={},
        timeout_seconds=timeout_seconds,
        retries=retries,
        sleep_seconds=sleep_seconds,
        verbose=verbose,
    )
    all_comparisons = [
        item
        for item in comparisons_payload.get("data", [])
        if isinstance(item, dict)
        and _clean_text(item.get("phenotype1"))
        and _clean_text(item.get("phenotype2"))
        and HEALTH_MESH_ID in {_clean_text(item.get("phenotype1")), _clean_text(item.get("phenotype2"))}
    ]
    if max_comparisons is not None and max_comparisons >= 0:
        all_comparisons = all_comparisons[:max_comparisons]

    records: list[dict[str, object]] = []
    failed_pairs: list[dict[str, str]] = []
    for index, comparison in enumerate(all_comparisons, start=1):
        mesh1 = _clean_text(comparison.get("phenotype1"))
        mesh2 = _clean_text(comparison.get("phenotype2"))
        disease_mesh = mesh2 if mesh1 == HEALTH_MESH_ID else mesh1
        disease_name = _canonical_disease_name(disease_mesh, term_by_mesh.get(disease_mesh, ""))
        if verbose:
            print(f"[info] [{index}/{len(all_comparisons)}] fetching pair {mesh1} vs {mesh2}", file=sys.stderr)
        try:
            detail_payload = _load_or_fetch_details(
                cache_dir=cache_dir,
                refresh_cache=refresh_cache,
                base_url=base_url,
                timeout_seconds=timeout_seconds,
                retries=retries,
                sleep_seconds=sleep_seconds,
                verbose=verbose,
                mesh1=mesh1,
                mesh2=mesh2,
            )
        except RuntimeError:
            failed_pairs.append({"mesh1": mesh1, "mesh2": mesh2})
            continue

        for row in detail_payload.get("alldata", []):
            mapped = _map_detail_row(
                row,
                disease_mesh_id=disease_mesh,
                disease_name=disease_name,
                mesh1=mesh1,
                mesh2=mesh2,
            )
            if mapped is None:
                continue
            if abs(float(mapped["lda_score"])) < min_abs_lda:
                continue
            if int(mapped["marker_nr_projects"]) < min_nrproj:
                continue
            if not include_conflict and int(mapped["marker_conflict_flag"]) != 0:
                continue
            records.append(mapped)

    if records:
        frame = pd.DataFrame.from_records(records)
        frame["confidence_rank"] = frame["relation_confidence"].map(CONFIDENCE_RANK).fillna(0)
        frame["abs_lda_score"] = pd.to_numeric(frame["lda_score"], errors="coerce").abs().fillna(0.0)
        frame.sort_values(
            ["confidence_rank", "marker_nr_projects", "abs_lda_score"],
            ascending=[False, False, False],
            inplace=True,
        )
        frame.drop_duplicates(
            subset=["disease_name", "microbe_key", "taxon_level", "desired_step1_effect"],
            keep="first",
            inplace=True,
        )
        frame.drop(columns=["confidence_rank", "abs_lda_score"], inplace=True)
        frame.insert(0, "reference_id", [f"GMR_{idx + 1:05d}" for idx in range(len(frame))])
    else:
        frame = pd.DataFrame(
            columns=[
                "reference_id",
                "source_sheet",
                "disease_name",
                "microbe_name_raw",
                "microbe_name_cn",
                "microbe_key",
                "genus_hint",
                "taxon_level",
                "abundance_change_raw",
                "abundance_shift",
                "microbe_role_raw",
                "microbe_role_in_disease",
                "disease_effect_on_microbe",
                "desired_step1_effect",
                "relation_confidence",
                "mechanism_note",
                "source_database",
                "source_url",
                "mesh_id_disease",
                "mesh_id_health",
                "mesh_id_comparison_1",
                "mesh_id_comparison_2",
                "lda_score",
                "marker_nr_projects",
                "marker_conflict_flag",
            ]
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_path, index=False)
    summary = {
        "source_database": "GMrepo v3",
        "source_url": SOURCE_URL,
        "base_url": base_url,
        "output_path": str(output_path),
        "n_rows": int(len(frame)),
        "n_diseases": int(frame["disease_name"].nunique()) if not frame.empty else 0,
        "n_comparisons_requested": int(len(all_comparisons)),
        "n_failed_pairs": int(len(failed_pairs)),
        "failed_pairs": failed_pairs[:20],
        "min_abs_lda": float(min_abs_lda),
        "min_nrproj": int(min_nrproj),
        "include_conflict": bool(include_conflict),
        "desired_step1_effect_counts": (
            {str(key): int(value) for key, value in frame["desired_step1_effect"].value_counts(dropna=False).to_dict().items()}
            if not frame.empty
            else {}
        ),
        "taxon_level_counts": (
            {str(key): int(value) for key, value in frame["taxon_level"].value_counts(dropna=False).to_dict().items()}
            if not frame.empty
            else {}
        ),
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Build disease-microbe supplement table from GMrepo phenotype comparisons.")
    parser.add_argument(
        "--output-path",
        default=ROOT / "data/reference/disease_microbe_gmrepo_supplement.csv",
        type=Path,
        help="Output CSV path for mapped disease-microbe relations.",
    )
    parser.add_argument(
        "--summary-path",
        default=ROOT / "data/reference/disease_microbe_gmrepo_supplement.summary.json",
        type=Path,
        help="Output JSON summary path.",
    )
    parser.add_argument(
        "--base-url",
        default="https://gmrepo.humangut.info/api",
        type=str,
        help="GMrepo API base URL.",
    )
    parser.add_argument("--min-abs-lda", default=2.0, type=float, help="Minimum absolute LDA to keep a marker.")
    parser.add_argument(
        "--min-nrproj",
        default=2,
        type=int,
        help="Minimum number of supporting projects (`nrproj`) to keep a marker.",
    )
    parser.add_argument(
        "--include-conflict",
        action="store_true",
        help="Keep rows marked as conflict by GMrepo (default drops conflict=1).",
    )
    parser.add_argument(
        "--max-comparisons",
        default=None,
        type=int,
        help="Optional cap on health-vs-disease comparisons to fetch (for quick dry-runs).",
    )
    parser.add_argument("--timeout-seconds", default=60.0, type=float, help="Per-request timeout.")
    parser.add_argument("--retries", default=3, type=int, help="Retry count per request.")
    parser.add_argument("--sleep-seconds", default=0.1, type=float, help="Sleep between API calls.")
    parser.add_argument(
        "--cache-dir",
        default=ROOT / "data/cache/gmrepo_phenotype_comparison",
        type=Path,
        help="Directory to cache per-comparison JSON payloads.",
    )
    parser.add_argument("--refresh-cache", action="store_true", help="Ignore cache files and refetch.")
    parser.add_argument("--verbose", action="store_true", help="Print progress and warnings to stderr.")
    args = parser.parse_args()

    summary = build_gmrepo_supplement(
        output_path=args.output_path,
        summary_path=args.summary_path,
        base_url=args.base_url,
        min_abs_lda=args.min_abs_lda,
        min_nrproj=args.min_nrproj,
        include_conflict=args.include_conflict,
        max_comparisons=args.max_comparisons,
        timeout_seconds=args.timeout_seconds,
        retries=args.retries,
        sleep_seconds=args.sleep_seconds,
        cache_dir=args.cache_dir,
        refresh_cache=args.refresh_cache,
        verbose=args.verbose,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
