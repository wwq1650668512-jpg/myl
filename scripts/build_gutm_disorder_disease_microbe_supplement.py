from __future__ import annotations

import argparse
import http.client
import json
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.utils.text import canonicalize_key as _canonicalize_key
from gut_drug_microbiome.utils.text import normalize_whitespace as _clean_text


BASE_URL = "https://bio-computing.hrbmu.edu.cn/gutMDisorder"
SOURCE_URL = f"{BASE_URL}/browse.dhtml"
LITERATURE_HUMAN_PHENOTYPE_NODE_ID = "Literature-based associations_Human_Condition_Phenotype"
RAW_HUMAN_PHENOTYPE_NODE_ID = "Raw data-based associations_Human_Condition_Phenotype"
LITERATURE_SOURCE_SHEET = "gutm_disorder_human_literature_phenotype"
RAW_SOURCE_SHEET = "gutm_disorder_human_raw_phenotype"
AGGREGATED_SOURCE_SHEET = "gutm_disorder_human_phenotype"
SOURCE_DATABASE = "gutMDisorder v2.0"
SOURCE_SPECS = (
    {
        "branch_name": "literature_human_phenotype",
        "node_parent_id": LITERATURE_HUMAN_PHENOTYPE_NODE_ID,
        "root_node": "Literature-based associations",
        "source_sheet": LITERATURE_SOURCE_SHEET,
        "cache_nodes_file": "literature_human_phenotype_nodes.json",
        "cache_table_dir": "table_data_literature",
        "mode": "literature",
    },
    {
        "branch_name": "raw_human_phenotype",
        "node_parent_id": RAW_HUMAN_PHENOTYPE_NODE_ID,
        "root_node": "Raw data-based associations",
        "source_sheet": RAW_SOURCE_SHEET,
        "cache_nodes_file": "raw_human_phenotype_nodes.json",
        "cache_table_dir": "table_data_raw",
        "mode": "raw",
    },
)

PROJECT_DISEASE_PATTERNS: list[tuple[str, tuple[str, ...]]] = [
    ("肠易激综合征-便秘型（IBS-C）", ("constipatedirritablebowelsyndrome", "irritablebowelsyndromeconstipationt", "irritablebowelsyndromeconstipation", "肠易激综合征便秘型", "ibsc")),
    ("肠易激综合征-腹泻型（IBS-D）", ("diarrheapredominantirritablebowelsyndrome", "irritablebowelsyndromediarrhea", "肠易激综合征腹泻型", "ibsd")),
    ("克罗恩病（CD）", ("crohndisease", "crohnsdisease", "克罗恩病")),
    ("溃疡性结肠炎（UC）", ("ulcerativecolitis", "colitisulcerative", "溃疡性结肠炎")),
    ("炎症性肠病（IBD）", ("inflammatoryboweldiseases", "inflammatoryboweldisease", "炎症性肠病")),
    ("肠易激综合征（IBS）", ("irritablebowelsyndrome", "肠易激综合征", "ibs")),
    ("便秘（Constipation）", ("constipation", "便秘")),
    ("腹泻（Diarrhea）", ("diarrhea", "腹泻")),
    ("结直肠癌（CRC）", ("colorectalneoplasms", "colorectalcancer", "结直肠癌", "crc")),
    ("多发性硬化（MS）", ("multiplesclerosis", "多发性硬化")),
    ("系统性红斑狼疮（SLE）", ("lupuserythematosussystemic", "systemiclupuserythematosus", "系统性红斑狼疮")),
    ("类风湿关节炎（RA）", ("arthritisrheumatoid", "rheumatoidarthritis", "类风湿关节炎")),
    ("自身免疫性肝炎（AIH）", ("hepatitisautoimmune", "autoimmuneliverdisease", "autoimmunehepatitis", "自身免疫性肝炎")),
    ("自身免疫性甲状腺疾病（AITD）", ("hashimotodisease", "autoimmunethyroiddisease", "自身免疫性甲状腺疾病")),
]

CONTROL_LIKE_PATTERNS = (
    "health",
    "healthy",
    "control",
    "controls",
    "non",
    "no",
    "without",
    "normal",
    "typicaldevelopment",
    "tissuedonors",
    "uninfected",
    "uninfection",
    "unaffected",
    "nondiarrheal",
    "nondiarrhea",
    "nonconstipated",
    "healthyrelatives",
    "neurotypical",
)


def _get_json(
    base_url: str,
    endpoint: str,
    params: dict[str, object],
    *,
    timeout_seconds: float,
    retries: int,
    sleep_seconds: float,
    verbose: bool,
) -> object:
    query = urllib.parse.urlencode({key: value for key, value in params.items() if value is not None})
    url = f"{base_url.rstrip('/')}/{endpoint.lstrip('/')}"
    if query:
        url = f"{url}?{query}"
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout_seconds) as response:
                payload = response.read().decode("utf-8")
            if sleep_seconds > 0:
                time.sleep(sleep_seconds)
            return json.loads(payload)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, http.client.HTTPException, OSError) as exc:
            last_error = exc
            if verbose:
                print(f"[warn] GET {endpoint} attempt {attempt}/{retries} failed: {exc}", file=sys.stderr)
            if attempt < retries:
                time.sleep(min(1.5 * attempt, 5.0))
    curl_path = shutil.which("curl")
    if curl_path is not None:
        if verbose:
            print(f"[warn] falling back to curl for {endpoint}", file=sys.stderr)
        try:
            completed = subprocess.run(
                [curl_path, "-LsS", "--max-time", str(int(max(timeout_seconds, 20))), url],
                check=True,
                capture_output=True,
                text=True,
            )
            if sleep_seconds > 0:
                time.sleep(sleep_seconds)
            return json.loads(completed.stdout)
        except (subprocess.CalledProcessError, json.JSONDecodeError) as exc:
            last_error = exc
    raise RuntimeError(f"Failed to fetch gutMDisorder endpoint: {endpoint}") from last_error


def _load_or_fetch_json(
    *,
    cache_path: Path | None,
    refresh_cache: bool,
    base_url: str,
    endpoint: str,
    params: dict[str, object],
    timeout_seconds: float,
    retries: int,
    sleep_seconds: float,
    verbose: bool,
) -> object:
    if cache_path is not None and cache_path.exists() and not refresh_cache:
        return json.loads(cache_path.read_text(encoding="utf-8"))
    payload = _get_json(
        base_url=base_url,
        endpoint=endpoint,
        params=params,
        timeout_seconds=timeout_seconds,
        retries=retries,
        sleep_seconds=sleep_seconds,
        verbose=verbose,
    )
    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return payload


def _canonical_condition_key(value: object) -> str:
    return _canonicalize_key(_clean_text(value), keep_cjk=True)


def _is_control_like_condition(value: object) -> bool:
    key = _canonical_condition_key(value)
    if not key:
        return False
    if key.startswith("health") or key.endswith("health"):
        return True
    if key.startswith("control") or key.endswith("control") or "controls" in key:
        return True
    if key.startswith("non") or key.startswith("no") or key.startswith("without"):
        return True
    return any(pattern in key for pattern in CONTROL_LIKE_PATTERNS)


def _canonical_project_disease_name(value: object) -> str:
    text = _clean_text(value)
    if not text or _is_control_like_condition(text):
        return ""
    key = _canonical_condition_key(text)
    for canonical_name, patterns in PROJECT_DISEASE_PATTERNS:
        if any(pattern in key for pattern in patterns):
            return canonical_name
    return ""


def _parse_node_pair(node_label: str) -> tuple[str, str] | None:
    if "/" not in node_label:
        return None
    left, right = node_label.split("/", 1)
    left = _clean_text(left)
    right = _clean_text(right)
    if not left or not right:
        return None
    return left, right


def _infer_disease_side(condition1: str, condition2: str) -> tuple[str, int] | None:
    disease1 = _canonical_project_disease_name(condition1)
    disease2 = _canonical_project_disease_name(condition2)
    if disease1 and disease2:
        return None
    if disease1:
        return disease1, 1
    if disease2:
        return disease2, 2
    return None


def _condition1_enriched(alteration: object) -> bool | None:
    text = _clean_text(alteration).lower()
    if text in {"increase", "present"}:
        return True
    if text in {"decrease", "absent"}:
        return False
    return None


def _condition1_enriched_from_raw_lda(value: object) -> bool | None:
    try:
        score = float(value)
    except (TypeError, ValueError):
        return None
    if score > 0:
        return True
    if score < 0:
        return False
    return None


def _infer_taxon_level(raw_value: object, microbe_name: object) -> str:
    text = _clean_text(raw_value).lower()
    if text in {"species", "genus", "family", "order", "class", "phylum"}:
        return text
    microbe = _clean_text(microbe_name)
    return "species" if " " in microbe else ("genus" if microbe else "unknown")


def _genus_hint(microbe_name: str) -> str:
    tokens = microbe_name.split()
    return tokens[0] if tokens else ""


def _map_literature_row(
    row: dict[str, object],
    *,
    node_label: str,
    node_id: str,
) -> dict[str, object] | None:
    condition1 = _clean_text(row.get("Condition1"))
    condition2 = _clean_text(row.get("Condition2"))
    disease_assignment = _infer_disease_side(condition1, condition2)
    if disease_assignment is None:
        return None
    disease_name, disease_side = disease_assignment

    condition1_is_enriched = _condition1_enriched(row.get("Alteration"))
    if condition1_is_enriched is None:
        return None

    disease_enriched = (condition1_is_enriched and disease_side == 1) or (not condition1_is_enriched and disease_side == 2)
    disease_effect = "increase" if disease_enriched else "decrease"
    desired_effect = "inhibit" if disease_enriched else "promote"
    microbe_role = "risk" if disease_enriched else "protective"
    comparator = condition2 if disease_side == 1 else condition1
    microbe_name = _clean_text(row.get("GutMicrobe"))
    if not microbe_name:
        return None

    return {
        "source_sheet": LITERATURE_SOURCE_SHEET,
        "disease_name": disease_name,
        "microbe_name_raw": microbe_name,
        "microbe_name_cn": "",
        "microbe_key": _canonicalize_key(microbe_name, keep_cjk=True),
        "genus_hint": _genus_hint(microbe_name),
        "taxon_level": _infer_taxon_level(row.get("Classification"), microbe_name),
        "abundance_change_raw": _clean_text(row.get("Alteration")).lower(),
        "abundance_shift": disease_effect,
        "microbe_role_raw": "disease_enriched" if disease_enriched else "control_enriched",
        "microbe_role_in_disease": microbe_role,
        "disease_effect_on_microbe": disease_effect,
        "desired_step1_effect": desired_effect,
        "relation_confidence": "low",
        "mechanism_note": "",
        "source_database": SOURCE_DATABASE,
        "source_url": SOURCE_URL,
        "source_branch": "literature_human_phenotype",
        "source_pmid": _clean_text(row.get("PMID")),
        "source_project_number": "",
        "source_record_id": _clean_text(row.get("PMID")),
        "condition_comparator": comparator,
        "sequencing_technology": _clean_text(row.get("SequencingTechnology")),
        "node_label": node_label,
        "node_id": node_id,
        "condition1_label": condition1,
        "condition2_label": condition2,
        "condition1_id": _clean_text(row.get("Condition1ID")),
        "condition2_id": _clean_text(row.get("Condition2ID")),
        "raw_lda_score": "",
    }


def _map_raw_row(
    row: dict[str, object],
    *,
    node_label: str,
    node_id: str,
) -> dict[str, object] | None:
    condition1 = _clean_text(row.get("Condition1"))
    condition2 = _clean_text(row.get("Condition2"))
    disease_assignment = _infer_disease_side(condition1, condition2)
    if disease_assignment is None:
        return None
    disease_name, disease_side = disease_assignment

    condition1_is_enriched = _condition1_enriched_from_raw_lda(row.get("LDAscore"))
    if condition1_is_enriched is None:
        return None

    disease_enriched = (condition1_is_enriched and disease_side == 1) or (not condition1_is_enriched and disease_side == 2)
    disease_effect = "increase" if disease_enriched else "decrease"
    desired_effect = "inhibit" if disease_enriched else "promote"
    microbe_role = "risk" if disease_enriched else "protective"
    comparator = condition2 if disease_side == 1 else condition1
    microbe_name = _clean_text(row.get("GutMicrobe"))
    if not microbe_name:
        return None

    project_number = _clean_text(row.get("ProjectNumber"))
    run_data = _clean_text(row.get("RunData"))
    source_record_id = run_data or project_number or _clean_text(node_id)
    return {
        "source_sheet": RAW_SOURCE_SHEET,
        "disease_name": disease_name,
        "microbe_name_raw": microbe_name,
        "microbe_name_cn": "",
        "microbe_key": _canonicalize_key(microbe_name, keep_cjk=True),
        "genus_hint": _genus_hint(microbe_name),
        "taxon_level": _infer_taxon_level(row.get("Classification"), microbe_name),
        "abundance_change_raw": f"LDA={_clean_text(row.get('LDAscore'))}",
        "abundance_shift": disease_effect,
        "microbe_role_raw": "disease_enriched" if disease_enriched else "control_enriched",
        "microbe_role_in_disease": microbe_role,
        "disease_effect_on_microbe": disease_effect,
        "desired_step1_effect": desired_effect,
        "relation_confidence": "low",
        "mechanism_note": "",
        "source_database": SOURCE_DATABASE,
        "source_url": SOURCE_URL,
        "source_branch": "raw_human_phenotype",
        "source_pmid": "",
        "source_project_number": project_number,
        "source_record_id": source_record_id,
        "condition_comparator": comparator,
        "sequencing_technology": _clean_text(row.get("SequencingTechnology")),
        "node_label": node_label,
        "node_id": node_id,
        "condition1_label": condition1,
        "condition2_label": condition2,
        "condition1_id": _clean_text(row.get("Condition1ID")),
        "condition2_id": _clean_text(row.get("Condition2ID")),
        "raw_lda_score": _clean_text(row.get("LDAscore")),
    }


def _confidence_from_group(
    *,
    dominant_evidence_count: int,
    conflict_evidence_count: int,
    dominant_project_count: int,
    dominant_mean_abs_lda: float,
) -> str:
    if conflict_evidence_count > 0:
        return "low"
    if dominant_project_count >= 3 or dominant_evidence_count >= 4 or dominant_mean_abs_lda >= 4.0:
        return "high"
    if dominant_project_count >= 2 or dominant_evidence_count >= 2 or dominant_mean_abs_lda >= 2.5:
        return "medium"
    return "low"


def _aggregate_records(records: list[dict[str, object]]) -> pd.DataFrame:
    if not records:
        return pd.DataFrame()

    frame = pd.DataFrame(records)
    output_rows: list[dict[str, object]] = []
    grouped = frame.groupby(["disease_name", "microbe_key", "taxon_level"], dropna=False, sort=True)
    for (_, _, _), group in grouped:
        evidence_ids_by_direction = {
            direction: sorted({_clean_text(value) for value in subset["source_record_id"].tolist() if _clean_text(value)})
            for direction, subset in group.groupby("disease_effect_on_microbe", dropna=False)
        }
        increase_evidence_ids = evidence_ids_by_direction.get("increase", [])
        decrease_evidence_ids = evidence_ids_by_direction.get("decrease", [])
        direction_score = len(increase_evidence_ids) - len(decrease_evidence_ids)
        if direction_score == 0:
            continue

        dominant_effect = "increase" if direction_score > 0 else "decrease"
        dominant_rows = group[group["disease_effect_on_microbe"] == dominant_effect].copy()
        conflict_rows = group[group["disease_effect_on_microbe"] != dominant_effect].copy()
        dominant_evidence_ids = sorted({_clean_text(value) for value in dominant_rows["source_record_id"].tolist() if _clean_text(value)})
        conflict_evidence_ids = sorted({_clean_text(value) for value in conflict_rows["source_record_id"].tolist() if _clean_text(value)})
        dominant_pmids = sorted({_clean_text(value) for value in dominant_rows["source_pmid"].tolist() if _clean_text(value)})
        dominant_projects = sorted({_clean_text(value) for value in dominant_rows["source_project_number"].tolist() if _clean_text(value)})
        conflict_projects = sorted({_clean_text(value) for value in conflict_rows["source_project_number"].tolist() if _clean_text(value)})
        source_branches = sorted({_clean_text(value) for value in dominant_rows["source_branch"].tolist() if _clean_text(value)})
        comparator_values = sorted({_clean_text(value) for value in dominant_rows["condition_comparator"].tolist() if _clean_text(value)})
        technology_values = sorted({_clean_text(value) for value in dominant_rows["sequencing_technology"].tolist() if _clean_text(value)})
        node_values = sorted({_clean_text(value) for value in dominant_rows["node_label"].tolist() if _clean_text(value)})
        condition1_ids = sorted({_clean_text(value) for value in dominant_rows["condition1_id"].tolist() if _clean_text(value)})
        condition2_ids = sorted({_clean_text(value) for value in dominant_rows["condition2_id"].tolist() if _clean_text(value)})
        dominant_lda_values = pd.to_numeric(dominant_rows["raw_lda_score"], errors="coerce").abs().dropna()

        preferred_name = Counter(dominant_rows["microbe_name_raw"].tolist()).most_common(1)[0][0]
        confidence = _confidence_from_group(
            dominant_evidence_count=len(dominant_evidence_ids),
            conflict_evidence_count=len(conflict_evidence_ids),
            dominant_project_count=len(dominant_projects),
            dominant_mean_abs_lda=float(dominant_lda_values.mean()) if not dominant_lda_values.empty else 0.0,
        )
        output_rows.append(
            {
                "source_sheet": AGGREGATED_SOURCE_SHEET,
                "disease_name": dominant_rows.iloc[0]["disease_name"],
                "microbe_name_raw": preferred_name,
                "microbe_name_cn": "",
                "microbe_key": dominant_rows.iloc[0]["microbe_key"],
                "genus_hint": _genus_hint(preferred_name),
                "taxon_level": dominant_rows.iloc[0]["taxon_level"],
                "abundance_change_raw": (
                    f"{dominant_effect}; dominant_evidence={len(dominant_evidence_ids)}; "
                    f"dominant_projects={len(dominant_projects)}; conflicting_evidence={len(conflict_evidence_ids)}"
                ),
                "abundance_shift": dominant_effect,
                "microbe_role_raw": "disease_enriched" if dominant_effect == "increase" else "control_enriched",
                "microbe_role_in_disease": "risk" if dominant_effect == "increase" else "protective",
                "disease_effect_on_microbe": dominant_effect,
                "desired_step1_effect": "inhibit" if dominant_effect == "increase" else "promote",
                "relation_confidence": confidence,
                "mechanism_note": (
                    "gutMDisorder human phenotype aggregate; "
                    f"branches={'; '.join(source_branches) or 'unknown'}; "
                    f"dominant_pmids={','.join(dominant_pmids[:8]) or 'none'}; "
                    f"dominant_projects={','.join(dominant_projects[:8]) or 'none'}; "
                    f"comparators={'; '.join(comparator_values[:4]) or 'unspecified'}; "
                    f"sequencing={'; '.join(technology_values[:3]) or 'unspecified'}; "
                    f"conflicting_evidence={','.join(conflict_evidence_ids[:6]) or 'none'}"
                ),
                "source_database": SOURCE_DATABASE,
                "source_url": SOURCE_URL,
                "source_pmid_list": ";".join(dominant_pmids),
                "source_project_list": ";".join(dominant_projects),
                "dominant_evidence_count": int(len(dominant_evidence_ids)),
                "conflicting_evidence_count": int(len(conflict_evidence_ids)),
                "dominant_pmid_count": int(len(dominant_pmids)),
                "conflicting_project_count": int(len(conflict_projects)),
                "condition_comparator_list": ";".join(comparator_values),
                "node_label_list": ";".join(node_values[:8]),
                "condition1_id_list": ";".join(condition1_ids),
                "condition2_id_list": ";".join(condition2_ids),
                "source_branch_list": ";".join(source_branches),
            }
        )

    output = pd.DataFrame(output_rows)
    if output.empty:
        return output
    output.sort_values(["disease_name", "microbe_name_raw", "taxon_level"], inplace=True)
    output.reset_index(drop=True, inplace=True)
    output.insert(0, "reference_id", [f"GMD_{index + 1:05d}" for index in range(len(output))])
    return output


def build_gutm_disorder_supplement(
    *,
    output_path: Path,
    summary_path: Path,
    base_url: str,
    timeout_seconds: float,
    retries: int,
    sleep_seconds: float,
    cache_dir: Path | None,
    refresh_cache: bool,
    max_nodes: int | None,
    verbose: bool,
) -> dict[str, object]:
    candidate_nodes: list[dict[str, str]] = []
    candidate_node_counts_by_branch: dict[str, int] = {}
    for spec in SOURCE_SPECS:
        node_cache_path = cache_dir / spec["cache_nodes_file"] if cache_dir is not None else None
        raw_nodes = _load_or_fetch_json(
            cache_path=node_cache_path,
            refresh_cache=refresh_cache,
            base_url=base_url,
            endpoint="browse/tree-nodes.dhtml",
            params={"id": spec["node_parent_id"]},
            timeout_seconds=timeout_seconds,
            retries=retries,
            sleep_seconds=sleep_seconds,
            verbose=verbose,
        )
        if not isinstance(raw_nodes, list):
            raise RuntimeError(f"Unexpected gutMDisorder node payload for {spec['branch_name']}.")

        branch_candidate_count = 0
        for item in raw_nodes:
            if not isinstance(item, dict):
                continue
            node_label = _clean_text(item.get("name"))
            node_id = _clean_text(item.get("id"))
            pair = _parse_node_pair(node_label)
            if not node_label or not node_id or pair is None:
                continue
            disease_assignment = _infer_disease_side(pair[0], pair[1])
            if disease_assignment is None:
                continue
            candidate_nodes.append(
                {
                    "node_label": node_label,
                    "node_id": node_id,
                    "root_node": str(spec["root_node"]),
                    "source_sheet": str(spec["source_sheet"]),
                    "branch_name": str(spec["branch_name"]),
                    "mode": str(spec["mode"]),
                    "cache_table_dir": str(spec["cache_table_dir"]),
                }
            )
            branch_candidate_count += 1
        candidate_node_counts_by_branch[str(spec["branch_name"])] = branch_candidate_count

    if max_nodes is not None and max_nodes >= 0:
        candidate_nodes = candidate_nodes[:max_nodes]

    mapped_rows: list[dict[str, object]] = []
    failures: list[dict[str, str]] = []
    mapped_rows_by_branch: Counter[str] = Counter()
    for index, item in enumerate(candidate_nodes, start=1):
        node_id = item["node_id"]
        node_label = item["node_label"]
        branch_name = item["branch_name"]
        if verbose:
            print(f"[info] fetching {index}/{len(candidate_nodes)} [{branch_name}] {node_label}", file=sys.stderr)
        cache_path = cache_dir / item["cache_table_dir"] / f"{node_id}.json" if cache_dir is not None else None
        try:
            payload = _load_or_fetch_json(
                cache_path=cache_path,
                refresh_cache=refresh_cache,
                base_url=base_url,
                endpoint="browse/table-data.dhtml",
                params={
                    "rootNode": item["root_node"],
                    "species": "Human",
                    "group": "Phenotype",
                    "value": node_label,
                    "index": node_id,
                },
                timeout_seconds=timeout_seconds,
                retries=retries,
                sleep_seconds=sleep_seconds,
                verbose=verbose,
            )
        except Exception as exc:  # noqa: BLE001
            failures.append({"node_id": node_id, "node_label": node_label, "branch_name": branch_name, "error": str(exc)})
            continue

        if not isinstance(payload, list):
            failures.append({"node_id": node_id, "node_label": node_label, "branch_name": branch_name, "error": "non-list payload"})
            continue
        for row in payload:
            if not isinstance(row, dict):
                continue
            if item["mode"] == "raw":
                mapped = _map_raw_row(row, node_label=node_label, node_id=node_id)
            else:
                mapped = _map_literature_row(row, node_label=node_label, node_id=node_id)
            if mapped is not None:
                mapped_rows.append(mapped)
                mapped_rows_by_branch[branch_name] += 1

    output = _aggregate_records(mapped_rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, index=False)

    summary = {
        "output_path": str(output_path),
        "n_candidate_nodes": int(len(candidate_nodes)),
        "candidate_node_counts_by_branch": {str(key): int(value) for key, value in candidate_node_counts_by_branch.items()},
        "n_mapped_rows_raw": int(len(mapped_rows)),
        "mapped_rows_by_branch": {str(key): int(value) for key, value in mapped_rows_by_branch.items()},
        "n_rows": int(len(output)),
        "n_diseases": int(output["disease_name"].nunique()) if not output.empty else 0,
        "disease_counts": (
            {str(key): int(value) for key, value in output["disease_name"].value_counts().to_dict().items()}
            if not output.empty
            else {}
        ),
        "confidence_counts": (
            {str(key): int(value) for key, value in output["relation_confidence"].value_counts().to_dict().items()}
            if not output.empty
            else {}
        ),
        "failed_nodes": failures,
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a project-aligned disease-microbe supplement from gutMDisorder v2.0 human phenotype records."
    )
    parser.add_argument(
        "--output-path",
        default=ROOT / "data/reference/disease_microbe_gutm_disorder_supplement.csv",
        type=Path,
    )
    parser.add_argument(
        "--summary-path",
        default=ROOT / "data/reference/disease_microbe_gutm_disorder_supplement.summary.json",
        type=Path,
    )
    parser.add_argument("--base-url", default=BASE_URL)
    parser.add_argument("--timeout-seconds", default=60.0, type=float)
    parser.add_argument("--retries", default=3, type=int)
    parser.add_argument("--sleep-seconds", default=0.05, type=float)
    parser.add_argument(
        "--cache-dir",
        default=ROOT / "data/cache/gutm_disorder_human_literature_phenotype",
        type=Path,
    )
    parser.add_argument("--refresh-cache", action="store_true")
    parser.add_argument("--max-nodes", default=None, type=int)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    summary = build_gutm_disorder_supplement(
        output_path=args.output_path,
        summary_path=args.summary_path,
        base_url=args.base_url,
        timeout_seconds=args.timeout_seconds,
        retries=args.retries,
        sleep_seconds=args.sleep_seconds,
        cache_dir=args.cache_dir,
        refresh_cache=args.refresh_cache,
        max_nodes=args.max_nodes,
        verbose=args.verbose,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
