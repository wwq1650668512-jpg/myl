from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from gut_drug_microbiome.utils.text import canonicalize_key as _canonicalize_key_base
from gut_drug_microbiome.utils.text import normalize_whitespace as _clean_text


def _canonicalize_key(value: object) -> str:
    """Build a lowercase alphanumeric key for loose dictionary matching."""
    return _canonicalize_key_base(value, keep_cjk=True)


def _parse_shift(value: object) -> str:
    """Map raw abundance/effect wording to increase, decrease, or unknown."""
    text = _clean_text(value).lower()
    if not text:
        return "unknown"
    if any(token in text for token in ["↑", "增加", "升高", "上升", "促进", "增殖", "bloom", "increase"]):
        return "increase"
    if any(token in text for token in ["↓", "减少", "降低", "下降", "抑制", "decrease", "drop"]):
        return "decrease"
    return "unknown"


def _parse_microbe_role(value: object) -> str:
    """Map disease-facing wording to protective, risk, or unknown."""
    text = _clean_text(value).lower()
    if not text:
        return "unknown"
    if any(token in text for token in ["保护", "beneficial", "protect"]):
        return "protective"
    if any(token in text for token in ["风险", "危险", "致病", "促病", "有害", "risk", "harm"]):
        return "risk"
    return "unknown"


def _infer_taxon_level(value: object) -> str:
    """Infer whether a disease reference points to species, genus, family, or phylum scope."""
    text = _clean_text(value)
    lower = text.lower()
    if not text:
        return "unknown"
    if "门" in text or "phylum" in lower:
        return "phylum"
    if "纲" in text or "class" in lower:
        return "class"
    if "目" in text or "order" in lower:
        return "order"
    if "科" in text or "aceae" in lower or "family" in lower:
        return "family"
    if "属" in text or " spp" in lower or " sp." in lower:
        return "genus"
    if " " not in text:
        return "genus"
    return "species"


def _extract_genus_hint(value: object) -> str:
    """Extract a genus-like token for broad matching when exact species names are unavailable."""
    text = _clean_text(value)
    if not text:
        return ""
    tokens = re.split(r"[\s/]+", text)
    return tokens[0] if tokens else ""


def _desired_effect_from_microbe_sheet(abundance_shift: str, microbe_role: str) -> str:
    """Infer the desirable drug->microbe effect from disease-associated abundance and role."""
    if abundance_shift == "decrease" and microbe_role == "protective":
        return "promote"
    if abundance_shift == "increase" and microbe_role == "risk":
        return "inhibit"
    return "unknown"


def _desired_effect_from_disease_sheet(disease_effect_on_microbe: str) -> str:
    """Invert disease->microbe direction into a therapeutic correction target."""
    if disease_effect_on_microbe == "increase":
        return "inhibit"
    if disease_effect_on_microbe == "decrease":
        return "promote"
    return "unknown"


def _normalize_disease_microbe_sheet(path: Path) -> pd.DataFrame:
    """Convert the two-sheet disease-microbe workbook into one normalized reference table."""
    workbook = pd.ExcelFile(path)
    frames: list[pd.DataFrame] = []

    if "微生物对疾病的影响" in workbook.sheet_names:
        raw = pd.read_excel(path, sheet_name="微生物对疾病的影响")
        raw.columns = [_clean_text(column) for column in raw.columns]
        microbe_column = "菌株/菌属（学名）"
        disease_column = "疾病"
        chinese_column = "中文名称"
        abundance_column = "丰度变化"
        role_column = "作用方向"
        mechanism_column = "核心机制/临床意义"

        work = raw.copy()
        work[disease_column] = work[disease_column].map(_clean_text).replace("", np.nan).ffill()
        work[microbe_column] = work[microbe_column].map(_clean_text)
        work = work[
            work[microbe_column].ne("")
            & ~work[microbe_column].str.contains("微生物对", na=False)
            & work[disease_column].fillna("").ne("")
        ].copy()
        work["source_sheet"] = "microbe_to_disease"
        work["disease_name"] = work[disease_column].map(_clean_text)
        work["microbe_name_raw"] = work[microbe_column]
        work["microbe_name_cn"] = work.get(chinese_column, pd.Series("", index=work.index)).map(_clean_text)
        work["microbe_key"] = work["microbe_name_raw"].map(_canonicalize_key)
        work["genus_hint"] = work["microbe_name_raw"].map(_extract_genus_hint)
        work["taxon_level"] = work["microbe_name_raw"].map(_infer_taxon_level)
        work["abundance_change_raw"] = work.get(abundance_column, pd.Series("", index=work.index)).map(_clean_text)
        work["abundance_shift"] = work["abundance_change_raw"].map(_parse_shift)
        work["microbe_role_raw"] = work.get(role_column, pd.Series("", index=work.index)).map(_clean_text)
        work["microbe_role_in_disease"] = work["microbe_role_raw"].map(_parse_microbe_role)
        work["disease_effect_on_microbe"] = work["abundance_shift"]
        work["desired_step1_effect"] = [
            _desired_effect_from_microbe_sheet(shift, role)
            for shift, role in zip(work["abundance_shift"], work["microbe_role_in_disease"])
        ]
        work["relation_confidence"] = np.where(work["desired_step1_effect"].eq("unknown"), "low", "high")
        work["mechanism_note"] = work.get(mechanism_column, pd.Series("", index=work.index)).map(_clean_text)
        frames.append(
            work[
                [
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
                ]
            ].copy()
        )

    if "疾病对微生物的影响" in workbook.sheet_names:
        raw = pd.read_excel(path, sheet_name="疾病对微生物的影响")
        raw.columns = [_clean_text(column) for column in raw.columns]
        disease_column = "疾病"
        microbe_column = "菌株/菌属（学名）"
        chinese_column = "中文名称"
        impact_column = "具体影响"
        direction_column = "影响方向"
        mechanism_column = "机制/临床意义"
        note_column = "备注"

        work = raw.copy()
        work[disease_column] = work[disease_column].map(_clean_text).replace("", np.nan).ffill()
        work[microbe_column] = work[microbe_column].map(_clean_text)
        work = work[
            work[microbe_column].ne("")
            & ~work[disease_column].fillna("").str.contains("疾病对微生物", na=False)
            & work[disease_column].fillna("").ne("")
        ].copy()
        work["source_sheet"] = "disease_to_microbe"
        work["disease_name"] = work[disease_column].map(_clean_text)
        work["microbe_name_raw"] = work[microbe_column]
        work["microbe_name_cn"] = work.get(chinese_column, pd.Series("", index=work.index)).map(_clean_text)
        work["microbe_key"] = work["microbe_name_raw"].map(_canonicalize_key)
        work["genus_hint"] = work["microbe_name_raw"].map(_extract_genus_hint)
        work["taxon_level"] = work["microbe_name_raw"].map(_infer_taxon_level)
        work["abundance_change_raw"] = work.get(impact_column, pd.Series("", index=work.index)).map(_clean_text)
        work["abundance_shift"] = work["abundance_change_raw"].map(_parse_shift)
        direction_raw = work.get(direction_column, pd.Series("", index=work.index)).map(_clean_text)
        work["microbe_role_raw"] = direction_raw
        work["microbe_role_in_disease"] = "unknown"
        work["disease_effect_on_microbe"] = np.where(
            work["abundance_shift"].eq("unknown"),
            direction_raw.map(_parse_shift),
            work["abundance_shift"],
        )
        work["desired_step1_effect"] = work["disease_effect_on_microbe"].map(_desired_effect_from_disease_sheet)
        work["relation_confidence"] = np.where(work["desired_step1_effect"].eq("unknown"), "low", "medium")
        mechanism = work.get(mechanism_column, pd.Series("", index=work.index)).map(_clean_text)
        notes = work.get(note_column, pd.Series("", index=work.index)).map(_clean_text)
        work["mechanism_note"] = (mechanism + " " + notes).str.strip()
        frames.append(
            work[
                [
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
                ]
            ].copy()
        )

    if not frames:
        return pd.DataFrame(
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
            ]
        )

    output = pd.concat(frames, ignore_index=True)
    output.insert(0, "reference_id", [f"DMR_{index+1:04d}" for index in range(len(output))])
    return output.drop_duplicates(subset=["source_sheet", "disease_name", "microbe_name_raw", "mechanism_note"]).reset_index(drop=True)


def _normalize_disease_drug_workbook(path: Path) -> pd.DataFrame:
    """Convert the disease->marketed-drug workbook into a tidy disease catalog."""
    workbook = pd.ExcelFile(path)
    rows: list[dict[str, object]] = []
    for sheet_name in workbook.sheet_names:
        frame = pd.read_excel(path, sheet_name=sheet_name)
        if frame.empty:
            continue
        first_column = frame.columns[0]
        for value in frame[first_column].tolist():
            drug_name = _clean_text(value)
            if not drug_name:
                continue
            rows.append(
                {
                    "disease_name": _clean_text(sheet_name),
                    "marketed_drug_name_raw": drug_name,
                    "marketed_drug_key": _canonicalize_key(drug_name),
                }
            )
    output = pd.DataFrame(rows)
    if output.empty:
        return pd.DataFrame(columns=["reference_id", "disease_name", "marketed_drug_name_raw", "marketed_drug_key"])
    output.insert(0, "reference_id", [f"DDR_{index+1:04d}" for index in range(len(output))])
    return output.drop_duplicates(subset=["disease_name", "marketed_drug_name_raw"]).reset_index(drop=True)


def build_disease_reference_tables(
    disease_microbe_xlsx_path: str | Path,
    disease_drug_xlsx_path: str | Path,
    disease_microbe_output_path: str | Path,
    disease_drug_output_path: str | Path,
    summary_path: str | Path | None = None,
) -> dict[str, object]:
    """Normalize manually curated disease-microbe and disease-drug workbooks into CSV references."""
    disease_microbe_xlsx_path = Path(disease_microbe_xlsx_path)
    disease_drug_xlsx_path = Path(disease_drug_xlsx_path)
    disease_microbe_output_path = Path(disease_microbe_output_path)
    disease_drug_output_path = Path(disease_drug_output_path)
    summary_path = (
        disease_microbe_output_path.with_name("disease_reference_summary.json")
        if summary_path is None
        else Path(summary_path)
    )

    disease_microbe = _normalize_disease_microbe_sheet(disease_microbe_xlsx_path)
    disease_drug = _normalize_disease_drug_workbook(disease_drug_xlsx_path)

    disease_microbe_output_path.parent.mkdir(parents=True, exist_ok=True)
    disease_drug_output_path.parent.mkdir(parents=True, exist_ok=True)
    disease_microbe.to_csv(disease_microbe_output_path, index=False)
    disease_drug.to_csv(disease_drug_output_path, index=False)

    summary = {
        "disease_microbe_input": str(disease_microbe_xlsx_path),
        "disease_drug_input": str(disease_drug_xlsx_path),
        "disease_microbe_output": str(disease_microbe_output_path),
        "disease_drug_output": str(disease_drug_output_path),
        "n_disease_microbe_rows": int(len(disease_microbe)),
        "n_diseases_with_microbe_relations": int(disease_microbe["disease_name"].nunique()) if not disease_microbe.empty else 0,
        "n_disease_drug_rows": int(len(disease_drug)),
        "n_diseases_with_marketed_drugs": int(disease_drug["disease_name"].nunique()) if not disease_drug.empty else 0,
        "desired_step1_effect_counts": {
            str(key): int(value)
            for key, value in disease_microbe["desired_step1_effect"].fillna("missing").value_counts().to_dict().items()
        }
        if not disease_microbe.empty
        else {},
        "taxon_level_counts": {
            str(key): int(value)
            for key, value in disease_microbe["taxon_level"].fillna("missing").value_counts().to_dict().items()
        }
        if not disease_microbe.empty
        else {},
    }
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    return summary


def build_disease_adjusted_community(
    microbe_metadata: pd.DataFrame,
    disease_name: str,
    disease_microbe_reference: pd.DataFrame,
) -> pd.DataFrame:
    """Create a simple disease-adjusted community table from curated disease-microbe effects."""
    disease_key = _canonicalize_key(disease_name)
    reference = disease_microbe_reference.copy()
    if reference.empty:
        raise ValueError("Disease microbe reference table is empty.")
    reference = reference[reference["disease_name"].map(_canonicalize_key).eq(disease_key)].copy()
    if reference.empty:
        raise ValueError(f"未找到疾病参考: {disease_name}")

    work = microbe_metadata.copy()
    work["weight"] = 1.0
    work["match_count"] = 0
    work["species_key"] = work.get("species_label", pd.Series("", index=work.index)).map(_canonicalize_key)
    work["microbe_label_key"] = work.get("microbe_label", pd.Series("", index=work.index)).map(_canonicalize_key)
    work["genus_key"] = work.get("genus", pd.Series("", index=work.index)).map(_canonicalize_key)
    work["family_key"] = work.get("family", pd.Series("", index=work.index)).map(_canonicalize_key)
    work["phylum_key"] = work.get("phylum", pd.Series("", index=work.index)).map(_canonicalize_key)

    level_multipliers = {
        "species": {"increase": 1.7, "decrease": 0.55},
        "genus": {"increase": 1.45, "decrease": 0.70},
        "family": {"increase": 1.25, "decrease": 0.82},
        "phylum": {"increase": 1.15, "decrease": 0.90},
        "class": {"increase": 1.10, "decrease": 0.92},
        "order": {"increase": 1.10, "decrease": 0.92},
        "unknown": {"increase": 1.0, "decrease": 1.0},
    }

    for _, row in reference.iterrows():
        shift = str(row.get("disease_effect_on_microbe", "unknown"))
        if shift not in {"increase", "decrease"}:
            continue
        taxon_level = str(row.get("taxon_level", "unknown"))
        key = _canonicalize_key(row.get("microbe_name_raw"))
        if not key:
            continue
        if taxon_level == "species":
            mask = work["species_key"].eq(key) | work["microbe_label_key"].eq(key)
        elif taxon_level == "genus":
            genus_hint = _canonicalize_key(row.get("genus_hint"))
            mask = work["genus_key"].eq(genus_hint or key)
        elif taxon_level == "family":
            mask = work["family_key"].eq(key)
        elif taxon_level == "phylum":
            mask = work["phylum_key"].eq(key)
        else:
            genus_hint = _canonicalize_key(row.get("genus_hint"))
            mask = work["genus_key"].eq(genus_hint) if genus_hint else pd.Series(False, index=work.index)
        if not mask.any():
            continue
        multiplier = level_multipliers.get(taxon_level, level_multipliers["unknown"])[shift]
        work.loc[mask, "weight"] = work.loc[mask, "weight"] * multiplier
        work.loc[mask, "match_count"] = work.loc[mask, "match_count"] + 1

    work["weight"] = work["weight"].clip(lower=1e-6)
    total = float(work["weight"].sum())
    work["abundance"] = work["weight"] / total if total > 0 else 1.0 / max(len(work), 1)
    return work.loc[:, ["nt_code", "abundance", "match_count"]].copy()
