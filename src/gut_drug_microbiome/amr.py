from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import pandas as pd

try:
    from rdkit import Chem
except ImportError:  # pragma: no cover - optional at import time
    Chem = None


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_AMR_REFERENCE_PATH = ROOT / "data/processed/amr/microbe_amr_reference.csv"
DEFAULT_AMR_RULES_PATH = ROOT / "data/processed/amr/drug_resistance_rules.csv"
DEFAULT_INHIBIT_THRESHOLD = 0.5
DEFAULT_PROMOTE_THRESHOLD = 0.2

_BETA_LACTAM_RING = Chem.MolFromSmarts("N1C(=O)CC1") if Chem is not None else None
_RULE_COLUMNS = [
    "rule_id",
    "drug_class",
    "drug_name",
    "species_label",
    "genus",
    "expected_phenotype",
    "rule_level",
    "mechanism_hint",
    "rule_strength",
    "source_name",
    "source_url",
    "source_version",
    "notes",
]
_STRENGTH_RANK = {"supporting": 1, "moderate": 2, "strong": 3}
_RESISTANT_FACTORS = {
    "supporting": (0.70, 0.80),
    "moderate": (0.45, 0.55),
    "strong": (0.20, 0.35),
}


def _normalize_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return "".join(character for character in str(value).strip().lower() if character.isalnum())


def _ordered_unique(items: list[str]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for item in items:
        if item and item not in seen:
            result.append(item)
            seen.add(item)
    return result


def _safe_float(value: object) -> float | None:
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(numeric):
        return None
    return float(numeric)


def _hybrid_effect_label(
    inhibit_probability: float | None,
    effect_score: float | None,
    inhibit_probability_threshold: float = DEFAULT_INHIBIT_THRESHOLD,
    promote_score_threshold: float = DEFAULT_PROMOTE_THRESHOLD,
) -> str:
    if inhibit_probability is not None and not math.isnan(inhibit_probability):
        if inhibit_probability >= inhibit_probability_threshold:
            return "inhibit"
    if effect_score is not None and not math.isnan(effect_score):
        if effect_score >= promote_score_threshold:
            return "promote"
    return "no_effect"


class AmrRuleEngine:
    def __init__(
        self,
        rules_path: str | Path = DEFAULT_AMR_RULES_PATH,
        microbe_reference_path: str | Path = DEFAULT_AMR_REFERENCE_PATH,
    ) -> None:
        self.rules_path = Path(rules_path)
        self.microbe_reference_path = Path(microbe_reference_path)
        self.rules = self._load_rules()
        self.microbe_reference = self._load_microbe_reference()
        self.reference_by_species_key = self._build_reference_lookup("species_label")
        self.reference_by_microbe_key = self._build_reference_lookup("microbe_label")

    def _load_rules(self) -> pd.DataFrame:
        if not self.rules_path.exists():
            frame = pd.DataFrame(columns=_RULE_COLUMNS)
        else:
            frame = pd.read_csv(self.rules_path, low_memory=False)
        for column in _RULE_COLUMNS:
            if column not in frame.columns:
                frame[column] = ""
        frame = frame.loc[:, _RULE_COLUMNS].copy()
        frame["species_key"] = frame["species_label"].map(_normalize_text)
        frame["genus_key"] = frame["genus"].map(_normalize_text)
        frame["drug_class_key"] = frame["drug_class"].map(_normalize_text)
        frame["drug_name_key"] = frame["drug_name"].map(_normalize_text)
        frame["expected_phenotype_key"] = frame["expected_phenotype"].map(_normalize_text)
        frame["rule_strength_key"] = frame["rule_strength"].map(_normalize_text)
        frame["strength_rank"] = frame["rule_strength_key"].map(_STRENGTH_RANK).fillna(0).astype(int)
        return frame

    def _load_microbe_reference(self) -> pd.DataFrame:
        if not self.microbe_reference_path.exists():
            return pd.DataFrame()
        frame = pd.read_csv(self.microbe_reference_path, low_memory=False)
        if "nt_code" not in frame.columns:
            return pd.DataFrame()
        return frame.drop_duplicates(subset=["nt_code"]).reset_index(drop=True)

    def _build_reference_lookup(self, column: str) -> dict[str, dict[str, object]]:
        if self.microbe_reference.empty or column not in self.microbe_reference.columns:
            return {}
        lookup: dict[str, dict[str, object]] = {}
        for _, row in self.microbe_reference.iterrows():
            key = _normalize_text(row.get(column))
            if not key:
                continue
            record = {
                "genus": row.get("genus"),
                "phylum": row.get("phylum"),
                "gram_stain": row.get("gram_stain"),
            }
            if key not in lookup or any(pd.isna(lookup[key].get(field)) or not lookup[key].get(field) for field in record):
                lookup[key] = {
                    field: (
                        record[field]
                        if record[field] is not None and not pd.isna(record[field]) and str(record[field]).strip()
                        else lookup.get(key, {}).get(field)
                    )
                    for field in record
                }
        return lookup

    def _fill_taxonomy_fields(self, frame: pd.DataFrame) -> pd.DataFrame:
        work = frame.copy()
        for index, row in work.iterrows():
            species_key = _normalize_text(row.get("species_label"))
            microbe_key = _normalize_text(row.get("microbe_label"))
            reference = self.reference_by_species_key.get(species_key) or self.reference_by_microbe_key.get(microbe_key) or {}

            genus = row.get("genus")
            if pd.isna(genus) or not str(genus).strip():
                inferred_genus = reference.get("genus")
                if not inferred_genus and isinstance(row.get("species_label"), str):
                    parts = str(row.get("species_label")).strip().split()
                    if len(parts) >= 2:
                        inferred_genus = parts[0]
                if inferred_genus:
                    work.at[index, "genus"] = inferred_genus

            phylum = row.get("phylum")
            if (pd.isna(phylum) or not str(phylum).strip()) and reference.get("phylum"):
                work.at[index, "phylum"] = reference.get("phylum")

            gram_stain = row.get("gram_stain")
            if (pd.isna(gram_stain) or not str(gram_stain).strip()) and reference.get("gram_stain"):
                work.at[index, "gram_stain"] = reference.get("gram_stain")
        return work

    def infer_drug_context(self, row: pd.Series) -> dict[str, object]:
        classes: list[str] = []
        text_fields = [
            row.get("chemical_name"),
            row.get("therapeutic_class"),
            row.get("therapeutic_effect"),
            row.get("atc_primary_l1"),
            row.get("atc_primary_l3"),
            row.get("atc_primary_l4"),
        ]
        normalized_texts = [_normalize_text(value) for value in text_fields]
        joined = " ".join(part for part in normalized_texts if part)

        penicillin_tokens = [
            "penicillin",
            "benzylpenicillin",
            "ampicillin",
            "amoxicillin",
            "oxacillin",
            "cloxacillin",
            "dicloxacillin",
            "flucloxacillin",
            "nafcillin",
            "piperacillin",
            "ticarcillin",
        ]
        cephalosporin_tokens = ["cephalospor", "cef", "ceft", "cefa"]
        carbapenem_tokens = ["carbapenem", "imipenem", "meropenem", "ertapenem", "doripenem"]
        monobactam_tokens = ["monobactam", "aztreonam"]
        glycopeptide_tokens = ["glycopeptide", "vancomycin", "teicoplanin", "dalbavancin", "oritavancin", "telavancin"]
        lipopeptide_tokens = ["lipopeptide", "daptomycin"]
        aminoglycoside_tokens = [
            "aminoglycoside",
            "gentamicin",
            "amikacin",
            "tobramycin",
            "streptomycin",
            "kanamycin",
            "neomycin",
            "plazomicin",
            "paromomycin",
        ]
        polymyxin_tokens = ["polymyxin", "polymyxinb", "colistin", "polymyxine"]
        low_anaerobe_fluoroquinolone_tokens = ["ciprofloxacin", "levofloxacin", "ofloxacin", "norfloxacin"]

        if any(token in joined for token in penicillin_tokens):
            classes.extend(["penicillin", "beta_lactam"])
        if any(token in joined for token in cephalosporin_tokens):
            classes.extend(["cephalosporin", "beta_lactam"])
        if any(token in joined for token in carbapenem_tokens):
            classes.extend(["carbapenem", "beta_lactam"])
        if any(token in joined for token in monobactam_tokens):
            classes.extend(["monobactam", "beta_lactam"])
        if any(token in joined for token in glycopeptide_tokens):
            if "vancomycin" in joined:
                classes.append("vancomycin")
            classes.append("glycopeptide")
        if any(token in joined for token in lipopeptide_tokens):
            if "daptomycin" in joined:
                classes.append("daptomycin")
            classes.append("lipopeptide")
        if any(token in joined for token in aminoglycoside_tokens):
            if "gentamicin" in joined:
                classes.append("gentamicin")
            if "amikacin" in joined:
                classes.append("amikacin")
            if "tobramycin" in joined:
                classes.append("tobramycin")
            if "streptomycin" in joined:
                classes.append("streptomycin")
            classes.append("aminoglycoside")
        if any(token in joined for token in polymyxin_tokens):
            if "colistin" in joined:
                classes.append("colistin")
            if "polymyxinb" in joined:
                classes.append("polymyxin_b")
            classes.append("polymyxin")
        if any(token in joined for token in low_anaerobe_fluoroquinolone_tokens):
            if "ciprofloxacin" in joined:
                classes.append("ciprofloxacin")
            if "levofloxacin" in joined:
                classes.append("levofloxacin")
            if "ofloxacin" in joined:
                classes.append("ofloxacin")
            if "norfloxacin" in joined:
                classes.append("norfloxacin")
            classes.append("fluoroquinolone_low_anaerobe")

        smiles = str(row.get("canonical_smiles_rdkit") or row.get("smiles") or row.get("main_component_smiles") or "").strip()
        if smiles and Chem is not None and _BETA_LACTAM_RING is not None:
            molecule = Chem.MolFromSmiles(smiles)
            if molecule is not None and molecule.HasSubstructMatch(_BETA_LACTAM_RING):
                classes.append("beta_lactam")

        ordered_classes = _ordered_unique(classes)
        primary = ordered_classes[0] if ordered_classes else ""
        return {
            "drug_name_key": _normalize_text(row.get("chemical_name")),
            "drug_classes": ordered_classes,
            "primary_drug_class": primary,
        }

    def _match_rule(self, row: pd.Series, drug_context: dict[str, object]) -> pd.Series | None:
        if self.rules.empty:
            return None
        species_key = _normalize_text(row.get("species_label"))
        genus_key = _normalize_text(row.get("genus"))
        drug_name_key = str(drug_context.get("drug_name_key", ""))
        drug_classes = set(str(value) for value in drug_context.get("drug_classes", []))
        primary_drug_class = str(drug_context.get("primary_drug_class", ""))

        best_row: pd.Series | None = None
        best_score = -1
        for _, candidate in self.rules.iterrows():
            species_match = bool(candidate["species_key"]) and candidate["species_key"] == species_key
            genus_match = bool(candidate["genus_key"]) and candidate["genus_key"] == genus_key
            if not species_match and not genus_match:
                continue

            if candidate["drug_name_key"]:
                if candidate["drug_name_key"] != drug_name_key:
                    continue
                drug_specificity = 25
            else:
                if candidate["drug_class_key"] and candidate["drug_class_key"] not in drug_classes:
                    continue
                drug_specificity = 15 if candidate["drug_class_key"] == primary_drug_class else 10

            taxon_specificity = 30 if species_match else 20
            score = taxon_specificity + drug_specificity + int(candidate["strength_rank"]) * 5
            if score > best_score:
                best_row = candidate
                best_score = score
        return best_row

    def annotate_frame(
        self,
        frame: pd.DataFrame,
        *,
        label_column: str | None,
        probability_column: str | None,
        score_column: str | None,
    ) -> pd.DataFrame:
        work = self._fill_taxonomy_fields(frame)

        raw_label_values = work[label_column] if label_column and label_column in work.columns else pd.Series([None] * len(work))
        if len(raw_label_values) != len(work):
            raw_label_values = pd.Series([None] * len(work), index=work.index)
        elif getattr(raw_label_values, "index", None) is None or not raw_label_values.index.equals(work.index):
            raw_label_values = pd.Series(list(raw_label_values), index=work.index)
        raw_probability_values = (
            pd.to_numeric(work[probability_column], errors="coerce")
            if probability_column and probability_column in work.columns
            else pd.Series([math.nan] * len(work), index=work.index)
        )
        raw_score_values = (
            pd.to_numeric(work[score_column], errors="coerce")
            if score_column and score_column in work.columns
            else pd.Series([math.nan] * len(work), index=work.index)
        )

        work["raw_step1_predicted_effect_label"] = raw_label_values
        work["raw_step1_predicted_inhibit_probability"] = raw_probability_values
        work["raw_step1_predicted_effect_score"] = raw_score_values
        work["display_step1_predicted_effect_label"] = raw_label_values
        work["display_step1_predicted_inhibit_probability"] = raw_probability_values
        work["display_step1_predicted_effect_score"] = raw_score_values
        work["amr_expected_phenotype"] = ""
        work["amr_rule_id"] = ""
        work["amr_rule_strength"] = ""
        work["amr_rule_level"] = ""
        work["amr_mechanism_hint"] = ""
        work["amr_source_name"] = ""
        work["amr_source_url"] = ""
        work["amr_conflict_flag"] = False
        work["amr_correction_applied"] = False
        work["amr_correction_reason"] = ""

        if work.empty:
            return work

        drug_context = self.infer_drug_context(work.iloc[0])
        work["amr_primary_drug_class"] = str(drug_context.get("primary_drug_class", ""))
        work["amr_drug_classes"] = "|".join(str(value) for value in drug_context.get("drug_classes", []))

        for index, row in work.iterrows():
            matched_rule = self._match_rule(row, drug_context)
            if matched_rule is None:
                continue

            expected = str(matched_rule.get("expected_phenotype_key", ""))
            strength = str(matched_rule.get("rule_strength_key", ""))
            work.at[index, "amr_expected_phenotype"] = str(matched_rule.get("expected_phenotype", "") or "")
            work.at[index, "amr_rule_id"] = str(matched_rule.get("rule_id", "") or "")
            work.at[index, "amr_rule_strength"] = str(matched_rule.get("rule_strength", "") or "")
            work.at[index, "amr_rule_level"] = str(matched_rule.get("rule_level", "") or "")
            work.at[index, "amr_mechanism_hint"] = str(matched_rule.get("mechanism_hint", "") or "")
            work.at[index, "amr_source_name"] = str(matched_rule.get("source_name", "") or "")
            work.at[index, "amr_source_url"] = str(matched_rule.get("source_url", "") or "")

            if expected != "resistant":
                continue

            raw_label_value = work.at[index, "raw_step1_predicted_effect_label"]
            raw_label = "" if pd.isna(raw_label_value) else str(raw_label_value)
            raw_probability = _safe_float(work.at[index, "raw_step1_predicted_inhibit_probability"])
            raw_score = _safe_float(work.at[index, "raw_step1_predicted_effect_score"])

            conflict_flag = raw_label == "inhibit"
            work.at[index, "amr_conflict_flag"] = bool(conflict_flag)

            probability_factor, score_factor = _RESISTANT_FACTORS.get(strength, _RESISTANT_FACTORS["supporting"])
            corrected_probability = raw_probability
            corrected_score = raw_score
            if raw_probability is not None:
                corrected_probability = float(raw_probability) * float(probability_factor)
            if raw_score is not None and raw_score < 0:
                corrected_score = float(raw_score) * float(score_factor)

            corrected_label = _hybrid_effect_label(corrected_probability, corrected_score)
            correction_applied = (
                corrected_label != raw_label
                or (raw_probability is not None and corrected_probability is not None and abs(corrected_probability - raw_probability) > 1e-12)
                or (raw_score is not None and corrected_score is not None and abs(corrected_score - raw_score) > 1e-12)
            )

            work.at[index, "display_step1_predicted_effect_label"] = corrected_label
            work.at[index, "display_step1_predicted_inhibit_probability"] = corrected_probability
            work.at[index, "display_step1_predicted_effect_score"] = corrected_score
            work.at[index, "amr_correction_applied"] = bool(correction_applied)
            if correction_applied:
                work.at[index, "amr_correction_reason"] = (
                    f"Applied {strength or 'supporting'} AMR resistant prior for "
                    f"{matched_rule.get('genus') or matched_rule.get('species_label') or row.get('microbe_label')}"
                )

        return work
