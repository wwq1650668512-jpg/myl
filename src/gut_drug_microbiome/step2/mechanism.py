from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import dump
from joblib import load


MECHANISM_GROUP_SPECS: tuple[dict[str, object], ...] = (
    {"name": "drug_nt", "drug_column": "prestwick_id", "microbe_column": "nt_code", "weight": 1.00},
    {"name": "drug_species", "drug_column": "prestwick_id", "microbe_column": "species_label", "weight": 0.95},
    {"name": "drug_genus", "drug_column": "prestwick_id", "microbe_column": "genus", "weight": 0.85},
    {"name": "scaffold_species", "drug_column": "murcko_scaffold", "microbe_column": "species_label", "weight": 0.75},
    {"name": "scaffold_genus", "drug_column": "murcko_scaffold", "microbe_column": "genus", "weight": 0.65},
    {"name": "class_genus", "drug_column": "therapeutic_class", "microbe_column": "genus", "weight": 0.55},
    {"name": "atc_l1_genus", "drug_column": "atc_primary_l1", "microbe_column": "genus", "weight": 0.45},
    {"name": "drug_global", "drug_column": "prestwick_id", "microbe_column": None, "weight": 0.50},
    {"name": "scaffold_global", "drug_column": "murcko_scaffold", "microbe_column": None, "weight": 0.40},
    {"name": "genus_global", "drug_column": None, "microbe_column": "genus", "weight": 0.30},
)


def _canonicalize_token(value: object) -> str:
    """Normalize a token used to build mechanism lookup keys."""
    if pd.isna(value):
        return ""
    text = str(value).strip().lower()
    return "".join(character for character in text if character.isalnum())


def _pick_column(frame: pd.DataFrame, candidates: list[str]) -> str | None:
    """Return the first available column name from a candidate list."""
    for candidate in candidates:
        if candidate in frame.columns:
            return candidate
    return None


def _split_multi_value(value: object) -> list[str]:
    """Split a semicolon-delimited field into a list of non-empty tokens."""
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    return [item.strip() for item in text.split(";") if item.strip()]


def _build_group_key(row: pd.Series, drug_column: str | None, microbe_column: str | None) -> str | None:
    """Build a normalized grouping key for one drug/microbe aggregation scope."""
    parts: list[str] = []
    if drug_column is not None:
        drug_key = _canonicalize_token(row.get(drug_column))
        if not drug_key:
            return None
        parts.append(drug_key)
    else:
        parts.append("__drug__")

    if microbe_column is not None:
        microbe_key = _canonicalize_token(row.get(microbe_column))
        if not microbe_key:
            return None
        parts.append(microbe_key)
    else:
        parts.append("__microbe__")
    return "::".join(parts)


def _empty_entry() -> dict[str, object]:
    """Create an empty aggregation entry for mechanism support counts."""
    return {
        "n_pairs": 0,
        "support_mass": 0.0,
        "reaction_counts": {},
        "product_counts": {},
        "gene_counts": {},
    }


def build_step2_mechanism_reference(
    modeling_table_path: str | Path | pd.DataFrame,
    output_path: str | Path | None = None,
) -> dict[str, object]:
    """Aggregate known metabolized pairs into a reusable mechanism lookup reference.

    Args:
        modeling_table_path: Step 2 modeling table path or an in-memory DataFrame.
        output_path: Optional joblib path where the reference is persisted.

    Returns:
        A dictionary containing grouped reaction, product, and gene support.
    """
    if isinstance(modeling_table_path, pd.DataFrame):
        frame = modeling_table_path.copy()
        modeling_table_repr = "<dataframe>"
    else:
        frame = pd.read_csv(modeling_table_path, low_memory=False)
        modeling_table_repr = str(modeling_table_path)

    label_column = _pick_column(frame, ["step2_metabolism_label", "metabolism_label"])
    reaction_column = _pick_column(frame, ["step2_reaction_class", "reaction_class"])
    product_column = _pick_column(frame, ["step2_product_ids", "product_ids"])
    gene_column = _pick_column(frame, ["step2_evidence_gene_ids", "evidence_gene_ids"])
    depletion_column = _pick_column(
        frame,
        ["step2_parent_depletion_fraction", "parent_depletion_fraction", "predicted_parent_depletion_fraction"],
    )
    if label_column is None:
        raise RuntimeError("Could not locate a Step 2 metabolism label column when building mechanism reference.")

    labeled = frame[frame[label_column].fillna("").astype(str).str.strip().eq("metabolized")].copy()
    if labeled.empty:
        raise RuntimeError("No metabolized rows were found in the Step 2 modeling table.")

    group_index: dict[str, dict[str, dict[str, object]]] = {str(spec["name"]): {} for spec in MECHANISM_GROUP_SPECS}

    for _, row in labeled.iterrows():
        depletion_value = pd.to_numeric(pd.Series([row.get(depletion_column)]), errors="coerce").iloc[0]
        support_mass = float(max(abs(float(depletion_value)) if not pd.isna(depletion_value) else 0.0, 0.25))
        reaction_class = str(row.get(reaction_column) or "").strip()
        product_ids = _split_multi_value(row.get(product_column))
        gene_ids = _split_multi_value(row.get(gene_column))

        for spec in MECHANISM_GROUP_SPECS:
            name = str(spec["name"])
            key = _build_group_key(
                row,
                drug_column=spec["drug_column"],  # type: ignore[arg-type]
                microbe_column=spec["microbe_column"],  # type: ignore[arg-type]
            )
            if key is None:
                continue
            entry = group_index[name].setdefault(key, _empty_entry())
            entry["n_pairs"] = int(entry["n_pairs"]) + 1
            entry["support_mass"] = float(entry["support_mass"]) + support_mass
            if reaction_class:
                reaction_counts = dict(entry["reaction_counts"])
                reaction_counts[reaction_class] = float(reaction_counts.get(reaction_class, 0.0)) + support_mass
                entry["reaction_counts"] = reaction_counts
            if product_ids:
                product_counts = dict(entry["product_counts"])
                for product_id in product_ids:
                    product_counts[product_id] = float(product_counts.get(product_id, 0.0)) + support_mass
                entry["product_counts"] = product_counts
            if gene_ids:
                gene_counts = dict(entry["gene_counts"])
                for gene_id in gene_ids:
                    gene_counts[gene_id] = float(gene_counts.get(gene_id, 0.0)) + support_mass
                entry["gene_counts"] = gene_counts

    reference = {
        "version": 1,
        "modeling_table_path": modeling_table_repr,
        "label_column": label_column,
        "reaction_column": reaction_column,
        "product_column": product_column,
        "gene_column": gene_column,
        "group_specs": [dict(spec) for spec in MECHANISM_GROUP_SPECS],
        "group_index": group_index,
        "n_metabolized_rows": int(len(labeled)),
    }
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        dump(reference, output_path)
    return reference


class Step2MechanismProjector:
    """Project reaction/product/gene hints onto new Step 2 predictions."""

    def __init__(self, reference: dict[str, object] | None = None) -> None:
        self.reference = reference or {}
        self.group_specs = [dict(item) for item in self.reference.get("group_specs", [])]
        self.group_index = {
            str(key): dict(value) for key, value in self.reference.get("group_index", {}).items()
        }

    @classmethod
    def from_joblib(cls, path: str | Path | None) -> "Step2MechanismProjector":
        """Load a mechanism projector from a joblib file if one is available."""
        if path is None:
            return cls(reference=None)
        candidate = Path(path)
        if not candidate.exists():
            return cls(reference=None)
        return cls(reference=load(candidate))

    def is_available(self) -> bool:
        """Return whether the projector has usable mechanism reference content."""
        return bool(self.group_specs and self.group_index)

    def annotate_frame(
        self,
        frame: pd.DataFrame,
        predicted_probability_column: str = "predicted_metabolized_probability",
        predicted_label_column: str = "predicted_metabolism_label",
        probability_floor: float = 0.35,
    ) -> pd.DataFrame:
        """Annotate predictions with mechanism projections derived from grouped evidence."""
        result = frame.copy()
        defaults = {
            "predicted_mechanism_projection_flag": False,
            "predicted_reaction_class": "",
            "predicted_reaction_confidence": np.nan,
            "predicted_reaction_support_pairs": 0,
            "predicted_mechanism_support_score": np.nan,
            "predicted_mechanism_support_scopes": "",
            "predicted_candidate_product_ids": "",
            "predicted_candidate_product_count": 0,
            "predicted_evidence_gene_ids": "",
            "predicted_evidence_gene_count": 0,
        }
        for column, default in defaults.items():
            result[column] = default

        if result.empty or not self.is_available():
            return result

        for row_index, row in result.iterrows():
            predicted_label = str(row.get(predicted_label_column) or "").strip()
            predicted_probability = pd.to_numeric(pd.Series([row.get(predicted_probability_column)]), errors="coerce").iloc[0]
            if predicted_label != "metabolized" and (pd.isna(predicted_probability) or float(predicted_probability) < probability_floor):
                continue

            reaction_support: Counter[str] = Counter()
            product_support: Counter[str] = Counter()
            gene_support: Counter[str] = Counter()
            support_pairs = 0
            support_score = 0.0
            matched_scopes: list[str] = []

            for spec in self.group_specs:
                name = str(spec["name"])
                key = _build_group_key(
                    row,
                    drug_column=spec.get("drug_column"),  # type: ignore[arg-type]
                    microbe_column=spec.get("microbe_column"),  # type: ignore[arg-type]
                )
                if key is None:
                    continue
                entry = self.group_index.get(name, {}).get(key)
                if not entry:
                    continue
                scope_weight = float(spec.get("weight", 0.0))
                support_pairs += int(entry.get("n_pairs", 0))
                support_score += scope_weight * float(entry.get("support_mass", 0.0))
                matched_scopes.append(name)

                for reaction_class, count in dict(entry.get("reaction_counts", {})).items():
                    reaction_support[str(reaction_class)] += float(count) * scope_weight
                for product_id, count in dict(entry.get("product_counts", {})).items():
                    product_support[str(product_id)] += float(count) * scope_weight
                for gene_id, count in dict(entry.get("gene_counts", {})).items():
                    gene_support[str(gene_id)] += float(count) * scope_weight

            if not matched_scopes:
                continue

            result.at[row_index, "predicted_mechanism_projection_flag"] = True
            result.at[row_index, "predicted_reaction_support_pairs"] = int(support_pairs)
            result.at[row_index, "predicted_mechanism_support_score"] = float(support_score)
            result.at[row_index, "predicted_mechanism_support_scopes"] = ";".join(matched_scopes)

            if reaction_support:
                top_reaction, top_value = reaction_support.most_common(1)[0]
                total_reaction = float(sum(reaction_support.values()))
                result.at[row_index, "predicted_reaction_class"] = top_reaction
                result.at[row_index, "predicted_reaction_confidence"] = (
                    float(top_value) / total_reaction if total_reaction > 0 else np.nan
                )

            top_products = [product_id for product_id, _ in product_support.most_common(3)]
            top_genes = [gene_id for gene_id, _ in gene_support.most_common(5)]
            result.at[row_index, "predicted_candidate_product_ids"] = ";".join(top_products)
            result.at[row_index, "predicted_candidate_product_count"] = int(len(top_products))
            result.at[row_index, "predicted_evidence_gene_ids"] = ";".join(top_genes)
            result.at[row_index, "predicted_evidence_gene_count"] = int(len(top_genes))
        return result
