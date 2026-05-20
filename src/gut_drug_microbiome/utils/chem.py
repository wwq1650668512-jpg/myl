from __future__ import annotations

import math

import pandas as pd


def count_smiles_token(smiles: object, token: str) -> float:
    """Count a token inside a SMILES string and return NaN for missing values."""
    if not isinstance(smiles, str) or not smiles:
        return math.nan
    return float(smiles.count(token))


def compute_smiles_descriptors(
    frame: pd.DataFrame,
    primary_smiles_column: str = "main_component_smiles",
    fallback_smiles_column: str = "smiles",
) -> pd.DataFrame:
    """Add lightweight text-derived SMILES descriptors to a table."""
    primary = frame.get(primary_smiles_column)
    fallback = frame.get(fallback_smiles_column)
    if primary is None and fallback is None:
        smiles = pd.Series([math.nan] * len(frame), index=frame.index, dtype=object)
    elif primary is None:
        smiles = fallback
    elif fallback is None:
        smiles = primary
    else:
        smiles = primary.fillna(fallback)

    result = frame.copy()
    result["smiles_length"] = smiles.apply(lambda value: float(len(value)) if isinstance(value, str) else math.nan)
    result["smiles_uppercase_count"] = smiles.apply(
        lambda value: float(sum(1 for char in value if char.isalpha() and char.isupper()))
        if isinstance(value, str)
        else math.nan
    )
    result["smiles_ring_index_count"] = smiles.apply(
        lambda value: float(sum(1 for char in value if char.isdigit()))
        if isinstance(value, str)
        else math.nan
    )
    result["smiles_branch_count"] = smiles.apply(lambda value: count_smiles_token(value, "("))
    result["smiles_double_bond_count"] = smiles.apply(lambda value: count_smiles_token(value, "="))
    result["smiles_halogen_count"] = smiles.apply(
        lambda value: (
            count_smiles_token(value, "Cl")
            + count_smiles_token(value, "Br")
            + count_smiles_token(value, "F")
            + count_smiles_token(value, "I")
        )
        if isinstance(value, str)
        else math.nan
    )
    return result
