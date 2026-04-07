from __future__ import annotations

import math
from typing import Iterable

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold


MORGAN_BITS = 256
MORGAN_RADIUS = 2
MORGAN_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=MORGAN_RADIUS, fpSize=MORGAN_BITS)


def _pick_smiles_series(frame: pd.DataFrame, smiles_columns: Iterable[str]) -> pd.Series:
    """Choose the first non-empty SMILES value across the preferred input columns."""
    series = pd.Series([math.nan] * len(frame), index=frame.index, dtype=object)
    for column in smiles_columns:
        if column not in frame.columns:
            continue
        values = frame[column]
        mask = series.isna() & values.notna()
        series.loc[mask] = values.loc[mask]
    return series


def _safe_float(value: float | int) -> float:
    """Convert a numeric RDKit result to float, falling back to NaN for None."""
    return float(value) if value is not None else math.nan


def _rdkit_row_from_smiles(smiles: str | float) -> dict[str, str | float]:
    """Compute RDKit descriptors and Morgan bits for one SMILES string."""
    result: dict[str, str | float] = {
        "canonical_smiles_rdkit": math.nan,
        "inchikey": math.nan,
        "murcko_scaffold": math.nan,
        "rdkit_formula": math.nan,
        "rdkit_valid_smiles": 0.0,
        "rdkit_exact_mol_wt": math.nan,
        "rdkit_logp": math.nan,
        "rdkit_tpsa": math.nan,
        "rdkit_molar_refractivity": math.nan,
        "rdkit_formal_charge": math.nan,
        "rdkit_heavy_atom_count": math.nan,
        "rdkit_hbond_donor_count": math.nan,
        "rdkit_hbond_acceptor_count": math.nan,
        "rdkit_rotatable_bond_count": math.nan,
        "rdkit_ring_count": math.nan,
        "rdkit_aromatic_ring_count": math.nan,
        "rdkit_aliphatic_ring_count": math.nan,
        "rdkit_hetero_atom_count": math.nan,
        "rdkit_fraction_csp3": math.nan,
    }
    for bit_index in range(MORGAN_BITS):
        result[f"morgan_fp_{bit_index:03d}"] = 0.0

    if not isinstance(smiles, str) or not smiles.strip():
        return result

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return result

    canonical_smiles = Chem.MolToSmiles(molecule, canonical=True)
    result.update(
        {
            "canonical_smiles_rdkit": canonical_smiles,
            "inchikey": Chem.MolToInchiKey(molecule),
            "murcko_scaffold": MurckoScaffold.MurckoScaffoldSmiles(mol=molecule) or math.nan,
            "rdkit_formula": rdMolDescriptors.CalcMolFormula(molecule),
            "rdkit_valid_smiles": 1.0,
            "rdkit_exact_mol_wt": _safe_float(Descriptors.ExactMolWt(molecule)),
            "rdkit_logp": _safe_float(Crippen.MolLogP(molecule)),
            "rdkit_tpsa": _safe_float(rdMolDescriptors.CalcTPSA(molecule)),
            "rdkit_molar_refractivity": _safe_float(Crippen.MolMR(molecule)),
            "rdkit_formal_charge": _safe_float(sum(atom.GetFormalCharge() for atom in molecule.GetAtoms())),
            "rdkit_heavy_atom_count": _safe_float(molecule.GetNumHeavyAtoms()),
            "rdkit_hbond_donor_count": _safe_float(Lipinski.NumHDonors(molecule)),
            "rdkit_hbond_acceptor_count": _safe_float(Lipinski.NumHAcceptors(molecule)),
            "rdkit_rotatable_bond_count": _safe_float(Lipinski.NumRotatableBonds(molecule)),
            "rdkit_ring_count": _safe_float(rdMolDescriptors.CalcNumRings(molecule)),
            "rdkit_aromatic_ring_count": _safe_float(rdMolDescriptors.CalcNumAromaticRings(molecule)),
            "rdkit_aliphatic_ring_count": _safe_float(rdMolDescriptors.CalcNumAliphaticRings(molecule)),
            "rdkit_hetero_atom_count": _safe_float(rdMolDescriptors.CalcNumHeteroatoms(molecule)),
            "rdkit_fraction_csp3": _safe_float(rdMolDescriptors.CalcFractionCSP3(molecule)),
        }
    )

    fingerprint = MORGAN_GENERATOR.GetFingerprint(molecule)
    bits = np.zeros((MORGAN_BITS,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fingerprint, bits)
    for bit_index, bit_value in enumerate(bits):
        result[f"morgan_fp_{bit_index:03d}"] = float(bit_value)
    return result


def enrich_drug_table_with_rdkit(frame: pd.DataFrame, smiles_columns: Iterable[str]) -> pd.DataFrame:
    """Append RDKit descriptors, scaffold fields, and fingerprints to a drug table.

    Args:
        frame: Input table containing one or more SMILES columns.
        smiles_columns: Ordered SMILES column names to try for each row.

    Returns:
        A copy of the input table with RDKit-derived feature columns added.
    """
    smiles_series = _pick_smiles_series(frame, smiles_columns)
    rdkit_rows = [_rdkit_row_from_smiles(smiles) for smiles in smiles_series]
    rdkit_frame = pd.DataFrame(rdkit_rows, index=frame.index)
    result = pd.concat([frame.reset_index(drop=True), rdkit_frame.reset_index(drop=True)], axis=1)

    if "molecular_formula" in result.columns:
        result["molecular_formula"] = result["molecular_formula"].fillna(result["rdkit_formula"])
    return result
