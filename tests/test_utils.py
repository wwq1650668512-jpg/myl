from __future__ import annotations

import math

import pandas as pd

from gut_drug_microbiome.utils.chem import compute_smiles_descriptors
from gut_drug_microbiome.utils.text import canonicalize_key
from gut_drug_microbiome.utils.text import normalize_whitespace


def test_normalize_whitespace_and_canonicalize_key() -> None:
    assert normalize_whitespace("  Foo\tbar  ") == "Foo bar"
    assert canonicalize_key(" Foo-bar 123 ") == "foobar123"
    assert canonicalize_key("肠 道 Microbe", keep_cjk=True) == "肠道microbe"


def test_compute_smiles_descriptors_uses_primary_then_fallback() -> None:
    frame = pd.DataFrame(
        [
            {"main_component_smiles": "CC(=O)Cl", "smiles": "ignored"},
            {"main_component_smiles": math.nan, "smiles": "C1=CC=CC=C1"},
        ]
    )
    result = compute_smiles_descriptors(frame)
    assert result.loc[0, "smiles_length"] == float(len("CC(=O)Cl"))
    assert result.loc[0, "smiles_branch_count"] == 1.0
    assert result.loc[0, "smiles_double_bond_count"] == 1.0
    assert result.loc[0, "smiles_halogen_count"] == 1.0
    assert result.loc[1, "smiles_ring_index_count"] == 2.0
