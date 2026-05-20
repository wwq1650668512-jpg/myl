from __future__ import annotations

import pandas as pd
import pytest

from gut_drug_microbiome.web.service import GutPredictionService


def test_build_custom_drug_table_rejects_invalid_smiles() -> None:
    service = GutPredictionService.__new__(GutPredictionService)

    with pytest.raises(ValueError, match="SMILES 格式无效"):
        service._build_custom_drug_table(
            drug_name="Invalid demo",
            smiles="C1=CC",
            drug_id="Custom-test",
        )


def test_build_custom_drug_table_accepts_valid_smiles() -> None:
    service = GutPredictionService.__new__(GutPredictionService)

    frame = service._build_custom_drug_table(
        drug_name="Valid demo",
        smiles="CCO",
        drug_id="Custom-test",
    )

    assert isinstance(frame, pd.DataFrame)
    assert len(frame) == 1
    assert float(frame.loc[0, "rdkit_valid_smiles"]) == 1.0
