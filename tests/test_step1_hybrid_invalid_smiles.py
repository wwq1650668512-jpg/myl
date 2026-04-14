from __future__ import annotations

import json

import pandas as pd
from joblib import dump

from gut_drug_microbiome.step1.hybrid import _predict_chemprop_inhibit_probability


def test_predict_chemprop_skips_invalid_smiles_without_subprocess(tmp_path, monkeypatch) -> None:
    prepare_dir = tmp_path / "prepare"
    prepare_dir.mkdir(parents=True, exist_ok=True)
    dump({"placeholder": 1}, prepare_dir / "descriptor_preprocessor.joblib")
    (prepare_dir / "descriptor_schema.json").write_text(
        json.dumps(
            {
                "numeric_descriptor_features": [],
                "categorical_descriptor_features": [],
            }
        ),
        encoding="utf-8",
    )

    frame = pd.DataFrame(
        [
            {
                "hybrid_row_id": 0,
                "pair_id": "p0",
                "prestwick_id": "d0",
                "nt_code": "m0",
                "smiles": "C1=CC",  # invalid ring closure
                "rdkit_valid_smiles": 0,
            }
        ]
    )

    monkeypatch.setattr(
        "gut_drug_microbiome.step1.hybrid.subprocess.run",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("subprocess.run should not be called")),
    )

    result = _predict_chemprop_inhibit_probability(
        frame=frame,
        classification_prepare_dir=prepare_dir,
        chemprop_model_path=tmp_path / "dummy.pt",
        output_dir=tmp_path,
    )

    assert len(result) == 1
    assert pd.isna(result.loc[0, "predicted_inhibit_probability"])
