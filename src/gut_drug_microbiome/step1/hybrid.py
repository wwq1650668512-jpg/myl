# hybrid是一个集成了Chemprop分类器和基于RDKit特征的回归模型的预测模块，用于生成Step 1的混合效果标签和相关预测结果。它首先确保输入数据具有必要的标识符列，然后分别调用Chemprop分类器和回归模型进行预测，最后将这些预测结果融合成一个综合的效果标签，并输出详细的预测文件和摘要信息。
from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import load

from .chemprop_scaffold import _pick_smiles
from .chemprop_scaffold import _prepare_descriptor_frame
from .chemprop_scaffold import _resolve_chemprop_bin
from .train_baseline import _prepare_feature_frame

# _ensure_identifier_columns：确保输入数据框中存在一个稳定的pair_id列，如果没有，则根据prestwick_id和nt_code生成一个组合列，或者使用行索引生成唯一标识符，以便后续的预测结果能够正确地与输入数据对应。
def _ensure_identifier_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Guarantee that inference rows have a stable pair_id column."""
    result = frame.copy()
    if "pair_id" not in result.columns:
        if {"prestwick_id", "nt_code"}.issubset(result.columns):
            result["pair_id"] = (
                result["prestwick_id"].fillna("unknown_drug").astype(str)
                + "::"
                + result["nt_code"].fillna("unknown_microbe").astype(str)
            )
        else:
            result["pair_id"] = [f"row_{index}" for index in range(len(result))]
    return result

# _ensure_columns：检查数据框中是否存在指定的列，如果缺失则添加这些列并填充NaN值，以确保后续的模型预测步骤能够顺利进行而不会因为缺少必要的特征列而出错。
def _ensure_columns(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """Add any missing columns needed by a downstream model as NaN placeholders."""
    result = frame.copy()
    for column in columns:
        if column not in result.columns:
            result[column] = np.nan
    return result

# _hybrid_effect_label：根据Chemprop的抑制概率、辅助促进概率和回归效果分数，结合预设的阈值，生成一个综合的效果标签。优先考虑抑制标签，如果满足抑制概率阈值则返回"inhibit"，否则如果满足促进概率或效果分数的阈值则返回"promote"，否则返回"no_effect"。
def _hybrid_effect_label(
    inhibit_probability: float | None,
    effect_score: float | None,
    inhibit_probability_threshold: float,
    promote_score_threshold: float,
) -> str:
    """Fuse classification probability and regression score into one hybrid label."""
    if inhibit_probability is not None and not np.isnan(inhibit_probability):
        if inhibit_probability >= inhibit_probability_threshold:
            return "inhibit"
    if effect_score is not None and not np.isnan(effect_score):
        if effect_score >= promote_score_threshold:
            return "promote"
    return "no_effect"

# _predict_chemprop_inhibit_probability：使用Chemprop分类器对输入数据进行预测，返回每行的抑制概率。它首先准备输入数据，包括选择合适的SMILES列、处理缺失值、计算化学描述符等，然后调用Chemprop的命令行接口进行预测，并将结果与输入数据合并返回。
def _predict_chemprop_inhibit_probability(
    frame: pd.DataFrame,
    classification_prepare_dir: str | Path,
    chemprop_model_path: str | Path,
    output_dir: str | Path,
    accelerator: str = "cpu",
    devices: str = "1",
) -> pd.DataFrame:
    """Run the Chemprop classifier and return inhibit probabilities for each row."""
    classification_prepare_dir = Path(classification_prepare_dir)
    chemprop_model_path = Path(chemprop_model_path)
    output_dir = Path(output_dir)

    descriptor_preprocessor_path = classification_prepare_dir / "descriptor_preprocessor.joblib"
    descriptor_schema_path = classification_prepare_dir / "descriptor_schema.json"
    if not descriptor_preprocessor_path.exists():
        raise FileNotFoundError(
            f"Expected Chemprop descriptor preprocessor at {descriptor_preprocessor_path}. "
            "Re-run scripts/prepare_step1_chemprop.py after updating the repository."
        )
    if not descriptor_schema_path.exists():
        raise FileNotFoundError(
            f"Expected Chemprop descriptor schema at {descriptor_schema_path}. "
            "Re-run scripts/prepare_step1_chemprop.py after updating the repository."
        )

    descriptor_preprocessor = load(descriptor_preprocessor_path)
    descriptor_schema = json.loads(descriptor_schema_path.read_text(encoding="utf-8"))
    numeric_features = [str(value) for value in descriptor_schema.get("numeric_descriptor_features", [])]
    categorical_features = [str(value) for value in descriptor_schema.get("categorical_descriptor_features", [])]

    inference = frame.copy()
    inference["smiles"] = _pick_smiles(inference)
    valid_mask = inference["smiles"].notna() & inference["smiles"].astype(str).str.strip().ne("")

    result = inference.loc[:, ["hybrid_row_id", "pair_id", "prestwick_id", "nt_code", "smiles"]].copy()
    result["predicted_inhibit_probability"] = np.nan

    if not valid_mask.any():
        return result

    valid = inference.loc[valid_mask].copy()
    valid = _ensure_columns(valid, numeric_features + categorical_features)
    valid = _prepare_descriptor_frame(valid, numeric_features=numeric_features, categorical_features=categorical_features)

    descriptors = descriptor_preprocessor.transform(valid.loc[:, numeric_features + categorical_features])
    predict_input = valid.loc[:, ["hybrid_row_id", "pair_id", "prestwick_id", "nt_code", "smiles"]].copy()

    with tempfile.TemporaryDirectory(prefix="chemprop_predict_", dir=output_dir) as tmp_dir_name:
        tmp_dir = Path(tmp_dir_name)
        predict_csv = tmp_dir / "predict_input.csv"
        descriptors_path = tmp_dir / "predict_descriptors.npz"
        predictions_path = tmp_dir / "predictions.csv"
        predict_input.to_csv(predict_csv, index=False)
        np.savez_compressed(descriptors_path, descriptors)

        command = [
            _resolve_chemprop_bin(require_installed=True),
            "predict",
            "--test-path",
            str(predict_csv),
            "--smiles-columns",
            "smiles",
            "--descriptors-path",
            str(descriptors_path),
            "--model-paths",
            str(chemprop_model_path),
            "--output",
            str(predictions_path),
            "--accelerator",
            accelerator,
            "--devices",
            devices,
            "--num-workers",
            "0",
        ]
        subprocess.run(command, check=True)
        predicted = pd.read_csv(predictions_path, low_memory=False)

    predicted = predicted.rename(columns={"target": "predicted_inhibit_probability"})
    predicted = predicted.loc[:, ["hybrid_row_id", "predicted_inhibit_probability"]].copy()
    result = result.merge(predicted, on="hybrid_row_id", how="left", suffixes=("", "_new"))
    if "predicted_inhibit_probability_new" in result.columns:
        result["predicted_inhibit_probability"] = result["predicted_inhibit_probability_new"].combine_first(
            result["predicted_inhibit_probability"]
        )
        result = result.drop(columns=["predicted_inhibit_probability_new"])
    return result

# _predict_regression_effect_score：使用基于RDKit特征的回归模型对输入数据进行预测，返回每行的效果分数。它首先加载训练好的回归模型和相关的特征列表，然后确保输入数据包含这些特征列，并进行必要的预处理，最后调用模型进行预测并将结果与输入数据合并返回。
def _predict_regression_effect_score(
    frame: pd.DataFrame,
    regressor_path: str | Path,
    regressor_metrics_path: str | Path,
) -> pd.DataFrame:
    """Run the RDKit baseline regressor and return predicted effect scores."""
    regressor_path = Path(regressor_path)
    regressor_metrics_path = Path(regressor_metrics_path)

    regressor = load(regressor_path)
    regressor_metrics = json.loads(regressor_metrics_path.read_text(encoding="utf-8"))
    numeric_features = [str(value) for value in regressor_metrics.get("numeric_features", [])]
    categorical_features = [str(value) for value in regressor_metrics.get("categorical_features", [])]
    feature_columns = numeric_features + categorical_features

    inference = _ensure_columns(frame, feature_columns)
    inference = _prepare_feature_frame(inference, numeric_features=numeric_features, categorical_features=categorical_features)
    predicted_score = regressor.predict(inference.loc[:, feature_columns])

    result = frame.loc[:, ["hybrid_row_id", "pair_id", "prestwick_id", "nt_code"]].copy()
    result["predicted_effect_score"] = predicted_score
    return result

# _predict_promote_probability：运行可选的辅助促进分类器并返回促进概率。它首先加载促进分类器和相关的特征列表，然后确保输入数据包含这些特征列，并进行必要的预处理，最后调用分类器进行预测并将结果与输入数据合并返回。
# predict_step1_hybrid：生成Step 1的混合预测结果，通过结合Chemprop的抑制概率和基于RDKit特征的回归效果分数，来生成一个综合的效果标签。它还会输出详细的预测文件和一个包含统计摘要的JSON文件。
def predict_step1_hybrid(
    input_table_path: str | Path,
    output_dir: str | Path,
    classification_prepare_dir: str | Path,
    chemprop_model_path: str | Path,
    regressor_path: str | Path,
    regressor_metrics_path: str | Path,
    inhibit_probability_threshold: float = 0.5,
    promote_score_threshold: float = 0.2,
    accelerator: str = "cpu",
    devices: str = "1",
) -> dict[str, object]:
    """Generate Step 1 hybrid predictions by combining Chemprop and baseline regression.

    Args:
        input_table_path: Input candidate pair table for inference.
        output_dir: Directory for prediction CSVs and summary JSON.
        classification_prepare_dir: Chemprop preparation directory with descriptor metadata.
        chemprop_model_path: Trained Chemprop classification checkpoint.
        regressor_path: Trained baseline regression model.
        regressor_metrics_path: Metrics JSON containing regressor feature lists.
        inhibit_probability_threshold: Probability cutoff for an inhibit call.
        promote_score_threshold: Regression cutoff for a promote call.
        accelerator: Chemprop inference accelerator name.
        devices: Device specification passed to Chemprop.

    Returns:
        A summary dictionary describing the generated prediction files.
    """
    input_table_path = Path(input_table_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(input_table_path, low_memory=False)
    raw = _ensure_identifier_columns(raw)
    raw = raw.reset_index(drop=True).copy()
    raw["hybrid_row_id"] = np.arange(len(raw))

    chemprop_predictions = _predict_chemprop_inhibit_probability(
        frame=raw,
        classification_prepare_dir=classification_prepare_dir,
        chemprop_model_path=chemprop_model_path,
        output_dir=output_dir,
        accelerator=accelerator,
        devices=devices,
    )
    regression_predictions = _predict_regression_effect_score(
        frame=raw,
        regressor_path=regressor_path,
        regressor_metrics_path=regressor_metrics_path,
    )

    predictions = raw.merge(
        chemprop_predictions.loc[:, ["hybrid_row_id", "predicted_inhibit_probability"]],
        on="hybrid_row_id",
        how="left",
    )
    predictions = predictions.merge(
        regression_predictions.loc[:, ["hybrid_row_id", "predicted_effect_score"]],
        on="hybrid_row_id",
        how="left",
    )
    predictions["predicted_binary_effect_label"] = np.where(
        predictions["predicted_inhibit_probability"] >= inhibit_probability_threshold,
        "inhibit",
        "no_effect",
    )
    predictions["predicted_effect_label_hybrid"] = [
        _hybrid_effect_label(
            inhibit_probability=inhibit_probability,
            effect_score=effect_score,
            inhibit_probability_threshold=inhibit_probability_threshold,
            promote_score_threshold=promote_score_threshold,
        )
        for inhibit_probability, effect_score in zip(
            predictions["predicted_inhibit_probability"],
            predictions["predicted_effect_score"],
        )
    ]
    predictions["predicted_effect_magnitude"] = predictions["predicted_effect_score"].abs()

    predictions_output = predictions.drop(columns=["hybrid_row_id"])
    predictions_output_path = output_dir / "predictions.csv"
    predictions_output.to_csv(predictions_output_path, index=False)

    slim_columns = [
        "pair_id",
        "prestwick_id",
        "nt_code",
        "smiles",
        "effect_label",
        "binary_effect_label",
        "effect_score",
        "predicted_inhibit_probability",
        "predicted_binary_effect_label",
        "predicted_effect_score",
        "predicted_effect_label_hybrid",
        "predicted_effect_magnitude",
    ]
    existing_slim_columns = [column for column in slim_columns if column in predictions_output.columns]
    predictions_slim_path = output_dir / "predictions_slim.csv"
    predictions_output.loc[:, existing_slim_columns].to_csv(predictions_slim_path, index=False)

    summary = {
        "input_table_path": str(input_table_path),
        "output_dir": str(output_dir),
        "classification_prepare_dir": str(classification_prepare_dir),
        "chemprop_model_path": str(chemprop_model_path),
        "regressor_path": str(regressor_path),
        "regressor_metrics_path": str(regressor_metrics_path),
        "inhibit_probability_threshold": float(inhibit_probability_threshold),
        "promote_score_threshold": float(promote_score_threshold),
        "n_rows": int(len(predictions_output)),
        "n_rows_with_smiles": int(predictions_output["predicted_inhibit_probability"].notna().sum()),
        "n_rows_with_regression_score": int(predictions_output["predicted_effect_score"].notna().sum()),
        "predicted_binary_effect_label_counts": {
            str(key): int(value)
            for key, value in predictions_output["predicted_binary_effect_label"].value_counts().to_dict().items()
        },
        "predicted_effect_label_hybrid_counts": {
            str(key): int(value)
            for key, value in predictions_output["predicted_effect_label_hybrid"].value_counts().to_dict().items()
        },
        "predictions_path": str(predictions_output_path),
        "predictions_slim_path": str(predictions_slim_path),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    return summary
