# Step 1 Hybrid Predictor

## 目标

把 Step 1 的两条最稳路线合成一个统一接口：

- 分类：`Chemprop + descriptor`
- 连续效应回归：`RDKit + ExtraTrees`

这样后续 Step 2/3 可以直接读取一份标准输出，而不需要分别调用两个模型。

## 默认模型

- 分类模型：
  - [best.pt](../models/step1/chemprop_scaffold_classification_v1/model_0/best.pt)
- 回归模型：
  - [regressor.joblib](../models/step1/gold_scaffold_split_rdkit_40/regressor.joblib)
  - [metrics.json](../models/step1/gold_scaffold_split_rdkit_40/metrics.json)

## 推理命令

推荐使用安装了 `Chemprop` 的独立环境运行：

```bash
/tmp/microbe_env/bin/python scripts/predict_step1_hybrid.py \
  --input-table data/processed/step1/step1_modeling_table.csv \
  --output-dir predictions/step1/hybrid_scaffold_v1
```

## 输出文件

- [predictions.csv](../predictions/step1/hybrid_scaffold_v1/predictions.csv)
  - 保留原始标准化特征列，并附加预测结果
- [predictions_slim.csv](../predictions/step1/hybrid_scaffold_v1/predictions_slim.csv)
  - 只保留 `pair_id / drug / microbe / 真实值 / 预测值`
- [summary.json](../predictions/step1/hybrid_scaffold_v1/summary.json)
  - 记录模型路径、阈值和预测标签分布

## 标签合成规则

- `predicted_binary_effect_label`
  - 当 `predicted_inhibit_probability >= 0.5` 时记为 `inhibit`
  - 否则记为 `no_effect`
- `predicted_effect_label_hybrid`
  - 如果 `Chemprop` 判定 `inhibit`，输出 `inhibit`
  - 否则如果 `predicted_effect_score >= 0.2`，输出 `promote`
  - 否则输出 `no_effect`

## 当前全表推理结果

- 输入：`43,109` 个药物-微生物 pair
- 其中带有效 `SMILES` 并完成 `Chemprop` 分类的有 `42,949` 个
- 混合标签分布：
  - `35,591 no_effect`
  - `7,361 inhibit`
  - `157 promote`

## 说明

- 这份接口当前更适合做 Step 1 到 Step 3 的工程串联，不代表新的独立 benchmark。
- 如果后续更换 `Chemprop` 或 `RDKit` 模型，只需要替换脚本参数中的模型路径，不需要改接口格式。
