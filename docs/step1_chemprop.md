# Step 1 Chemprop Scaffold

## 目标

在现有 `RDKit + ExtraTrees` baseline 之外，新增一条以 `SMILES` 图表示为核心的 `Chemprop` 训练路线。

这条路线当前遵循两个原则：

- 金标准测试集仍然来自 `Maier 2018`
- `MDIPID` 银标只并入分类训练，不进入回归训练和主测试评估

## 当前产物

通过 [prepare_step1_chemprop.py](../scripts/prepare_step1_chemprop.py) 会生成：

- `classification/dataset.csv`
- `classification/descriptors.npz`
- `classification/descriptor_feature_names.json`
- `classification/descriptor_preprocessor.joblib`
- `classification/descriptor_schema.json`
- `regression/dataset.csv`
- `regression/descriptors.npz`
- `regression/descriptor_feature_names.json`
- `regression/descriptor_preprocessor.joblib`
- `regression/descriptor_schema.json`
- `chemprop_prepare_summary.json`

通过 [train_step1_chemprop.py](../scripts/train_step1_chemprop.py) 和
[summarize_step1_chemprop.py](../scripts/summarize_step1_chemprop.py) 会额外生成：

- `chemprop_train_request.json`
- `config.toml`
- `model_0/best.pt`
- `model_0/test_predictions.csv`
- `metrics_summary.json`

## 任务定义

### 分类

- 目标列：`target`
- 标签定义：`effect_label == "inhibit"` 记为 `1`，其余记为 `0`
- 默认切分：`scaffold`

### 回归

- 目标列：`target`
- 目标值：`effect_score`
- 仅使用 `gold` 数据

## Descriptor 设计

Chemprop 学习药物分子图本身，额外的 `descriptor` 用于补充非图信息：

- 药物暴露和理化字段：`dose_umol`, `estimated_colon_concentration_um`, `molecular_weight`, `xlogp`, `tpsa` 等
- 微生物 taxonomy 和培养信息：`species_label`, `phylum`, `genus`, `gram_stain`, `medium_preference` 等
- 粗粒度药理类别：`therapeutic_class`, `atc_primary_l1/l3/l4`

## 准备命令

```bash
python scripts/prepare_step1_chemprop.py \
  --split-mode scaffold \
  --output-dir data/processed/step1/chemprop_scaffold
```

## 训练命令

推荐先在单独环境中安装：

```bash
python3 -m venv /tmp/microbe_env
/tmp/microbe_env/bin/python -m pip install -r requirements-step1-chemprop.txt
```

然后运行分类：

```bash
/tmp/microbe_env/bin/python scripts/train_step1_chemprop.py \
  --dataset-csv data/processed/step1/chemprop_scaffold/classification/dataset.csv \
  --descriptors-path data/processed/step1/chemprop_scaffold/classification/descriptors.npz \
  --output-dir models/step1/chemprop_scaffold_classification_v1 \
  --task-type classification \
  --epochs 10 \
  --extra-arg=--metrics \
  --extra-arg=roc \
  --extra-arg=prc \
  --extra-arg=accuracy \
  --extra-arg=f1 \
  --extra-arg=--class-balance \
  --extra-arg=--accelerator \
  --extra-arg=cpu \
  --extra-arg=--devices \
  --extra-arg=1 \
  --extra-arg=--save-data-splits
```

回归：

```bash
/tmp/microbe_env/bin/python scripts/train_step1_chemprop.py \
  --dataset-csv data/processed/step1/chemprop_scaffold/regression/dataset.csv \
  --descriptors-path data/processed/step1/chemprop_scaffold/regression/descriptors.npz \
  --output-dir models/step1/chemprop_scaffold_regression_v1 \
  --task-type regression \
  --epochs 10 \
  --extra-arg=--metrics \
  --extra-arg=rmse \
  --extra-arg=mae \
  --extra-arg=r2 \
  --extra-arg=--accelerator \
  --extra-arg=cpu \
  --extra-arg=--devices \
  --extra-arg=1 \
  --extra-arg=--save-data-splits
```

训练结束后可统一汇总：

- `descriptor_preprocessor.joblib` 和 `descriptor_schema.json` 会被后续 hybrid 推理接口直接复用

```bash
/tmp/microbe_env/bin/python scripts/summarize_step1_chemprop.py \
  --output-dir models/step1/chemprop_scaffold_classification_v1 \
  --task-type classification
```

```bash
/tmp/microbe_env/bin/python scripts/summarize_step1_chemprop.py \
  --output-dir models/step1/chemprop_scaffold_regression_v1 \
  --task-type regression
```

## 当前结果

- 分类结果在 [metrics_summary.json](../models/step1/chemprop_scaffold_classification_v1/metrics_summary.json)
  - `ROC-AUC`: `0.9007`
  - `PR-AUC`: `0.7545`
  - `accuracy`: `0.8880`
  - `balanced accuracy`: `0.8513`
  - `F1`: `0.6591`
- 回归结果在 [metrics_summary.json](../models/step1/chemprop_scaffold_regression_v1/metrics_summary.json)
  - `RMSE`: `0.2040`
  - `MAE`: `0.1275`
  - `R²`: `0.4893`
  - `Spearman`: `0.3416`

## 结论

- `Chemprop + descriptor` 的二分类 `scaffold` 外推已经跑通，且当前 `balanced accuracy` 高于 `RDKit + ExtraTrees + MDIPID` 的 `0.7996`。
- 回归端目前仍落后于 `RDKit` scaffold baseline，后者 `R² = 0.6169`、`Spearman = 0.4261`，说明 `effect_score` 强度建模暂时还是树模型更稳。
- 因此当前最合理的用法是：把 `Chemprop` 作为 Step 1 分类增强路线保留，把连续效应回归继续交给 `RDKit` baseline。
