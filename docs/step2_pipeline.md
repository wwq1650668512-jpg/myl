# Step 2 Pipeline

## 目标

完成 Step 2 的两层数据底座：

1. 把不同来源的药物代谢标签统一到一个标准 schema。
2. 把 Step 1 hybrid 结果接成 Step 2 的候选 pair / modeling table。

当前仓库先做的是“数据层”和“输入层”，还没有开始正式训练 Step 2 模型。

## 当前代码

- 标准化模块：
  - [normalize.py](../src/gut_drug_microbiome/step2/normalize.py)
- `Zimmermann 2019` source-specific 模块：
  - [zimmermann_2019.py](../src/gut_drug_microbiome/step2/zimmermann_2019.py)
- Step 2 baseline 训练模块：
  - [train_baseline.py](../src/gut_drug_microbiome/step2/train_baseline.py)
- Step 2 预测模块：
  - [predict.py](../src/gut_drug_microbiome/step2/predict.py)
- 建模表组装模块：
  - [assemble.py](../src/gut_drug_microbiome/step2/assemble.py)
- 通用标准化脚本：
  - [normalize_step2_generic.py](../scripts/normalize_step2_generic.py)
- `Zimmermann 2019` 标准化脚本：
  - [normalize_step2_zimmermann.py](../scripts/normalize_step2_zimmermann.py)
- Step 2 baseline 训练脚本：
  - [train_step2_baseline.py](../scripts/train_step2_baseline.py)
- Step 2 baseline 预测脚本：
  - [predict_step2_baseline.py](../scripts/predict_step2_baseline.py)
- 候选/建模表组装脚本：
  - [assemble_step2_inputs.py](../scripts/assemble_step2_inputs.py)

## 原始数据放置目录

- [zimmermann_2019](../data/raw/step2/zimmermann_2019)
- [javdan_2020](../data/raw/step2/javdan_2020)
- [agora2](../data/raw/step2/agora2)
- [gutmgene_v2](../data/raw/step2/gutmgene_v2)
- 说明文件见 [README.md](../data/raw/step2/README.md)

## 标准标签 schema

标准化后的 Step 2 标签表至少包含这些字段：

- `pair_id`
- `prestwick_id`
- `nt_code`
- `metabolism_label`
- `reaction_class`
- `parent_depletion_fraction`
- `product_ids`
- `evidence_gene_ids`
- `source_dataset`
- `label_tier`
- `source_scope`

标签约定与 [labeling_rules.yaml](../configs/labeling_rules.yaml) 保持一致：

- `metabolism_label`
  - `metabolized`
  - `not_metabolized`
  - `uncertain`
- `reaction_class`
  - `reduction`
  - `hydrolysis`
  - `deacetylation`
  - `dehydroxylation`
  - `demethylation`
  - `deconjugation`
  - `ring_cleavage`
  - `bioaccumulation_or_unresolved_depletion`
  - `other`

## 当前产物

已经基于 Step 1 hybrid 输出生成：

- [step2_candidate_pairs_full.csv](../data/processed/step2/step2_candidate_pairs_full.csv)
- [step2_candidate_pairs_slim.csv](../data/processed/step2/step2_candidate_pairs_slim.csv)
- [step2_modeling_table.csv](../data/processed/step2/step2_modeling_table.csv)
- [step2_summary.json](../data/processed/step2/step2_summary.json)

当前统计：

- 候选 pair：`43,109`
- 其中带有效 `SMILES`：`42,949`
- 当前还没有并入真实 Step 2 标签，所以 `n_labeled_modeling_rows = 0`
- 候选的 Step 1 hybrid 标签分布：
  - `35,591 no_effect`
  - `7,361 inhibit`
  - `157 promote`

`Zimmermann 2019` 真实 isolate-level `gold` 标签已经标准化到：

- [zimmermann_2019_label_table.csv](../data/processed/step2/zimmermann_2019/zimmermann_2019_label_table.csv)
- [zimmermann_2019_drug_table.csv](../data/processed/step2/zimmermann_2019/zimmermann_2019_drug_table.csv)
- [zimmermann_2019_microbe_table.csv](../data/processed/step2/zimmermann_2019/zimmermann_2019_microbe_table.csv)
- [zimmermann_2019_metabolite_candidates.csv](../data/processed/step2/zimmermann_2019/zimmermann_2019_metabolite_candidates.csv)
- [zimmermann_2019_metabolite_long.csv](../data/processed/step2/zimmermann_2019/zimmermann_2019_metabolite_long.csv)
- [zimmermann_2019_gene_links.csv](../data/processed/step2/zimmermann_2019/zimmermann_2019_gene_links.csv)
- [zimmermann_2019_modeling_table.csv](../data/processed/step2/zimmermann_2019/zimmermann_2019_modeling_table.csv)
- [zimmermann_2019_summary.json](../data/processed/step2/zimmermann_2019/zimmermann_2019_summary.json)

当前 `Zimmermann 2019` 统计：

- `20,596` 个 isolate-level drug-microbe pairs
- `271` 个药物
- `76` 个菌株
- `2,627 metabolized`
- `17,969 not_metabolized`
- `2,911` 个 pair 带候选代谢产物
- `6,183` 个唯一候选代谢产物
- `89` 条 drug-gene 证据链接
- `1,520` 个 pair 带 gene evidence

## 当前 Step 2 Baseline

已经在 `Zimmermann 2019` 上完成第一版真实 baseline 训练：

- `scaffold split`：
  - [metrics.json](../models/step2/zimmermann_scaffold_split/metrics.json)
- `drug split`：
  - [metrics.json](../models/step2/zimmermann_drug_split/metrics.json)
- `microbe split`：
  - [metrics.json](../models/step2/zimmermann_microbe_split/metrics.json)

当前默认部署模型是 `scaffold split + full-fit`，因为它更接近“新药骨架外推”场景：

- 分类：
  - `balanced accuracy = 0.6982`
  - `macro-F1 = 0.6367`
  - `ROC-AUC = 0.7740`
  - `PR-AUC = 0.3160`
- 回归：
  - `RMSE = 0.2409`
  - `MAE = 0.1278`
  - `R² = -0.1098`
  - `Spearman = 0.3897`

两个对照切分：

- `drug split`
  - 分类：`balanced accuracy = 0.7123`，`macro-F1 = 0.6522`，`ROC-AUC = 0.8171`
  - 回归：`R² = 0.0842`，`Spearman = 0.4914`
- `microbe split`
  - 分类：`balanced accuracy = 0.8686`，`macro-F1 = 0.7485`，`ROC-AUC = 0.9486`
  - 回归：`R² = 0.7355`，`Spearman = 0.5745`

说明：

- 当前 Step 2 已经能稳定做两件事：
  - 预测 `metabolized / not_metabolized`
  - 预测连续 `parent_depletion_fraction`
- 当前还不能把“会代谢成什么精确产物”泛化成一个稳定的新药预测器。
  - 原因是 `Zimmermann 2019` 的 `reaction_class` 大多缺失，已解析到标准 reaction class 的只有 `bioaccumulation_or_unresolved_depletion`
  - `product_ids` 目前是来源内候选代谢物注释，适合做 source annotation，不适合直接当跨数据集的新药产物标签

## 当前命令

先生成无标签 Step 2 输入表：

```bash
/tmp/microbe_env/bin/python scripts/assemble_step2_inputs.py \
  --output-dir data/processed/step2
```

如果你已经手工整理出一张 `Zimmermann` 或 `Javdan` 的核心表，可以先做通用标准化：

```bash
/tmp/microbe_env/bin/python scripts/normalize_step2_generic.py \
  --input-path path/to/source_table.csv \
  --output-path data/processed/step2/zimmermann_2019_normalized.csv \
  --source-dataset zimmermann_2019 \
  --label-tier gold \
  --source-scope isolate
```

如果你已经放好了 `Zimmermann 2019` 的补充 workbook，可以直接跑 source-specific 标准化：

```bash
/tmp/microbe_env/bin/python scripts/normalize_step2_zimmermann.py \
  --input-path data/raw/step2/zimmermann_2019/NIHMS1530152-supplement-Supplementary_Tables_1-21.xlsx \
  --output-dir data/processed/step2/zimmermann_2019
```

训练 Step 2 baseline：

```bash
/tmp/microbe_env/bin/python scripts/train_step2_baseline.py \
  --modeling-table data/processed/step2/zimmermann_2019/zimmermann_2019_modeling_table.csv \
  --split-mode scaffold \
  --output-dir models/step2/zimmermann_scaffold_split
```

用默认 `scaffold full-fit` 模型给全量候选 pair 打分：

```bash
/tmp/microbe_env/bin/python scripts/predict_step2_baseline.py \
  --input-table data/processed/step2/step2_candidate_pairs_full.csv \
  --output-dir predictions/step2/baseline_scaffold_v1 \
  --classifier-path models/step2/zimmermann_scaffold_split/classifier_full.joblib \
  --regressor-path models/step2/zimmermann_scaffold_split/regressor_full.joblib \
  --metrics-path models/step2/zimmermann_scaffold_split/metrics.json \
  --applicability-reference-path models/step2/zimmermann_scaffold_split/applicability_reference.joblib
```

然后把标签并进 modeling table：

```bash
/tmp/microbe_env/bin/python scripts/assemble_step2_inputs.py \
  --output-dir data/processed/step2 \
  --step2-label-table data/processed/step2/zimmermann_2019_normalized.csv
```

## 当前预测产物

默认 `scaffold full-fit` 模型已经跑在 Step 1 产出的 Step 2 候选对上：

- [predictions.csv](../predictions/step2/baseline_scaffold_v1/predictions.csv)
- [predictions_slim.csv](../predictions/step2/baseline_scaffold_v1/predictions_slim.csv)
- [summary.json](../predictions/step2/baseline_scaffold_v1/summary.json)

当前统计：

- `43,109` 个候选 pair
- `5,189 predicted metabolized`
- `37,920 predicted not_metabolized`
- `39,075` 个 pair 落在当前 applicability 范围内

## Smoke Test

通用标准化脚本已经用一个小样例跑通过，能正确把：

- `Depleted` 映射成 `metabolized`
- `not metabolized` 映射成 `not_metabolized`
- `demethylation` 保留成标准 reaction class
- `-35` 转成 `-0.35`

## 下一步

- 接入 `Javdan 2020` 的 community-level 标签表，作为外部验证或银标补充。
- 接入 `MagMD / gutMGene / AGORA2` 的酶-反应先验，补强 `reaction_class / product` 预测。
- 把当前 [predictions.csv](../predictions/step2/baseline_scaffold_v1/predictions.csv) 直接作为 Step 3 推演模型输入之一。
