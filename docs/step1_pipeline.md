# Step 1 Pipeline

## 目标

完成 Step 1 的三件核心工作：

1. 下载 Maier et al. 2018 的公开原始数据与补充表。
2. 将原始文件标准化为药物表、微生物表和药物-微生物交互表。
3. 训练一个可复用的 baseline，用于预测药物对肠道微生物的影响。
4. 接入弱监督银标数据，为后续扩大训练信号做准备。

## 当前使用的数据

- Figshare dataset: `10.6084/m9.figshare.4813882.v1`
- Springer supplementary tables:
  - Supplementary Table 1
  - Supplementary Table 2
  - Supplementary Table 3
  - Source data to Fig. 1
- MDIPID 2025:
  - drug / substance metadata
  - microbiota metadata
  - DEIM: drug or other exogenous substances impact on microbiota
- MASI v2.0:
  - microbe / substance interaction records
  - microbe metadata
  - substance metadata

说明：

- `Maier 2018` 是当前 Step 1 的 `gold` 数据。
- `MDIPID DEIM` 当前作为 `silver` 数据，仅并入分类训练，不参与回归训练和主测试集评估。
- `MASI` 当前已接入下载与标准化骨架，但受外部站点稳定性影响，可能需要手工放置原始文件后再做标准化。

## 生成文件

- [step1_drug_table.csv](../data/processed/step1/step1_drug_table.csv)
- [step1_microbe_table.csv](../data/processed/step1/step1_microbe_table.csv)
- [step1_interactions.csv](../data/processed/step1/step1_interactions.csv)
- [step1_modeling_table.csv](../data/processed/step1/step1_modeling_table.csv)
- [step1_summary.json](../data/processed/step1/step1_summary.json)
- [step1_silver_mdipid.csv](../data/processed/step1/step1_silver_mdipid.csv)
- [step1_silver_mdipid_summary.json](../data/processed/step1/step1_silver_mdipid_summary.json)
- [step1_silver_masi.csv](../data/processed/step1/step1_silver_masi.csv)
- [step1_silver_masi_summary.json](../data/processed/step1/step1_silver_masi_summary.json)
- [step1_silver_masi_curated.csv](../data/processed/step1/step1_silver_masi_curated.csv)
- [step1_silver_masi_curated_summary.json](../data/processed/step1/step1_silver_masi_curated_summary.json)
- [step1_silver_mdipid_masi.csv](../data/processed/step1/step1_silver_mdipid_masi.csv)
- [step1_silver_mdipid_masi_summary.json](../data/processed/step1/step1_silver_mdipid_masi_summary.json)
- [step1_silver_mdipid_masi_curated.csv](../data/processed/step1/step1_silver_mdipid_masi_curated.csv)
- [step1_silver_mdipid_masi_curated_summary.json](../data/processed/step1/step1_silver_mdipid_masi_curated_summary.json)

## 标签逻辑

- 连续效应值：`effect_score = reported_norm_auc - 1.0`
- 分类标签：
  - `inhibit`: `effect_score <= -0.20` 且 `q < 0.05`
  - `promote`: `effect_score >= 0.20` 且 `q < 0.05`
  - `no_effect`: 其他情况

说明：

- 当前这份公开数据中几乎没有显著 `promote` 样本，因此默认 baseline 会退化为以 `inhibit / no_effect` 为主的分类任务。
- 三分类接口保留不变，后续接入新数据时无需重写训练代码。
- `MDIPID` 银标通过文献描述词规则推断 `promote / inhibit`，属于弱监督标签，适合做增量训练和外部验证，不适合替代统一体外筛选金标准。

## Baseline 特征

### 药物

- 理化性质：`molecular_weight`, `xlogp`, `tpsa`, `complexity`, `volume3d`
- 暴露相关字段：`dose_umol`, `estimated_intestine_concentration_um`, `estimated_colon_concentration_um`
- 简单双字符结构代理：`smiles_length`, `smiles_ring_index_count` 等
- 粗粒度药理分类：`therapeutic_class`, `therapeutic_effect`, `atc_primary_l1/l3/l4`
- RDKit 描述符：`exact_mol_wt`, `MolLogP`, `TPSA`, `MR`, `formal_charge`, `ring_count`, `fraction_csp3` 等
- `Morgan fingerprint`：`radius=2`, `nBits=256`
- `Murcko scaffold` 与 `rdkit_formula`

### 微生物

- `species_label`
- taxonomy: `phylum`, `class`, `order`, `family`, `genus`
- `gram_stain`, `medium_preference`, `biosafety`

## 评估设计

- `random split`: sanity check
- `drug split`: 模拟见过菌、没见过药物的泛化
- `scaffold split`: 模拟没见过化学骨架的新药泛化

说明：

- 当前 `silver` 只加入分类训练集，测试集始终保持 `gold`。
- 当前回归任务只在 `gold` 上训练和评估，用于保留连续效应强度建模。

## 当前结果摘要

- `gold drug split + RDKit (40 trees)`:
  - balanced accuracy: `0.7886`
  - macro-F1: `0.8229`
  - regression R²: `0.5250`
  - Spearman: `0.3453`
- `gold + MDIPID silver, drug split + RDKit (40 trees)`:
  - balanced accuracy: `0.7938`
  - macro-F1: `0.8265`
  - 对 `gold-only` 有小幅分类增益
- `gold scaffold split + RDKit (40 trees)`:
  - balanced accuracy: `0.8063`
  - macro-F1: `0.8385`
  - regression R²: `0.6169`
  - Spearman: `0.4261`
- `gold + MDIPID silver, scaffold split + RDKit (40 trees)`:
  - balanced accuracy: `0.7996`
  - macro-F1: `0.8299`
  - 银标在更严格的新骨架外推上暂未带来稳定增益
- `MASI silver` 实测标准化结果：
  - 全量 `drug/non-drug mixed`: `8251` 条
  - `drug-like only`: `7165` 条
  - `MASI curated` (`drug-like + Gut + (Human or In vitro)`): `504` 条
- `gold + MASI drug-like, drug split + RDKit (40 trees)`:
  - balanced accuracy: `0.7654`
  - macro-F1: `0.8090`
  - 相比 `gold-only` 明显下降，说明全量 MASI 银标噪声较大
- `gold + MASI drug-like, scaffold split + RDKit (40 trees)`:
  - balanced accuracy: `0.7616`
  - macro-F1: `0.8056`
- `gold + MDIPID + MASI, drug split + RDKit (40 trees)`:
  - balanced accuracy: `0.7483`
  - macro-F1: `0.7970`
  - 说明简单拼接银标源并不可取
- `gold + MASI curated, drug split + RDKit (40 trees)`:
  - balanced accuracy: `0.7852`
  - macro-F1: `0.8178`
  - 比全量 MASI 明显更稳，但仍不如 `MDIPID only`
- `gold + MDIPID + MASI curated, source-aware weighted, drug split + RDKit (40 trees)`:
  - weights: `mdipid_deim=1.0`, `masi_v2=0.25`
  - balanced accuracy: `0.7803`
  - macro-F1: `0.8217`
  - 加权后比未筛选 MASI 稳定，但仍未超过 `MDIPID only`
- `gold microbe split + RDKit (40 trees)`:
  - balanced accuracy: `0.9321`
  - macro-F1: `0.8426`
  - regression R²: `0.6960`
  - Spearman: `0.6556`
- `gold + MDIPID, microbe split + RDKit (40 trees)`:
  - balanced accuracy: `0.9300`
  - macro-F1: `0.8452`
- `gold + MDIPID + MASI curated, source-aware weighted, microbe split + RDKit (40 trees)`:
  - weights: `mdipid_deim=1.0`, `masi_v2=0.25`
  - balanced accuracy: `0.9297`
  - macro-F1: `0.8495`
- `Chemprop + descriptor, scaffold classification v1`:
  - output: [metrics_summary.json](../models/step1/chemprop_scaffold_classification_v1/metrics_summary.json)
  - ROC-AUC: `0.9007`
  - PR-AUC: `0.7545`
  - accuracy: `0.8880`
  - balanced accuracy: `0.8513`
  - F1: `0.6591`
  - 当前二分类 scaffold 外推优于 `RDKit + MDIPID` 的 balanced accuracy
- `Chemprop + descriptor, scaffold regression v1`:
  - output: [metrics_summary.json](../models/step1/chemprop_scaffold_regression_v1/metrics_summary.json)
  - RMSE: `0.2040`
  - MAE: `0.1275`
  - R²: `0.4893`
  - Spearman: `0.3416`
  - 当前回归端仍弱于 `RDKit` scaffold baseline

## 当前推荐配方

- 主训练增强：`gold + MDIPID`
- `MASI` 的更合适用途：`curated` 子集作为补充分析或低权重辅助，不建议直接全量并入
- 更严格泛化评估：至少同时保留 `drug split`、`scaffold split`、`microbe split`

## 运行命令

在已安装 `requirements-step1.txt` 依赖的环境中：

```bash
python scripts/download_step1_data.py
python scripts/normalize_step1_data.py
python scripts/train_step1_baseline.py --split-mode drug
```

如果需要构建 `MDIPID silver`：

```bash
python scripts/download_step1_weak_supervision.py
python scripts/normalize_step1_mdipid.py
```

如果已经手工放入 `MASI` 原始文件，可以继续：

```bash
python scripts/normalize_step1_masi.py
```

也可以做不同切分策略的评估：

```bash
python scripts/train_step1_baseline.py \
  --split-mode random \
  --output-dir models/step1/baseline_random_split
```

```bash
python scripts/train_step1_baseline.py \
  --split-mode scaffold \
  --n-estimators 40 \
  --output-dir models/step1/gold_scaffold_split_rdkit_40
```

```bash
python scripts/train_step1_baseline.py \
  --split-mode drug \
  --n-estimators 40 \
  --silver-table data/processed/step1/step1_silver_mdipid.csv \
  --output-dir models/step1/gold_plus_mdipid_drug_split_rdkit_40
```

```bash
python scripts/train_step1_baseline.py \
  --split-mode microbe \
  --n-estimators 40 \
  --output-dir models/step1/gold_microbe_split_rdkit_40
```

```bash
python scripts/train_step1_baseline.py \
  --split-mode drug \
  --n-estimators 40 \
  --silver-table data/processed/step1/step1_silver_mdipid_masi_curated.csv \
  --source-weight mdipid_deim=1.0 \
  --source-weight masi_v2=0.25 \
  --output-dir models/step1/gold_plus_mdipid_masi_curated_weighted_drug_split_rdkit_40
```

## Chemprop 路线

- `RDKit + ExtraTrees` 仍然是当前主回归 baseline
- `Chemprop` 的数据准备、训练和结果汇总都已经接入，详见 [step1_chemprop.md](step1_chemprop.md)
- 当前推荐：
  - 分类任务可保留 `Chemprop + descriptor` 作为更强的 scaffold 外推路线
  - 连续效应回归仍优先使用 `RDKit + ExtraTrees`

## Hybrid 接口

- 已新增统一推理接口，详见 [step1_hybrid.md](step1_hybrid.md)
- 默认组合：
  - 分类：`models/step1/chemprop_scaffold_classification_v1/model_0/best.pt`
  - 回归：`models/step1/gold_scaffold_split_rdkit_40/regressor.joblib`
- 默认输出：
  - [predictions.csv](../predictions/step1/hybrid_scaffold_v1/predictions.csv)
  - [predictions_slim.csv](../predictions/step1/hybrid_scaffold_v1/predictions_slim.csv)
  - [summary.json](../predictions/step1/hybrid_scaffold_v1/summary.json)
