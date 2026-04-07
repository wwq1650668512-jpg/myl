# Step 3 Real Cohort Integration Plan

## 目标

把当前 `panel-proxy` 的 Step 3 起始场景，逐步升级成基于真实 cohort abundance 初始化的模拟流程。

当前仓库其实已经具备一半基础：

- Step 3 模拟器已经支持 `--community-table`
- 自定义群落表只需要能映射到现有 Step 1 / Step 3 的微生物面板

这份文档的目标，是把“真实 cohort 怎么接进来”拆成一个可执行的最小方案。

## 当前最小可运行链路

1. 准备 cohort abundance 表
2. 映射到当前 Step 3 面板
3. 导出 `community_table.csv`
4. 用 `scripts/run_step3_simulation.py --community-table ...` 直接模拟

对应新增脚本：

- [prepare_step3_cohort_community.py](../scripts/prepare_step3_cohort_community.py)

## 推荐目录结构

建议后续真实 cohort 数据按下面结构组织：

```text
data/
  raw/
    step3/
      cohorts/
        <dataset_name>/
          abundance_table.csv
          sample_metadata.csv
  processed/
    step3/
      cohorts/
        <dataset_name>/
          <sample_id>_community.csv
          <sample_id>_mapping_report.csv
          <sample_id>_unmapped.csv
          <sample_id>_summary.json
```

## 输入表最低要求

原始 abundance 表至少需要两类信息：

- 一个 taxon 列
- 一个 abundance 列

如果是一张多样本表，还需要：

- 一个 sample id 列

当前预处理脚本会自动尝试识别这些别名：

- taxon 列：
  - `nt_code`
  - `microbe_label`
  - `species_label`
  - `species_name`
  - `species`
  - `taxon`
  - `taxon_name`
  - `genus`
- abundance 列：
  - `abundance`
  - `relative_abundance`
  - `relative_frequency`
  - `weight`
  - `biomass`
- sample 列：
  - `sample_id`
  - `sample`
  - `sample_name`
  - `subject_id`

## 映射规则

当前版本采用一个偏保守但可复现的映射策略：

1. 先尝试 exact mapping
   - `nt_code`
   - `microbe_label`
   - `species_label`
   - `species_name`
2. 如果 exact mapping 失败，再尝试 genus-level mapping
3. genus 级别命中多个面板菌时，按等比例拆分 abundance

因此输出的 mapping report 会标明：

- `mapping_mode = exact`
- `mapping_mode = genus_split`

没有映射到面板的 taxon 会进入 `*_unmapped.csv`，方便后续补 alias 或扩面板。

## 输出格式

预处理脚本的主输出是 `community_table.csv`，当前最关键字段是：

- `nt_code`
- `abundance`

它已经可以直接作为：

```bash
/tmp/microbe_env/bin/python scripts/run_step3_simulation.py \
  --drug-query Digoxin \
  --community-table data/processed/step3/cohorts/example/sample_01_community.csv
```

的输入。

同时会额外导出：

- `*_mapping_report.csv`
- `*_unmapped.csv`
- `*_summary.json`

## 当前解释边界

这条真实 cohort 接入链路现在还只是第一版，边界需要写清楚：

1. 当前只是在“初始化群落”层面接入真实 cohort，还没有把 Step 3 校准到真实纵向结局。
2. genus-level split 是工程上的折中，不是严格的菌株分辨。
3. 如果原始 cohort taxon 命名风格与当前面板差异很大，还需要补 alias 表或更强的 taxonomy resolver。
4. 当前面板本身仍受 Step 1 / Step 2 可覆盖菌范围限制。

## 建议的下一批迭代

### 1. 先接样本级真实群落

优先把这些数据源接成 `community_table.csv`：

- `curatedMetagenomicData`
- `GMrepo`
- `MGnify`

目标是先让 Step 3 从“手工代理场景”变成“真实样本初始化”。

### 2. 再接 cohort metadata

为每个 sample 增加：

- dataset
- disease status
- diet / treatment 标签
- collection timepoint

这样网页和脚本后面就能直接展示：

- `healthy donor`
- `T2D baseline`
- `post-antibiotic`

这类真实场景名。

### 3. 最后接验证

建议把 Step 3 的外部验证拆成三块：

- 群落健康方向验证
- 药物母药保留/代谢方向验证
- 已知菌群依赖药物案例验证

## 推荐命令

把真实 cohort abundance 表整理成 Step 3 community table：

```bash
/tmp/microbe_env/bin/python scripts/prepare_step3_cohort_community.py \
  --input-table data/raw/step3/cohorts/example/abundance_table.csv \
  --sample-id sample_01 \
  --output-path data/processed/step3/cohorts/example/sample_01_community.csv
```

在该群落上运行 Step 3：

```bash
/tmp/microbe_env/bin/python scripts/run_step3_simulation.py \
  --drug-query Digoxin \
  --community-table data/processed/step3/cohorts/example/sample_01_community.csv
```
