# Step 3 Pipeline

## 目标

在现有 Step 1 和 Step 2 结果之上，完成一个可运行的离散时间肠道群落推演模型：

1. 输入起始肠道群落与给药方案；
2. 用 Step 1 预测药物对菌群的影响；
3. 用 Step 2 预测菌群对药物代谢的影响；
4. 输出时间轨迹、肠道健康指数和启发式开发评分。

当前版本已经是可运行的第一版，但仍然是 `panel-proxy simulation`，还不是最终的真实 cohort 版本。

## 当前代码

- Step 3 模块：
  - [simulation.py](../src/gut_drug_microbiome/step3/simulation.py)
- Step 3 脚本：
  - [run_step3_simulation.py](../scripts/run_step3_simulation.py)
  - [prepare_step3_cohort_community.py](../scripts/prepare_step3_cohort_community.py)
- Step 3 候选筛选脚本：
  - [screen_step3_candidates.py](../scripts/screen_step3_candidates.py)
- 导出入口：
  - [__init__.py](../src/gut_drug_microbiome/step3/__init__.py)
- 真实 cohort 设计：
  - [step3_real_cohort_plan.md](step3_real_cohort_plan.md)

## 当前输入

Step 3 当前直接消费集成后的 Step 2 预测表：

- [predictions.csv](../predictions/step2/baseline_scaffold_v1/predictions.csv)

这张表里已经包含：

- Step 1 预测：
  - `step1_predicted_inhibit_probability`
  - `step1_predicted_effect_score`
  - `step1_predicted_effect_label_hybrid`
- Step 2 预测：
  - `predicted_metabolized_probability`
  - `predicted_parent_depletion_fraction`
  - `applicability_flag`
  - `drug_max_fingerprint_jaccard`

## 当前模拟框架

当前 Step 3 使用离散时间更新：

- 菌群更新：
  - 用 `Step 1 predicted_effect_score` 作为药物压力项
  - 用 `ecology_strength` 把群落缓慢拉回起始背景
- 药物更新：
  - 用 `Step 2 predicted_metabolized_probability` 和 `predicted_parent_depletion_fraction` 近似代谢消耗
  - 同时加一个简单的母药 clearance
- 代谢物更新：
  - 当前是 `aggregate metabolite pool`
  - 还不是精确产物级模拟
- 健康指数：
  - 当前是 `interaction-aware heuristic`
  - 保留 `diversity、beneficial genera、risk genera、stability`
  - 同时新增一个受文献启发的 `ENBI-like interaction balance` 层，用 `positive vs negative interactions` 评估群落是更接近竞争主导还是交叉喂养主导
- 开发评分：
  - 当前会同时保留 `legacy score`
  - 新版综合分会把 `母药保留`、`群落保真度`、`dysbiosis penalty`、`interaction dysbiosis penalty`、`uncertainty penalty`、`metabolite burden` 一起纳入
  - 当前仍用于项目内部排序，不应解释成临床药效

## 当前内置场景

已经实现 4 个内置起始群落场景：

- `healthy_reference`
- `high_fiber`
- `high_fat`
- `antibiotic_perturbed`

这些场景当前已经可以挂到现成的 83 菌扩展面板上运行，但本质上仍然是 `panel proxy`，不是从 `curatedMetagenomicData / GMrepo` 直接回填出来的真实 cohort。

## 当前命令

单场景模拟：

```bash
/tmp/microbe_env/bin/python scripts/run_step3_simulation.py \
  --drug-query "Metformin hydrochloride" \
  --scenario healthy_reference
```

批量跑全部内置场景：

```bash
/tmp/microbe_env/bin/python scripts/run_step3_simulation.py \
  --drug-query "Metformin hydrochloride" \
  --all-scenarios \
  --output-dir predictions/step3/metformin_hydrochloride
```

批量比较多个候选药物：

```bash
/tmp/microbe_env/bin/python scripts/screen_step3_candidates.py \
  --output-dir predictions/step3/candidate_screen_demo \
  --scenario healthy_reference \
  --drug-query "Metformin hydrochloride" \
  --drug-query Digoxin \
  --drug-query Sulfasalazine \
  --drug-query "Diltiazem hydrochloride"
```

## 当前产物

示例药物 `Metformin hydrochloride` 已经跑通，输出在：

- [scenario_grid_summary.csv](../predictions/step3/metformin_hydrochloride/scenario_grid_summary.csv)
- [scenario_grid_summary.json](../predictions/step3/metformin_hydrochloride/scenario_grid_summary.json)

每个场景目录下都有：

- `trajectory_metrics.csv`
- `trajectory_abundances.csv`
- `trajectory_abundances_wide.csv`
- `top_microbe_changes.csv`
- `summary.json`

当前 `Metformin hydrochloride` 的 4 个场景结果概览：

- `healthy_reference`
  - final health index: `49.38`
  - final parent retention ratio: `0.4820`
  - development score: `0.00`
- `high_fiber`
  - final health index: `50.74`
  - final parent retention ratio: `0.4824`
  - development score: `0.00`
- `high_fat`
  - final health index: `46.37`
  - final parent retention ratio: `0.4779`
  - development score: `3.53`
- `antibiotic_perturbed`
  - final health index: `38.21`
  - final parent retention ratio: `0.4814`
  - development score: `12.25`

候选药物比较示例已经输出到：

- [candidate_ranking.csv](../predictions/step3/candidate_screen_demo/candidate_ranking.csv)
- [candidate_ranking.json](../predictions/step3/candidate_screen_demo/candidate_ranking.json)

当前 `healthy_reference` 下的示例排序：

- `Sulfasalazine`
  - development score: `38.87`
  - final health index: `84.33`
- `Diltiazem hydrochloride`
  - development score: `38.60`
  - final health index: `85.53`
- `Digoxin`
  - development score: `36.21`
  - final health index: `82.20`
- `Metformin hydrochloride`
  - development score: `0.00`
  - final health index: `49.38`

## 当前解释边界

当前 Step 3 已经能回答：

- 在不同起始菌群背景里，给药后群落大致朝什么方向变化；
- 哪些菌最可能被压低或放大；
- 母药保留和微生物代谢负担大致如何变化；
- 在当前模型假设下，哪个场景更可能出现较大的 dysbiosis penalty。

当前 Step 3 还不能严格回答：

- 真实人群 cohort 上的定量健康结局；
- 精确产物级、代谢物级药效增减；
- 真实 GMWI2 终值；
- 临床层面的“值不值得开发新药”。

## 最近升级

- 新增 `interaction-aware dynamics`
  - 会基于 `cross_feeding_edges.csv`、`compound_semantic_family`、Step 2 酶先验、代谢概率等信号，构建一个轻量的正负交互网络
  - 在每个时间步输出：
    - `positive_interaction_strength`
    - `negative_interaction_strength`
    - `interaction_balance_rho`
    - `interaction_balance_shift`
- 这个 `interaction_balance_rho` 不是原论文的全文复现，而是按同样思想实现的 `ENBI-like proxy`
  - `rho < 0` 更接近竞争主导
  - `rho > 0` 更接近交叉喂养主导
  - `shift > 0` 表示相对起始状态更往 dysbiosis-like 的合作网络偏移

## 下一步

- 先用 [prepare_step3_cohort_community.py](../scripts/prepare_step3_cohort_community.py) 把真实 abundance 表转换成 `community_table.csv`。
- 再接入 `curatedMetagenomicData / GMrepo / GMWI2`，把当前 `panel proxy` 替换成真实 cohort 初始化。
- 接入 `Javdan 2020` 的 community-level 药物代谢数据，对 Step 3 的药物转化轨迹做外部验证。
- 接入 `AGORA2 / gutMGene / MagMD`，把 `aggregate metabolite pool` 升级成反应类和候选产物层。
- 在 Step 3 上补 `scenario batch runner` 和 `drug portfolio ranking`，做真正的候选药物比较。
