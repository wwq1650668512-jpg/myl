# Web App

## 目标

给当前项目补一个本地网页入口，把：

- Step 1：药物对微生物影响
- Step 2：微生物对药物代谢
- Step 3：群落时间推演

统一到一个页面里可视化展示，并支持直接交互预测。当前主入口已经切换成“新药 SMILES 预测优先”。

## 当前实现

网页应用采用“轻量本地服务 + 原生前端”的方式实现，不依赖额外的前端框架。

后端文件：

- [service.py](../src/gut_drug_microbiome/web/service.py)
- [server.py](../src/gut_drug_microbiome/web/server.py)
- [run_web_app.py](../scripts/run_web_app.py)

前端文件：

- [index.html](../webapp/static/index.html)
- [app.js](../webapp/static/app.js)
- [styles.css](../webapp/static/styles.css)

## 页面能力

### 新药输入

- 只输入 `SMILES` 就可以触发整条预测链；
- 药物名字是可选项，仅用于页面展示；
- 后端会自动生成内部 `Custom-<session_id>` 药物 ID；
- 默认会在当前 83 个微生物扩展面板上完成 Step 1 和 Step 2 预测；
- 支持在同一页面继续运行 Step 3 场景模拟。

### Step 1

- 选择一个微生物；
- 查看该 pair 的 `predicted_effect_label`、`inhibit probability`、`effect_score`；
- 查看该药物在全部微生物上的影响强度条形图；
- 查看该药物下 Top 影响微生物表。
- Step 1 结果会附带现实性约束元数据（例如 `step1_drug_profile`、`step1_constraint_applied`），用于解释为何某些药物不会被判成“无差别强抑制”。

### Step 2

- 查看该 pair 的 `predicted_metabolism_label`、`predicted_metabolized_probability`、`predicted_parent_depletion_fraction`；
- 查看 applicability 相关信息；
- 查看该药物在全部微生物上的代谢概率条形图；
- 查看 Top 代谢微生物表。

### Step 3

- 选择内置场景：`healthy_reference`, `high_fiber`, `high_fat`, `antibiotic_perturbed`
- 调整时间步数、剂量、清除率、代谢强度等参数；
- 运行单场景模拟；
- 查看健康指数、母药保留、开发评分轨迹；
- 对比全部内置场景；
- 查看 Top 微生物变化表。

### 机制层 + 疾病候选

新药预测会返回 `candidate_diseases`，每个疾病条目同时包含：

- 原始菌层分：`disease_score_raw_only`
- 机制融合分：`disease_score_mechanism`
- 机制细分分量：`mechanism_scores`（如 `butyrate_support_score`、`barrier_protection_score`、`toxin_risk_score` 等）
- 融合模式标记：`fusion_mode`（服务默认 `weighted_0.65_0.35`）
- 机制证据：`mechanism_top_contributors`、`evidence_examples`

说明：当前 Web 默认返回“服务内置融合模式”的结果。  
固定 case 的多融合模式对比（`raw_only / mechanism_only / weighted_0.3_0.7`）由离线脚本产出，详见下方“评估脚本”。

## 数据来源

网页同时支持两种数据来源：

- 库内药物模式：直接读取 [predictions.csv](../predictions/step2/baseline_scaffold_v1/predictions.csv)
- 新药 SMILES 模式：基于输入 SMILES 动态调用 Step 1 hybrid、Step 2 baseline 和 Step 3 simulation

库内预测表已经包含：

- Step 1 hybrid 输出
- Step 2 baseline 输出
- 药物和微生物元数据

因此：

- 库内药物不需要重新训练模型，就能直接查询和模拟；
- 新药 SMILES 也不需要网页端训练，只是调用已经训练好的模型做实时推理。

## 疾病目录与候选空间更新

为避免功能性肠病候选缺失，服务启动时会对疾病参考表执行标准化与扩展：

- 统一疾病命名（含 IBS 别名归一）
- 强制保证以下候选在目录中可见：
  - `肠易激综合征（IBS）`
  - `肠易激综合征-腹泻型（IBS-D）`
  - `肠易激综合征-便秘型（IBS-C）`
- 当参考表中仅有 IBS 基线条目时，自动扩展生成 IBS-D / IBS-C 候选覆盖

实现位置：[service.py](../src/gut_drug_microbiome/web/service.py)

## 运行方式

推荐使用已有的 Python 环境运行，因为默认系统 `python3` 当前缺少 `pandas`。

```bash
/tmp/microbe_env/bin/python scripts/run_web_app.py --host 127.0.0.1 --port 8080
```

然后在浏览器中打开：

```text
http://127.0.0.1:8080
```

如果你已经在自己的环境里装好了 [requirements-step1.txt](../requirements-step1.txt) 的依赖，也可以直接用自己的 Python 解释器启动。

## API 概览

当前内置接口：

- `GET /api/bootstrap`
- `GET /api/drug-profile?drug=<prestwick_id_or_name>`
- `GET /api/pair-prediction?drug=<prestwick_id_or_name>&microbe=<nt_code_or_name>`
- `POST /api/custom-drug/predict`
- `GET /api/custom-drug/pair?session_id=<id>&microbe=<nt_code_or_name>`
- `POST /api/step3/simulate`
- `POST /api/step3/scenario-grid`
- `POST /api/custom-drug/step3/simulate`
- `POST /api/custom-drug/step3/scenario-grid`

`POST /api/custom-drug/predict` 的最小请求体只需要：

```json
{
  "smiles": "CC(=O)Oc1ccccc1C(=O)O"
}
```

返回体中的 `profile.aggregated` 还会包含跨菌群聚合指标，例如：

- `step1_counts`
- `mean_predicted_effect_score`
- `mean_predicted_inhibit_probability`
- `candidate_diseases`

## 说明边界

- Step 1 与 Step 2 的模型仍然是离线训练好的；网页做的是实时推理，不是网页端再训练。
- Step 3 是实时调用现有模拟器，因此是动态的。
- 库内药物查询与新药 SMILES 预测现在都统一到现成 83 菌扩展面板。
- 新药输入支持任意合法 SMILES，但预测可信度仍然受当前 applicability 范围限制。

## 评估脚本（离线）

机制层行为对比与 case-based benchmark 在离线脚本中运行：

- [run_fusion_comparison.py](../scripts/run_fusion_comparison.py)
- [evaluate_case_based_sanity_benchmark.py](../scripts/evaluate_case_based_sanity_benchmark.py)
- [build_revised_case_outputs.py](../scripts/build_revised_case_outputs.py)

主要输出目录：

- `predictions/evaluation/fusion_comparison/`
  - `rifaximin.csv`, `vancomycin.csv`, `lubiprostone.csv`, `metronidazole.csv`
  - `sanity_summary.csv`
  - `ranking_benchmark_summary.csv`
  - `ecology_benchmark_summary.csv`
  - `revised_case_based_results.csv`
  - `revised_case_based_summary.md`
