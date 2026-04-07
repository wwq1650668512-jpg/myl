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

## 说明边界

- Step 1 与 Step 2 的模型仍然是离线训练好的；网页做的是实时推理，不是网页端再训练。
- Step 3 是实时调用现有模拟器，因此是动态的。
- 库内药物查询与新药 SMILES 预测现在都统一到现成 83 菌扩展面板。
- 新药输入支持任意合法 SMILES，但预测可信度仍然受当前 applicability 范围限制。
