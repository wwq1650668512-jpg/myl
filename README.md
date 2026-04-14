# Gut Drug-Microbiome Interaction Project

这个仓库用于搭建一个基于真实开源数据库和论文数据集的肠道药物-微生物相互作用建模框架，目标分为三步：

1. 预测药物对微生物的影响程度，输出抑制、促进、无影响，并尽量保留连续效应值。
2. 预测微生物是否影响药物代谢，以及可能生成什么代谢产物。
3. 在给定肠道菌群环境下模拟加药后的群落变化，并结合肠道健康指数评估药物开发价值。

## 当前设计原则

- 只优先采用真实开源数据库、公开论文补充数据、公开补充材料和可复用的标准化资源。
- 先做可复现 baseline，再做更复杂的图模型或机制模型。
- 先以细菌为主，后续再扩展到真菌、古菌和噬菌体。
- 先以 species/strain level 为核心粒度，必要时向 genus level 回退。

## 首选数据来源

- Step 1 主金标准：Maier et al., 2018，1197 drugs x 40 gut strains 的体外生长影响实验。
- Step 1 公开 benchmark：McCoubrey et al., 2021 和 Wang et al., 2023，均基于 Maier 数据做机器学习建模。
- Step 2 主金标准：Zimmermann et al., 2019，76 gut bacteria x 271 oral drugs 的药物代谢测定。
- Step 2 社区层补充：Javdan et al., 2020，Microbiome-Derived Metabolism screen。
- Step 2 机制补充：AGORA2、Rhea、gutMGene v2.0、ChEBI、PubChem、NCBI Taxonomy。
- Step 3 人群与环境：curatedMetagenomicData、GMrepo v3、MGnify、iHMP、paired microbiome-metabolome collection。
- Step 3 健康评分：GMWI2。

详细清单见 [configs/data_sources.yaml](/mnt/e/毕业/configs/data_sources.yaml)。

## 建模总路线

### Step 1: drug -> microbe effect

- 样本单元：`drug x microbe x assay`
- 目标：
  - 三分类：`inhibit / promote / no_effect`
  - 连续值：标准化 growth effect score
- 推荐起步模型：
  - baseline：XGBoost / ExtraTrees / CatBoost
  - 进阶：bipartite graph model 或 multi-task deep tabular model

### Step 2: microbe -> drug metabolism

- 样本单元：`drug x microbe` 或 `drug x community`
- 目标：
  - 是否代谢
  - 代谢反应类型
  - 可能产物
- 推荐起步模型：
  - baseline：二分类 + 反应类型多分类
  - 进阶：候选反应模板检索 + 产物排序模型

### Step 3: community simulation

- 状态变量：
  - 菌群丰度
  - 药物母体浓度
  - 药物代谢产物浓度
  - 关键宿主相关代谢物
  - 健康指数
- 推荐起步模型：
  - 离散时间更新模型
  - 后续接 AGORA2 / MICOM / mgPipe 做机制增强

## 仓库内容

- [docs/technical_route_detailed.md](/mnt/e/毕业/docs/technical_route_detailed.md)：项目技术路线与实现过程总说明。
- [docs/project_blueprint.md](/mnt/e/毕业/docs/project_blueprint.md)：三步模型的详细设计蓝图。
- [docs/step1_pipeline.md](/mnt/e/毕业/docs/step1_pipeline.md)：Step 1 的下载、标准化和 baseline 训练说明。
- [docs/step2_pipeline.md](/mnt/e/毕业/docs/step2_pipeline.md)：Step 2 的标签标准化和建模表组装说明。
- [docs/web_app.md](/mnt/e/毕业/docs/web_app.md)：三步模型可视化网页说明。
- [configs/data_sources.yaml](/mnt/e/毕业/configs/data_sources.yaml)：真实数据源注册表。
- [configs/labeling_rules.yaml](/mnt/e/毕业/configs/labeling_rules.yaml)：标签定义和建模假设。
- [src/gut_drug_microbiome/schemas.py](/mnt/e/毕业/src/gut_drug_microbiome/schemas.py)：核心数据结构。

## Step 1 快速开始

## 环境与命令

仓库现在以 `pyproject.toml` 作为统一的工程入口，建议在虚拟环境中安装：

```bash
python -m pip install -e .
```

如果需要测试工具：

```bash
python -m pip install -e ".[dev]"
```

如果需要 Chemprop 相关能力：

```bash
python -m pip install -e ".[chemprop]"
```

统一 CLI 入口：

```bash
gut-drug-microbiome step1 normalize
gut-drug-microbiome step1 train-baseline --split-mode drug
gut-drug-microbiome step2 assemble
gut-drug-microbiome step3 simulate --drug-query metformin
gut-drug-microbiome web serve --host 127.0.0.1 --port 8080
```

自动化测试：

```bash
pytest
```

安装 Step 1 最小依赖后，依次运行：

```bash
python scripts/download_step1_data.py
python scripts/normalize_step1_data.py
python scripts/train_step1_baseline.py --split-mode drug
```

如果想先看一个最基础的 sanity check，可以再跑：

```bash
python scripts/train_step1_baseline.py \
  --split-mode random \
  --output-dir models/step1/baseline_random_split
```

## 建议下一步

1. 先把 Step 1 的原始数据落盘并统一 drug 与 taxon 标识。
2. 用 Maier 数据做首个三分类 + 回归 baseline。
3. 接 Zimmermann + Javdan 建 Step 2 的代谢标签表。
4. 选一组真实 cohort 做 Step 3 的初始群落环境。
5. 最后把 GMWI2 接成终端评估模块。
