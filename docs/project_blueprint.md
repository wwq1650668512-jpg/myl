# 三步模型实施蓝图

## 1. 任务定义

### Step 1

输入一个药物和一个微生物，预测药物是否抑制、促进或不影响该微生物，并尽可能输出连续效应强度。

### Step 2

输入一个药物和一个微生物，预测该微生物是否会代谢该药物；若会，则进一步预测代谢反应类型和候选产物。

### Step 3

输入一个起始肠道群落、给药方案和可选饮食/宿主背景，模拟菌群与药物在时间维度上的共同变化，并输出药效相关指标与肠道健康指数变化。

## 2. 数据源分工

| 数据源 | 作用 | 对应步骤 |
| --- | --- | --- |
| Maier et al., 2018 | 药物对 40 个肠道菌株的生长影响主金标准 | Step 1 |
| McCoubrey et al., 2021 | Step 1 公开 baseline 对照 | Step 1 |
| Wang et al., 2023 | Step 1 公开 baseline 对照 | Step 1 |
| Zimmermann et al., 2019 | 单菌株药物代谢主金标准 | Step 2 |
| Javdan et al., 2020 | 社区层药物代谢补充数据 | Step 2, Step 3 |
| AGORA2 | 菌株级代谢网络、药物转化能力、基因上下文 | Step 2, Step 3 |
| Rhea + ChEBI + PubChem | 反应模板、产物标准化、结构标准化 | Step 2 |
| gutMGene v2.0 | 微生物-代谢物-宿主基因知识图谱 | Step 2, Step 3 |
| curatedMetagenomicData | 标准化 cohort 丰度和通路数据 | Step 3 |
| GMrepo v3 | 大规模公开人群 gut microbiome 数据库 | Step 3 |
| MGnify / iHMP | 公共多组学和纵向环境数据 | Step 3 |
| microbiome-metabolome collection | 微生物-代谢物联合验证数据 | Step 3 |
| GMWI2 | 健康指数与终端评估 | Step 3 |

## 3. 实体标准化策略

### 药物

- 主 ID：PubChem CID
- 辅助 ID：ChEBI、InChIKey、SMILES
- 统一规则：
  - 优先保留 canonical SMILES
  - 合并盐型和配伍体
  - 对代谢产物保留 parent-child 关系

### 微生物

- 主 ID：NCBI Taxonomy ID
- 粒度：
  - 有菌株数据时保留 strain
  - 跨 cohort 聚合时回退到 species 或 genus
- 与 AGORA2 的 strain/model ID 建映射表

## 4. Step 1 设计

### 4.1 标签

- 主标签：`inhibit / promote / no_effect`
- 辅标签：连续效应值 `effect_score`
- 说明：
  - Maier 原始实验天然更偏向抑制端点。
  - “促进”通常需要从连续生长值中二次构建，而不是直接来自论文原始三分类标签。

### 4.2 特征

#### 药物特征

- 结构指纹：Morgan / MACCS / PubChem fingerprints
- 理化描述符：MW、logP、TPSA、HBA、HBD、rotatable bonds
- 药理上下文：ATC class 或 target family，后续再接入

#### 微生物特征

- taxonomy embedding
- AGORA2 或公开基因组导出的代谢能力特征
- drug-modifying enzyme presence/absence
- SCFA、胆汁酸、糖代谢等 pathway 模块

#### 配对特征

- 药物结构特征与菌属/菌种 embedding 的交互项
- 是否存在已知耐药或转化相关酶
- 已知培养条件和暴露浓度

### 4.3 模型

- baseline：
  - XGBoost
  - ExtraTrees
  - CatBoost
- 进阶：
  - multi-task neural network
  - drug-microbe heterogeneous graph model

### 4.4 评估

- Macro-F1
- Balanced accuracy
- AUROC
- Spearman/Pearson for continuous effect
- Quadratic weighted kappa for ordinal interpretation

### 4.5 数据切分

- random split 只做最初 sanity check
- scaffold split 用于检验新药泛化
- leave-one-strain-out 用于检验跨微生物泛化
- leave-one-class-of-drug-out 用于检验真实开发场景

## 5. Step 2 设计

### 5.1 任务拆分

Step 2 不建议一开始直接端到端做“产物生成”，更稳妥的顺序是：

1. 二分类：会不会被代谢
2. 多分类：属于哪类反应
3. 候选产物排序：在给定反应模板下产物是什么

### 5.2 标签层级

- Level 1：`metabolized / not_metabolized`
- Level 2：reaction class
  - reduction
  - hydrolysis
  - deacetylation
  - dehydroxylation
  - demethylation
  - deconjugation
  - ring cleavage
  - bioaccumulation or depletion without resolved product
- Level 3：具体产物结构

### 5.3 特征

#### 药物侧

- Step 1 中的所有结构与理化特征
- 可反应位点
- 是否属于前药、糖苷化产物、葡萄糖醛酸结合物等

#### 微生物侧

- AGORA2 中与药物转化相关的反应和基因
- Rhea/UniProt/ChEBI 支持下的酶反应注释
- gutMGene 中与宿主代谢和微生物代谢物相关的上下文

#### 配对侧

- 药物可反应基团 x 酶能力的匹配分数
- 已知相似药物的被代谢邻域

### 5.4 模型

- baseline：
  - 二分类模型预测是否代谢
  - 多分类模型预测反应类型
- 进阶：
  - 反应模板检索模型
  - 基于候选 SMARTS/Rhea 反应的产物重排与排序

### 5.5 评估

- 二分类：AUROC、AUPRC、balanced accuracy
- 反应类型：macro-F1
- 产物预测：
  - exact match
  - top-k hit rate
  - Tanimoto similarity

## 6. Step 3 设计

### 6.1 状态变量

- `A_i(t)`：第 i 个菌的相对丰度或 biomass
- `D_j(t)`：第 j 个药物母体浓度
- `P_k(t)`：第 k 个代谢产物浓度
- `M_l(t)`：关键宿主相关代谢物浓度
- `H(t)`：肠道健康指数

### 6.2 离散时间更新框架

建议先做轻量离散时间模型，而不是一上来做全约束代谢网络：

`A_i(t+1) = A_i(t) * exp(base_growth_i + drug_effect_i + ecology_i + interaction_i)`

`D_j(t+1) = D_j(t) - sum_i metabolism_rate_ij * A_i(t)`

`P_k(t+1) = P_k(t) + sum_(i,j) product_flux_ijk - clearance_k`

`H(t+1) = f(taxa, metabolites, diversity, health_signature_proxy, interaction_balance_proxy)`

这里建议明确一层新的 Step 3 设计原则：

- 不只看“哪些菌变多/变少”
- 还要看群落更接近：
  - `competition-dominated`
  - 还是 `cross-feeding-dominated`

这条思路受《Imbalance in gut microbial interactions as a marker of health and disease》启发，当前仓库里适合先落成 `ENBI-like interaction proxy`，再逐步向真实 cohort interaction inference 过渡。

### 6.3 输入环境

- 起始菌群丰度：来自 curatedMetagenomicData / GMrepo / iHMP
- 功能背景：HUMAnN pathway abundance 或 AGORA2 映射
- 饮食或培养基：先做几个固定场景
  - balanced diet
  - high fiber
  - high fat
  - antibiotic perturbed

### 6.4 输出指标

- 母药保留比例
- 关键活性代谢物比例
- 受影响菌的数量与方向
- alpha diversity 变化
- GMWI2 或其简化替代值
- interaction balance 指标
  - `positive_interaction_strength`
  - `negative_interaction_strength`
  - `interaction_balance_rho`
  - `interaction_balance_shift`
- 综合开发评分

## 7. 药物开发价值评分建议

可以先做一个启发式综合分，而不是直接声称临床药效预测：

`development_score = benefit_subscore - risk_subscore`

其中：

- `efficacy_proxy`：母药保留、活性代谢物生成、目标暴露 proxy
- `dysbiosis_penalty`：群落失衡、健康指数下降、关键益生菌受损
- `interaction_dysbiosis_penalty`：交互网络向 dysbiosis-like 正互作放大状态偏移
- `uncertainty_penalty`：超出训练分布、菌群背景缺失、产物不确定

推荐更接近当前实现的写法是：

`benefit_subscore = efficacy_proxy + community_preservation_score`

`risk_subscore = dysbiosis_penalty + interaction_dysbiosis_penalty + uncertainty_penalty + metabolite_burden_penalty`

## 8. 验证路线

### Step 1

- 对照已公开 benchmark
- 做 scaffold split 和 leave-one-strain-out

### Step 2

- 先复现 Zimmermann / Javdan 已知现象
- 再做新药或新菌泛化测试

### Step 3

- 用已知扰动场景做回顾性验证：
  - 抗生素暴露
  - metformin
  - FMT
  - 高脂饮食相关 cohort

## 9. 当前最合理的落地顺序

### Sprint 1

- 下载 Step 1 主数据
- 完成药物与 taxon 标准化
- 训练首个三分类 + 回归 baseline

### Sprint 2

- 整理 Step 2 标签表
- 建立药物反应模板与产物标准化流程
- 完成二分类代谢模型

### Sprint 3

- 接入一个真实 cohort
- 实现离散时间模拟器
- 接入健康指数

### Sprint 4

- 融合 AGORA2 / MICOM / mgPipe
- 做更强机制模拟
- 开始形成可发表结果

## 10. 关键风险

- Step 1 的“促进”标签样本通常远少于抑制和无影响。
- Step 2 的“产物是什么”远难于“会不会代谢”。
- Step 3 容易被过度解释，必须把“机制推演”与“临床结论”严格区分。
- 不同 cohort 的批次效应、测序平台差异和菌种分辨率差异会显著影响模拟稳定性。

## 11. 当前建议

如果我们现在就开始动手，最稳的起点不是 Step 3，而是：

1. 先把 Step 1 做到可复现、可解释、可外推。
2. 再做 Step 2 的二分类和反应类型预测。
3. 最后在真实 cohort 上接入 Step 3 模拟器。

这样成功率最高，也最容易逐步产出论文级结果。
