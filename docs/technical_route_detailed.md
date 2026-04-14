# 技术路线与实现过程详解

更新日期：`2026-04-14`

本文用于系统说明当前仓库的总体技术路线、数据工程流程、建模策略、评估设计、已有结果与后续扩展方向。它面向三类场景：

- 项目开题或中期汇报；
- 论文方法部分撰写；
- 代码仓库交接和后续迭代。

## 0. 阅读指南与常见名词解释

这份文档默认读者已经知道不少机器学习和微生物组领域术语。为了让第一次接触项目的人也能顺着读下来，先把最常见的概念集中解释一下。

建议阅读顺序：

1. 先看第 `1` 节和第 `2` 节，理解为什么一定要拆成三步；
2. 再看第 `3` 节，理解全局数据工程原则；
3. 然后顺着读 `Step 1 -> Step 2 -> Step 3`；
4. 最后看第 `8` 节和第 `11` 节，把系统边界和后续路线对齐。

### 0.1 pair-level 是什么

`pair-level` 指的是：一个训练样本不是单独的“药物”或“微生物”，而是一对实体。

例如：

- Step 1 学的是 `drug x microbe` 的影响关系；
- Step 2 学的也是 `drug x microbe` 的代谢关系。

所以 Step 1 和 Step 2 都是“成对关系预测器”，而 Step 3 才是“群落系统级模拟器”。

### 0.2 schema 是什么

`schema` 可以理解为“统一的数据结构规范”。

不同论文和数据库里，药物 id、菌株名字、效应标签、显著性字段、代谢字段的叫法经常不同。  
如果不先统一 schema，后面的训练脚本就会写出大量特判逻辑，既难维护也容易出错。

### 0.3 gold / silver 是什么

- `gold`：实验体系比较统一、可信度更高、适合作为主训练和主评估依据的数据。
- `silver`：覆盖更广，但异质性更大、噪声更多的数据。

更直白地说：

- `gold` 更像主教材；
- `silver` 更像辅助习题集。

`silver` 的价值不在于替代 `gold`，而在于帮助模型看到更多样本空间，但必须低权重、分来源处理。

### 0.4 split 是什么

`split` 指训练集和测试集怎么切分。  
这不是单纯的工程细节，而是在定义“模型到底在测什么能力”。

- `random split`：随机打散样本，最容易，但最可能高估效果。
- `drug split`：测试集中是训练时没见过的药。
- `scaffold split`：测试集中是训练时没见过的化学骨架，通常更接近新药外推。
- `microbe split`：测试集中是训练时没见过的菌或菌背景。

### 0.5 scaffold 是什么

这里的 `scaffold` 主要指 `Murcko scaffold`，也就是分子的核心骨架。

它不是完整分子本身，而是把外围取代基简化后保留下来的“结构主干”。  
做 `scaffold split` 的目的，是防止训练集和测试集虽然药名不同，但骨架几乎一样，从而让模型分数看起来不真实地偏高。

### 0.6 applicability 是什么

`applicability` 可以理解为“这个预测落不落在模型比较熟悉的区域里”。

例如：

- 这个药和训练集药物在指纹上很像；
- 这个骨架模型见过；
- 这个菌的门/属训练里也见过；

那么这个预测通常更可信。反过来，如果药物骨架很新、菌的分类层级也陌生，那么模型虽然仍然会输出结果，但不确定性会更高。

### 0.7 heuristic 是什么

`heuristic` 常译为“启发式”。  
它表示这个指标是为了排序、比较、原型系统运行而设计的，不等同于严格的机制真值或临床终点。

例如 Step 3 里的：

- `GMWI2-like heuristic`
- `development_score`

都应理解为“内部比较指标”，而不是临床结论。

### 0.8 panel-proxy 是什么

`panel-proxy` 表示当前 Step 3 使用的是“已有菌面板上的代理群落”，而不是直接来自真实 cohort 样本的群落丰度。

所以 Step 3 当前的定位是：

- 已经能把 Step 1 / Step 2 接进动态系统；
- 默认演示仍以 `panel-proxy` 为主；
- 但代码层已经支持用真实 `community_table.csv` 初始化群落。

### 0.9 AMR 修正是什么

`AMR` 是 `antimicrobial resistance`，也就是抗微生物药物抗性。

在这个项目里，AMR 修正不是重新训练一个独立模型，而是先引入一层“与已知耐药规律对照的知识层”，专门压住那些统计模型容易犯、但生物学上明显不合理的高风险假阳性。

例如：

- `Bacteroides × penicillin / beta-lactam`
- `gram-negative × vancomycin`
- `strict anaerobes × aminoglycoside`
- `gram-positive × polymyxin`

当前实现里，这层 AMR 修正主要服务于网页与接口预测解释，不是说所有离线预测表都已经永久写死了修正后的值。

## 1. 项目目标与为什么要拆成三步

本项目希望回答一个最终问题：在不同肠道微生物环境中，给定一个药物后，菌群会如何变化、药物会不会被代谢、最终这种变化对肠道健康和药物开发价值意味着什么。

这个问题不能直接用一个端到端模型一次性解决，主要有 4 个原因：

1. 数据来源的粒度不一致。  
   药物对微生物影响的数据通常是 `drug x strain` 级别，药物代谢数据也是 `drug x strain` 或 `drug x community` 级别，而健康结局和人群队列数据又是 `community x subject` 级别。
2. 标签类型不一致。  
   Step 1 更偏生长抑制/促进效应，Step 2 更偏代谢转化，Step 3 则是时间动态和系统级指标。
3. 真实公开数据并不支持“药物结构 -> 最终临床健康收益”这种直接监督。
4. 从工程上，先把 pair-level 预测器做稳，再把它们接进群落模拟器，比直接做超大一体化模型更可复现，也更容易解释。

因此，当前系统采用三步拆解：

1. `Step 1: drug -> microbe effect`  
   预测药物是否抑制、促进或不影响某个微生物，并尽量输出连续效应值。
2. `Step 2: microbe -> drug metabolism`  
   预测某个微生物是否会代谢该药物，以及代谢强度和候选产物线索。
3. `Step 3: community simulation`  
   在给定起始群落和给药条件下，利用 Step 1 和 Step 2 的 pair-level 预测结果进行时间推演，输出群落变化、健康指数和启发式开发评分。

换句话说，这三步分别在回答三个层级不同的问题：

- Step 1：药物会不会直接影响某个菌；
- Step 2：某个菌会不会反过来代谢这个药；
- Step 3：如果把很多菌放在一起并持续给药，整个系统最后会走向哪里。

## 2. 当前总体架构

当前仓库的主流程可以概括为：

```text
真实开源数据 / 论文补充表
    -> source-specific 标准化
    -> 药物与微生物统一 schema
    -> Step 1 药物影响模型
    -> Step 1 hybrid 输出
    -> Step 1 drug-profile-aware realism constraints
    -> AMR / 机制先验修正与解释层
    -> 疾病目录标准化与 IBS / IBS-D / IBS-C 候选补齐
    -> Step 2 候选 pair 组装 + 代谢模型
    -> reaction class / product / gene evidence 投影
    -> 集成 pair-level 预测表
    -> Step 3 群落时间推演
    -> 真实 cohort 初始化 / TCG proxy 健康签名
    -> 健康指数与开发评分
```

当前对应的主文档与入口如下：

- 总蓝图：[project_blueprint.md](project_blueprint.md)
- Step 1 细化说明：[step1_pipeline.md](step1_pipeline.md)
- Step 1 Chemprop 路线：[step1_chemprop.md](step1_chemprop.md)
- Step 1 hybrid 接口：[step1_hybrid.md](step1_hybrid.md)
- Step 2 细化说明：[step2_pipeline.md](step2_pipeline.md)
- Step 3 细化说明：[step3_pipeline.md](step3_pipeline.md)
- 统一数据结构：[schemas.py](../src/gut_drug_microbiome/schemas.py)

## 3. 全流程共用的数据工程原则

### 3.1 只使用真实公开来源

当前系统优先采用公开论文补充数据、开放数据库和可追溯的标准化资源，不使用人工捏造数据作为主训练来源。核心来源包括：

- Step 1：`Maier 2018`, `MDIPID`, `MASI`
- Step 2：`Zimmermann 2019`
- Step 3：当前先用内置 `panel-proxy` 场景，后续计划接 `curatedMetagenomicData`, `GMrepo`, `GMWI2`, `Javdan 2020`, `AGORA2`

### 3.2 gold / silver 分层

这是整个项目最关键的工程原则之一。

- `gold`：统一实验体系、标签定义稳定、适合主训练和主评估。
- `silver`：文献整理型或弱监督型标签，适合增广训练、外部验证或低权重补充。

这样做的原因是，不同数据库之间实验条件、浓度、菌株命名、方向定义和证据强度差异很大。如果直接混成同权训练集，通常会降低泛化性能。

所以 `gold / silver` 不是一个好听的命名，而是一种显式承认“数据来源不平等”的工程做法。

### 3.3 先保留连续信号，再派生离散标签

本项目尽量保留原始连续端点，例如：

- Step 1 的 `effect_score`
- Step 2 的 `parent_depletion_fraction`

原因是 Step 3 的动态模拟更依赖连续强度，而不是只有硬标签。

### 3.4 评估优先模拟真实外推场景

随机切分只用于 sanity check，真正重要的是：

- `drug split`
- `scaffold split`
- `microbe split`

这是因为项目最终关心的是“新药是否可用”和“换一个微生物背景还能不能推得动”，而不是只在随机打散后的同分布测试集上拿高分。

### 3.5 统一 schema

项目的核心实体和记录已经统一定义在 [schemas.py](../src/gut_drug_microbiome/schemas.py) 中，包括：

- `DrugEntity`
- `MicrobeEntity`
- `DrugMicrobeEffectRecord`
- `MicrobeDrugMetabolismRecord`
- `CommunityState`
- `SimulationResult`

这层 schema 的作用是保证不同来源的数据在进入建模前先被压到统一结构，减少后续脚本之间的隐式耦合。

这里的“隐式耦合”指的是：某个脚本偷偷依赖另一个脚本生成的特殊列名或特殊格式，但这种依赖没有写成稳定接口。  
统一 schema 的意义，就是尽量把这种脆弱依赖收束到显式的数据层里。

### 3.6 知识修正尽量写成显式数据层

这条原则是最近一轮迭代里新增的。

原因很简单：如果某个纠偏逻辑只藏在代码里，后面很难解释“为什么修了、修了哪些、依据是什么”；但如果它被拆成显式表和独立脚本，就可以持续扩展、复核和引用。

当前已经落地的两类知识层就是：

- `AMR` 修正层
  - [microbe_amr_reference.csv](../data/processed/amr/microbe_amr_reference.csv)
  - [drug_resistance_rules.csv](../data/processed/amr/drug_resistance_rules.csv)
  - [species_drugclass_phenotype_prior.csv](../data/processed/amr/species_drugclass_phenotype_prior.csv)
- `health signature / TCG proxy` 层
  - [microbe_reference_genome_mapping.csv](../data/processed/health_signature/microbe_reference_genome_mapping.csv)
  - [microbe_tcg_proxy_mapping.csv](../data/processed/health_signature/microbe_tcg_proxy_mapping.csv)
  - [health_signature_target_genomes.csv](../data/processed/health_signature/health_signature_target_genomes.csv)

这样做的好处是：

1. 文档、代码和数据表可以互相追溯；
2. 后续如果要换规则、补来源或升级机制层，不需要把逻辑重新散落到多个脚本里；
3. 网页展示时也能把“原始预测”和“知识修正后结果”分开呈现。

## 4. Step 1 技术路线：药物对微生物影响建模

### 4.1 Step 1 要解决什么问题

输入一个药物和一个微生物，预测：

- `inhibit`
- `promote`
- `no_effect`

同时尽量输出连续效应值 `effect_score`，供 Step 3 使用。

### 4.2 为什么选 Maier 2018 做主金标准

Step 1 当前主金标准是 `Maier et al., 2018`，原因有三点：

1. 数据规模足够大，覆盖 `1197 drugs x 40 gut strains` 的系统筛选。
2. 实验体系相对统一，适合作为 first-principles baseline。
3. 已经有公开 benchmark 论文可作为参考对照，例如 `McCoubrey 2021` 和 `Wang 2023`。

在当前实现中，Step 1 的主处理脚本为：

- [download_step1_data.py](../scripts/download_step1_data.py)
- [normalize_step1_data.py](../scripts/normalize_step1_data.py)
- [train_step1_baseline.py](../scripts/train_step1_baseline.py)

生成的核心表包括：

- [step1_drug_table.csv](../data/processed/step1/step1_drug_table.csv)
- [step1_microbe_table.csv](../data/processed/step1/step1_microbe_table.csv)
- [step1_interactions.csv](../data/processed/step1/step1_interactions.csv)
- [step1_modeling_table.csv](../data/processed/step1/step1_modeling_table.csv)

### 4.3 Step 1 标签构建逻辑

当前标签规则来自 [labeling_rules.yaml](../configs/labeling_rules.yaml) 和 [step1_pipeline.md](step1_pipeline.md)，核心定义为：

- 连续值：`effect_score = reported_norm_auc - 1.0`
- `inhibit`: `effect_score <= -0.20` 且 `q < 0.05`
- `promote`: `effect_score >= 0.20` 且 `q < 0.05`
- `no_effect`: 其他情况

这里要特别强调一个方法学边界：

- `promote` 不是 Maier 原文直接给出的天然三分类，而是根据连续生长值和显著性阈值派生出来的工作标签。

这也是为什么当前公开数据几乎没有显著 `promote` 样本，Step 1 实际上在很多设置下会退化为 `inhibit / no_effect` 的主任务。

这里顺手解释两个统计名词：

- `reported_norm_auc` 可以粗略理解成“标准化后的生长曲线面积”，减去 `1.0` 后更方便表示偏抑制还是偏促进；
- `q < 0.05` 可以理解成“多重检验校正后的显著性阈值”，目的是在大量 pair 同时检验时控制假阳性。

### 4.4 Step 1 特征工程

Step 1 当前的 baseline 主要是表格学习，但已经接入比较完整的药物化学特征与菌株侧元数据。

### 药物侧特征

基础药物特征包括：

- 分子量、`xlogp`、`tpsa`、`complexity`、`volume3d`
- 给药与暴露相关字段：`dose_umol`、肠道估计浓度、血浆浓度、粪便/尿液排泄比例
- SMILES 简单统计特征：长度、支链数、双键数、卤素数、环索引数
- 药理分类：`therapeutic_class`, `therapeutic_effect`, `ATC`

进一步，项目通过 [chem_features.py](../src/gut_drug_microbiome/step1/chem_features.py) 接入了 RDKit 化学特征：

- `canonical_smiles_rdkit`
- `InChIKey`
- `Murcko scaffold`
- `rdkit_formula`
- `ExactMolWt`
- `MolLogP`
- `TPSA`
- `MolMR`
- `formal_charge`
- `heavy_atom_count`
- `HBD / HBA`
- `rotatable_bond_count`
- `ring_count`
- `aromatic_ring_count`
- `aliphatic_ring_count`
- `hetero_atom_count`
- `fraction_csp3`
- `256-bit Morgan fingerprint`

几个容易让非化学背景读者卡住的名词解释如下：

- `RDKit`：常用的开源化学信息学工具包，用来从 SMILES 计算结构特征。
- `TPSA`：拓扑极性表面积，常用来粗略反映分子的极性和跨膜能力。
- `MolLogP / xlogp`：常用来描述脂溶性。
- `Murcko scaffold`：分子的核心骨架。
- `Morgan fingerprint`：把分子的局部结构编码成一串二进制位向量，便于模型比较两个分子在结构上是否相似。

### 微生物侧特征

当前主用的是 Step 1 可直接获得的菌株元数据：

- `species_label`
- `phylum / class / order / family / genus`
- `gram_stain`
- `medium_preference`
- `biosafety`

这是一套偏“可复现、低假设”的特征方案。更机制化的基因组、代谢网络和酶特征计划在后续与 `AGORA2`、`gutMGene` 结合后加入。

### 4.5 Step 1 baseline 模型为什么选 ExtraTrees

当前主 baseline 实现在 [train_baseline.py](../src/gut_drug_microbiome/step1/train_baseline.py)，核心模型为：

- 分类：`ExtraTreesClassifier`
- 回归：`ExtraTreesRegressor`

选择它的原因是：

1. 对数值 + 类别拼接后的宽表特征适配度高；
2. 对非线性和高维稀疏指纹特征比较稳；
3. 对小样本多类别菌株数据不需要太重的调参；
4. 训练速度快，适合大量 split 对比实验。

预处理使用：

- 数值列：`median imputation`
- 类别列：`constant imputation + one-hot`

这里的：

- `median imputation` 指缺失的数值列用中位数补；
- `one-hot` 指把一个类别变量拆成多个 0/1 指示列。

ExtraTrees 这类树模型对这种“数值列 + one-hot 类别列”的混合宽表输入比较友好。

### 4.6 为什么要引入 MDIPID 和 MASI

Step 1 只有 Maier gold 数据时，最大问题不是不能训练，而是对新场景的覆盖有限。因此我们接入了两类银标：

- `MDIPID DEIM`
- `MASI v2.0`

但这两类银标的用途并不一样。

### MDIPID

`MDIPID DEIM` 的方向性比较适合 Step 1，当前作为主要 `silver` 增强源。它只并入分类训练，不参与回归训练和主测试集评估。

### MASI

`MASI` 的覆盖很大，但异质性更强。实际实验后发现：

- 全量 `MASI drug-like` 并入会明显拉低表现；
- `MASI curated` 子集，也就是 `drug-like + Gut + (Human or In vitro)`，更稳，但仍不如 `MDIPID only`。

因此，当前推荐策略是：

- `MDIPID` 作为主银标；
- `MASI curated` 作为低权重补充或外部验证；
- 不建议把全量 `MASI` 直接并入主训练。

### 4.7 silver 的 source-aware weighting

为了降低不同银标源之间的噪声冲突，Step 1 训练支持按来源加权。当前代码允许：

- 对 `gold` 和 `silver` 设置默认权重；
- 对特定来源单独设置权重，例如 `mdipid_deim=1.0`, `masi_v2=0.25`

实验结果说明一个重要结论：

- 银标不是越多越好；
- 更重要的是“分层”和“加权”；
- 目前最稳的主配方依然是 `gold + MDIPID`。

### 4.8 Step 1 评估设计

当前 Step 1 保留三类切分：

- `drug split`
- `scaffold split`
- `microbe split`

它们分别对应的真实问题是：

- `drug split`：见过这些菌，但没见过这批药；
- `scaffold split`：没见过这类化学骨架的新药；
- `microbe split`：换一个微生物背景是否还能泛化。

评估指标也建议在这里统一理解一下：

- `balanced accuracy`：对类别不平衡更稳，不会因为多数类太多而虚高；
- `macro-F1`：先分别算每一类的 F1，再做平均，能更公平地反映少数类；
- `ROC-AUC`：衡量分类模型整体排序能力；
- `PR-AUC`：在正类较少时通常更敏感；
- `R²`：回归拟合优度，越接近 1 越好，小于 0 往往说明效果很弱；
- `Spearman`：看排序相关性，更强调“大小顺序对不对”。

### 4.9 Step 1 当前结果与默认路线

截至当前版本，Step 1 的结论已经比较稳定：

- `gold + MDIPID, drug split + RDKit`
  - balanced accuracy: `0.7938`
  - macro-F1: `0.8265`
- `gold-only, scaffold split + RDKit`
  - balanced accuracy: `0.8063`
  - macro-F1: `0.8385`
  - regression `R² = 0.6169`
  - regression `Spearman = 0.4261`
- `gold + MDIPID, scaffold split + RDKit`
  - balanced accuracy: `0.7996`
  - macro-F1: `0.8299`
- `gold-only, microbe split + RDKit`
  - balanced accuracy: `0.9321`
  - macro-F1: `0.8426`

这说明：

1. RDKit 化学特征已经显著提升了 Step 1 的基线能力；
2. `MDIPID` 对分类有小幅稳定增益；
3. `MASI` 更适合做谨慎补充，而不是无脑扩库。

这部分结果其实说明了一个很重要的现实：  
在生物医药建模里，“更多数据”并不自动等于“更好模型”。如果新增数据带来的主要是体系差异和标签噪声，效果反而可能下降。

### 4.10 为什么又上 Chemprop

Step 1 的 RDKit 路线已经够用，但对真正的新骨架外推，还需要更强的分子表示。因此新增了 Chemprop 路线：

- 数据准备：[prepare_step1_chemprop.py](../scripts/prepare_step1_chemprop.py)
- 训练：[train_step1_chemprop.py](../scripts/train_step1_chemprop.py)
- 说明：[step1_chemprop.md](step1_chemprop.md)

当前 Chemprop 路线采用：

- 输入：SMILES 图表示
- 辅助输入：描述符矩阵
- 切分：`scaffold split`

`Chemprop` 可以把它简单理解成“直接在分子图上学习表示”的模型。  
相比只吃表格描述符的模型，它通常更擅长捕捉新骨架上的结构模式，因此很适合作为 Step 1 hybrid 里的分类部分。

当前结果：

- 分类 `ROC-AUC = 0.9007`
- 分类 `PR-AUC = 0.7545`
- 分类 `balanced accuracy = 0.8513`

这条线已经在分类端超过了当前的 RDKit 表格 baseline，但回归端仍不如 `RDKit + ExtraTrees`。

### 4.11 Step 1 hybrid：为什么分类和回归分开走

当前 Step 1 的最优工程方案并不是单一模型，而是组合路线：

- 分类：`Chemprop`
- 回归：`RDKit + ExtraTrees`

这一套逻辑实现在 [hybrid.py](../src/gut_drug_microbiome/step1/hybrid.py)。

Step 1 hybrid 的规则是：

1. 先由 Chemprop 给出 `predicted_inhibit_probability`
2. 再由 RDKit 回归器给出 `predicted_effect_score`
3. 根据阈值组合为最终的 `predicted_effect_label_hybrid`

这样做的原因是：

- Chemprop 在新骨架分类外推上更强；
- RDKit 回归器对连续强度更稳；
- Step 3 同时需要方向和强度，因此 hybrid 输出更适合做下游输入。

当前 hybrid 全量预测输出在：

- [predictions.csv](../predictions/step1/hybrid_scaffold_v1/predictions.csv)
- [predictions_slim.csv](../predictions/step1/hybrid_scaffold_v1/predictions_slim.csv)
- [summary.json](../predictions/step1/hybrid_scaffold_v1/summary.json)

当前共覆盖 `43,109` 个 drug-microbe pair，其中：

- `7,361 inhibit`
- `157 promote`
- `35,591 no_effect`

这张表就是 Step 2 和 Step 3 的直接上游输入之一。

### 4.12 Step 1 新增的 AMR 修正与解释层

Step 1 的统计模型本身仍然是“药物结构 + 微生物元数据 -> 作用方向/强度”的 pair-level 预测器。  
但在实际网页使用中，我们已经发现一个很现实的问题：某些抗生素-菌组合即使统计上被打成 `inhibit`，也可能与已知内在耐药规律明显冲突。

因此，当前系统在 Step 1 展示层之外，额外加了一层服务时 AMR 修正引擎，核心实现位于：

- [amr.py](../src/gut_drug_microbiome/amr.py)

它做的不是“替代模型输出”，而是：

1. 先识别药物上下文，例如：
   - `penicillin / beta_lactam`
   - `vancomycin / glycopeptide`
   - `daptomycin / lipopeptide`
   - `aminoglycoside`
   - `polymyxin`
   - 低厌氧覆盖 `fluoroquinolone`
2. 再按 `species / genus / gram stain / strict anaerobe` 规则，判断这个 pair 是否属于已知高风险假阳性区域；
3. 如果命中规则，就：
   - 下调 `predicted_inhibit_probability`
   - 回退 `predicted_effect_label`
   - 同时保留原始值和修正值，方便前端对照展示。

当前第一批已经落地的规则包括：

- `Bacteroides × penicillin / beta-lactam`
- `gram-negative × vancomycin / glycopeptide`
- `gram-negative × daptomycin / lipopeptide`
- `strict anaerobes × aminoglycoside`
- `gram-positive × polymyxin`
- `strict anaerobes × low-anaerobe-coverage fluoroquinolone`

相关脚本和数据层包括：

- [prepare_amr_reference_tables.py](../scripts/prepare_amr_reference_tables.py)
- [seed_beta_lactam_amr_rules.py](../scripts/seed_beta_lactam_amr_rules.py)
- [seed_high_risk_antibiotic_amr_rules.py](../scripts/seed_high_risk_antibiotic_amr_rules.py)
- [seed_aminoglycoside_anaerobe_rules.py](../scripts/seed_aminoglycoside_anaerobe_rules.py)
- [seed_polymyxin_and_low_anaerobe_fq_rules.py](../scripts/seed_polymyxin_and_low_anaerobe_fq_rules.py)

当前 `drug_resistance_rules.csv` 已经有 `263` 条规则，`microbe_amr_reference.csv` 已覆盖当前 `83` 菌面板。  
需要明确的是：这层修正现在主要体现在网页 / API 返回与交互式解释中，默认离线预测表不一定已经把修正结果永久写回。

### 4.13 Step 1 的 drug-profile-aware realism constraints（新增）

为避免“统计可行但机制不合理”的 Step 1 输出，当前在 `hybrid` 推理末端新增了药物类型感知约束，代码位于 [hybrid.py](../src/gut_drug_microbiome/step1/hybrid.py)。

新增两类 profile：

- `eubiotic_modulator`
- `host_pathway_agent`

核心约束逻辑：

1. `eubiotic_modulator`（代表：Rifaximin）  
   对核心产丁酸菌（`Faecalibacterium prausnitzii`, `Roseburia spp.`, `Eubacterium rectale`）的“强抑制”进行惩罚/裁剪，避免出现“核心产丁酸菌几乎全部强抑制”的非现实模式。
2. `host_pathway_agent`（代表：Lubiprostone）  
   全局下调直接微生物效应，默认偏向 `no_effect`，仅保留高证据尾部效应。

这层约束的工程定位是“推理后处理现实性约束”，它不会改训练权重，但会直接影响下游机制层聚合和疾病排序。

对应输出字段也已补齐：

- `step1_drug_profile`
- `step1_constraint_applied`
- `step1_constraint_reason`
- `summary.json` 中的 `step1_drug_profile_counts` 与 `step1_constraint_summary`

## 5. Step 2 技术路线：微生物对药物代谢建模

### 5.1 Step 2 要解决什么问题

Step 2 的目标不是一开始就做“精确产物生成器”，而是按难度分层推进：

1. 先预测 `会不会被代谢`
2. 再预测 `代谢强度`
3. 再补 `反应类型`
4. 再补 `候选产物` 和 `gene evidence`

这样拆分的原因是，公开数据里最稳定的监督信号是“母药是否被消耗”，而不是“精确产物结构”。

### 5.2 为什么选 Zimmermann 2019 做 Step 2 gold

Step 2 当前主金标准是 `Zimmermann et al., 2019` 的补充 workbook。原因是：

1. 它提供了 `271` 个药物和 `76` 个菌株的系统代谢筛选；
2. 粒度是标准的 `drug x isolate`，很适合 pair-level 建模；
3. 同时给出了候选代谢物和部分基因线索，方便后续机制增强。

核心标准化逻辑在 [zimmermann_2019.py](../src/gut_drug_microbiome/step2/zimmermann_2019.py) 中，主要解析：

- `Supplementary Table 1`: 菌株信息
- `Supplementary Table 2`: 药物信息
- `Supplementary Table 3`: 药物在不同菌株中的 parent depletion 主筛结果
- 候选代谢物相关表
- 药物-基因证据表

### 5.3 Step 2 source-specific 标准化过程

Step 2 不是直接拿原始 Excel 做建模，而是先压成标准标签表。当前核心脚本为：

- [normalize_step2_zimmermann.py](../scripts/normalize_step2_zimmermann.py)
- [normalize.py](../src/gut_drug_microbiome/step2/normalize.py)

标准 schema 至少包括：

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

当前标准化后的主要文件为：

- [zimmermann_2019_label_table.csv](../data/processed/step2/zimmermann_2019/zimmermann_2019_label_table.csv)
- [zimmermann_2019_drug_table.csv](../data/processed/step2/zimmermann_2019/zimmermann_2019_drug_table.csv)
- [zimmermann_2019_microbe_table.csv](../data/processed/step2/zimmermann_2019/zimmermann_2019_microbe_table.csv)
- [zimmermann_2019_metabolite_candidates.csv](../data/processed/step2/zimmermann_2019/zimmermann_2019_metabolite_candidates.csv)
- [zimmermann_2019_gene_links.csv](../data/processed/step2/zimmermann_2019/zimmermann_2019_gene_links.csv)
- [zimmermann_2019_modeling_table.csv](../data/processed/step2/zimmermann_2019/zimmermann_2019_modeling_table.csv)

当前统计结果：

- `20,596` 个 isolate-level pair
- `271` 个药物
- `76` 个菌株
- `2,627 metabolized`
- `17,969 not_metabolized`
- `2,911` 个 pair 带候选代谢物
- `6,183` 个唯一候选代谢产物
- `89` 条 drug-gene links
- `1,520` 个 pair 带 gene evidence

### 5.4 Step 2 标签定义

当前 Step 2 主标签为：

- `metabolized`
- `not_metabolized`
- `uncertain`

当前连续端点为：

- `parent_depletion_fraction`

这里要把两个概念分开理解：

- `metabolized / not_metabolized`：离散事件，回答“有没有发生代谢”；
- `parent_depletion_fraction`：连续强度，回答“母药被消耗了多少”。

Step 3 更依赖后者，因为动态模拟不仅要知道“会不会被代谢”，还要知道“代谢得多快、多强”。

在 `Zimmermann 2019` 的标准化中，主代谢标签依赖：

- `% consumed`
- adaptive threshold
- `p(FDR)`

连续值当前按“母药消耗比例”近似保留，用于回归建模和 Step 3 的药物动态更新。

### 5.5 Step 2 候选输入表如何来自 Step 1

Step 2 并不是只在 Gold 表里做闭环训练，还会把 Step 1 的全量 hybrid 预测接成“候选推理空间”。这一部分由 [assemble.py](../src/gut_drug_microbiome/step2/assemble.py) 负责。

处理逻辑是：

1. 读取 Step 1 hybrid 输出；
2. 统一 `pair_id`；
3. 将 Step 1 预测字段重命名为 `step1_*` 前缀；
4. 生成 `step2_candidate_pairs_full.csv`
5. 如果有标准化后的 Step 2 标签，再合并成 `step2_modeling_table.csv`

当前生成的输入表包括：

- [step2_candidate_pairs_full.csv](../data/processed/step2/step2_candidate_pairs_full.csv)
- [step2_candidate_pairs_slim.csv](../data/processed/step2/step2_candidate_pairs_slim.csv)
- [step2_modeling_table.csv](../data/processed/step2/step2_modeling_table.csv)

这样做的好处是，Step 2 和 Step 3 可以直接继承 Step 1 的 pair-level 结果，而不需要重新做药物和微生物匹配。

从工程角度看，这一步相当于把“上游预测结果”变成了“下游模型可直接消费的特征表”，是全链路能否真正串起来的关键桥接层。

### 5.6 Step 2 baseline 特征与模型

Step 2 当前 baseline 与 Step 1 的思路类似，仍然使用表格学习，但特征更偏代谢相关的结构 + 分类信息。

当前 baseline 代码在 [train_baseline.py](../src/gut_drug_microbiome/step2/train_baseline.py) 中，核心模型为：

- 分类：`ExtraTreesClassifier`
- 回归：`ExtraTreesRegressor`

当前主要特征包括：

- 药物基础理化特征
- RDKit 描述符
- `Morgan fingerprint`
- `Murcko scaffold`
- 药理分类和 ATC
- 菌株 `species / genus / phylum`
- 从 `Zimmermann` 解析出的物种名称与门水平描述

### 5.7 Step 2 评估和 applicability 设计

Step 2 当前保留：

- `random split`
- `drug split`
- `scaffold split`
- `microbe split`

其中当前默认部署的是 `scaffold split + full-fit`，因为它最接近“面对未见化学骨架新药”时的外推场景。

此外，Step 2 已经引入了可用性边界估计，也就是 `applicability`。当前参考信号包括：

- 药物与训练集的指纹相似度
- scaffold 是否在训练集中出现过
- 微生物 genus / phylum 是否在训练集中出现过

这层设计很重要，因为 Step 3 的开发评分里，已经把低可用性区域作为 `uncertainty penalty` 的一部分。

可以把这些可用性信号理解成三种不同层次的“熟悉度”：

- 指纹相似度：这个药在局部结构上和训练药像不像；
- `scaffold_seen_in_training`：这个核心骨架以前见没见过；
- `microbe_genus/phylum_seen_in_training`：这个菌的分类层级训练里熟不熟。

### 5.8 Step 2 当前结果与默认路线

截至当前版本，Step 2 的主结果如下。

### scaffold split

- balanced accuracy: `0.6982`
- macro-F1: `0.6367`
- ROC-AUC: `0.7740`
- PR-AUC: `0.3160`
- regression `RMSE = 0.2409`
- regression `MAE = 0.1278`
- regression `R² = -0.1098`
- regression `Spearman = 0.3897`

### drug split

- balanced accuracy: `0.7123`
- macro-F1: `0.6522`
- ROC-AUC: `0.8171`
- regression `R² = 0.0842`
- regression `Spearman = 0.4914`

### microbe split

- balanced accuracy: `0.8686`
- macro-F1: `0.7485`
- ROC-AUC: `0.9486`
- regression `R² = 0.7355`
- regression `Spearman = 0.5745`

这些结果说明：

1. 对新药骨架的代谢外推仍然比 Step 1 更难；
2. 跨微生物泛化反而相对稳定，说明当前微生物侧特征已经足以抓住一部分门/属层级信号；
3. Step 2 目前最稳的是“会不会被代谢”和“代谢强度”，而不是精确产物结构。

这组结果的结构也很有代表性：

- `microbe split` 相对较高，说明当前门/属层面的菌特征已经能抓住一部分代谢模式；
- `scaffold split` 更难，说明真正的新化学骨架外推仍然是 Step 2 的主要瓶颈；
- 回归 `R²` 在最难切分上偏低，说明“代谢强度”的绝对数值还没有学得足够稳。

### 5.9 Step 2 新增的半机制投影层

在当前版本里，Step 2 已经不再是“只有发生与否 + 强度”的纯统计输出。  
在 baseline 预测之外，我们新增了一层半机制投影模块，核心实现位于：

- [mechanism.py](../src/gut_drug_microbiome/step2/mechanism.py)

它的思路不是重新训练一个产物生成模型，而是利用 `Zimmermann 2019` 已解析出的：

- `reaction_class`
- `product_ids`
- `evidence_gene_ids`

去构建一个分层支持库，然后把新的代谢阳性 pair 投影到最相近的机制证据上。

当前投影支持的层级包括：

- `drug_nt`
- `drug_species`
- `drug_genus`
- `scaffold_species`
- `scaffold_genus`
- `class_genus`
- `atc_l1_genus`
- `drug_global`
- `scaffold_global`
- `genus_global`

可以把这套东西理解成一种“加权邻域投影”：

- 如果新 pair 和历史上同药、同菌或同骨架、同属的代谢阳性记录很接近，就更容易继承它们的 `reaction class / product / gene evidence`；
- 如果只能命中更宽泛的 `class_genus` 或 `global` 层级，那就仍然能给出线索，但可信度会更保守。

这层输出当前主要包括：

- `predicted_mechanism_projection_flag`
- `predicted_reaction_class`
- `predicted_reaction_confidence`
- `predicted_candidate_product_ids`
- `predicted_evidence_gene_ids`

需要强调的是，这仍然不是“精确产物生成器”，而是一层半机制线索投影器。  
它的定位是：

- 比“只有会不会代谢”更可解释；
- 能为 Step 2 和网页解释提供 `reaction / product / gene` 证据线索；
- 但还不能替代真正基于反应规则或代谢网络的产物预测系统。

所以当前 Step 2 的工程定位可以更新为：

- 已完成 `代谢发生与否` 和 `代谢强度` 的稳定 baseline；
- 已新增 `reaction class / product / gene evidence` 的半机制投影层；
- 但“会精确代谢成什么结构”仍处于待增强阶段。

### 5.10 Step 2 全量候选预测输出

当前默认部署模型已经跑在全量候选对上，输出在：

- [predictions.csv](../predictions/step2/baseline_scaffold_v1_83/predictions.csv)
- [predictions_slim.csv](../predictions/step2/baseline_scaffold_v1_83/predictions_slim.csv)
- [summary.json](../predictions/step2/baseline_scaffold_v1_83/summary.json)

当前统计：

- `89,557` 个候选 pair
- `1,079` 个药物
- `83` 个微生物
- `9,134 predicted metabolized`
- `80,423 predicted not_metabolized`
- `41,123` 个 pair 位于当前 applicability 范围内

这张 83 菌扩展集成表就是当前 Step 3 的主输入。  
需要补充说明的是：

- 离线 `predictions.csv` 主要保存 Step 1/Step 2 的基础集成结果；
- `AMR` 修正与 Step 2 机制投影目前更多是在服务层按需注释，用于网页/API 解释和交互展示；
- 因此“离线预测表”和“网页展示结果”在解释字段上可能会比基础列更丰富。

## 6. Step 3 技术路线：肠道群落推演与开发评分

### 6.1 Step 3 要解决什么问题

Step 1 和 Step 2 都是 pair-level 模型，但项目最终关心的是“在一个真实或近似真实的菌群环境里，加药后群落会怎么变”。这一步不能只看单个 `drug x microbe` 关系，而必须把它们在时间轴上耦合起来。

因此 Step 3 的核心任务是：

1. 设定一个起始菌群状态；
2. 加入给药方案；
3. 利用 Step 1 影响和 Step 2 代谢预测更新系统状态；
4. 输出时间轨迹和健康相关指标。

### 6.2 为什么当前先做离散时间模拟器

当前 Step 3 并没有直接上全约束代谢网络或全微分方程系统，而是先实现一个轻量离散时间模型。原因是：

1. 可以先把 Step 1 / Step 2 的 pair-level 结果真正跑起来；
2. 更容易解释每个模块在系统中的作用；
3. 以后接 `AGORA2`, `MICOM`, `mgPipe` 时可以平滑升级，而不是推倒重来。

核心实现位于 [simulation.py](../src/gut_drug_microbiome/step3/simulation.py)。

### 6.3 Step 3 当前输入

当前 Step 3 直接消费 Step 2 的集成预测表：

- [predictions.csv](../predictions/step2/baseline_scaffold_v1_83/predictions.csv)

这张表已经同时包含：

- Step 1 输出：
  - `step1_predicted_inhibit_probability`
  - `step1_predicted_effect_score`
  - `step1_predicted_effect_label_hybrid`
- Step 2 输出：
  - `predicted_metabolized_probability`
  - `predicted_parent_depletion_fraction`
  - `applicability_flag`
  - `drug_max_fingerprint_jaccard`

也就是说，Step 3 不需要重新分别读取 Step 1 和 Step 2，而是直接建立在统一集成表之上。

和早期版本相比，这里有两个重要更新：

1. 当前网页和默认服务已经统一切到 `83` 菌扩展面板，而不是最初的 `40` 菌库内模式。
2. 如果用户提供真实 `community_table.csv`，Step 3 可以跳过内置场景，直接用真实 cohort abundance 初始化群落。

### 6.4 Step 3 的状态变量与更新逻辑

当前 Step 3 使用的核心状态包括：

- 微生物丰度
- 母药浓度
- 聚合代谢物池
- 健康指数

可以把这 4 个状态理解成两类对象：

- “系统本身的组成”：
  - 微生物丰度
- “药物相关状态”：
  - 母药浓度
  - 聚合代谢物池
- “系统的评价结果”：
  - 健康指数

也就是说，Step 3 一边在更新系统状态，一边在更新系统的评价指标。

直观上，它做的是下面几件事：

1. 用 Step 1 的 `predicted_effect_score` 调整每个菌的增长/抑制趋势；
2. 用 `ecology_strength` 把群落缓慢拉回起始背景，避免系统发散；
3. 用一个轻量的 `interaction-aware` 层，把群落内的正向交互与竞争交互接进动力学；
4. 用 Step 2 的 `predicted_metabolized_probability` 和 `predicted_parent_depletion_fraction` 近似药物代谢消耗；
5. 用一个简化的 `aggregate metabolite pool` 累积代谢负荷；
6. 每个时间点计算健康指数和开发评分。

最近一轮升级后，这里的核心变化是：

- Step 3 不再只把群落看成“很多独立菌 + 一个药物压力项”
- 现在还会额外近似计算：
  - 哪些菌更像 `producer / cross-feeding source`
  - 哪些菌更像 `consumer / beneficiary`
  - 哪些菌之间更可能因为酶功能、反应类和分类重叠而产生竞争
- 这个思路主要受《Imbalance in gut microbial interactions as a marker of health and disease》启发
  - 健康态更接近 `competition-dominated`
  - dysbiosis-like 状态更接近 `positive-interaction / cross-feeding dominated`
  - 当前仓库实现的是 `ENBI-like proxy`，不是原论文 interaction inference 的全文复现

从公式上可以理解为：

```text
microbe_abundance(t+1)
    = normalize(
        microbe_abundance(t)
        * exp(
            effect_scale * step1_effect * drug_exposure
            + ecology_pull
            + interaction_delta
          )
      )

parent_drug(t+1)
    = parent_drug(t) + dose_input - metabolism_loss - clearance
    where metabolism_loss is mildly boosted in cross-feeding-dominant states

metabolite_pool(t+1)
    = metabolite_pool(t) + metabolism_gain - metabolite_clearance
```

这里几个关键参数的直观含义如下：

- `effect_scale`：药物对微生物抑制/促进作用被放大的程度；
- `ecology_strength`：群落被拉回初始生态背景的力度；
- `interaction_scale`：群落内部正向/负向交互反馈被放大的程度；
- `drug_clearance_rate`：母药自然清除速度；
- `product_clearance_rate`：代谢产物自然清除速度。

其中 `normalize(...)` 很重要，因为当前丰度是相对丰度，更新完以后通常要重新归一化，使所有菌丰度之和回到 1。

当前 `interaction_delta` 的数据来源不是单一表，而是多源近似拼出来的：

- 手工整理的 `cross_feeding_edges.csv`
- 化合物语义层 `compound_semantic_family / aliases / keywords`
- Step 2 的 `predicted_enzyme_ids / predicted_enzyme_step1_promote_support_score`
- Step 2 的 `predicted_metabolized_probability / predicted_parent_depletion_fraction`
- 反应类和酶集重叠，用来近似 resource overlap / competition

### 6.5 Step 3 当前内置场景

当前版本仍然内置了 4 个 `panel-proxy` 起始场景：

- `healthy_reference`
- `high_fiber`
- `high_fat`
- `antibiotic_perturbed`

这些场景不是直接从真实 cohort 回填出来的，而是基于当前可覆盖的 `83` 菌扩展面板构造的代理群落。实现方式是：

1. 先读取该药物可覆盖的菌株集合；
2. 根据场景为不同 genus / phylum 赋予乘数；
3. 归一化成起始相对丰度分布；
4. 将其作为模拟初始态。

这一步的优势是，现在就可以跑通系统；不足是，它仍然不是基于真实人群的初始群落。

这里的“按 genus / phylum 乘数构造场景”可以理解为一种工程近似：  
它不是在说真实世界里高脂饮食群落就一定等于某几个属乘 `1.4`、另几个属乘 `0.7`，而是在当前缺乏真实 cohort 初始化的前提下，先人为构造出方向不同的起始生态背景。

不过和旧版不同的是，当前 Step 3 已经支持真实 cohort 初始化。  
也就是说，`panel-proxy` 仍然是默认演示入口，但不再是唯一入口。

当前真实 cohort 接入链路包括：

- [prepare_step3_cohort_community.py](../scripts/prepare_step3_cohort_community.py)
- 网页/API 的 `community_table_path`
- `screen_step3_candidates.py --community-table ...`

实际流程是：

1. 把原始 abundance 表映射到现有 `83` 菌面板；
2. 导出标准 `community_table.csv`；
3. Step 3 用该文件直接初始化群落；
4. 此时网页中的场景比较会退化成“当前真实样本群落”的单场景模拟，而不是强制再跑 4 个内置场景。

### 6.6 健康指数如何计算

当前健康指数已经从单纯 `GMWI2-like heuristic`，升级成了 `interaction-aware heuristic`。  
它仍然不是官方 GMWI2 终值，但现在同时结合了“群落组成”和“群落交互状态”两层信息。

组成层分量包括：

- 多样性 `diversity`
- 有益菌比例 `beneficial_fraction`
- 风险菌比例 `risk_fraction`
- 与起始群落的稳定性 `stability`

交互层分量包括：

- `positive_interaction_strength`
- `negative_interaction_strength`
- `interaction_balance_rho`
- `interaction_balance_shift`
- `interaction_component`

在代码中，可以把它近似理解为：

```text
balance = clip(0.5 + beneficial_fraction - risk_fraction, 0, 1)
interaction_component
    = f(interaction_balance_rho, interaction_balance_shift, interaction_coverage)

health_index
    = 100 * (
        0.30 * diversity
        + 0.20 * balance
        + 0.20 * stability
        + 0.30 * interaction_component
      )
```

几个分量的白话解释如下：

- `diversity`：是不是被少数几个菌“垄断”了。越多样，通常越稳。
- `beneficial_fraction`：当前启发式定义的有益菌占比。
- `risk_fraction`：当前启发式定义的风险菌占比。
- `stability`：群落和起始状态相比偏离了多少。
- `balance`：把“有益菌多不多”和“风险菌高不高”压成一个平衡项。
- `interaction_balance_rho`：当前交互网络更偏正向还是更偏竞争。
  - `rho < 0` 更接近竞争主导
  - `rho > 0` 更接近交叉喂养主导
- `interaction_balance_shift`：相对起始群落，是否朝更 `dysbiosis-like` 的交互状态偏移。

这里再次强调：

- `health_index` 现在已经不是纯旧版 `GMWI2-like`
- 但 `interaction_balance_rho / shift` 也不是原论文 ENBI 的严格复现
- 更准确的说法是：当前 Step 3 用一个可运行的 `ENBI-like interaction proxy` 去增强原有健康分

当前有益菌和风险菌集合是手工定义的启发式名单，用于支持第一版排序和场景比较。

最近一轮迭代里，这一层又新增了一条并行路线：`TCG-inspired health signature proxy`。  
核心思路来自《A core microbiome signature as an indicator of health》，但当前实现仍是保守代理版，而不是原论文的严格 genome-level 复现。

已经落地的内容包括：

- `83菌 -> 参考基因组` 模板
- `83菌 -> guild / TCG proxy` 模板
- Step 3 读取 [microbe_tcg_proxy_mapping.csv](../data/processed/health_signature/microbe_tcg_proxy_mapping.csv) 后，额外输出：
  - `tcg_health_index`
  - `tcg_guild_1_fraction`
  - `tcg_guild_2_fraction`
  - `tcg_mapped_fraction`

当前这条线的定位是：

- 作为旧 `beneficial / risk genera` 启发式健康分旁边的 secondary readout；
- 在 `guild` 映射尚未补齐时，可以老老实实显示 `N/A` 或低覆盖率；
- 等目标论文的 genome/guild 成员表与 83 菌映射补齐后，再逐步升级成主健康分候选。

因此，当前 Step 3 的健康分其实已经进入“双增强状态”：

- 一条线是 `TCG-inspired health signature proxy`
- 另一条线是 `ENBI-like interaction balance proxy`

前者更偏“谁属于健康签名成员”，后者更偏“成员之间现在形成了什么生态关系”。

### 6.7 开发评分如何计算

当前的开发评分不是临床药效分，而是项目内部排序用的启发式综合分。

这部分近期已经从“单一线性扣分”升级成了“收益分 - 风险分，再做 sigmoid 拉伸”的双层结构。

当前核心分量包括：

- 收益项
  - `efficacy_proxy`：母药保留比例越高越好
  - `community_preservation_score`：群落健康和稳定性保留得越好越好
- 风险项
  - `dysbiosis_penalty`：对肠道健康破坏越大惩罚越高
  - `interaction_dysbiosis_penalty`：如果交互网络明显向 `cross-feeding-dominated dysbiosis-like` 偏移，则额外惩罚
  - `uncertainty_penalty`：落在 applicability 之外越多惩罚越高
  - `metabolite_burden_penalty`：聚合代谢负担越高惩罚越高

代码中的近似形式为：

```text
benefit_subscore
    = 0.55 * efficacy_proxy
    + 0.45 * community_preservation_score

risk_subscore
    = 0.38 * dysbiosis_penalty
    + 0.22 * interaction_dysbiosis_penalty
    + 0.22 * uncertainty_penalty
    + 0.18 * metabolite_burden_penalty

development_score_balance
    = benefit_subscore - risk_subscore

development_score
    = sigmoid(development_score_balance)
```

同时，旧版线性 `legacy score` 仍然保留，用于前后对照和网页解释。

这里要特别注意一层解释边界：

- `interaction_dysbiosis_penalty` 适合做系统内部排序和场景比较
- 但它不能被直接解释成真实临床 dysbiosis 风险
- 它反映的是“在当前 Step 3 代理交互网络里，系统是否更偏向正互作放大和竞争减弱”

这一步的意义不是直接判断药物有没有临床价值，而是为候选药物筛选提供一个统一的、可解释的内部排序指标。

因此它更适合回答：

- “在当前系统里，哪个候选药更值得优先往下看？”

而不适合直接回答：

- “这个药临床上一定更好”
- “这个药一定更容易开发成功”

### 6.8 Step 3 当前产物与演示结果

Step 3 当前脚本包括：

- [run_step3_simulation.py](../scripts/run_step3_simulation.py)
- [screen_step3_candidates.py](../scripts/screen_step3_candidates.py)

目前它已经能稳定输出：

- 单药单场景或多场景的时间轨迹
- `summary.json / trajectory_metrics.csv / top_microbe_changes.csv`
- 场景对比表
- 候选药物排序表
- 真实 cohort 单样本初始化后的单场景结果

其中最近新增、最值得看的 Step 3 字段包括：

- `health_index_legacy`
- `positive_interaction_strength`
- `negative_interaction_strength`
- `interaction_balance_rho`
- `interaction_balance_shift`
- `interaction_dysbiosis_penalty`

网页和服务层已经额外支持展示：

- 新版 `development score`
- 旧版 `legacy score`
- `benefit / risk` 拆解
- `TCG` 签名分和覆盖率

因此 Step 3 当前已经可以承担两类任务：

1. 同一药物在不同菌群背景下的相对比较；
2. 同一场景内多个候选药物的相对排序。

但需要注意，当前最适合解读的是“相对比较”和“分量拆解”，而不是把某个绝对数值当成最终生物学定论。

## 7. 当前推荐的默认系统配方

如果今天要把这个项目作为一个可运行基线提交，当前推荐配方如下。

### Step 1

- 主 gold：`Maier 2018`
- 主 silver：`MDIPID DEIM`
- `MASI curated`：低权重补充或外部验证
- 分类：`Chemprop`
- 回归：`RDKit + ExtraTrees`
- 输出：`hybrid`

### Step 2

- 主 gold：`Zimmermann 2019`
- 模型：`ExtraTrees classifier + regressor`
- 默认部署：`scaffold split + full-fit`
- 输出：全量 pair 的代谢概率、代谢强度和 applicability
- 增强层：`reaction class / product / gene evidence` 半机制投影

### Step 3

- 输入：Step 2 集成预测表
- 模拟器：离散时间群落模拟器
- 默认初始化：`panel-proxy`
- 可选初始化：真实 `community_table.csv`
- 输出：时间轨迹、健康指数、TCG proxy 签名分、开发评分、候选药物排序

## 8. 当前系统的解释边界

为了保证方法学表述准确，当前系统有 7 个边界需要明确写进任何对外说明中。

1. Step 1 的 `promote` 标签是派生标签，不是 Maier 原文天然三分类。
2. Step 2 当前已经能较稳定地预测“是否代谢”和“代谢强度”，但还不能稳健地做“新药精确产物预测”。
3. Step 2 的 `product_ids` 目前更适合做候选线索和来源内注释。
4. Step 2 的 `reaction / product / gene evidence` 当前是半机制投影，不是精确产物生成模型。
5. AMR 修正当前主要是服务层知识纠偏，用来压住高风险假阳性，不等同于菌株级 AST 真值。
6. Step 3 默认仍是 `panel-proxy simulation`，虽然已经支持真实 cohort 初始化，但还没有校准到纵向真实结局。
7. `development_score` 和 `tcg_health_index` 都是内部比较指标，不应解释成临床药效或注册开发结论。

## 9. 当前代码与数据资产分层

从仓库组织上看，当前系统已经形成了比较清晰的四层结构：

### 原始数据层

- [data/raw/step1](../data/raw/step1)
- [data/raw/step2](../data/raw/step2)

### 标准化与建模表层

- [data/processed/step1](../data/processed/step1)
- [data/processed/step2](../data/processed/step2)
- [data/processed/amr](../data/processed/amr)
- [data/processed/health_signature](../data/processed/health_signature)
- [data/processed/step3/cohorts](../data/processed/step3/cohorts)

### 模型资产层

- [models/step1](../models/step1)
- [models/step2](../models/step2)

### 预测与模拟输出层

- [predictions/step1](../predictions/step1)
- [predictions/step2](../predictions/step2)
- [predictions/step3](../predictions/step3)

这种分层的好处是：后续接入新数据源时，只需要在 `raw -> processed` 之间增加新的 source-specific 标准化模块，而不必改动整个系统。

理想情况下，后面如果接入一篇新论文，应该只需要：

1. 新增一个 `normalize_xxx.py`；
2. 让它输出仍然符合统一 schema 的表；
3. 下游训练和推理脚本几乎不用改。

## 10. 推荐复现实验顺序

如果要按当前仓库已经跑通的方式复现，推荐按照下面顺序执行。

### 10.1 Step 1 数据下载、标准化与 baseline

在已安装 [requirements-step1.txt](../requirements-step1.txt) 的环境中：

```bash
python scripts/download_step1_data.py
python scripts/normalize_step1_data.py
python scripts/train_step1_baseline.py --split-mode scaffold
```

如果需要补银标：

```bash
python scripts/normalize_step1_mdipid.py
python scripts/normalize_step1_masi.py
```

### 10.2 Step 1 Chemprop 与 hybrid

在已安装 [requirements-step1-chemprop.txt](../requirements-step1-chemprop.txt) 的环境中：

```bash
/tmp/microbe_env/bin/python scripts/prepare_step1_chemprop.py
/tmp/microbe_env/bin/python scripts/train_step1_chemprop.py \
  --dataset-csv data/processed/step1/chemprop_scaffold/classification/dataset.csv \
  --descriptors-path data/processed/step1/chemprop_scaffold/classification/descriptors.npz \
  --output-dir models/step1/chemprop_scaffold_classification_v1 \
  --task-type classification

/tmp/microbe_env/bin/python scripts/train_step1_chemprop.py \
  --dataset-csv data/processed/step1/chemprop_scaffold/regression/dataset.csv \
  --descriptors-path data/processed/step1/chemprop_scaffold/regression/descriptors.npz \
  --output-dir models/step1/chemprop_scaffold_regression_v1 \
  --task-type regression
```

然后生成 hybrid 输出：

```bash
/tmp/microbe_env/bin/python scripts/predict_step1_hybrid.py
```

### 10.3 Step 2 标准化、训练与全量预测

先标准化 `Zimmermann 2019`：

```bash
/tmp/microbe_env/bin/python scripts/normalize_step2_zimmermann.py \
  --input-path data/raw/step2/zimmermann_2019/NIHMS1530152-supplement-Supplementary_Tables_1-21.xlsx \
  --output-dir data/processed/step2/zimmermann_2019
```

再组装 Step 2 输入表：

```bash
/tmp/microbe_env/bin/python scripts/assemble_step2_inputs.py \
  --output-dir data/processed/step2 \
  --step2-label-table data/processed/step2/zimmermann_2019/zimmermann_2019_label_table.csv
```

训练 baseline：

```bash
/tmp/microbe_env/bin/python scripts/train_step2_baseline.py \
  --modeling-table data/processed/step2/zimmermann_2019/zimmermann_2019_modeling_table.csv \
  --split-mode scaffold \
  --output-dir models/step2/zimmermann_scaffold_split
```

给全量候选 pair 打分：

```bash
/tmp/microbe_env/bin/python scripts/predict_step2_baseline.py \
  --input-table data/processed/step2/step2_candidate_pairs_full.csv \
  --output-dir predictions/step2/baseline_scaffold_v1_83 \
  --classifier-path models/step2/zimmermann_scaffold_split/classifier_full.joblib \
  --regressor-path models/step2/zimmermann_scaffold_split/regressor_full.joblib \
  --metrics-path models/step2/zimmermann_scaffold_split/metrics.json \
  --applicability-reference-path models/step2/zimmermann_scaffold_split/applicability_reference.joblib
```

如果要重建网页默认使用的 83 菌库内预测表，当前推荐直接运行：

```bash
/tmp/microbe_env/bin/python scripts/build_library_83_panel_predictions.py
```

如果要同时生成 Step 2 半机制参考库：

```bash
/tmp/microbe_env/bin/python scripts/build_step2_mechanism_reference.py
```

### 10.4 Step 3 推演与候选筛选

单药多场景推演：

```bash
/tmp/microbe_env/bin/python scripts/run_step3_simulation.py \
  --drug-query "Metformin hydrochloride" \
  --all-scenarios \
  --output-dir predictions/step3/metformin_hydrochloride
```

多药比较：

```bash
/tmp/microbe_env/bin/python scripts/screen_step3_candidates.py \
  --output-dir predictions/step3/candidate_screen_demo \
  --scenario healthy_reference \
  --drug-query "Metformin hydrochloride" \
  --drug-query Digoxin \
  --drug-query Sulfasalazine \
  --drug-query "Diltiazem hydrochloride"
```

如果要在真实群落上运行：

```bash
/tmp/microbe_env/bin/python scripts/prepare_step3_cohort_community.py \
  --input-table data/raw/step3/cohorts/example/abundance_table.csv \
  --sample-id sample_01 \
  --output-path data/processed/step3/cohorts/example/sample_01_community.csv

/tmp/microbe_env/bin/python scripts/run_step3_simulation.py \
  --drug-query Digoxin \
  --community-table data/processed/step3/cohorts/example/sample_01_community.csv
```

如果要初始化 AMR 规则与 health signature 模板：

```bash
/tmp/microbe_env/bin/python scripts/prepare_amr_reference_tables.py --overwrite
/tmp/microbe_env/bin/python scripts/seed_beta_lactam_amr_rules.py
/tmp/microbe_env/bin/python scripts/seed_high_risk_antibiotic_amr_rules.py
/tmp/microbe_env/bin/python scripts/seed_aminoglycoside_anaerobe_rules.py
/tmp/microbe_env/bin/python scripts/seed_polymyxin_and_low_anaerobe_fq_rules.py

/tmp/microbe_env/bin/python scripts/prepare_health_signature_reference_tables.py --overwrite
/tmp/microbe_env/bin/python scripts/prepare_health_signature_target_tables.py --overwrite
```

这套顺序对应的是当前仓库已经跑通的主链路：`Step 1 -> Step 2 -> Step 3`。

### 10.5 机制层融合与 case-based benchmark（新增）

固定药物 case 的融合对比与 sanity benchmark 脚本：

```bash
/tmp/microbe_env/bin/python scripts/run_fusion_comparison.py
/tmp/microbe_env/bin/python scripts/evaluate_case_based_sanity_benchmark.py
/tmp/microbe_env/bin/python scripts/build_revised_case_outputs.py
```

输出目录：`predictions/evaluation/fusion_comparison/`

关键文件包括：

- `rifaximin.csv`
- `vancomycin.csv`
- `lubiprostone.csv`
- `metronidazole.csv`
- `sanity_summary.csv`
- `ranking_benchmark_summary.csv`
- `ecology_benchmark_summary.csv`
- `revised_case_based_results.csv`
- `revised_case_based_summary.md`

其中 `ranking_benchmark_summary.csv` 与 `ecology_benchmark_summary.csv` 是拆分后的双视角评估：

- ranking benchmark：检查疾病排序行为是否向目标方向移动
- ecology benchmark：单独检查生态现实性（如产丁酸菌保护、广谱抑制风险、生态风险平衡）

## 11. 建议的后续完善路线

在当前基础上，最值得继续推进的是下面几件事。

### 11.1 Step 1

- 把 AMR 规则从当前高风险组合扩成更系统的 `drug class × species/genus` 先验
- 将 `microbe_amr_reference.csv` 中的机制字段并入训练特征，而不只是服务层后修正
- 加 `domain of applicability` 与置信度校准
- 增加 `leave-drug-class-out` 评估
- 接入更多微生物机制特征，例如 `AGORA2` 或菌株代谢能力 embedding

### 11.2 Step 2

- 接入 `Javdan 2020` 做 community-level 外部验证
- 接入 `MagMD`, `gutMGene`, `AGORA2`, `Rhea` 做 `reaction class` 和 `product` 增强
- 将 Step 2 从“发生与否 + 强度 + 半机制投影”继续升级为“发生与否 + 反应类 + 候选产物排序 + 机制证据校准”

### 11.3 Step 3

- 用 `curatedMetagenomicData / GMrepo / MGnify` 扩大真实 cohort 初始化样本库
- 用更完整的 `health signature / TCG` 映射替换当前手工 `beneficial / risk genera`
- 用真实 `GMWI2` 或更严格 genome-level signature 替换当前 `GMWI2-like heuristic`
- 加入 community-level 外部验证数据，校准模拟器参数
- 加入显式微生物种间互作项，而不只依赖共享药物池和归一化产生的间接竞争
- 增加批量药物筛选和组合给药场景

## 12. 一句话总结当前项目状态

截至当前版本，这个项目已经不是单纯的规划，而是形成了一个可以真实运行的三段式系统：

- Step 1 已有稳定的药物影响预测器、hybrid 输出和第一批 AMR 服务层修正；
- Step 2 已有真实 gold 数据支撑的代谢 baseline、83 菌全量候选打分和半机制投影层；
- Step 3 已能在多场景或真实 `community_table` 初始化下完成群落推演，并开始接入 TCG proxy 健康签名。

下一阶段的重点，不再是“从 0 到 1 能不能跑起来”，而是“如何把机制信息、真实 cohort 和更严格外部验证接进来，把系统从可运行原型推进到更可信的研究级框架”。
