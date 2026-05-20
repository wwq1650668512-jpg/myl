# Reference Master List

这份文档汇总当前项目里我已经实际参考过、并且对最近 `promote / cross-feeding / Step 2 机制解释 / Step 3 interaction-aware simulation` 改动产生影响的文献与数据库。

阅读方式：

- `链接`：优先给 PubMed / PMC / 官方站点
- `参考内容`：我实际借用了什么信息
- `落地位置`：这些信息最终进入了哪个文档、种子表或代码逻辑

## 一. Step 1 主金标准与 Step 2 主金标准

### 1. Maier et al., 2018

- 标题：Extensive impact of non-antibiotic drugs on human gut bacteria
- 链接：
  - Nature: https://www.nature.com/articles/nature25979
- 参考内容：
  - Step 1 主金标准来源
  - `1197 drugs x 40 strains` 的基础实验框架
  - `auc` / `q` 统计量对应的主标签构建思路
- 落地位置：
  - `src/gut_drug_microbiome/step1/normalize.py`
  - `configs/labeling_rules.yaml`
  - `docs/step1_pipeline.md`

### 2. Zimmermann et al., 2019

- 标题：Mapping human microbiome drug metabolism by gut bacteria and their genes
- 链接：
  - Nature: https://www.nature.com/articles/s41586-019-1291-3
- 参考内容：
  - Step 2 主金标准来源
  - drug-microbe metabolism / parent depletion / mechanism projection 的主数据框架
- 落地位置：
  - `src/gut_drug_microbiome/step2/zimmermann_2019.py`
  - `src/gut_drug_microbiome/step2/train_baseline.py`
  - `docs/step2_pipeline.md`

## 二. 已实际用于 promote 种子表的文献

这些文献已经进入：

- `data/reference/promote_literature_seed_template.csv`
- `data/reference/promote_literature_seed_table.csv`
- `data/processed/step1/step1_silver_promote_literature.csv`

### 3. PMID: 39436683

- 标题：Akkermansia muciniphila Growth Promoted by Lychee Major Flavonoid through Bacteroides uniformis Metabolism
- 链接：
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/39436683/
- 参考内容：
  - `QRR -> Akkermansia muciniphila`
  - 关键不是“直接促进”，而是 `Bacteroides uniformis` 代谢后支持 `Akkermansia`
  - 这是当前 `metabolism_supported_promote` / `cross-feeding` 最关键的一条高质量证据
- 落地位置：
  - `data/reference/promote_literature_seed_template.csv`
  - `data/reference/cross_feeding_seed_template.csv`
  - `data/reference/cross_feeding_edges.csv`

### 4. PMID: 25845659

- 标题：Dietary Polyphenols Promote Growth of the Gut Bacterium Akkermansia muciniphila and Attenuate High-Fat Diet-Induced Metabolic Syndrome
- 链接：
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/25845659/
- 参考内容：
  - `dietary polyphenols -> Akkermansia muciniphila`
  - 作为 `direct_promote` / 宿主相关 promote 线索
- 落地位置：
  - `data/reference/promote_literature_seed_template.csv`
  - `docs/promote_data_integration_notes.md`

### 5. PMID: 29571008

- 标题：Grape proanthocyanidin-induced intestinal bloom of Akkermansia muciniphila is dependent on its baseline abundance and precedes activation of host genes related to metabolic health
- 链接：
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/29571008/
- 参考内容：
  - `grape proanthocyanidins -> Akkermansia muciniphila`
  - 作为 `direct_promote` / animal-level promote 证据
- 落地位置：
  - `data/reference/promote_literature_seed_template.csv`
  - `docs/promote_data_integration_notes.md`

### 6. PMID: 32598202

- 标题：Green Tea Encourages Growth of Akkermansia muciniphila
- 链接：
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/32598202/
- 参考内容：
  - `green tea extract / EGCG related cues -> Akkermansia`
  - 作为茶多酚方向的 promote 文献入口
- 落地位置：
  - `data/reference/promote_literature_seed_template.csv`
  - `docs/promote_data_integration_notes.md`

### 7. PMID: 34855864

- 标题：In vitro co-metabolism of epigallocatechin-3-gallate (EGCG) by the mucin-degrading bacterium Akkermansia muciniphila
- 链接：
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/34855864/
- 参考内容：
  - `EGCG -> Akkermansia`
  - 支持“自身代谢/代谢相关 promote”这条机制线
  - 也用于补 `EGCG` 的结构化种子记录
- 落地位置：
  - `data/reference/promote_literature_seed_template.csv`
  - `data/reference/promote_literature_seed_table.csv`

### 8. PMID: 28400010

- 标题：Exploitation of grape marc as functional substrate for lactic acid bacteria and bifidobacteria growth and enhanced antioxidant activity
- 链接：
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/28400010/
- 参考内容：
  - `grape marc -> Bifidobacterium / Lactobacillus`
  - 作为益生菌促进的功能型证据
- 落地位置：
  - `data/reference/promote_literature_seed_template.csv`

### 9. PMID: 36192946

- 标题：Nutraceutical formulations combining Limosilactobacillus fermentum quercetin and or resveratrol with beneficial impacts on the abundance of intestinal bacterial populations metabolite production and antioxidant capacity during colonic fermentation
- 链接：
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/36192946/
- 参考内容：
  - `quercetin/resveratrol + fermentation context`
  - 作为 `functional_promote` / colonic fermentation 证据
- 落地位置：
  - `data/reference/promote_literature_seed_template.csv`

### 10. PMID: 31608042

- 标题：Effects of Quercetin and Resveratrol on in vitro Properties Related to the Functionality of Potentially Probiotic Lactobacillus Strains
- 链接：
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/31608042/
- 参考内容：
  - `quercetin/resveratrol -> Lactobacillus`
  - 更偏功能增强，不完全等同于纯生长促进
- 落地位置：
  - `data/reference/promote_literature_seed_template.csv`

### 11. PMID: 23554103

- 标题：Flavonols enhanced production of anti-inflammatory substance(s) by Bifidobacterium adolescentis: prebiotic effects of phytochemicals?
- 链接：
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/23554103/
- 参考内容：
  - `quercetin -> Bifidobacterium adolescentis`
  - 归为 `functional_promote`
  - 用来补 Bifidobacterium 上的正向功能增强监督
- 落地位置：
  - `data/reference/promote_literature_seed_template.csv`
  - `docs/promote_data_integration_notes.md`

### 12. PMID: 25721815

- 标题：Effects of phytochemicals on in vitro anti-inflammatory activity of Bifidobacterium adolescentis
- 链接：
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/25721815/
- 参考内容：
  - `EGCG -> Bifidobacterium adolescentis`
  - 归为 `functional_promote`
  - 用来补 EGCG 在 Bifidobacterium 上的结构化正向监督
- 落地位置：
  - `data/reference/promote_literature_seed_template.csv`
  - `docs/promote_data_integration_notes.md`

### 13. PMID: 32751457

- 标题：Resveratrol Favors Adhesion and Biofilm Formation of Lacticaseibacillus paracasei subsp. paracasei Strain ATCC334
- 链接：
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/32751457/
- 参考内容：
  - `resveratrol -> Lactobacillus/Lacticaseibacillus paracasei`
  - 更偏定植/黏附/生物膜相关的 `functional_promote`
- 落地位置：
  - `data/reference/promote_literature_seed_template.csv`
  - `docs/promote_data_integration_notes.md`

### 14. PMID: 30013359

- 标题：Analysis of Temporal Changes in Growth and Gene Expression for Commensal Gut Microbes in Response to the Polyphenol Naringenin
- 链接：
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/30013359/
- 参考内容：
  - `naringenin -> Bifidobacterium catenulatum`
  - 当前作为一条较弱但明确的 promote 记录
- 落地位置：
  - `data/reference/promote_literature_seed_template.csv`
  - `data/reference/promote_literature_seed_table.csv`

### 15. PMID: 16701572

- 标题：Feruloyl oligosaccharides stimulate the growth of Bifidobacterium bifidum
- 链接：
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/16701572/
- 参考内容：
  - `feruloyl oligosaccharides -> Bifidobacterium bifidum`
  - 作为寡糖/功能底物促进的直接证据
- 落地位置：
  - `data/reference/promote_literature_seed_template.csv`

## 三. 用于 cross-feeding 参考边的文献

这些文献已经进入：

- `data/reference/cross_feeding_seed_template.csv`
- `data/reference/cross_feeding_edges.csv`

### 16. PMID: 36461198

- 标题：Co-culture fermentations suggest cross-feeding among Bacteroides ovatus DSMZ 1896 Lactiplantibacillus plantarum WCFS1 and Bifidobacterium adolescentis DSMZ 20083 for utilizing dietary galactomannans
- 链接：
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/36461198/
- 参考内容：
  - `Bacteroides ovatus -> Bifidobacterium adolescentis`
  - `Bacteroides ovatus -> Lactiplantibacillus plantarum`
  - 关键点是 `galactomannans / beta-manno-oligosaccharides` 的 producer-consumer 关系
- 落地位置：
  - `data/reference/cross_feeding_seed_template.csv`
  - `data/reference/cross_feeding_edges.csv`
  - `src/gut_drug_microbiome/step1/hybrid.py`

### 17. PMID: 16887514

- 标题：Arabinogalactan utilization in continuous cultures of Bifidobacterium longum effect of co-culture with Bacteroides thetaiotaomicron
- 链接：
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/16887514/
- 参考内容：
  - `Bacteroides thetaiotaomicron -> Bifidobacterium longum`
  - `arabinogalactan` 连续培养 / 共培养利用证据
- 落地位置：
  - `data/reference/cross_feeding_seed_template.csv`
  - `data/reference/cross_feeding_edges.csv`

## 四. 用于化合物标准化或候选整理的数据库

### 18. PhytoHub

- 链接：
  - 官方站点: https://phytohub.eu/
- 参考内容：
  - 化合物别名
  - phytochemical 分类
  - 作为后续补 `SMILES` / 规范名的参考来源
- 落地位置：
  - `docs/promote_data_integration_notes.md`
  - 当前主要作为策略参考，未直接自动接库

### 19. Phenol-Explorer

- 链接：
  - 数据库论文: https://pubmed.ncbi.nlm.nih.gov/20428313/
  - 代谢扩展版: https://pubmed.ncbi.nlm.nih.gov/22879444/
- 参考内容：
  - polyphenol / flavonoid 候选清单
  - 名称标准化、食物来源、代谢物信息
- 落地位置：
  - `docs/promote_data_integration_notes.md`

## 五. 用于银标与机制扩展判断的数据库或资源

### 20. MDIPID DEIM

- 官方下载入口：
  - 站点: https://mdipid.idrblab.net/
- 参考内容：
  - Step 1 `silver` 主来源
  - 借用了 `increase / decrease / promote / inhibit` 这类方向性描述
- 落地位置：
  - `src/gut_drug_microbiome/step1/weak_supervision.py`
  - `docs/step1_pipeline.md`
  - `docs/technical_route_detailed.md`

### 21. MASI v2.0

- 推断下载入口：
  - https://www.aiddlab.com/MASI2025/download/MASI_v2.0_download_microbeInfo.xlsx
  - https://www.aiddlab.com/MASI2025/download/MASI_v2.0_download_substanceInfo.xlsx
  - https://www.aiddlab.com/MASI2025/download/MASI_v2.0_download_microbeSubstanceInteractionRecords.xlsx
- 参考内容：
  - Step 1 辅助 `silver` 来源
  - 更适合低权重补充，不适合直接替代主金标准
- 落地位置：
  - `src/gut_drug_microbiome/step1/weak_supervision.py`
  - `docs/step1_pipeline.md`

### 22. GutCP

- 链接：
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/33637740/
- 参考内容：
  - 我参考它来判断“它更像 cross-feeding interaction 方法，而不是 compound -> strain growth effect 数据库”
  - 所以最终没有直接接入主训练
- 落地位置：
  - `docs/promote_data_integration_notes.md`

## 六. 已实际用于 Step 3 健康签名与交互平衡的文献

### 23. Corral Lopez et al., 2026

- 标题：Imbalance in gut microbial interactions as a marker of health and disease
- 链接：
  - Science: https://www.science.org/doi/10.1126/science.ady1729
  - bioRxiv: https://www.biorxiv.org/content/10.1101/2025.04.30.651474v1
- 参考内容：
  - 这篇文献对当前 Step 3 最重要的启发，不是“再加一组有益菌/风险菌名单”，而是把健康和失衡理解成两种不同的生态交互状态
  - 我实际借用了它的核心判断：
    - 健康态更接近 `competition-dominated`
    - dysbiosis-like 状态更接近 `cross-feeding / positive-interaction dominated`
  - 由此推动 Step 3 从只看组成的 `GMWI2-like heuristic`，升级成加入 `ENBI-like interaction balance proxy` 的 `interaction-aware` 版本
  - 当前仓库没有直接复现原文的 cohort-level interaction inference，而是基于：
    - `cross_feeding_edges.csv`
    - `compound_semantic_family`
    - Step 2 酶先验与代谢概率
    - 反应类/酶重叠近似竞争
    来构造一个可运行的 `positive vs negative interaction` 代理层
- 落地位置：
  - `src/gut_drug_microbiome/step3/simulation.py`
  - `docs/step3_pipeline.md`
  - `docs/technical_route_detailed.md`

### 24. A core microbiome signature as an indicator of health

- 标题：A core microbiome signature as an indicator of health
- 链接：
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/39378879/
- 参考内容：
  - 这篇文献继续作为 Step 3 `TCG-inspired health signature proxy` 的核心来源
  - 主要用于：
    - `83菌 -> guild / TCG proxy` 模板设计
    - `tcg_health_index / tcg_guild_1_fraction / tcg_guild_2_fraction` 这些 secondary readout 的定义
  - 在最近这轮 Step 3 升级里，它和上面的 `ady1729` 论文是互补关系：
    - 前者更偏“核心健康签名 / guild”
    - 后者更偏“交互网络平衡 / dysbiosis regime”
- 落地位置：
  - `data/processed/health_signature/microbe_tcg_proxy_mapping.csv`
  - `docs/health_signature_integration_plan.md`
  - `src/gut_drug_microbiome/step3/simulation.py`

## 七. 已查阅但暂不直接写成 promote seed 的混合证据

### 25. PMID: 26619254

- 标题：The impact of polyphenols on Bifidobacterium growth
- 链接：
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/26619254/
- 参考内容：
  - 这篇被用来判断 `quercetin / naringenin / rutin` 不能在摘要层面直接粗暴标成 promote
  - 其中存在明显的浓度依赖与菌株依赖双向效应
- 落地位置：
  - `docs/promote_data_integration_notes.md`
  - 当前作为“已查阅但暂不并入 promote seed”的保守依据

### 26. gutMGene v2.0

- 链接：
  - PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC11701569/
- 参考内容：
  - 用于 Step 2 / Step 3 后续机制增强路线判断
  - 不是这次 promote 直接入库的主来源，但被列为后续机制增强的重要资源
- 落地位置：
  - `docs/project_blueprint.md`
  - `docs/food_tcm_microbe_relation_plan.md`
  - `docs/step2_pipeline.md`

### 27. MagMD

- 链接：
  - PubMed: https://pubmed.ncbi.nlm.nih.gov/36467581/
- 参考内容：
  - 用于 Step 2 后续机制增强路线判断
  - 主要帮助考虑 drug metabolism mechanism / product knowledge 的补强方向
- 落地位置：
  - `docs/step2_pipeline.md`
  - `docs/technical_route_detailed.md`

## 八. 我实际从这些参考里拿了什么

最近这几轮改动里，最主要借用的是三类信息：

### 1. 方向性证据

- 某化合物是否让某菌 `increase / promote / bloom`
- 用来构造：
  - `promote_literature_seed_template.csv`
  - `step1_silver_promote_literature.csv`

### 2. 机制类型

- 是 `direct_promote`
- 还是 `metabolism_supported_promote`
- 还是 `cross_feeding_supported_promote`
- 用来构造：
  - `effect_type`
  - `promote_evidence_type`
  - `predicted_promote_evidence_type`

### 3. producer-consumer 边

- 哪个 producer 微生物先处理底物
- 哪个 consumer 微生物后获益
- 底物或家族关键词是什么
- 用来构造：
  - `cross_feeding_seed_template.csv`
  - `cross_feeding_edges.csv`
  - `predicted_cross_feeding_*`

### 4. interaction regime

- 群落当前更像是：
  - `competition-dominated`
  - 还是 `cross-feeding-dominated`
- 起始状态到末态是往哪个方向偏移
- 用来构造：
  - `positive_interaction_strength`
  - `negative_interaction_strength`
  - `interaction_balance_rho`
  - `interaction_balance_shift`
  - `interaction_dysbiosis_penalty`

## 九. 你现在优先去哪里看

如果你想快速查：

- 总体研究判断：
  - `docs/promote_data_integration_notes.md`
- 已转成可用 promote 记录的文献：
  - `data/reference/promote_literature_seed_template.csv`
- 已转成 cross-feeding 边的文献：
  - `data/reference/cross_feeding_seed_template.csv`
- Step 3 的交互平衡与健康签名逻辑：
  - `docs/technical_route_detailed.md`
  - `docs/step3_pipeline.md`
- 这份总表：
  - `docs/reference_master_list.md`
