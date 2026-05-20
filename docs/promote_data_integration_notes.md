# Promote Data Integration Notes

## 目标

补充 Step 1 中稀缺的 `promote` 标签，优先接入：

- 有 `compound -> microbe` 明确方向的数据
- 能补 `SMILES` 和标准化化合物名的数据
- 能区分 `direct_promote` 和 `metabolism_supported_promote` 的数据

## 已核实、值得优先接入

### 1. 论文级 promote 数据

这些来源最适合直接做人审抽取或半自动抽取，因为它们给出了明确的化合物和菌响应关系。

- Akkermansia muciniphila Growth Promoted by Lychee Major Flavonoid through Bacteroides uniformis Metabolism
  - PMID: 39436683
  - 链接: https://pubmed.ncbi.nlm.nih.gov/39436683/
  - 价值:
    - 直接支持 `metabolism_supported_promote`
    - 说明不是所有菌一起获益，而是特定代谢菌驱动促进

- Dietary Polyphenols Promote Growth of the Gut Bacterium Akkermansia muciniphila and Attenuate High-Fat Diet-Induced Metabolic Syndrome
  - PMID: 25845659
  - 链接: https://pubmed.ncbi.nlm.nih.gov/25845659/
  - 价值:
    - `Akkermansia` 经典 promote 来源
    - 适合做第一批 promote 种子样本

- Grape proanthocyanidin-induced intestinal bloom of Akkermansia muciniphila...
  - PMID: 29571008
  - 链接: https://pubmed.ncbi.nlm.nih.gov/29571008/
  - 价值:
    - 可补充葡萄多酚 / 原花青素方向

- Green Tea Encourages Growth of Akkermansia muciniphila
  - PMID: 32598202
  - 链接: https://pubmed.ncbi.nlm.nih.gov/32598202/
  - 价值:
    - 可补茶多酚 / EGCG 相关 promote 证据

- In vitro co-metabolism of epigallocatechin-3-gallate (EGCG) by the mucin-degrading bacterium Akkermansia muciniphila
  - PMID: 34855864
  - 链接: https://pubmed.ncbi.nlm.nih.gov/34855864/
  - 价值:
    - 支持 `代谢相关促进` 这条机制线

- The impact of polyphenols on Bifidobacterium growth
  - PMID: 26619254
  - 链接: https://pubmed.ncbi.nlm.nih.gov/26619254/
  - 价值:
    - 直接做了 `B. adolescentis` / `B. bifidum` 的微孔板生长实验
    - 测试了 naringenin、hesperidin、rutin、quercetin 及多种 phenolic acids
    - 很适合提取 `compound + strain + concentration + direction`

- Flavonols enhanced production of anti-inflammatory substance(s) by Bifidobacterium adolescentis
  - PMID: 23554103
  - 链接: https://pubmed.ncbi.nlm.nih.gov/23554103/
  - 价值:
    - 不一定是纯生长促进，但说明 quercetin / galangin / fisetin 对 `B. adolescentis` 有正向调节
    - 适合做 `functional_promote` 子类

- Effects of phytochemicals on in vitro anti-inflammatory activity of Bifidobacterium adolescentis
  - PMID: 25721815
  - 链接: https://pubmed.ncbi.nlm.nih.gov/25721815/
  - 价值:
    - 可补 `B. adolescentis + phytochemicals` 的功能增强型记录

- Effects of Quercetin and Resveratrol on in vitro Properties Related to the Functionality of Potentially Probiotic Lactobacillus Strains
  - PMID: 31608042
  - 链接: https://pubmed.ncbi.nlm.nih.gov/31608042/
  - 价值:
    - 可补 `Lactobacillus + quercetin/resveratrol`
    - 更偏功能增强，不一定直接映射为增长促进

- Flavonols enhanced production of anti-inflammatory substance(s) by Bifidobacterium adolescentis: prebiotic effects of phytochemicals?
  - PMID: 23554103
  - 链接: https://pubmed.ncbi.nlm.nih.gov/23554103/
  - 价值:
    - 可补 `quercetin -> Bifidobacterium adolescentis`
    - 更适合归到 `functional_promote`

- Effects of phytochemicals on in vitro anti-inflammatory activity of Bifidobacterium adolescentis
  - PMID: 25721815
  - 链接: https://pubmed.ncbi.nlm.nih.gov/25721815/
  - 价值:
    - 可补 `EGCG -> Bifidobacterium adolescentis`
    - 适合做 `functional_promote`，不要硬当成纯生长促进

- Resveratrol Favors Adhesion and Biofilm Formation of Lacticaseibacillus paracasei subsp. paracasei Strain ATCC334
  - PMID: 32751457
  - 链接: https://pubmed.ncbi.nlm.nih.gov/32751457/
  - 价值:
    - 可补 `resveratrol -> Lactobacillus/Lacticaseibacillus paracasei`
    - 更像益生功能增强和定植促进

- Nutraceutical formulations combining Limosilactobacillus fermentum, quercetin, and/or resveratrol...
  - PMID: 36192946
  - 链接: https://pubmed.ncbi.nlm.nih.gov/36192946/
  - 价值:
    - 可补 colonic fermentation 条件下的益菌上升证据

### 2. 可用于化合物标准化的数据库

- PhytoHub
  - 链接: https://phytohub.eu/
  - 价值:
    - 补 `SMILES`
    - 补化合物别名
    - 补 phytochemical 分类
    - 补已知代谢物信息

- Phenol-Explorer
  - 内容库论文: https://pubmed.ncbi.nlm.nih.gov/20428313/
  - 代谢扩展版: https://pubmed.ncbi.nlm.nih.gov/22879444/
  - 价值:
    - 适合整理 polyphenol / flavonoid 候选物列表
    - 适合标准化名称、食物来源、代谢物
    - 不适合直接当 `compound -> strain -> effect_score` 训练表

## 暂不建议直接接入主训练

以下名称这轮没有核实到足够可信的官方数据库或主论文，不建议直接写进主流程：

- `GMBD / gmbd.bmicc.cn`
- `Flamingo Database`
- `MDP – Microbiome Drug Pharmacophore Database`
- `PMP DB`

说明:

- `GutCP` 已核实存在，但它是 `cross-feeding interaction` 预测方法，不是 `compound -> strain growth effect` 数据库
  - PMID: 33637740
  - 链接: https://pubmed.ncbi.nlm.nih.gov/33637740/

## 已查阅但暂不直接写成 promote 标签的混合证据

- The impact of polyphenols on Bifidobacterium growth
  - PMID: 26619254
  - 链接: https://pubmed.ncbi.nlm.nih.gov/26619254/
  - 原因:
    - 这篇很有价值，但摘要层面已经明确提示存在浓度依赖和菌株依赖的双向效应
    - 其中 quercetin 更偏持续抑制，naringenin/rutin 也不是稳定单方向促进
    - 所以更适合留作“后续人工精读后再拆浓度级记录”，而不是现在直接并成 promote silver

## 建议的数据模式

建议新增一张 `promote_literature_seed_table.csv`，字段至少包括：

- `record_id`
- `compound_name_raw`
- `compound_name_normalized`
- `smiles`
- `inchikey`
- `microbe_name_raw`
- `microbe_label_normalized`
- `strain_label`
- `effect_direction`
- `effect_type`
  - `direct_promote`
  - `metabolism_supported_promote`
  - `functional_promote`
- `effect_score_proxy`
- `effect_score_proxy_type`
  - `od_ratio`
  - `cfu_change`
  - `relative_abundance_change`
  - `qualitative_increase`
- `dose_value`
- `dose_unit`
- `culture_context`
  - `mono_culture`
  - `co_culture`
  - `fecal_fermentation`
  - `animal`
- `supporting_microbe`
- `source_pmid`
- `source_title`
- `evidence_level`

## 建议的接入策略

### 第一层：安全接入

- 先做人审或半自动抽取
- 只保留明确 `increase / promote / bloom / viable count increase` 的记录
- 暂时不强行映射成主 `effect_score`
- 先做 `promote_literature_silver`

### 第二层：分类型 promote

- `direct_promote`
- `metabolism_supported_promote`
- `functional_promote`

这样可以和当前 Step 2 修正层对齐。

### 第三层：弱定量

把文献里不同类型的正向证据先映射成 `effect_score_proxy`：

- 强促进: `0.20`
- 中促进: `0.12`
- 弱促进: `0.05`

这里只能先做 proxy，不能冒充和 `auc - 1.0` 完全同量纲。

## 模型改进建议

### 1. 不要把 promote 和 inhibit 用同一标签逻辑

当前已经改为：

- `inhibit`: `effect_score <= -0.20`, `q < 0.05`
- `promote`: `effect_score >= 0.12`, `q < 0.30`

### 2. 保持 promote 是单独头

推荐保留：

- 主分类器: `inhibit / no_effect`
- promote 辅助头: `promote vs not_promote`

### 3. Step 2 只修正“代谢该药物的菌”

当前 promote 重打分已经按这个方向收紧：

- 非代谢菌: 不 uplift
- 弱代谢菌: 轻微 uplift
- 强代谢 + 机制支持: 明显 uplift

### 4. 最终需要做 OOF Step 2 特征

如果后面把 Step 2 结果作为 Step 1 promote 头训练特征，最好构建：

- `out-of-fold Step 2 feature table`

避免训练信息泄漏。

## 推荐下一步

优先抽取以下化合物的 promote 文献数据：

- quercetin
- rutin
- EGCG
- resveratrol
- grape proanthocyanidins
- lychee flavonoids

优先关注以下菌：

- `Akkermansia muciniphila`
- `Bifidobacterium adolescentis`
- `Bifidobacterium bifidum`
- `Lactobacillus` / `Limosilactobacillus` strains

## 当前未纳入 promote 种子表的说明

- `rutin`
  - 已补结构来源，但当前核实到的 `PMID:26619254` 摘要只说明总体上既有促进也有抑制，未在摘要里明确指出 `rutin` 对目标双歧杆菌是稳定 promote
- `hesperidin`
  - 当前核实到的 `PMID:26619254` 摘要明确提示其抑制效应更稳定，因此暂不纳入 promote 种子表

后续如果拿到正文图表或补充材料里更细的菌株级增长曲线，可以再把 `rutin / hesperidin` 精细拆分后补入。
