# Core Microbiome Health Signature 接入方案

## 目标

把当前 Step 3 里基于手工 `beneficial / risk genera` 的启发式健康分，逐步升级成更接近
《A core microbiome signature as an indicator of health》思路的健康签名层。

这里要先明确一件事：

- 当前仓库的 83 菌面板主要是 `species / strain-hint / nt_code` 层级；
- 目标论文强调的是更偏 `genome-specific` 的核心菌群签名；
- 所以仓库里更现实的第一步不是“直接复现原版全文分数”，而是先建立：
  - `83菌 -> 参考基因组`
  - `83菌 -> TCG / guild proxy`

## 已新增的数据模板

这次先把数据骨架落到：

- [microbe_reference_genome_mapping.csv](../data/processed/health_signature/microbe_reference_genome_mapping.csv)
- [microbe_tcg_proxy_mapping.csv](../data/processed/health_signature/microbe_tcg_proxy_mapping.csv)
- [microbe_reference_sequence_sets.csv](../data/processed/health_signature/microbe_reference_sequence_sets.csv)
- [health_signature_target_genomes.csv](../data/processed/health_signature/health_signature_target_genomes.csv)
- [health_signature_source_registry.csv](../data/processed/health_signature/health_signature_source_registry.csv)
- [template_summary.json](../data/processed/health_signature/template_summary.json)
- [target_template_summary.json](../data/processed/health_signature/target_template_summary.json)

生成脚本是：

- [prepare_health_signature_reference_tables.py](../scripts/prepare_health_signature_reference_tables.py)
- [prepare_health_signature_target_tables.py](../scripts/prepare_health_signature_target_tables.py)

## 三张表分别做什么

### 1. microbe_reference_genome_mapping.csv

这张表回答：

- 每个 `nt_code` 对应哪个参考基因组；
- 当前是 `exact strain`、`same species proxy`，还是 `same genus proxy`；
- 后续要去哪里下载参考序列。

关键字段：

- `nt_code`
- `species_label`
- `canonical_species_name`
- `strain_hint`
- `culture_collection_hint`
- `ncbi_taxid`
- `ncbi_assembly_accession`
- `refseq_genome_accession`
- `gtdb_genome_id`
- `reference_genome_label`
- `mapping_level`
- `mapping_confidence`

建议的填写顺序：

1. 先补 `canonical_species_name`
2. 再补 `culture_collection_hint / DSM / ATCC`
3. 再定 `TaxID`
4. 再定 `Assembly accession`
5. 最后判断 `mapping_level`

### 2. microbe_tcg_proxy_mapping.csv

这张表回答：

- 当前 83 菌能否映射到目标论文里的核心健康签名成员；
- 是精确命中，还是物种级 / 属级代理；
- 是否已经可以进入 Step 3 健康分计算。

关键字段：

- `nt_code`
- `reference_genome_label`
- `tcg_membership`
- `tcg_mapping_level`
- `tcg_confidence`
- `tcg_target_genome_id`
- `tcg_target_label`
- `tcg_support_evidence`
- `ready_for_step3`

推荐枚举值：

- `tcg_membership`:
  - `guild_1`
  - `guild_2`
  - `unmapped`
- `tcg_mapping_level`:
  - `exact_genome`
  - `same_species_proxy`
  - `same_genus_proxy`
  - `unmapped`

### 3. microbe_reference_sequence_sets.csv

这张表回答：

- 真正要落到本地或外部对象存储里的序列集是什么；
- 每个 `nt_code` 后续用哪个 `FASTA / assembly / genome accession`。

这张表更偏“可运行资产索引”，方便后面接：

- NCBI Assembly
- RefSeq
- GTDB representative genomes
- 自己下载到本地的 fasta 文件

## 推荐工作流

### 第一步：把 83 菌先映射到参考基因组

最小目标不是一上来做到 83/83 全部 exact strain，而是先分层：

- `exact_strain`
- `same_species_proxy`
- `same_genus_proxy`

这样 Step 3 后续在解释健康签名时，就能区分：

- “这是强对应”
- “这是保守代理”

### 第二步：整理目标论文签名成员

把目标论文里的 genome / guild / cluster 成员整理成独立参考表。  
这一步现在已经有单独模板：

- [health_signature_target_genomes.csv](../data/processed/health_signature/health_signature_target_genomes.csv)
- [health_signature_source_registry.csv](../data/processed/health_signature/health_signature_source_registry.csv)

字段建议：

- `target_genome_id`
- `guild_membership`
- `species_name`
- `taxid`
- `assembly_accession`
- `source_name`
- `source_url`

建议这一步的填表顺序是：

1. 先把论文中出现的 guild / cluster member genome 全部抄入 `health_signature_target_genomes.csv`
2. 再补 `taxid / assembly accession / refseq accession`
3. 最后把 `ready_for_matching` 标成 `yes`

### 第三步：建立 83 菌 -> guild proxy

映射原则建议非常保守：

1. 先 `exact genome`
2. 再 `same species`
3. 最后才 `same genus`

不要反过来。否则 Step 3 健康分会被过度放大解释。

### 第四步：替换 Step 3 的健康分中间项

当前健康分实现还在：

- [simulation.py](../src/gut_drug_microbiome/step3/simulation.py#L12)
- [simulation.py](../src/gut_drug_microbiome/step3/simulation.py#L349)

也就是：

- 先定义 `BENEFICIAL_GENERA`
- 再定义 `RISK_GENERA`
- 再算 `beneficial_fraction / risk_fraction`

后面可以升级成：

- `guild_1_fraction`
- `guild_2_fraction`
- `guild_balance`

再和：

- `diversity`
- `stability`

一起组成新的 `health_signature_score`

## 当前边界

这一步现在已经进入“双轨状态”：

- 数据层已经有 `83菌 -> 参考基因组 / TCG proxy / 序列集` 的模板；
- Step 3 代码已经能读取 [microbe_tcg_proxy_mapping.csv](../data/processed/health_signature/microbe_tcg_proxy_mapping.csv)，并输出：
  - `tcg_health_index`
  - `tcg_guild_1_fraction`
  - `tcg_guild_2_fraction`
  - `tcg_mapped_fraction`
- 但当前默认模板还是空白映射，所以网页里这部分大多会显示为 `N/A` 或覆盖率接近 `0`。

还没有完成的部分仍然是：

- 自动联网下载参考基因组；
- 自动抓取目标论文的 genome member 列表；
- 把 `TCG-inspired` 指标升级成主健康分，而不是当前的 secondary readout。
所以当前状态应该理解为：

- 已经具备稳定落表位置；
- 已经具备稳定的 Step 3 读取接口；
- 已经能开始人工或半自动补 `TaxID / Assembly / guild`；
- 但还没有完成原论文级别的 full reproduction。

## 推荐下一步

如果继续往前推，优先顺序建议是：

1. 先补 83 菌的 `TaxID + Assembly accession`
2. 再整理目标论文的核心签名 genome 列表
3. 再做 `microbe_tcg_proxy_mapping.csv`
4. 最后把 Step 3 健康分替换成新的 `TCG-inspired` 指标
