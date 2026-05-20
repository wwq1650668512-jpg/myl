# Step 2 Enzyme Prior Layer

## 目标

给当前 Step 2 再补一层“可人工整理、可逐步证据化”的机制先验：

1. `83菌 -> 候选酶能力`
2. `酶 -> 可能作用的键 / 底物 / 反应类型`
3. 把这层先验转成 Step 2 的软特征，再反向补给 Step 1 promote 解释

这样做的定位不是替代真实实验标签，而是把“gutMGene / AGORA2 / 文献知识 / 人工整理”的知识层，
变成现有管线可以直接消费的表和列。

## 当前表

运行：

```bash
gut-drug-microbiome step2 build-enzyme-priors
```

或者：

```bash
/tmp/microbe_env/bin/python scripts/build_step2_enzyme_reference.py
```

会生成 5 个参考文件：

- `data/reference/step2_enzyme_function_catalog.csv`
- `data/reference/step2_microbe_enzyme_evidence_ledger.csv`
- `data/reference/step2_microbe_enzyme_prior_long.csv`
- `data/reference/step2_microbe_enzyme_prior_matrix.csv`
- `data/reference/step2_microbe_enzyme_curation_template.csv`

以及一个摘要：

- `data/reference/step2_enzyme_prior_summary.json`

## 表含义

### 1. `step2_microbe_enzyme_curation_template.csv`

这是给 species / strain 级文献整理准备的工作表。

特点：

- 是 `83菌 x 16酶类` 的全量网格
- 保留 `starter_*` 列，方便你看到当前 genus-level 初始建议
- 预留 `curated_*` 列，方便你逐条填写文献/基因组/人工校正结果

建议重点填写这些列：

- `curated_presence_call`
- `curated_evidence_scope`
- `curated_evidence_source`
- `curated_source_database`
- `curated_literature_citation`
- `curated_pmid`
- `curated_doi`
- `curated_genome_accession`
- `curated_strain_match_level`
- `curated_evidence_note`
- `curated_curation_status`

这张表本身也可以直接作为 `--literature-evidence-path` 输入，只要你把 `curated_*` 列填起来即可。

### 2. `step2_microbe_enzyme_evidence_ledger.csv`

这是“starter prior + 你补充的 species/strain 文献证据”的合并台账。

用途：

- 保留所有证据行，不丢历史痕迹
- 便于排查某个 `nt_code x enzyme` 是被哪条证据覆盖的
- 为后续做版本化和 reviewer 审核留底

### 3. `step2_microbe_enzyme_prior_long.csv`

一行代表一个 `nt_code x enzyme_id` 关系。

核心字段：

- `nt_code`
- `microbe_label`
- `species_label`
- `genus / family / phylum`
- `enzyme_id`
- `enzyme_name`
- `presence_call`
- `presence_weight`
- `evidence_scope`
- `evidence_source`
- `evidence_note`
- `curation_status`

现在这张表已经是 `resolved` 结果：

- 默认来自 genus-level starter prior
- 一旦你补了更强的 species / strain 证据，会自动覆盖 starter 行

当前 `presence_call` 可以来自 starter seed，也可以来自文献校正。建议后续按证据逐步改成：

- `curated_present`
- `likely_present`
- `genus_prior`
- `weak_prior`
- `absent`

### 4. `step2_microbe_enzyme_prior_matrix.csv`

是同一层信息的宽表版，方便快速看 83 菌每个菌有哪些酶能力先验。

### 5. `step2_enzyme_function_catalog.csv`

一行代表一个酶类及其机制解释。

核心字段：

- `enzyme_id`
- `enzyme_name`
- `enzyme_family`
- `ec_number`
- `reaction_class`
- `bond_target`
- `substrate_scope`
- `compound_semantic_families`
- `substrate_keywords`
- `likely_products_or_outcomes`
- `step2_mechanistic_role`
- `step1_feedback_role`
- `metabolism_weight`
- `step1_promote_weight`
- `step1_inhibit_weight`
- `notes`

## 当前如何接入管线

### Step 2

`predict_step2_baseline.py` 现在会在 baseline 预测和 mechanism projection 之后，再补一层 enzyme prior 注释。

新增的主要输出列：

- `predicted_enzyme_prior_flag`
- `predicted_enzyme_match_count`
- `predicted_enzyme_ids`
- `predicted_enzyme_names`
- `predicted_enzyme_reaction_classes`
- `predicted_enzyme_bond_targets`
- `predicted_enzyme_presence_score`
- `predicted_enzyme_support_score`
- `predicted_enzyme_step1_promote_support_score`
- `predicted_enzyme_step1_inhibit_risk_score`

它们的用途是：

- 给 Step 2 提供“这个 pair 从酶层面是否合理”的软支持
- 给网页和分析层补更可解释的机制提示
- 给 Step 1 的 promote 重评分提供下游反馈

### Step 1

`refine_step1_promote_with_step2(...)` 现在会把 enzyme prior support 一起纳入 support score。

这意味着：

- 如果某个菌本身带有更匹配当前化学家族/底物线索的酶能力
- 它在 Step 1 的 promote 倾向会得到更有机制解释的加权

同时：

- `predicted_enzyme_step1_inhibit_risk_score`
  会作为轻量风险项，避免把可能产生活性/毒性中间体的代谢简单当成 promote。

## 推荐后续人工整理顺序

1. 先校正 `83菌 -> genus / family / phylum` 的缺失项
2. 先人工确认高价值酶类：
   - `beta_glucuronidase`
   - `sulfatase`
   - `azoreductase`
   - `nitroreductase`
   - `carboxylesterase`
   - `amidase`
   - `beta_glucosidase`
   - `alpha_rhamnosidase`
   - `o_demethylase`
   - `dehydroxylase`
3. 对真实关心的药物家族补 `substrate_keywords / compound_semantic_families`
4. 先在 `step2_microbe_enzyme_curation_template.csv` 里填写 `curated_*` 列
5. 将这张表作为 `--literature-evidence-path` 重新跑构表脚本
6. 后续再把 gutMGene / AGORA2 / 文献证据持续写回 `evidence_source` / `pmid` / `doi` / `genome_accession`

## species / strain 证据覆盖规则

当同一个 `nt_code x enzyme_id` 同时存在多条证据时，当前会优先保留更强的证据：

1. `strain_literature`
2. `strain_genome`
3. `species_literature`
4. `species_genome`
5. `genus_literature`
6. `manual_curation`
7. `genus_prior_seed`

在同一层级内，还会参考：

- `curation_status`
- `strain_match_level`
- `presence_weight`
- 是否带 PMID / DOI / citation

这意味着：

- 你补的 species/strain 文献证据可以真正覆盖 starter prior
- 你补的 `absent` 证据也会覆盖旧的 `genus_prior`
- Step 2 / Step 1 后续消费到的是“解析后的最佳证据”，不是简单拼接

## 边界

- 当前这层是 `starter prior`，不是 strain-resolved genome truth。
- `presence_weight` 是软先验权重，不应等同于实验验证。
- 如果后续你补到了更准确的基因组或文献证据，优先更新长表，再重新跑 Step 2 预测。
