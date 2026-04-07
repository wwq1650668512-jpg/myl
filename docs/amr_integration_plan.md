# AMR Integration Plan

## 目标

把“微生物-药物抗性修正”从零散规则，整理成一条可以持续扩展的数据层。

这版先落的是最小骨架：

- `microbe_amr_reference.csv`
- `drug_resistance_rules.csv`
- `species_drugclass_phenotype_prior.csv`

它们的作用分别是：

- 面板菌的 AMR 机制与参考基因组画像
- 物种 / 属级的药类抗性规则
- 表型层的抗性先验与频率校准

## 推荐目录结构

```text
data/
  processed/
    amr/
      microbe_amr_reference.csv
      drug_resistance_rules.csv
      species_drugclass_phenotype_prior.csv
      template_summary.json
```

对应初始化脚本：

- [prepare_amr_reference_tables.py](../scripts/prepare_amr_reference_tables.py)
- [seed_beta_lactam_amr_rules.py](../scripts/seed_beta_lactam_amr_rules.py)
- [seed_high_risk_antibiotic_amr_rules.py](../scripts/seed_high_risk_antibiotic_amr_rules.py)
- [seed_aminoglycoside_anaerobe_rules.py](../scripts/seed_aminoglycoside_anaerobe_rules.py)
- [seed_polymyxin_and_low_anaerobe_fq_rules.py](../scripts/seed_polymyxin_and_low_anaerobe_fq_rules.py)

参考来源与引用点清单：

- [amr_reference_notes.md](amr_reference_notes.md)

## 三张表分别存什么

### 1. `microbe_amr_reference.csv`

面向当前 83 菌面板，按 `nt_code` 一行一个菌。

建议优先填这些字段：

- `reference_genome_id`
- `reference_genome_source`
- `has_beta_lactamase`
- `beta_lactamase_family`
- `has_efflux_amr`
- `has_target_alteration`
- `has_intrinsic_amr_evidence`
- `expected_beta_lactam_resistant`
- `amr_gene_summary`
- `amr_gene_source`
- `curation_status`

这张表更偏“机制层”与“菌面板基础画像”。

### 2. `drug_resistance_rules.csv`

面向规则修正层，一行表示一个“药类-菌”的规则。

核心字段：

- `drug_class`
- `drug_name`
- `species_label`
- `genus`
- `expected_phenotype`
- `rule_level`
- `mechanism_hint`
- `rule_strength`
- `source_name`
- `source_url`

推荐把 `expected_phenotype` 限定在：

- `resistant`
- `susceptible`
- `intermediate`
- `unknown`

推荐把 `rule_level` 限定在：

- `species`
- `genus`
- `drug`
- `drug_class`

### 3. `species_drugclass_phenotype_prior.csv`

面向表型校准层，一行表示一个“物种 / 属 - 药类”的统计先验。

核心字段：

- `species_label`
- `genus`
- `drug_class`
- `drug_name`
- `n_tested`
- `resistant_fraction`
- `intermediate_fraction`
- `susceptible_fraction`
- `mic50`
- `mic90`
- `source_name`
- `source_url`

这张表后面可以用来做：

- raw prediction 概率的再校准
- “与已知耐药先验冲突”的提示
- 场景内候选药物的风险排序

## 推荐数据接入顺序

### 第一步：先做规则层

优先把这些资源整理到 `drug_resistance_rules.csv`：

- EUCAST expected phenotypes
- EUCAST expert rules

这一步最适合先修正明显不合理的结果，比如：

- `Bacteroides` 对 `penicillin / beta-lactam` 的已知耐药先验

### 第二步：再做机制层

优先把这些资源整理到 `microbe_amr_reference.csv`：

- CARD
- NCBI AMRFinderPlus / NDARO

建议先集中补这些机制：

- `beta_lactamase`
- `efflux`
- `target alteration`

### 第三步：最后做表型校准层

优先把这些资源整理到 `species_drugclass_phenotype_prior.csv`：

- NCBI AST Browser
- BV-BRC phenotype / AMR tables

这一步不是替代模型，而是给模型补一个“现实世界的回拉力”。

## 建议接到模型里的方式

### 方案 A：先做结果修正

在 Step 1 / 自定义药物预测输出后新增字段：

- `drug_class`
- `expected_phenotype`
- `amr_conflict_flag`
- `amr_rule_source`

如果出现：

- `expected_phenotype = resistant`
- 但模型输出 `effect_label = inhibit`

就：

- 下调抑制概率
- 或回退成 `no_effect`
- 同时在前端打出“与 AMR 先验冲突”的标记

### 方案 B：再做特征增强

把 `microbe_amr_reference.csv` 中的机制字段并入 Step 1 / Step 2 特征：

- `has_beta_lactamase`
- `has_efflux_amr`
- `has_target_alteration`
- `expected_beta_lactam_resistant`

这一步能让模型从“结果后修正”逐步升级到“训练时知晓 AMR 机制”。

### 方案 C：最后做概率校准

使用 `species_drugclass_phenotype_prior.csv` 的统计先验，对：

- `inhibit_probability`
- `effect_score`

做一次校准，生成：

- `corrected_inhibit_probability`
- `corrected_effect_label`

## 初始化命令

首次生成模板文件：

```bash
/tmp/microbe_env/bin/python scripts/prepare_amr_reference_tables.py
```

如果已经存在旧模板，想重建：

```bash
/tmp/microbe_env/bin/python scripts/prepare_amr_reference_tables.py --overwrite
```

写入第一批 `beta-lactam / penicillin` 修正规则：

```bash
/tmp/microbe_env/bin/python scripts/seed_beta_lactam_amr_rules.py
```

写入扩展版“高风险假阳性”抗生素规则：

```bash
/tmp/microbe_env/bin/python scripts/seed_high_risk_antibiotic_amr_rules.py
```

写入 `aminoglycoside × strict anaerobes` 规则：

```bash
/tmp/microbe_env/bin/python scripts/seed_aminoglycoside_anaerobe_rules.py
```

写入 `polymyxin × gram-positive` 与“低厌氧覆盖 fluoroquinolone × strict anaerobes”规则：

```bash
/tmp/microbe_env/bin/python scripts/seed_polymyxin_and_low_anaerobe_fq_rules.py
```

## 下一步最值得做的事

如果继续推进，我建议优先做这两件事：

1. 用 EUCAST 先生成第一版 `drug_resistance_rules.csv`
2. 给 83 菌面板补一批 `beta_lactamase` 与 `expected_beta_lactam_resistant`

这样可以最快修正当前像 `Penicillin G -> Bacteroides inhibit` 这类明显不合理的结果。

## 当前已落地的第一批规则

目前仓库已经先种下了一版偏保守的 `Bacteroides × beta-lactam / penicillin` 规则，定位是：

- `beta_lactam` 只做支持性下调
- `penicillin` 做更强的耐药先验修正
- 预测输出同时保留原始值与修正值，方便对照

在这基础上，当前还额外扩到了：

- `gram-negative × vancomycin / glycopeptide`
- `gram-negative × daptomycin / lipopeptide`
- `Lactobacillus × vancomycin / glycopeptide`
- `strict anaerobes × gentamicin / amikacin / tobramycin / streptomycin / aminoglycoside`
- `gram-positive × colistin / polymyxin_b / polymyxin`
- `strict anaerobes × ciprofloxacin / levofloxacin / ofloxacin / norfloxacin / fluoroquinolone_low_anaerobe`

这版主要是为了先压住最明显的假阳性，不代表已经完成了完整的 AMR 校准。
