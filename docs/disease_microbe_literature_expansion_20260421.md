# 疾病-菌群关系扩充建议（2026-04-21）

## 现状（项目内）
- 当前 `data/reference/disease_microbe_dictionary.csv` 共 `235` 条关系，覆盖 `16` 个疾病。
- 字段以方向性为核心（`desired_step1_effect`），但来源可追溯信息不足（缺少统一 `source_database/PMID/source_url`）。
- 疾病名仍有跨来源同义问题（中英文/缩写/MeSH 名称并存），会导致“同病不同名”的重复。

## 可立即接入的高价值来源

### 1) GMrepo v3（优先级最高）
- 论文：NAR 2026，GMrepo v3（`320,208` runs、`1,107` phenotypes；含 marker taxa 资源）。
- 价值点：
  - 有“健康 vs 疾病”的跨项目 marker 对比；
  - 可得到方向（疾病富集 or 健康富集），可直接映射为 `inhibit/promote`；
  - 可使用跨项目一致性（`nrproj`、conflict）作为置信度。
- 项目内已新增脚本：
  - `scripts/build_gmrepo_disease_microbe_supplement.py`
  - 输出到 `data/reference/disease_microbe_gmrepo_supplement.csv`
  - 默认使用更稳健阈值：`min_abs_lda=2.0`、`min_nrproj=2`。

### 2) gutMDisorder v2.0（第二优先级）
- 论文：NAR 2023（数据库更新）。
- 价值点：
  - 包含疾病-菌群方向关系；
  - 同时有饮食/药物干预与菌群变化的关联，可补 `mechanism_note` 与干预背景。
- 建议：作为 GMrepo 的“外部验证层”，优先补充 GMrepo 未覆盖疾病。

### 3) gutMGene 2.0（机制增强）
- 论文：NAR 2025（`3,323` associations，覆盖肠道菌/宿主基因/代谢物关系）。
- 价值点：
  - 不仅是“疾病-菌”，还能补“菌-基因/代谢物”机制证据；
  - 适合增强你 Step3 的 `mechanism_note`、疾病目标收益解释项。
- 建议：不直接替代疾病 marker，而是增强关系解释与机制可信度。

## 已落地的自动化脚本

### 脚本
- `scripts/build_gmrepo_disease_microbe_supplement.py`

### 功能
- 调用 GMrepo API：
  - `get_all_phenotypes`
  - `get_all_phenotype_comparisons`
  - `getPhenotypeComparisonsDetails`
- 自动筛选“健康 vs 疾病”对比；
- 将 LDA 方向映射到：
  - `disease_effect_on_microbe`（increase/decrease）
  - `desired_step1_effect`（inhibit/promote）
- 依据 `nrproj` + `|LDA|` + conflict 生成 `relation_confidence`；
- 去重并输出可直接并入现有词典的结构化表。

### 运行示例
```bash
PYTHONPATH=src /tmp/microbe_env/bin/python scripts/build_gmrepo_disease_microbe_supplement.py \
  --output-path data/reference/disease_microbe_gmrepo_supplement.csv \
  --summary-path data/reference/disease_microbe_gmrepo_supplement.summary.json \
  --min-abs-lda 2.0 \
  --min-nrproj 2 \
  --verbose
```

### 与主词典合并（去重）
```bash
PYTHONPATH=src /tmp/microbe_env/bin/python scripts/merge_disease_microbe_references.py \
  --primary-path data/reference/disease_microbe_dictionary.csv \
  --supplement-path data/reference/disease_microbe_gmrepo_supplement.csv \
  --output-path data/reference/disease_microbe_dictionary_merged.csv \
  --summary-path data/reference/disease_microbe_dictionary_merged.summary.json
```

## 与你当前诉求的对应关系
- “预测未见过 SMILES”：关系库扩充后，疾病候选评分不再高度依赖少量手工关系，能降低 OOD 情况下的候选抖动。
- “疾病信息重复”：脚本内已对常见疾病术语（CRC/CD/UC/IBD/IBS/便秘/腹泻）做了规范映射，减少同病不同名重复。
- “疾病目标收益项 + 无药对照归一化”：扩充后的关系会直接增强 `disease_target_profile` 的覆盖面和稳定性。

## 主要来源
- GMrepo v3（NAR 2026）  
  https://academic.oup.com/nar/article/54/D1/D734/8340991
- GMrepo 文档：Cross-project comparisons of disease markers  
  https://evolgeniusteam.github.io/gmrepodocumentation/usage/crossprojectcomprisons/
- gutMDisorder v2.0（NAR 2023）  
  https://academic.oup.com/nar/article/51/D1/D717/6754909
- gutMGene 2.0（NAR 2025）  
  https://academic.oup.com/nar/article/53/D1/D783/7850954
