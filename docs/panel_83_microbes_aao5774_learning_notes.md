# 83菌面板与 10.1126/science.aao5774 学习表

生成文件：`docs/panel_83_microbes_aao5774_learning_table.csv`

说明：
- 本表的 `english_name / chinese_name / 门纲目科属` 来自当前仓库 83 菌面板及其补充映射。
- `aao5774_role` 不是论文原文逐株给出的官方“好菌/坏菌名单”。
- 论文摘要明确强调的是两类方向：
  - 被高纤维促进的 `SCFA-producing strains`
  - 被抑制的、与 `indole / hydrogen sulfide` 等代谢不利产物相关的菌
- 因此本表把当前 83 菌按“论文启发的学习用途”分成：
  - `beneficial_like`
  - `detrimental_like`
  - `unclear_or_mixed`
- 如果你后面想做更严格版本，建议再去补论文正文 / supplement 中的 strain-level 清单，而不是把本表直接当最终金标准。

主要来源：
- 当前 83 菌面板：`data/processed/step1/step1_microbe_table.csv`
- 参考映射：`data/processed/health_signature/microbe_reference_genome_mapping.csv`
- 论文：Gut bacteria selectively promoted by dietary fibers alleviate type 2 diabetes
  DOI: 10.1126/science.aao5774
