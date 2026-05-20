# Custom SMILES 去重后全菌预测结果（更新版）

- 生成时间: 2026-04-17T14:34:09.373734
- 输入 SMILES: `CC1=C(C#CC2=CN=C3C=CC(N4CCC(O)CC4)=NN23)C=CC=C1C(=O)NC1=CC(C(F)(F)F)=CC(Cl)=C1`
- session_id: `4f396f45b0a2`
- 去重规则: `species_label first, then microbe_label, then nt_code`
- 原始菌行数: **83**
- 去重后菌行数: **69**

## Confidence

- confidence_score: **0.35**
- confidence_tier: **low**
- warning_flags: `['core-butyrate-suppression', 'ecology-risk', 'over-suppression']`
- confidence_explanation: 当前预测置信度低（0.35），主要风险来自：核心产丁酸菌出现强抑制（5/5）。另检测到 2 项风险信号。

## 关键变化

- 去重后表中 `enzyme_prior_support_rate > 0` 的菌数: **63 / 69**
- 结构语义关键词已由 SMILES 自动提取（本分子触发 `amide|lactam`）。
- taxonomy 缺失已在服务侧补齐，显著减少“同值分组”现象。

完整 CSV: `/mnt/e/毕业/docs/test_results/custom_smiles_test_20260417_143343_dedup_all_microbes.csv`
Latest CSV: `/mnt/e/毕业/docs/test_results/custom_smiles_test_latest_dedup_all_microbes.csv`
JSON 摘要: `/mnt/e/毕业/docs/test_results/custom_smiles_test_20260417_143343_dedup_all_microbes.json`