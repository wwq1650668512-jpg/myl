# Step1 当前准确率评估（重点 promote）

- 预测文件: `/mnt/e/毕业/predictions/step1/hybrid_scaffold_v1_current_eval/predictions.csv`
- 预测列: `predicted_effect_label_hybrid`
- 真值列: `effect_label`

## 全量标注集

- n_rows: **43109**
- accuracy: **0.9218**
- balanced_accuracy: **0.6931**
- macro_f1: **0.5697**
- promote_precision: **0.0064**
- promote_recall: **0.2500**
- promote_f1: **0.0124**
- promote_support(true): **4**

混淆矩阵:

|                |   pred_inhibit |   pred_no_effect |   pred_promote |
|:---------------|---------------:|-----------------:|---------------:|
| true_inhibit   |           4631 |              484 |              0 |
| true_no_effect |           2729 |            35105 |            156 |
| true_promote   |              1 |                2 |              1 |

## Gold 子集

- n_rows: **43109**
- accuracy: **0.9218**
- balanced_accuracy: **0.6931**
- macro_f1: **0.5697**
- promote_precision: **0.0064**
- promote_recall: **0.2500**
- promote_f1: **0.0124**
- promote_support(true): **4**

混淆矩阵:

|                |   pred_inhibit |   pred_no_effect |   pred_promote |
|:---------------|---------------:|-----------------:|---------------:|
| true_inhibit   |           4631 |              484 |              0 |
| true_no_effect |           2729 |            35105 |            156 |
| true_promote   |              1 |                2 |              1 |

## Gold Promote（二分类）

- promote_binary_precision: **0.006369426751592357**
- promote_binary_recall: **0.25**
- promote_binary_f1: **0.012422360248447204**
- promote_binary_support: **4**
- promote_pred_positive: **157**

完整 JSON: `/mnt/e/毕业/docs/test_results/step1_current_accuracy_eval.json`