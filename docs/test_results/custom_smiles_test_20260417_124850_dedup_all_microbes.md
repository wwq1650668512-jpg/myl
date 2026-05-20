# Custom SMILES 去重后全菌预测结果

- 生成时间: 2026-04-17T12:49:17.704946
- 输入 SMILES: `CC1=C(C#CC2=CN=C3C=CC(N4CCC(O)CC4)=NN23)C=CC=C1C(=O)NC1=CC(C(F)(F)F)=CC(Cl)=C1`
- session_id: `0df3b3ca8269`

## 结果预览（Top 20）

| dedup_microbe            | n_rows | effect_label_majority | inhibit_prob_mean |    effect_score_mean | metabolism_label_majority | metabolized_prob_mean | enzyme_prior_support_rate | representative_nt_code |
| ------------------------ | -----: | --------------------- | ----------------: | -------------------: | ------------------------- | --------------------: | ------------------------: | ---------------------- |
| Roseburia intestinalis   |      1 | inhibit               |        0.94622815 |  -0.2793047261621862 | not_metabolized           |    0.2436558617142993 |                       0.0 | NT5011                 |
| Eubacterium rectale      |      1 | inhibit               |        0.93809116 |  -0.3686121686430975 | not_metabolized           |    0.2372026978207002 |                       0.0 | NT5009                 |
| Roseburia hominis        |      1 | inhibit               |          0.920772 |  -0.1553107813101785 | not_metabolized           |    0.2390697674886751 |                       0.0 | NT5079                 |
| Prevotella copri         |      1 | inhibit               |         0.9161046 |  -0.2219468195247023 | not_metabolized           |    0.2752028547363897 |                       0.0 | NT5019                 |
| Clostridium perfringens  |      2 | inhibit               |       0.899926665 | -0.16599459949092643 | not_metabolized           |   0.23393267309148524 |                       0.0 | NT5032                 |
| Bacteroides vulgatus     |      2 | inhibit               |        0.89814595 |   -0.137230803739515 | metabolized               |    0.2900953950117812 |                       0.0 | NT5001                 |
| Blautia obeum            |      1 | inhibit               |         0.8960232 |  -0.1331059165504944 | not_metabolized           |    0.2541932291545318 |                       0.0 | NT5069                 |
| Ruminococcus torques     |      1 | inhibit               |         0.8687094 |  -0.0539715217905881 | not_metabolized           |    0.2541932291545318 |                       0.0 | NT5047                 |
| Bifidobacterium animalis |      2 | inhibit               |        0.86793953 |  -0.0749623330608815 | not_metabolized           |    0.2336158658056326 |                       0.0 | NT5043                 |
| Collinsella aerofaciens  |      1 | inhibit               |         0.8535251 |  -0.1247909991248548 | not_metabolized           |    0.2171145859607738 |                       0.0 | NT5073                 |
| Alistipes shahii         |      1 | inhibit               |          0.839642 |  -0.0698092811450895 | not_metabolized           |    0.2336158658056326 |                       0.0 | NT5070                 |
| Bacteroides clarus       |      1 | inhibit               |          0.839642 |  -0.0698092811450895 | not_metabolized           |    0.2336158658056326 |                       0.0 | NT5052                 |
| Bacteroides coprocola    |      1 | inhibit               |          0.839642 |  -0.0698092811450895 | not_metabolized           |    0.2336158658056326 |                       0.0 | NT5053                 |
| Bacteroides dorei        |      1 | inhibit               |          0.839642 |  -0.0698092811450895 | not_metabolized           |    0.2400214698599232 |                       0.0 | NT5049                 |
| Blautia hansenii         |      1 | inhibit               |          0.839642 |  -0.0698092811450895 | not_metabolized           |    0.2470855708900345 |                       0.0 | NT5005                 |
| Butyrivibrio crossotus   |      1 | inhibit               |          0.839642 |  -0.0698092811450895 | not_metabolized           |    0.2336158658056326 |                       0.0 | NT5013                 |
| Coprococcus catus        |      1 | inhibit               |          0.839642 |  -0.0698092811450895 | not_metabolized           |    0.2336158658056326 |                       0.0 | NT5080                 |
| Desulfovibrio piger      |      1 | inhibit               |          0.839642 |  -0.0698092811450895 | not_metabolized           |    0.2336158658056326 |                       0.0 | NT5035                 |
| Dorea longicatena        |      1 | inhibit               |          0.839642 |  -0.0698092811450895 | not_metabolized           |    0.2336158658056326 |                       0.0 | NT5027                 |
| Eubacterium halii        |      1 | inhibit               |          0.839642 |  -0.0698092811450895 | not_metabolized           |    0.2336158658056326 |                       0.0 | NT5067                 |

## CSV 字段说明

| 字段                                                         | 含义                    | 备注                                                   |
| ------------------------------------------------------------ | ----------------------- | ------------------------------------------------------ |
| dedup_microbe                                                | 菌名键                  | 去重主键（优先 species_label，再用 microbe_label）     |
| n_rows                                                       | 该菌对应的原始行数      | >1 说明有重复条目被合并                                |
| nt_codes                                                     | 被合并的 nt_code 列表   | `;` 分隔                                             |
| species_labels / microbe_labels                              | 合并后的物种/菌名文本   | 当前这次 species_labels 为空                           |
| effect_label_majority                                        | Step1 主导效应标签      | 多数投票结果（inhibit/no_effect）                      |
| inhibit_prob_mean / max                                      | 抑制概率均值/最大值     | 越高越偏抑制                                           |
| promote_prob_refined_mean                                    | 修正后 promote 概率均值 | 结合 Step2 证据后的 promote 侧概率                     |
| effect_score_mean / min / max                                | 连续效应分数统计        | 负值偏抑制，正值偏促进                                 |
| metabolism_label_majority                                    | Step2 主导代谢标签      | 多数投票结果（metabolized/not_metabolized）            |
| metabolized_prob_mean / max                                  | 代谢概率均值/最大值     | 越高越偏 metabolized                                   |
| depletion_fraction_mean                                      | 母体消耗比例均值        | 当前输出多为负值，需结合模型定义解读                   |
| enzyme_prior_support_rate                                    | 酶先验覆盖率            | 0~1，表示该菌条目中酶先验命中占比                      |
| enzyme_support_score_mean                                    | 酶支持分均值            | 当前多数为空                                           |
| enzyme_names / enzyme_reaction_classes / enzyme_bond_targets | 酶机制信息汇总          | 当前这次为空                                           |
| reaction_classes                                             | 反应类型汇总            | 当前仅出现 `bioaccumulation_or_unresolved_depletion` |
| representative_nt_code                                       | 代表行 nt_code          | 由\|effect_score\| 最大的行选出                        |

## 表格分析

### 1) 去重效果与数据结构

- 原始 83 行，去重后 69 行，减少 14 行重复。
- `n_rows` 分布：62 个菌为 1 行，4 个菌为 2 行，2 个菌为 3 行，1 个菌为 7 行。
- 说明重复主要集中在少数菌，绝大多数菌条目是唯一的。

### 2) Step1 效应整体偏“抑制”

- `effect_label_majority`：`inhibit=59`，`no_effect=10`，`promote=0`。
- `inhibit_prob_mean` 范围：`0.0041 ~ 0.9462`，整体均值 `0.7227`。
- `effect_score_mean` 范围：`-0.3686 ~ -0.0076`，整体均值 `-0.0874`（整体偏负，和抑制方向一致）。

抑制最强 Top5（按 `inhibit_prob_mean`）：

- Roseburia intestinalis (`0.9462`)
- Eubacterium rectale (`0.9381`)
- Roseburia hominis (`0.9208`)
- Prevotella copri (`0.9161`)
- Clostridium perfringens (`0.8999`)

### 3) Step2 代谢判断以 not_metabolized 为主

- `metabolism_label_majority`：`not_metabolized=61`，`metabolized=8`。
- `metabolized_prob_mean` 范围：`0.2171 ~ 0.3269`，整体均值 `0.2472`，整体偏低。
- 解读：该分子在多数菌中更偏向“不被明显代谢”，仅少部分菌显示 metabolized 倾向。
