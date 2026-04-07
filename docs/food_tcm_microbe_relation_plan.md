# 83菌与食品 / 中医药关系整理方案

## 目标

把当前 `83` 个 Step 1 / Step 3 面板微生物，与：

- 食品 / 膳食模式 / 膳食成分
- 中医药 / 中药单味药 / 方剂 / 相关活性成分

之间的关系，整理成一套可持续补充的结构化表，而不是一次性聊天记录。

这里默认把用户说的“中医院”按“中医药 / 中药”理解。  
如果后续实际要查的是医院科室或临床机构协作关系，需要另外建表。

## 当前最值得用的外部资源

### 1. 食品方向

#### FGMDI

- 来源：Food Bioscience, 2024
- DOI: https://doi.org/10.1016/j.fbio.2024.104091
- 参考点：
  - `1806` 条 food-gut microbe associations
  - `495` 个 gut microbes
  - `313` 个 foods
- 用途：
  - 最适合做“食物 / 膳食成分 -> 微生物增减”初筛
  - 也适合查饮食模式和某些植物来源成分

#### gutMDisorder

- 来源：Nucleic Acids Research, 2020
- 链接：https://pmc.ncbi.nlm.nih.gov/articles/PMC6943049/
- 参考点：
  - `2263` 条 human curated associations
  - `579` 个 gut microbes
  - `77` 个 intervention measures
- 用途：
  - 适合补“食品 / 干预措施 -> 微生物”的宿主背景、疾病背景和实验上下文
  - 也能补一些 FGMDI 没有覆盖到的 intervention 语境

### 2. 中医药方向

#### MicrobeTCM

- 来源：Pharmacological Research, 2024
- 网站：https://www.microbetcm.com
- 数据库介绍页：https://ngdc.cncb.ac.cn/databasecommons/database/id/9463
- 参考点：
  - `725` 个 microbes
  - `1032` 个 herbs
  - `1468` 个 herb-formulas
  - `15780` 个 chemical compositions
- 用途：
  - 最适合做“中药 / 方剂 / 成分 -> 微生物”的主检索库
  - 也是当前最应该优先查的中医药资源

### 3. 补上下游机制

#### gutMGene v2.0

- 来源：Nucleic Acids Research, 2025
- 链接：https://pmc.ncbi.nlm.nih.gov/articles/PMC11701569/
- 用途：
  - 当某个食物或中药不是直接改变菌丰度，而是通过代谢物或宿主基因轴起作用时，可以用来补机制语境
  - 更适合做“解释增强”，不是第一张关系表

## 已新增的本地模板

这次已经把模板落到：

- [microbe_food_relation_table.csv](../data/processed/food_tcm/microbe_food_relation_table.csv)
- [microbe_tcm_relation_table.csv](../data/processed/food_tcm/microbe_tcm_relation_table.csv)
- [source_registry.csv](../data/processed/food_tcm/source_registry.csv)
- [template_summary.json](../data/processed/food_tcm/template_summary.json)

生成脚本：

- [prepare_food_tcm_microbe_reference_tables.py](../scripts/prepare_food_tcm_microbe_reference_tables.py)

## 本次新增的第一版结果

这轮已经把可直接访问的官方资源进一步落成了两类结果文件：

- [microbetcm_first_pass_hits.csv](../data/processed/food_tcm/microbetcm_first_pass_hits.csv)
- [microbetcm_relation_first_pass.csv](../data/processed/food_tcm/microbetcm_relation_first_pass.csv)
- [microbetcm_first_pass_summary.json](../data/processed/food_tcm/microbetcm_first_pass_summary.json)
- [fgmdi_high_frequency_food_first_pass.csv](../data/processed/food_tcm/fgmdi_high_frequency_food_first_pass.csv)
- [food_tcm_first_pass_build_summary.json](../data/processed/food_tcm/food_tcm_first_pass_build_summary.json)

生成脚本：

- [build_food_tcm_first_pass_tables.py](../scripts/build_food_tcm_first_pass_tables.py)

## 这次实际用到的公开入口

### 1. MicrobeTCM

当前最有用的是它的公开静态 JSON：

- 微生物目录：
  `https://www.microbetcm.com/microbetcm/static/BROWSE_JSON/Microbe.json`
- 关系主表：
  `https://www.microbetcm.com/microbetcm/static/AllDataForFigure.json`

这次第一轮命中清单和关系子表，就是从这两个入口解析出来的。

### 2. FGMDI

当前能稳定访问到的是文章官方页面和 DOI：

- DOI：
  `https://doi.org/10.1016/j.fbio.2024.104091`

但截至这次整理，公开可访问页面没有直接提供一个可程序化下载的 row-level 关系表。
所以当前 [fgmdi_high_frequency_food_first_pass.csv](../data/processed/food_tcm/fgmdi_high_frequency_food_first_pass.csv) 的定位是：

- 对高频菌先落一张标准化结果表；
- 能从官方文章正文示例直接确认的先填入；
- 其余条目标记为 `await_public_row_level_export`，避免误当成已经命中的数据库行。

## 三张表分别做什么

### 1. microbe_food_relation_table.csv

按 `nt_code` 先铺一行一个面板菌，后面逐步补：

- `food_name`
- `food_category`
- `food_component`
- `relation_direction`
- `relation_scope`
- `evidence_type`
- `host_context`
- `disease_context`
- `source_name`
- `source_url`
- `source_record_id`
- `pmid_or_doi`

推荐 `relation_direction` 用这些值：

- `increase`
- `decrease`
- `bidirectional`
- `associated`
- `unclear`

### 2. microbe_tcm_relation_table.csv

结构和 food 表类似，只是主实体换成：

- `tcm_name`
- `tcm_type`
- `tcm_component`

推荐 `tcm_type` 用这些值：

- `single_herb`
- `formula`
- `component`
- `acupoint`
- `acupoint_formula`

### 3. source_registry.csv

记录哪些数据库或文献是当前整理流程里真正用到的来源，便于后面追溯。

## 推荐整理顺序

### 第一步：先跑数据库覆盖

优先顺序建议：

1. `FGMDI`
2. `gutMDisorder`
3. `MicrobeTCM`
4. `gutMGene v2.0`

原因是前 3 个最适合补“关系表”，第 4 个更偏“机制解释”。

### 第二步：先做 species 精确匹配，再做 genus 代理

推荐匹配顺序：

1. `species_label` 完全匹配
2. `species_name` 近似匹配
3. `genus` 级代理

不要一开始就只做 genus 级别，否则后面会把很多关系放大解释。

### 第三步：把关系和上下文分开写

建议每条记录最少分清：

- 这是 `increase` 还是 `decrease`
- 这是在 `健康人 / 疾病人群 / 动物 / 体外`
- 这是 `food` 本体、`food component`，还是 `herb / formula / component`

否则后面很容易把不同研究条件下的结果混在一起。

## 当前最现实的工作目标

如果不是要完整复现所有文献，而是先为当前模型服务，我建议优先整理以下高频菌：

- `Akkermansia muciniphila`
- `Bacteroides vulgatus`
- `Bacteroides uniformis`
- `Bacteroides thetaiotaomicron`
- `Bifidobacterium adolescentis`
- `Bilophila wadsworthia`
- `Blautia obeum`
- `Collinsella aerofaciens`
- `Eggerthella lenta`
- `Faecalibacterium prausnitzii`
- `Fusobacterium nucleatum`
- `Lactobacillus acidophilus`
- `Lactobacillus gasseri`
- `Parabacteroides distasonis`
- `Prevotella copri`
- `Roseburia hominis`
- `Roseburia intestinalis`
- `Ruminococcus bromii`
- `Ruminococcus gnavus`

这些菌在饮食、宿主代谢和中药相关文献里通常更常见，也更有机会在现成数据库里直接命中。

## 推荐命令

初始化模板：

```bash
/tmp/microbe_env/bin/python scripts/prepare_food_tcm_microbe_reference_tables.py --overwrite
```

## 当前边界

这套模板现在解决的是“落点”和“字段规范”问题，还没有自动抓取远程数据库条目。

所以当前状态应该理解为：

- 已经有固定的数据层可写入；
- 已经把最值得优先用的外部资源整理出来；
- 但还没有把 `83` 个菌在外部库里的命中关系批量导入。
