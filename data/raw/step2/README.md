# Step 2 Raw Data Layout

把 Step 2 原始文件按来源放到这些目录：

- `data/raw/step2/zimmermann_2019`
- `data/raw/step2/javdan_2020`
- `data/raw/step2/agora2`
- `data/raw/step2/gutmgene_v2`

当前仓库已经提供：

- `scripts/normalize_step2_generic.py`
- `scripts/assemble_step2_inputs.py`

推荐顺序：

1. 先把原始表下载到对应目录。
2. 手工挑出最核心的一张或几张 `CSV/TSV/XLSX` 表。
3. 用 `normalize_step2_generic.py` 先转成统一 schema。
4. 再用 `assemble_step2_inputs.py` 把 Step 1 hybrid 结果和 Step 2 标签拼成建模表。
