# Step 1 / Step 2 / Step 3 结构图

这份文档把仓库当前的三步主链路整理成可直接复用的结构图，优先对应现有实现，而不是只画概念蓝图。

## 总览结构图

```mermaid
flowchart LR
    A[原始数据源<br/>Maier / MDIPID / MASI<br/>Zimmermann / Javdan<br/>AGORA2 / gutMGene / cohort] --> B[标准化与特征准备<br/>drug / microbe / pair schema]
    B --> C[Step 1<br/>drug -> microbe effect]
    C --> D[Step 1 预测输出<br/>effect_score / inhibit_probability / label]
    D --> E[Step 2 候选组装<br/>drug x microbe pairs]
    E --> F[Step 2<br/>microbe -> drug metabolism]
    F --> G[Step 2 预测输出<br/>metabolized_probability / depletion / mechanism evidence]
    G --> H[Step 3<br/>community simulation]
    H --> I[场景结果输出<br/>trajectory / health index / development score]

    J[真实 cohort 或内置场景<br/>healthy_reference / high_fiber<br/>high_fat / antibiotic_perturbed] --> H
    K[交叉喂养与健康签名先验<br/>cross_feeding / health signature / enzyme prior] --> H
    L[Step 2 机制证据] --> H
```

## Step 1 结构图

```mermaid
flowchart LR
    A1[输入<br/>drug + microbe + assay] --> B1[数据标准化<br/>药物表 / 微生物表 / 交互表]
    B1 --> C1[特征构建<br/>SMILES / RDKit / fingerprint<br/>taxonomy / phenotype]
    C1 --> D1[Baseline 与 Hybrid 推理<br/>分类头 + 回归头]
    D1 --> E1[输出<br/>inhibit / promote / no_effect<br/>effect_score]

    F1[Gold<br/>Maier 2018] --> B1
    G1[Silver<br/>MDIPID / MASI] --> B1
```

## Step 2 结构图

```mermaid
flowchart LR
    A2[输入一<br/>Step 1 预测结果] --> B2[候选 pair 组装]
    A3[输入二<br/>Zimmermann / Javdan 标签表<br/>AGORA2 / gutMGene / Rhea 先验] --> C2[标签标准化与机制整理]
    B2 --> D2[Step 2 建模表]
    C2 --> D2
    D2 --> E2[Baseline 预测<br/>是否代谢 + depletion]
    E2 --> F2[机制层补充<br/>reaction / enzyme / product evidence]
    F2 --> G2[输出<br/>metabolized_probability<br/>parent_depletion_fraction<br/>mechanism evidence]
```

## Step 3 结构图9

```mermaid
flowchart LR
    A4[输入一<br/>Step 2 集成预测表<br/>内含 Step 1 + Step 2 结果] --> B4[药物切片与目标菌群映射]
    A5[输入二<br/>内置场景或真实 cohort] --> C4[初始化群落状态]
    A6[输入三<br/>cross_feeding / health signature / enzyme prior] --> D4[交互与评分先验]
    B4 --> E4[离散时间模拟器]
    C4 --> E4
    D4 --> E4
    E4 --> F4[每个时间步更新<br/>菌群丰度 / 母药浓度 / 代谢池 / 健康指数]
    F4 --> G4[输出<br/>trajectory tables<br/>top microbe changes<br/>summary / development score]
```

## 当前代码映射

- Step 1：
  - `scripts/download_step1_data.py`
  - `scripts/normalize_step1_data.py`
  - `scripts/train_step1_baseline.py`
  - `scripts/predict_step1_hybrid.py`
- Step 2：
  - `scripts/normalize_step2_zimmermann.py`
  - `scripts/assemble_step2_inputs.py`
  - `scripts/train_step2_baseline.py`
  - `scripts/predict_step2_baseline.py`
- Step 3：
  - `scripts/run_step3_simulation.py`
  - `scripts/prepare_step3_cohort_community.py`
  - `scripts/screen_step3_candidates.py`

## 一句话理解

- Step 1：先判断“药物怎么影响菌”。
- Step 2：再判断“菌怎么处理药物”。
- Step 3：最后把前两步结果放进群落场景里，模拟“给药后整个系统会怎么变”。
