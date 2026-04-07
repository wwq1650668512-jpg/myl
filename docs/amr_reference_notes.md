# AMR References and Usage Notes

## 目的

这份文档整理当前 AMR 修正相关的外部文献与数据资源，并注明每条来源在项目里支持了什么判断、落在哪个文件或规则上。

阅读建议：

- 如果想看“当前已经真正写进规则里的证据”，看“已实际用于当前实现”
- 如果想看“后面准备接入的数据源”，看“已讨论、待后续接入”

## 已实际用于当前实现

### 1. EUCAST Expected Phenotypes

- 类型：官方方法学 / 预期表型原则
- 链接：https://www.eucast.org/bacteria/important-additional-information/expected-phenotypes/
- 参考点：
  - `expected phenotype` 的定义
  - 对“与预期表型相反的实验或预测结果应谨慎看待”的方法学依据
  - 当前 `expected_phenotype` / `amr_conflict_flag` 这套字段设计的概念来源
- 当前落地：
  - [amr.py](../src/gut_drug_microbiome/amr.py#L206)
  - [amr_integration_plan.md](amr_integration_plan.md#L93)
- 说明：
  - 这份来源主要提供“怎么理解预期耐药/敏感”的原则
  - 当前并没有直接从 EUCAST 页面自动抽取 `Bacteroides × penicillin` 的结构化规则表

### 2. NCBI Bookshelf: Bacteroides Fragilis - StatPearls

- 类型：综述 / 背景资料
- 链接：https://www.ncbi.nlm.nih.gov/sites/books/NBK553032/
- 参考点：
  - `Bacteroides fragilis` 与临床耐药背景
  - `Bacteroides` 相关的 `beta-lactamase` 与耐药认知
  - 用作当前属级 `Bacteroides × beta_lactam` 支持性先验的背景来源
- 当前落地：
  - [seed_beta_lactam_amr_rules.py](../scripts/seed_beta_lactam_amr_rules.py#L58)
  - [drug_resistance_rules.csv](../data/processed/amr/drug_resistance_rules.csv#L2)
  - [microbe_amr_reference.csv](../data/processed/amr/microbe_amr_reference.csv#L3)
- 对应规则：
  - `amr_bacteroides_beta_lactam_supporting`
  - `amr_buniformis_penicillin_moderate`
  - `amr_btheta_penicillin_moderate`
  - `amr_bfragilis_penicillin_strong`
- 说明：
  - 当前更多是“属级/物种级保守先验”的依据，不是替代菌株级 AST 真值

### 3. PMC352644: Factors Contributing to Resistance to Beta-Lactam Antibiotics in Bacteroides fragilis

- 类型：机制文献
- 链接：https://pmc.ncbi.nlm.nih.gov/articles/PMC352644/
- 参考点：
  - `Bacteroides fragilis` 对 `beta-lactam / penicillin` 耐药与 `beta-lactamase` 机制
  - 支撑当前 `penicillin` 规则强于一般 `beta_lactam` 规则
- 当前落地：
  - [seed_beta_lactam_amr_rules.py](../scripts/seed_beta_lactam_amr_rules.py#L74)
  - [drug_resistance_rules.csv](../data/processed/amr/drug_resistance_rules.csv#L3)
- 对应规则：
  - `amr_bacteroides_penicillin_strong`
- 说明：
  - 这条文献主要支撑“对明确 penicillin-like 药物，应该施加更强耐药先验”

### 4. PMC187887 / PubMed 8517690: Genetic and biochemical analysis of a novel Ambler class A beta-lactamase responsible for cefoxitin resistance in Bacteroides species

- 类型：原始研究
- PMC 链接：https://pmc.ncbi.nlm.nih.gov/articles/PMC187887/
- PubMed 链接：https://pubmed.ncbi.nlm.nih.gov/8517690/
- 参考点：
  - `Bacteroides vulgatus` 中 `beta-lactamase` 与 `beta-lactam` 耐药的直接物种级证据
  - 支撑把 `B. vulgatus` 的规则强度设为高于一般属级规则
- 当前落地：
  - [seed_beta_lactam_amr_rules.py](../scripts/seed_beta_lactam_amr_rules.py#L89)
  - [drug_resistance_rules.csv](../data/processed/amr/drug_resistance_rules.csv#L4)
- 对应规则：
  - `amr_bvulgatus_penicillin_strong`
- 说明：
  - 这是当前种子规则里最直接的物种级来源之一

### 5. PubChem Compound 5904: Penicillin G

- 类型：化合物数据库
- 链接：https://pubchem.ncbi.nlm.nih.gov/compound/5904
- 参考点：
  - 用于确认用户输入的自定义 SMILES 对应的是 `Penicillin G / Benzylpenicillin`
  - 支撑“当前修正场景确实是 penicillin-like 药物，不是泛指抗生素类别”
- 当前落地：
  - 没有直接写入代码规则
  - 主要用于解释和复核自定义 SMILES 的药物身份
- 说明：
  - 这条来源服务于“药物身份确认”，不是直接的 AMR 规则来源

### 6. NCBI Bookshelf: Vancomycin - StatPearls

- 类型：综述 / 药物背景资料
- 链接：https://www.ncbi.nlm.nih.gov/books/NBK459263/
- 参考点：
  - `vancomycin / glycopeptide` 的作用范围主要偏向革兰阳性菌
  - 作为当前 `gram-negative × vancomycin / glycopeptide` 高风险假阳性修正的来源
- 当前落地：
  - [seed_high_risk_antibiotic_amr_rules.py](../scripts/seed_high_risk_antibiotic_amr_rules.py#L47)
  - [drug_resistance_rules.csv](../data/processed/amr/drug_resistance_rules.csv#L8)
- 对应规则：
  - 当前 9 个革兰阴性属的：
    - `*_vancomycin_strong`
    - `*_glycopeptide_strong`
- 说明：
  - 这是“类边界非常明确”的规则，优先用来压住明显不合理的假阳性

### 7. PMC4846043: Daptomycin review

- 类型：综述
- 链接：https://pmc.ncbi.nlm.nih.gov/articles/PMC4846043/
- 参考点：
  - `daptomycin` 的活性范围偏向革兰阳性菌
  - 作为当前 `gram-negative × daptomycin / lipopeptide` 修正的来源
- 当前落地：
  - [seed_high_risk_antibiotic_amr_rules.py](../scripts/seed_high_risk_antibiotic_amr_rules.py#L47)
  - [drug_resistance_rules.csv](../data/processed/amr/drug_resistance_rules.csv#L10)
- 对应规则：
  - 当前 9 个革兰阴性属的：
    - `*_daptomycin_strong`
    - `*_lipopeptide_moderate`
- 说明：
  - 当前把显式 `daptomycin` 命名设成强规则，把泛 `lipopeptide` 设成一档更保守的中等强度

### 8. PubMed 10467540: Gram-positive organisms intrinsically resistant to vancomycin

- 类型：综述 / 提醒性资料
- 链接：https://pubmed.ncbi.nlm.nih.gov/10467540/
- 参考点：
  - `Lactobacillus` 的内在 `vancomycin` 耐药
  - 作为当前 `Lactobacillus × vancomycin / glycopeptide` 修正的来源
- 当前落地：
  - [seed_high_risk_antibiotic_amr_rules.py](../scripts/seed_high_risk_antibiotic_amr_rules.py#L113)
  - [drug_resistance_rules.csv](../data/processed/amr/drug_resistance_rules.csv#L44)
- 对应规则：
  - `amr_lactobacillus_vancomycin_strong`
  - `amr_lactobacillus_glycopeptide_moderate`
- 说明：
  - 这条是本轮扩展里唯一新增的革兰阳性高风险假阳性修正

### 9. Gentamicin - StatPearls

- 类型：综述 / 药物背景资料
- 链接：https://www.ncbi.nlm.nih.gov/books/NBK557550/
- 参考点：
  - `gentamicin` 对厌氧菌通常不活跃
  - `aminoglycoside` 进入依赖氧化依赖转运，严格厌氧环境下难以高效摄取
  - 作为当前 `strict anaerobes × gentamicin / aminoglycoside` 修正的主要来源
- 当前落地：
  - [seed_aminoglycoside_anaerobe_rules.py](../scripts/seed_aminoglycoside_anaerobe_rules.py#L51)
  - [drug_resistance_rules.csv](../data/processed/amr/drug_resistance_rules.csv#L46)
- 对应规则：
  - 当前 18 个面板内严格厌氧属的：
    - `*_gentamicin_strong`
    - `*_aminoglycoside_moderate`
- 说明：
  - 这是本轮扩展里最主要的机制来源

### 10. Aminoglycosides - StatPearls

- 类型：综述
- 链接：https://www.ncbi.nlm.nih.gov/books/NBK541105/
- 参考点：
  - `amikacin / tobramycin / streptomycin` 属于 `aminoglycoside` 类
  - 支撑把显式药名规则和泛类规则分开处理
- 当前落地：
  - [seed_aminoglycoside_anaerobe_rules.py](../scripts/seed_aminoglycoside_anaerobe_rules.py#L51)
  - [drug_resistance_rules.csv](../data/processed/amr/drug_resistance_rules.csv#L47)
- 对应规则：
  - 当前 18 个面板内严格厌氧属的：
    - `*_amikacin_strong`
    - `*_tobramycin_strong`
    - `*_streptomycin_moderate`
- 说明：
  - 当前把 `streptomycin` 设成略弱一档，是为了先保守一些

### 11. Anaerobic Infections - StatPearls

- 类型：综述 / 背景资料
- 链接：https://www.ncbi.nlm.nih.gov/books/NBK482349/
- 参考点：
  - 作为“严格厌氧菌”这一层筛选的背景来源
  - 支撑当前把面板内部分常见肠道厌氧属单独拉出做规则先验
- 当前落地：
  - 没有直接写成单条 CSV 来源
  - 主要用于本轮规则的“panel-specific strict anaerobe genera”选择边界
- 说明：
  - 当前严格厌氧属名单仍是工程上的保守筛选，不代表所有属内所有菌株都一概等价

### 12. Polymyxin review

- 类型：综述
- 链接：https://pmc.ncbi.nlm.nih.gov/articles/PMC5766840/
- 参考点：
  - `polymyxin / colistin / polymyxin B` 的作用依赖革兰阴性外膜与 LPS
  - 作为当前 `gram-positive × polymyxin` 高风险假阳性修正的来源
- 当前落地：
  - [seed_polymyxin_and_low_anaerobe_fq_rules.py](../scripts/seed_polymyxin_and_low_anaerobe_fq_rules.py#L69)
  - [drug_resistance_rules.csv](../data/processed/amr/drug_resistance_rules.csv#L136)
- 对应规则：
  - 当前所有革兰阳性属的：
    - `*_colistin_strong`
    - `*_polymyxin_b_strong`
    - `*_polymyxin_moderate`
- 说明：
  - 这是又一类“边界很清楚”的服务层修正规则

### 13. Fluoroquinolones and Anaerobes review

- 类型：综述
- 链接：https://pubmed.ncbi.nlm.nih.gov/10428911/
- 参考点：
  - 并非所有 `fluoroquinolone` 对厌氧菌都同样不活跃
  - `ciprofloxacin / ofloxacin / norfloxacin` 等较早期代表药对厌氧覆盖较弱
  - 作为当前“低厌氧覆盖 fluoroquinolone 子集 × strict anaerobes”规则的来源
- 当前落地：
  - [seed_polymyxin_and_low_anaerobe_fq_rules.py](../scripts/seed_polymyxin_and_low_anaerobe_fq_rules.py#L120)
  - [drug_resistance_rules.csv](../data/processed/amr/drug_resistance_rules.csv#L175)
- 对应规则：
  - 当前 18 个严格厌氧属的：
    - `*_ciprofloxacin_strong`
    - `*_levofloxacin_strong`
    - `*_ofloxacin_strong`
    - `*_norfloxacin_strong`
    - `*_fluoroquinolone_low_anaerobe_moderate`
- 说明：
  - 当前没有做“整个 fluoroquinolone 类”的统一耐药规则
  - 这是为了避免误伤 `moxifloxacin` 一类对厌氧菌覆盖更好的成员

## 当前实现与来源的对应关系

### 规则文件

- [drug_resistance_rules.csv](../data/processed/amr/drug_resistance_rules.csv#L1)
  - `L2`：`Bacteroides × beta_lactam` 支持性规则，主要参考 NCBI Bookshelf
  - `L3`：`Bacteroides × penicillin` 强规则，主要参考 PMC352644
  - `L4`：`Bacteroides vulgatus × penicillin` 强规则，主要参考 PMC187887
  - `L5-L7`：`Bacteroides uniformis / thetaiotaomicron / fragilis × penicillin`，当前属于保守先验，强度低于直接物种证据
  - `L8-L43`：当前 9 个革兰阴性属的 `vancomycin / glycopeptide / daptomycin / lipopeptide` 规则
  - `L44-L45`：`Lactobacillus × vancomycin / glycopeptide` 规则
  - `L46-L135`：当前 18 个严格厌氧属的 `gentamicin / amikacin / tobramycin / streptomycin / aminoglycoside` 规则
  - `L136-L174`：当前所有革兰阳性属的 `colistin / polymyxin_b / polymyxin` 规则
  - `L175-L264`：当前 18 个严格厌氧属的“低厌氧覆盖 fluoroquinolone 子集”规则

### 参考表

- [microbe_amr_reference.csv](../data/processed/amr/microbe_amr_reference.csv#L3)
  - 当前对 `Bacteroides` 属先统一打上：
    - `has_beta_lactamase = likely`
    - `has_intrinsic_amr_evidence = reported`
    - `expected_beta_lactam_resistant = penicillin-like beta-lactams`
  - 这部分是“面板菌参考画像”，不是菌株级最终真值

### 服务层修正逻辑

- [amr.py](../src/gut_drug_microbiome/amr.py#L120)
  - 负责识别 `penicillin / beta_lactam / vancomycin / glycopeptide / daptomycin / lipopeptide / aminoglycoside / polymyxin / low-anaerobe fluoroquinolone`
- [amr.py](../src/gut_drug_microbiome/amr.py#L173)
  - 负责按物种/属级规则匹配最佳来源
- [amr.py](../src/gut_drug_microbiome/amr.py#L206)
  - 负责生成：
    - `display_step1_predicted_effect_label`
    - `amr_conflict_flag`
    - `amr_correction_applied`
- [service.py](../src/gut_drug_microbiome/web/service.py#L433)
  - 负责把修正后的结果接回网页接口输出

## 已讨论、待后续接入的数据资源

这些来源已经在方案讨论里提过，但当前还没有直接写成结构化规则或自动化数据管线。

### 1. CARD

- 类型：AMR 基因与机制数据库
- 链接：https://card.mcmaster.ca/about
- 参考点：
  - 后续给 `microbe_amr_reference.csv` 补：
    - `has_beta_lactamase`
    - `has_efflux_amr`
    - `has_target_alteration`
    - `amr_gene_summary`
- 当前状态：
  - 已在方案中指定为机制层优先来源
  - 尚未自动接入到当前 83 菌面板
- 对应方案位置：
  - [amr_integration_plan.md](amr_integration_plan.md#L100)

### 2. NCBI NDARO

- 类型：官方 AMR 参考数据平台
- 链接：https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/
- 参考点：
  - 作为 AMR gene catalog、reference hierarchy、reference HMM catalog 的官方来源
  - 后续可用于参考基因与规则的标准化
- 当前状态：
  - 已讨论，未接入自动化流程

### 3. NCBI AMRFinderPlus

- 类型：AMR 基因识别工具
- 链接：https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/
- 参考点：
  - 后续可对 83 菌面板参考基因组做 `AMR gene / mutation` 检测
  - 用于把机制信息从“手工先验”升级为“基因组支持”
- 当前状态：
  - 已讨论，未开始跑参考基因组

### 4. NCBI AST Browser

- 类型：表型数据浏览与下载入口
- 链接：https://www.ncbi.nlm.nih.gov/pathogens/ast
- 文档链接：https://www.ncbi.nlm.nih.gov/pathogens/docs/ast/
- 参考点：
  - 后续为 `species_drugclass_phenotype_prior.csv` 提供：
    - `n_tested`
    - `resistant_fraction`
    - `susceptible_fraction`
    - `mic50 / mic90`
  - 文档还特别提醒 AST 数据为提交者提供，需要谨慎解释
- 当前状态：
  - 已在方案中列为表型校准层优先来源
  - 尚未接入当前规则

### 5. BV-BRC AMR

- 类型：AMR 数据说明与整合平台
- 链接：https://www.bv-brc.org/docs/data_protocols/antimicrobial_resistance.html
- 参考点：
  - 后续可补充 phenotype / genotype 关联数据
  - 适合作为 `species_drugclass_phenotype_prior.csv` 的补充来源
- 当前状态：
  - 已讨论，未接入

### 6. Center for Genomic Epidemiology / ResFinder

- 类型：在线 AMR 检测资源集合
- 链接：https://www.genomicepidemiology.org/services/
- 参考点：
  - `ResFinder` 适合 acquired resistance gene 检出
  - `PointFinder` 路线适合后续补靶点突变规则
- 当前状态：
  - 已讨论，未接入

## 使用边界

这批文献和资源当前主要支撑的是“第一版保守修正层”，不是完整的临床 AST 判定系统。

需要特别注意：

- 现在的 `Bacteroides × penicillin` 修正，大多仍是属级或保守物种级先验
- `gram-negative × vancomycin / daptomycin` 这类规则虽然边界更清楚，但仍是“优先压假阳性”的服务层修正
- `strict anaerobes × aminoglycoside` 这类规则依赖“属级厌氧先验”，适合先压明显假阳性，不适合替代菌株级 MIC
- `fluoroquinolone × strict anaerobes` 当前只实现了“低厌氧覆盖子集”，没有推广到整个类
- 它更适合压住明显假阳性，不适合被解释成菌株级精准耐药结论
- 真正往下一步走，还是要补：
  - 参考基因组
  - AMRFinderPlus / CARD 机制注释
  - AST 表型校准

## 建议维护方式

后面每次新增规则时，建议至少同时更新这 3 个位置：

1. [drug_resistance_rules.csv](../data/processed/amr/drug_resistance_rules.csv#L1)
2. [microbe_amr_reference.csv](../data/processed/amr/microbe_amr_reference.csv#L1)
3. 本文档 [amr_reference_notes.md](amr_reference_notes.md#L1)

这样后面回头看每条规则时，能立刻知道：

- 规则来自哪条文献
- 参考的是机制、综述还是表型数据
- 当前只是支持性先验，还是已经有更强的物种级证据
