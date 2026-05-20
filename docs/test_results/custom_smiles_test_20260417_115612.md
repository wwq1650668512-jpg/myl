# Custom SMILES Prediction Test Result

- 生成时间: 2026-04-17T11:56:46.615358
- 输入 SMILES: `CC1=C(C#CC2=CN=C3C=CC(N4CCC(O)CC4)=NN23)C=CC=C1C(=O)NC1=CC(C(F)(F)F)=CC(Cl)=C1`
- session_id: `140c5916c057`

## Top Effect Microbes (Top 12)

| nt_code | microbe_label              | effect_label | inhibit_prob |        effect_score | promote_prob_refined |
| ------- | -------------------------- | ------------ | -----------: | ------------------: | -------------------: |
| NT5009  | Eubacterium rectale        | inhibit      |   0.93809116 | -0.3686121686430976 |  0.08084375574333102 |
| NT5011  | Roseburia intestinalis     | inhibit      |   0.94622815 | -0.2793047261621862 |   0.0926712756935332 |
| NT5019  | Prevotella copri           | inhibit      |    0.9161046 | -0.2219468195247024 |  0.09217189236938585 |
| NT5074  | Parabacteroides distasonis | inhibit      |    0.8215716 | -0.2209203605641761 |  0.10534269357477798 |
| NT5071  | Parabacteroides merdae     | inhibit      |    0.6658631 | -0.1940053197901607 |  0.10537902884614621 |
| NT5032  | Clostridium perfringens    | inhibit      |   0.90427387 | -0.1835528338477287 |  0.10604981292783548 |
| NT5001  | Bacteroides vulgatus       | inhibit      |    0.8850944 | -0.1758968979137377 |  0.10792661123307803 |
| NT5079  | Roseburia hominis          | inhibit      |     0.920772 | -0.1553107813101785 |  0.09279204113707144 |
| NT5031  | Clostridium perfringens    | inhibit      |   0.89557946 | -0.1484363651341242 |  0.20190675657489657 |
| NT5003  | Bacteroides fragilis (NT)  | inhibit      |   0.67193365 | -0.1461086807086701 |   0.1079866086175289 |
| NT5028  | Bifidobacterium longum     | inhibit      |      0.73942 | -0.1430579850276739 |   0.0801957662631507 |
| NT5033  | Bacteroides fragilis (ET)  | inhibit      |   0.67284924 | -0.1412880560748628 |  0.10800956237919157 |

## Top Metabolism Microbes (Top 12)

| nt_code | microbe_label                | metabolism_label |   metabolized_prob |  depletion_fraction | enzyme_prior | enzyme_names |
| ------- | ---------------------------- | ---------------- | -----------------: | ------------------: | ------------ | ------------ |
| NT5002  | Bacteroides uniformis        | metabolized      |  0.340581579956787 | -0.1375727435391459 | False        | None         |
| NT5065  | Bacteroides uniformis        | metabolized      |  0.340581579956787 | -0.1375727435391459 | False        | None         |
| NT5057  | Bacteroides fragilis         | metabolized      | 0.3389619988528928 |  -0.106990013797591 | False        | None         |
| NT5001  | Bacteroides vulgatus         | metabolized      | 0.3354028941418317 | -0.1291324590051613 | False        | None         |
| NT5064  | Bacteroides xylanisolvens    | metabolized      | 0.3269079063683568 | -0.1103991264654787 | False        | None         |
| NT5033  | Bacteroides fragilis (ET)    | metabolized      | 0.3239800542868161 | -0.1084113050369072 | False        | None         |
| NT5003  | Bacteroides fragilis (NT)    | metabolized      | 0.3239800542868161 | -0.1084113050369072 | False        | None         |
| NT5004  | Bacteroides thetaiotaomicron | metabolized      | 0.3234675606015362 | -0.0982455157692305 | False        | None         |
| NT5054  | Bacteroides ovatus           | metabolized      | 0.3205594121365858 | -0.1069547717035739 | False        | None         |
| NT5050  | Bacteroides caccae           | metabolized      | 0.3176762732738805 | -0.0947484629734151 | False        | None         |
| NT5081  | Odoribacter splanchnicus     | not_metabolized  | 0.2823168689934499 | -0.1215339745726494 | False        | None         |
| NT5019  | Prevotella copri             | not_metabolized  | 0.2752028547363897 | -0.1083254696581194 | False        | None         |

## Candidate Diseases (Top 10)

| disease_name                  | support_score |
| ----------------------------- | ------------: |
| 痔疮                          |        0.3922 |
| 克罗恩病（CD）                |         0.353 |
| 类风湿关节炎（RA）            |        0.3274 |
| 结直肠癌（CRC）               |        0.3093 |
| 系统性红斑狼疮（SLE）         |        0.2627 |
| 肛周脓肿（Anorectal Abscess） |        0.2547 |
| 便秘（Constipation）          |        0.2352 |
| 肛瘘（Anal Fistula）          |        0.2284 |
| 肠易激综合征（IBS）           |        0.1362 |
| 肠易激综合征-腹泻型（IBS-D）  |        0.1362 |

完整原始结果见: `/mnt/e/毕业/docs/test_results/custom_smiles_test_20260417_115612.json`
