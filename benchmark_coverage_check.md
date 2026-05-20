# Benchmark Coverage Check

## 1. Disease Reference Coverage
| disease_or_scenario     |   reference_relation_count | in_reference   |
|:------------------------|---------------------------:|:---------------|
| 克罗恩病（CD）                |                          6 | YES            |
| 溃疡性结肠炎（UC）              |                          5 | YES            |
| 肠易激综合征（IBS）             |                          6 | YES            |
| 肠易激综合征-腹泻型（IBS-D）       |                          0 | NO             |
| 肠易激综合征-便秘型（IBS-C）       |                          0 | NO             |
| 便秘（Constipation）        |                         35 | YES            |
| 肛周脓肿（Anorectal Abscess） |                         15 | YES            |

## 2. Weighted Evaluation Panel Coverage (per drug)
| drug          | target                  | present_in_panel   |   rank |   final_score |
|:--------------|:------------------------|:-------------------|-------:|--------------:|
| Rifaximin     | 克罗恩病（CD）                | YES                |      8 |     0.126313  |
| Rifaximin     | 溃疡性结肠炎（UC）              | YES                |     12 |     0         |
| Rifaximin     | 肠易激综合征（IBS）             | YES                |      5 |     0.165405  |
| Rifaximin     | 肠易激综合征-腹泻型（IBS-D）       | YES                |      6 |     0.165405  |
| Rifaximin     | 肠易激综合征-便秘型（IBS-C）       | YES                |      7 |     0.165405  |
| Rifaximin     | 便秘（Constipation）        | YES                |      4 |     0.191486  |
| Rifaximin     | 肛周脓肿（Anorectal Abscess） | YES                |      3 |     0.194537  |
| Vancomycin    | 克罗恩病（CD）                | YES                |     10 |     0.0560752 |
| Vancomycin    | 溃疡性结肠炎（UC）              | YES                |     11 |     0         |
| Vancomycin    | 肠易激综合征（IBS）             | YES                |      6 |     0.109625  |
| Vancomycin    | 肠易激综合征-腹泻型（IBS-D）       | YES                |      5 |     0.109625  |
| Vancomycin    | 肠易激综合征-便秘型（IBS-C）       | YES                |      4 |     0.109625  |
| Vancomycin    | 便秘（Constipation）        | YES                |      3 |     0.112028  |
| Vancomycin    | 肛周脓肿（Anorectal Abscess） | YES                |      8 |     0.101855  |
| Lubiprostone  | 克罗恩病（CD）                | YES                |      8 |     0.143283  |
| Lubiprostone  | 溃疡性结肠炎（UC）              | YES                |     10 |     0         |
| Lubiprostone  | 肠易激综合征（IBS）             | YES                |      3 |     0.180582  |
| Lubiprostone  | 肠易激综合征-腹泻型（IBS-D）       | YES                |      2 |     0.180582  |
| Lubiprostone  | 肠易激综合征-便秘型（IBS-C）       | YES                |      1 |     0.180582  |
| Lubiprostone  | 便秘（Constipation）        | YES                |      5 |     0.162817  |
| Lubiprostone  | 肛周脓肿（Anorectal Abscess） | YES                |      9 |     0.135722  |
| Metronidazole | 克罗恩病（CD）                | YES                |      1 |     0.211572  |
| Metronidazole | 溃疡性结肠炎（UC）              | YES                |     12 |     0         |
| Metronidazole | 肠易激综合征（IBS）             | YES                |      7 |     0.156748  |
| Metronidazole | 肠易激综合征-腹泻型（IBS-D）       | YES                |      8 |     0.156748  |
| Metronidazole | 肠易激综合征-便秘型（IBS-C）       | YES                |      9 |     0.156748  |
| Metronidazole | 便秘（Constipation）        | YES                |     10 |     0.151152  |
| Metronidazole | 肛周脓肿（Anorectal Abscess） | YES                |      4 |     0.182202  |

## 3. Scenario Coverage
| drug          | constipation_related_present   |   constipation_best_rank | infection_abscess_related_present   |   infection_abscess_best_rank |
|:--------------|:-------------------------------|-------------------------:|:------------------------------------|------------------------------:|
| Rifaximin     | YES                            |                        4 | YES                                 |                             2 |
| Vancomycin    | YES                            |                        3 | YES                                 |                             7 |
| Lubiprostone  | YES                            |                        1 | YES                                 |                             8 |
| Metronidazole | YES                            |                        9 | YES                                 |                             1 |

## 4. Key Conclusions
- RFX_01 comparator panel is now complete for IBS/IBS-D/CD/UC scope (UC may be zero-evidence placeholder when no matched relations).
- VAN_02 now has explicit CD entry in panel: rank=10, score=0.0561.
- Vancomycin UC entry is present: rank=11, score=0.0000.