# Confidence Calibration Analysis

## Scope

- Benchmark assertions: 9
- Fixed test drugs: 4
- Correctness mapping: PASS=1.0, PARTIAL=0.5, FAIL=0.0
- Severe warning set: over-suppression, core-butyrate-suppression, ecology-risk, drug-profile-conflict

## Tier Distribution (PASS/PARTIAL/FAIL)

- high: PASS=0, PARTIAL=0, FAIL=0, pass_rate=N/A, fail_rate=N/A
- medium: PASS=3, PARTIAL=1, FAIL=0, pass_rate=0.750, fail_rate=0.000
- low: PASS=3, PARTIAL=2, FAIL=0, pass_rate=0.600, fail_rate=0.000

## Correlation With Benchmark Correctness

- pearson_confidence_vs_correctness: 0.135 (n=9)
- spearman_confidence_vs_correctness: 0.188 (n=9)

## Benchmark Alignment Checks

- FAIL 对齐（低置信或有 warning）: 0/0 = N/A
- PASS 对齐（无严重 warning）: 2/6 = 0.333

## Quick Read

- 如果 low tier 的 FAIL 占比高、high tier 的 PASS 占比高，说明分层有区分度。
- 如果相关系数为正，说明 confidence_score 趋势与 benchmark correctness 同向。
- 对齐检查可直接监控“FAIL 是否被提醒、PASS 是否不过度报警”。
