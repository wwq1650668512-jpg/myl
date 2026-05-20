# Research Notes: Drug -> Microbiome -> Mechanism -> Disease

Date: 2026-04-14  
Scope: IBD / CD / CRC / IBS, plus mechanism-layer feature design for this repository.

## 1) High-confidence disease-microbe patterns (for mechanism priors)

### 1.1 IBD / CD
- Conclusion: `Faecalibacterium prausnitzii` is consistently depleted in IBD and lower in active disease.
  - Evidence:
    - Systematic review/meta-analysis (16 studies): CD and UC both showed reduced abundance; active disease lower than remission.
    - Source: https://pubmed.ncbi.nlm.nih.gov/32815163/
  - Modeling rationale:
    - Supports `anti_inflammatory_score`, `butyrate_support_score`, `barrier_protection_score` as beneficial axes.

- Conclusion: UC dysbiosis includes reduced butyrate producers `Roseburia hominis` and `F. prausnitzii`.
  - Evidence:
    - Gut 2014: both species significantly reduced, inversely correlated with disease activity.
    - Source: https://pubmed.ncbi.nlm.nih.gov/24021287/
  - Modeling rationale:
    - Directly justifies a butyrate-centered mechanism feature instead of only generic beneficial/risk genus counts.

- Conclusion: `Ruminococcus gnavus` has mechanistic pro-inflammatory potential in CD context.
  - Evidence:
    - PNAS 2019: identified inflammatory glucorhamnan inducing TNF-alpha via TLR4.
    - Source: https://pubmed.ncbi.nlm.nih.gov/31182571/
  - Modeling rationale:
    - Supports explicit `pro_inflammatory_score` and `mucus_degradation_score` risk channels.

- Conclusion: AIEC (`Escherichia coli` pathovar) remains a strong CD-associated pathobiont hypothesis.
  - Evidence:
    - Gut 2026 review (25-year synthesis): AIEC repeatedly associated with CD, host interaction and virulence evidence summarized.
    - Source: https://pubmed.ncbi.nlm.nih.gov/40473402/
  - Modeling rationale:
    - Supports `pathobiont_load` and `toxin_risk_score` risk channels.

### 1.2 CRC
- Conclusion: `Fusobacterium nucleatum` is repeatedly linked to CRC and can promote tumorigenesis.
  - Evidence:
    - Cell Host & Microbe 2013 mechanistic tumorigenesis/immune microenvironment study.
    - Source: https://pubmed.ncbi.nlm.nih.gov/23954159/
    - Genome Research 2012: enrichment in human CRC.
    - Source: https://pubmed.ncbi.nlm.nih.gov/22009989/
  - Modeling rationale:
    - Supports `pro_inflammatory_score`, `toxin_risk_score`, `pathobiont_load`.

- Conclusion: ETBF is more prevalent in CRC and may rise with advanced stages.
  - Evidence:
    - 2025 systematic review/meta-analysis: higher ETBF prevalence in CRC vs controls.
    - Source: https://pubmed.ncbi.nlm.nih.gov/40125515/
  - Modeling rationale:
    - Supports toxin-centric risk feature (`toxin_risk_score`).

- Conclusion: `pks+ E. coli` context is relevant for CRC risk stratification.
  - Evidence:
    - Prospective cohort molecular-pathological analysis: stronger western-diet association in tumors with higher pks+ E. coli.
    - Source: https://pubmed.ncbi.nlm.nih.gov/35760086/
  - Modeling rationale:
    - Supports `toxin_risk_score` + `pathobiont_load` coupled risk representation.

### 1.3 IBS
- Conclusion: IBS has reproducible microbiota alterations at cohort level, but heterogeneity is high.
  - Evidence:
    - Systematic review/meta-analysis on IBS microbiota alterations.
    - Source: https://pubmed.ncbi.nlm.nih.gov/27300149/
  - Modeling rationale:
    - Motivates mechanism-layer aggregation (more robust than relying on single taxon direction).

- Conclusion: Guideline-level support exists for rifaximin in IBS-D.
  - Evidence:
    - ACG Clinical Guideline (GRADE): recommends rifaximin for global IBS-D symptoms.
    - Source: https://pubmed.ncbi.nlm.nih.gov/33315591/
  - Modeling rationale:
    - Supports clinical plausibility of “microbiota modulation” drug archetype in our case studies.

## 2) Drug-type signal: strong bactericidal vs microbiota-modulating

- Conclusion: Oral vancomycin behaves as a strong-disruption antibiotic for gut ecology.
  - Evidence:
    - Human longitudinal sequencing: broad depletion during therapy, expansion of opportunists (e.g., Klebsiella/Escherichia-Shigella), incomplete recovery in some subjects.
    - Source: https://pubmed.ncbi.nlm.nih.gov/27707993/
  - Modeling rationale:
    - Justifies that some drugs should score high on risk channels (`pathobiont_load`, `pro_inflammatory_score`) despite target efficacy.

- Conclusion: Rifaximin often behaves as a low-absorption, microbiota-modulating/eubiotic antibiotic rather than classic broad-kill systemic antibiotic.
  - Evidence:
    - Review on rifaximin eubiotic properties and microbiota modulation.
    - Source: https://pubmed.ncbi.nlm.nih.gov/28740337/
  - Modeling rationale:
    - Supports mechanism-layer differentiation between “kill-heavy” and “ecology-modulating” drugs.

## 3) Mechanism axes supported by literature

- Conclusion: SCFA/butyrate is a mechanistically meaningful intermediate axis for inflammation/barrier.
  - Evidence:
    - Nature Reviews Immunology 2024: SCFAs regulate epithelial barrier and mucosal/systemic immunity.
    - Source: https://pubmed.ncbi.nlm.nih.gov/38565643/
  - Modeling rationale:
    - Supports `butyrate_support_score` and contributes to `barrier_protection_score`.

- Conclusion: Ecological interaction balance (competition vs cross-feeding dominance) is a plausible dysbiosis marker.
  - Evidence:
    - Science 2026 ENBI: disease-associated states trend toward positive-interaction-dominant regimes.
    - Source: https://pubmed.ncbi.nlm.nih.gov/41747050/
    - Nat Rev Microbiology review: clear conceptual taxonomy for competition vs cooperation/cross-feeding.
    - Source: https://pubmed.ncbi.nlm.nih.gov/31163167/
  - Modeling rationale:
    - Supports `competition_vs_crossfeeding_proxy`.

- Conclusion: Diet-bile acid-pathobiont axis (e.g., `Bilophila wadsworthia`) links ecology to inflammation/barrier dysfunction.
  - Evidence:
    - Nature 2012: saturated-fat/taurocholate promotes Bilophila bloom and colitis in susceptible mice.
    - Source: https://pubmed.ncbi.nlm.nih.gov/22722865/
  - Modeling rationale:
    - Supports `pathobiont_load`, `pro_inflammatory_score`, and barrier-related risk interpretation.

## 4) Mechanism features selected for minimal implementable layer

Selected (implemented target):
- `anti_inflammatory_score`
- `pro_inflammatory_score`
- `butyrate_support_score`
- `barrier_protection_score`
- `toxin_risk_score`
- `mucus_degradation_score`
- `pathobiont_load`
- `competition_vs_crossfeeding_proxy`

Why this subset:
- These map directly to recurrent evidence themes (IBD/CD/CRC/IBS + antibiotic ecology effects).
- They are computable from existing repository columns (`step1` effect/probability, `step2` mechanism/enzyme/cross-feeding fields, taxonomic metadata).
- They are interpretable at microbe-contribution level and support ablation.

## 5) Notes on evidence quality and caution

- Strongest mechanistic anchors in this note:
  - Science ENBI (interaction balance): https://pubmed.ncbi.nlm.nih.gov/41747050/
  - Cell Host & Microbe Fusobacterium CRC mechanism: https://pubmed.ncbi.nlm.nih.gov/23954159/
  - Gut / meta-analysis IBD butyrate-producer depletion: https://pubmed.ncbi.nlm.nih.gov/24021287/ and https://pubmed.ncbi.nlm.nih.gov/32815163/
  - ACG guideline for rifaximin in IBS-D: https://pubmed.ncbi.nlm.nih.gov/33315591/
  - Human oral vancomycin disruption study: https://pubmed.ncbi.nlm.nih.gov/27707993/

- Some organism-level evidence remains context/strain-dependent (especially mucin degraders such as Akkermansia and Ruminococcus species), so mechanism scores should be treated as priors for ranking, not causal clinical proof.
