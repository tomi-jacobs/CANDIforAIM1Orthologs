# CAnDI Analysis of AIM1 Orthologs — Aizoaceae Phylogenomics

This repository contains a gene tree concordance and conflict analysis using **CAnDI (Conflict and Duplication Identifier)** on 4,845 orthologous gene trees from *Aizoaceae* (ice plants), representing 39 taxa. The analysis compares two ASTRAL species trees — one built from all orthologs and one built from a filtered high-quality subset — to assess phylogenetic signal and gene tree discordance across the species tree.

---

## Biological Context

This work is part of AIM 1 of a broader phylogenomics project on *Aizoaceae* (39 taxa). Orthologs were extracted from homologous gene families using a Maximum Inclusion (MI) paralog pruning approach. The central question addressed here is: **do high-quality orthologs (selected using information-theoretic scores and Pythia difficulty) recover the same species tree topology and concordance patterns as the full ortholog dataset?**

---

## Methods

### 1. Ortholog Filtering
Before running the concordance analysis, orthologs were pre-screened for quality using two independent metrics:

- **Pythia difficulty** (Haag & Stamatakis, 2025): A machine-learning predictor of phylogenetic analysis difficulty on a scale of 0 (easy) to 1 (difficult). Orthologs with Pythia difficulty < 0.5 were retained.
- **abs(TCA)**: The absolute Transfer Concordance with Adjustments score, computed using RAxML. Orthologs with abs(TCA) in the top 25% were retained.

Applying both filters yielded **1,095 high-quality orthologs** (23% of 4,845 total). Spearman correlations confirmed a strong negative relationship between all four information-theoretic scores (TC, RTC, TCA, RTCA) and Pythia difficulty (r = -0.66, p < 0.001), validating the use of these scores as proxies for phylogenetic signal quality.

### 2. ASTRAL Species Tree Inference
Two ASTRAL (v5.7.8) species trees were inferred using coalescent-based species tree estimation:
```
java -jar astral.5.7.8.jar -i all_orthologs_trees.tre -o ASTRAL_all_orthologs.tre
java -jar astral.5.7.8.jar -i filtered_orthologs_trees.tre -o ASTRAL_filtered_orthologs.tre
```

The two trees were compared using RF distance (RAxML-NG v1.2.2), yielding **RF distance = 0** — the two topologies are identical.

### 3. Gene Tree Rooting
All 4,845 ortholog gene trees were rooted using RootDigger v1.7.0 prior to CAnDI analysis.

### 4. CAnDI Concordance Analysis
CAnDI was run in normal mode on the rooted ortholog gene trees against each species tree, with and without a bootstrap cutoff of 95.

### 5. Pie Chart Visualization
Pie charts showing concordance, conflict, and uninformative proportions per node were generated using CAnDI Pie.py with -at 4841.

---

## Key Results

- Both species trees are topologically identical (RF distance = 0)
- The filtered subset recovers the same phylogeny using only 23% of the data
- Pie charts reveal patterns of concordance and conflict across all 37 internal nodes

---

## Tools and Versions

| Tool | Version | Purpose |
|------|---------|---------|
| ASTRAL | 5.7.8 | Coalescent species tree inference |
| RootDigger | 1.7.0 | Gene tree rooting |
| CAnDI | — | Concordance/conflict analysis |
| RAxML-NG | 1.2.2 | RF distance computation |
| Pythia | 2.0.0 | Phylogenetic difficulty prediction |

---

## Citations

> Haag, J. & Stamatakis, A. (2025). Pythia 2.0: New Data, New Prediction Model, New Features. bioRxiv. https://doi.org/10.1101/2025.03.25.645182

> Kozlov, A.M. et al. (2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. Bioinformatics, 35(21):4453-4455.

> Zhang, C. et al. (2018). ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees. BMC Bioinformatics, 19(S6):153.

> Robertson, H.M., Walker, J.F. & Moyroud, E. (2025). CAnDI: a new tool to investigate conflict in homologous gene trees and explain convergent trait evolution. Systematic Biology, syaf028. https://doi.org/10.1093/sysbio/syaf028

---

## Repository Structure
```
CANDIforAIM1Orthologs/
├── ASTRAL_all_orthologs.tre
├── ASTRAL_filtered_orthologs.tre
├── Correlation_TreeStats.csv
├── filtered_orthologs.csv
├── RootedTrees/
├── RootDigger_checkpoints/
├── CAnDI_all_95cutoff/
├── CAnDI_filtered_95cutoff/
├── CAnDI_all_nocutoff/
├── CAnDI_filtered_nocutoff/
├── PieCharts_all/
├── PieCharts_filtered/
└── README.md
```

---

## Author

**Tomi Jacobs** — PhD Candidate, Computational Biology / Phylogenomics
