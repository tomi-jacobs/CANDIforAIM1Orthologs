# 🌿 CAnDI Analysis of AIM1 Orthologs — *Aizoaceae* Phylogenomics

This repository contains a gene tree concordance and conflict analysis using **CAnDI (Conflict and Duplication Identifier)** on 4,845 orthologous gene trees from *Aizoaceae* (ice plants), representing 39 taxa. The analysis compares two ASTRAL species trees — one built from all orthologs and one built from a filtered high-quality subset — to assess phylogenetic signal and gene tree discordance across the species tree.

---

## 🧬 Biological Context

This work is part of AIM 1 of a broader phylogenomics project on *Aizoaceae* (39 taxa). Orthologs were extracted from homologous gene families using a Maximum Inclusion (MI) paralog pruning approach. The central question addressed here is: **do high-quality orthologs (selected using information-theoretic scores and Pythia difficulty) recover the same species tree topology and concordance patterns as the full ortholog dataset?**

---

## 📂 Repository Structure

```
CANDIforAIM1Orthologs/
│
├── ASTRAL_all_orthologs.tre          # ASTRAL species tree from all 4,845 orthologs
├── ASTRAL_filtered_orthologs.tre     # ASTRAL species tree from 1,095 filtered orthologs
├── Correlation_TreeStats.csv         # TC, RTC, TCA, RTCA, Taxa, Pythia scores per ortholog
├── filtered_orthologs.csv            # List of orthologs passing the quality filter
│
├── RootedTrees/                      # 4,845 ortholog gene trees rooted with RootDigger v1.7.0
├── RootDigger_checkpoints/           # RootDigger checkpoint files
│
├── CAnDI_all_95cutoff/               # CAnDI results (all orthologs, bootstrap cutoff = 95)
├── CAnDI_filtered_95cutoff/          # CAnDI results (filtered orthologs, bootstrap cutoff = 95)
├── CAnDI_all_nocutoff/               # CAnDI results (all orthologs, no cutoff) → total_analyzed.tre
├── CAnDI_filtered_nocutoff/          # CAnDI results (filtered orthologs, no cutoff) → total_analyzed.tre
│
├── PieCharts_all/                    # Pie chart SVGs for all orthologs (Node_0 to Node_36)
├── PieCharts_filtered/               # Pie chart SVGs for filtered orthologs (Node_0 to Node_36)
│
└── README.md
```

---

## 🔬 Methods

### 1. Ortholog Filtering
Before running the concordance analysis, orthologs were pre-screened for quality using two independent metrics:

- **Pythia difficulty** (Haag & Stamatakis, 2025): A machine-learning predictor of phylogenetic analysis difficulty on a scale of 0 (easy) to 1 (difficult). Orthologs with Pythia difficulty < 0.5 were retained.
- **abs(TCA)**: The absolute Transfer Concordance with Adjustments score, computed using RAxML. Orthologs with abs(TCA) in the top 25% were retained.

Applying both filters yielded **1,095 high-quality orthologs** (23% of 4,845 total). Spearman correlations confirmed a strong negative relationship between all four information-theoretic scores (TC, RTC, TCA, RTCA) and Pythia difficulty (r ≈ -0.66, p < 0.001), validating the use of these scores as proxies for phylogenetic signal quality.

### 2. ASTRAL Species Tree Inference
Two ASTRAL (v5.7.8) species trees were inferred using coalescent-based species tree estimation:

```bash
# All orthologs
java -jar astral.5.7.8.jar \
  -i all_orthologs_trees.tre \
  -o ASTRAL_all_orthologs.tre

# Filtered orthologs
java -jar astral.5.7.8.jar \
  -i filtered_orthologs_trees.tre \
  -o ASTRAL_filtered_orthologs.tre
```

The two trees were compared using RF distance (RAxML-NG v1.2.2), yielding **RF distance = 0** — the two topologies are identical. This confirms that the full ortholog dataset is phylogenetically robust and that filtering does not change the inferred topology.

### 3. Gene Tree Rooting
All 4,845 ortholog gene trees were rooted using **RootDigger v1.7.0** prior to CAnDI analysis:

```bash
cd InfoToCalc/

for aln in *-cln; do
    prefix=~/CANDIforAIM1Orthologs/${aln}.treefile
    if [ ! -f "${prefix}.rooted.tree" ]; then
        ./rd \
            --msa ${aln} \
            --tree ${aln}.treefile \
            --threads 2 \
            --prefix ${prefix}
    fi
done
```

### 4. CAnDI Concordance Analysis
CAnDI (Conflict And Duplication Identifier) was run in normal mode (`--mode n`) on the rooted ortholog gene trees against each species tree. Two runs were performed per species tree:

**With bootstrap cutoff (95)** — nodes below 95% bootstrap support are excluded from analysis:
```bash
python3 CAnDI.py \
  --mode n \
  --species_tree ASTRAL_all_orthologs.tre \
  --gene_folder RootedTrees/ \
  --outfile_prefix CAnDI_all_95cutoff/CAnDI_all \
  --cutoff 95
```

**Without cutoff** — required to generate `_total_analyzed.tre` for Pie.py:
```bash
python3 CAnDI.py \
  --mode n \
  --species_tree ASTRAL_all_orthologs.tre \
  --gene_folder RootedTrees/ \
  --outfile_prefix CAnDI_all_nocutoff/CAnDI_all_nocutoff
```

The same was repeated for the filtered orthologs species tree.

### 5. Pie Chart Visualization
Pie charts showing concordance, conflict, and uninformative proportions per node were generated using CAnDI's `Pie.py`:

```bash
python3 Pie.py \
  -f CAnDI_all_95cutoff/ \
  -p CAnDI_all \
  -t CAnDI_all_nocutoff/CAnDI_all_nocutoff_total_analyzed.tre \
  -o le \
  -a 3 \
  -at 4841
```

The `-at 4841` flag specifies the total number of orthologs analyzed (the maximum node count from `_total_analyzed.tre`), ensuring uninformative gene trees are correctly represented in the pie charts.

Each `Node_X.svg` file corresponds to a node in the species tree as labelled in `CAnDI_all_labels.tre`.

---

## 📊 Key Results

| Analysis | Species Tree | RF Distance |
|----------|-------------|-------------|
| All orthologs (n=4,845) | ASTRAL_all_orthologs.tre | — |
| Filtered orthologs (n=1,095) | ASTRAL_filtered_orthologs.tre | 0 |

- Both species trees are **topologically identical** (RF distance = 0)
- The filtered subset recovers the same phylogeny using only 23% of the data
- Pie charts reveal patterns of concordance and conflict across all 37 internal nodes

---

## 🛠️ Tools & Versions

| Tool | Version | Purpose |
|------|---------|---------|
| ASTRAL | 5.7.8 | Coalescent species tree inference |
| RootDigger | 1.7.0 | Gene tree rooting |
| CAnDI | — | Concordance/conflict analysis |
| RAxML-NG | 1.2.2 | RF distance computation |
| Pythia | 2.0.0 | Phylogenetic difficulty prediction |

---

## 📚 Citations

> Haag, J. & Stamatakis, A. (2025). *Pythia 2.0: New Data, New Prediction Model, New Features*. bioRxiv. https://doi.org/10.1101/2025.03.25.645182

> Kozlov, A.M. et al. (2019). *RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference*. Bioinformatics, 35(21):4453–4455.

> Zhang, C. et al. (2018). *ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees*. BMC Bioinformatics, 19(S6):153.

---

## 👤 Author

**Tomi Jacobs** — PhD Candidate, Computational Biology / Phylogenomics

> Robertson, H.M., Walker, J.F. & Moyroud, E. (2025). *CAnDI: a new tool to investigate conflict in homologous gene trees and explain convergent trait evolution*. Systematic Biology, syaf028. https://doi.org/10.1093/sysbio/syaf028
