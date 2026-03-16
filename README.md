# CAnDI Analysis of AIM1 Orthologs — Aizoaceae Phylogenomics

This repository contains a gene tree concordance and conflict analysis using **CAnDI (Conflict and Duplication Identifier)** on 4,845 orthologous gene trees from *Aizoaceae* (ice plants), representing 39 taxa. The analysis compares two ASTRAL species trees — one built from all orthologs and one built from a filtered high-quality subset — to assess phylogenetic signal and gene tree discordance across the species tree.

---

## Biological Context

This work is part of AIM 1 of a broader phylogenomics project on *Aizoaceae* (39 taxa). Orthologs used in this analysis were previously extracted from homologous gene families during the AIM1 pipeline using a Maximum Inclusion (MI) paralog pruning approach (Yang & Smith, 2014; Morales-Briones et al., 2021), re-aligned with PRANK, filtered with pxclsq (30% occupancy), and gene trees inferred with IQ-TREE2 using ModelFinder automatic model selection and 1000 ultrafast bootstraps. The central question addressed here is: **do high-quality orthologs (selected using information-theoretic scores and Pythia difficulty) recover the same species tree topology and concordance patterns as the full ortholog dataset?**

---

## Methods

### 1. Ortholog Filtering
Before running the concordance analysis, orthologs were pre-screened for quality using two independent metrics:

- **Pythia difficulty** (Haag & Stamatakis, 2025): A machine-learning predictor of phylogenetic analysis difficulty on a scale of 0 (easy) to 1 (difficult). Orthologs with Pythia difficulty < 0.5 were retained.
- **abs(TCA)**: The absolute Transfer Concordance with Adjustments score, computed using RAxML. Orthologs with abs(TCA) in the top 25% were retained.

Applying both filters yielded **1,095 high-quality orthologs** (23% of 4,845 total). Spearman correlations confirmed a strong negative relationship between all four information-theoretic scores (TC, RTC, TCA, RTCA) and Pythia difficulty (r = -0.66, p < 0.001), validating the use of these scores as proxies for phylogenetic signal quality.

### 2. ASTRAL Species Tree Inference
Two ASTRAL (v5.7.8) species trees were inferred using coalescent-based species tree estimation:

    # All orthologs
    java -jar astral.5.7.8.jar \
      -i all_orthologs_trees.tre \
      -o ASTRAL_all_orthologs.tre

    # Filtered orthologs
    java -jar astral.5.7.8.jar \
      -i filtered_orthologs_trees.tre \
      -o ASTRAL_filtered_orthologs.tre

The two trees were compared using RF distance (RAxML-NG v1.2.2), yielding **RF distance = 0** — the two topologies are identical. This confirms that the full ortholog dataset is phylogenetically robust and that filtering does not change the inferred topology.

### 3. Gene Tree Rooting
All 4,845 ortholog gene trees were rooted using RootDigger v1.7.0 prior to CAnDI analysis:

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

### 4. CAnDI Concordance Analysis
CAnDI was run in normal mode (--mode n) on the rooted ortholog gene trees against each species tree. Two runs were performed per species tree:

**With bootstrap cutoff (95)** — nodes below 95% bootstrap support are excluded from analysis:

    python3 CAnDI.py \
      --mode n \
      --species_tree ASTRAL_all_orthologs.tre \
      --gene_folder RootedTrees/ \
      --outfile_prefix CAnDI_all_95cutoff/CAnDI_all \
      --cutoff 95

**Without cutoff** — required to generate _total_analyzed.tre for Pie.py:

    python3 CAnDI.py \
      --mode n \
      --species_tree ASTRAL_all_orthologs.tre \
      --gene_folder RootedTrees/ \
      --outfile_prefix CAnDI_all_nocutoff/CAnDI_all_nocutoff

The same was repeated for the filtered orthologs species tree.

### 5. Pie Chart Visualization
Pie charts showing concordance, conflict, and uninformative proportions per node were generated using CAnDI's Pie.py:

    python3 Pie.py \
      -f CAnDI_all_95cutoff/ \
      -p CAnDI_all \
      -t CAnDI_all_nocutoff/CAnDI_all_nocutoff_total_analyzed.tre \
      -o le \
      -a 3 \
      -at 4841

The -at 4841 flag specifies the total number of orthologs analyzed (the maximum node count from _total_analyzed.tre), ensuring uninformative gene trees are correctly represented in the pie charts. Each Node_X.svg file corresponds to a node in the species tree as labelled in CAnDI_all_labels.tre.

### 6. BES (Branch Estimation Synthesizer) Analysis

#### What is BES?
BES computes concordance-informed branch lengths for a species tree. Unlike conventional branch lengths which reflect substitutions per site, BES branch lengths are derived exclusively from gene trees that are concordant with the species tree at each node. For each internal branch, BES identifies all gene trees that agree with that bipartition and computes summary statistics (mean, median, min, max, 95% CI, and count) of their branch lengths. This provides a powerful way to quantify how much evolutionary time separates successive speciation events, using only the most reliable phylogenetic signal.

#### Why run BES?
For a rapidly radiating group like *Aizoaceae* — the fastest radiating plant family on Earth — BES branch lengths are particularly informative. If speciation events occurred in rapid succession, there would be very little time for gene trees to sort between splits, resulting in near-zero concordance-informed branch lengths at internal nodes. This is a direct molecular signature of rapid radiation.

#### How BES was run
Both ASTRAL species trees were first unrooted using pxrr (Phyx v1.2):

    pxrr -t ASTRAL_all_orthologs.tre -u -o ASTRAL_all_orthologs.ur
    pxrr -t ASTRAL_filtered_orthologs.tre -u -o ASTRAL_filtered_orthologs.ur

BES was then run on both unrooted species trees using Python 2.7:

    # All orthologs
    python2.7 ~/data/BES/src/BES.py \
      -t all_orthologs_trees.tre \
      -m ASTRAL_all_orthologs.ur \
      -o BES_all/BES_all_orthologs

    # Filtered orthologs
    python2.7 ~/data/BES/src/BES.py \
      -t filtered_orthologs_trees.tre \
      -m ASTRAL_filtered_orthologs.ur \
      -o BES_filtered/BES_filtered_orthologs

Each output .tre file contains 7 trees: mean, median, minimum, maximum, lower 95% CI, upper 95% CI, and number of concordant edges per node. The median tree (Tree 2) was used for interpretation, consistent with the approach used for the homolog analysis.

---

## Key Results

| Analysis | Species Tree | RF Distance |
|----------|-------------|-------------|
| All orthologs (n=4,845) | ASTRAL_all_orthologs.tre | — |
| Filtered orthologs (n=1,095) | ASTRAL_filtered_orthologs.tre | 0 |

- Both species trees are **topologically identical** (RF distance = 0)
- The filtered subset recovers the same phylogeny using only 23% of the data
- Pie charts reveal patterns of concordance and conflict across all 37 internal nodes
- BES analysis confirms near-zero concordance-informed branch lengths at internal nodes across homologs, all orthologs, and filtered orthologs, providing molecular evidence of rapid radiation in *Aizoaceae*

---

## Tools and Versions

| Tool | Version | Purpose |
|------|---------|---------|
| ASTRAL | 5.7.8 | Coalescent species tree inference |
| RootDigger | 1.7.0 | Gene tree rooting |
| CAnDI | — | Concordance/conflict analysis |
| RAxML-NG | 1.2.2 | RF distance computation |
| Pythia | 2.0.0 | Phylogenetic difficulty prediction |
| BES | — | Concordance-informed branch lengths |
| Phyx | 1.2 | Tree manipulation |

---

## Citations

> Haag, J. & Stamatakis, A. (2025). Pythia 2.0: New Data, New Prediction Model, New Features. bioRxiv. https://doi.org/10.1101/2025.03.25.645182

> Kozlov, A.M. et al. (2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. Bioinformatics, 35(21):4453-4455.

> Zhang, C. et al. (2018). ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees. BMC Bioinformatics, 19(S6):153.

> Robertson, H.M., Walker, J.F. & Moyroud, E. (2025). CAnDI: a new tool to investigate conflict in homologous gene trees and explain convergent trait evolution. Systematic Biology, syaf028. https://doi.org/10.1093/sysbio/syaf028

> Yang, Y. & Smith, S.A. (2014). Orthology inference in nonmodel organisms using transcriptomes and low-coverage genomes. Molecular Biology and Evolution, 31(11):3081-3092.

> Morales-Briones, D.F. et al. (2021). Disentangling sources of gene tree discordance in phylogenomic data sets. Systematic Biology, 70(3):591-609.

---

## Repository Structure

    CANDIforAIM1Orthologs/
    |
    |-- ASTRAL_all_orthologs.tre          # ASTRAL species tree from all 4,845 orthologs
    |-- ASTRAL_all_orthologs.ur           # Unrooted version for BES
    |-- ASTRAL_filtered_orthologs.tre     # ASTRAL species tree from 1,095 filtered orthologs
    |-- ASTRAL_filtered_orthologs.ur      # Unrooted version for BES
    |-- Correlation_TreeStats.csv         # TC, RTC, TCA, RTCA, Taxa, Pythia scores per ortholog
    |-- filtered_orthologs.csv            # List of orthologs passing the quality filter
    |
    |-- RootedTrees/                      # 4,845 ortholog gene trees rooted with RootDigger v1.7.0
    |-- RootDigger_checkpoints/           # RootDigger checkpoint files
    |
    |-- CAnDI_all_95cutoff/               # CAnDI results (all orthologs, bootstrap cutoff = 95)
    |-- CAnDI_filtered_95cutoff/          # CAnDI results (filtered orthologs, bootstrap cutoff = 95)
    |-- CAnDI_all_nocutoff/               # CAnDI results (all orthologs, no cutoff)
    |-- CAnDI_filtered_nocutoff/          # CAnDI results (filtered orthologs, no cutoff)
    |
    |-- PieCharts_all/                    # Pie chart SVGs for all orthologs (Node_0 to Node_36)
    |-- PieCharts_filtered/               # Pie chart SVGs for filtered orthologs (Node_0 to Node_36)
    |
    |-- BES_all/                          # BES results for all orthologs
    |-- BES_filtered/                     # BES results for filtered orthologs
    |
    |-- README.md

---

## Author

**Tomi Jacobs** — PhD Candidate, Computational Biology / Phylogenomics
