# 02.3wpi RNA-seq: Medicago truncatula Root Nodule Transcriptome

Differential gene expression analysis of *Medicago truncatula* root nodules at 3 weeks post-inoculation (3 wpi), comparing 6 mutant/genotype lines against wild-type (WT). The study investigates the transcriptional basis of symbiotic nitrogen fixation (SNF) using combinatorial mutants of RSP (Root-Specific Peroxidase) genes and the Fix⁻ control mutant *dnf2*.

## Table of Contents

- [Experimental Design](#experimental-design)
- [Pipeline Overview](#pipeline-overview)
- [1. Quality Control](#1-quality-control)
- [2. Differential Expression Results](#2-differential-expression-results)
- [3. Cross-Genotype Comparisons](#3-cross-genotype-comparisons)
  - [3.1 UpSet Plots](#31-upset-plots)
  - [3.2 Combinatorial-Triple Venn Diagrams](#32-combinatorial-triple-venn-diagrams)
  - [3.3 Pairwise Venn Diagrams](#33-pairwise-venn-diagrams)
  - [3.4 A86485 vs dnf2 — Shared Symbiotic-Failure Signature](#34-a86485-vs-dnf2--shared-symbiotic-failure-signature)
- [4. Summary Heatmaps](#4-summary-heatmaps)
- [5. Curated Nodulation-Gene Panel (Data11)](#5-curated-nodulation-gene-panel-data11)
- [6. GO/KEGG Enrichment Analysis](#6-go-kegg-enrichment-analysis)
- [7. Key Biological Findings](#7-key-biological-findings)
- [Directory Structure](#directory-structure)
- [Dependencies](#dependencies)
- [Usage](#usage)

---

## Experimental Design

| Genotype | Description | Replicates |
|:---------|:------------|:----------:|
| **WT** | Wild-type (A17) | 3 |
| **A86** | *rsp86* single mutant | 3 |
| **A256** | *rsp256* single mutant | 3 |
| **A485** | *rsp485* single mutant | 3 |
| **A86256** | *rsp86* × *rsp256* double mutant | 3 |
| **A86485** | *rsp86* × *rsp485* double mutant | 3 |
| **dnf2** | *dnf2* Fix⁻ control mutant (loss-of-function) | 3 |

**Total:** 21 samples (7 genotypes × 3 biological replicates)

The combinatorial mutant design (A86256 = A86 + A256; A86485 = A86 + A485) enables dissection of epistatic and additive effects between RSP loci. The *dnf2* mutant serves as a known Fix⁻ (nitrogen fixation-deficient) control.

**Reference genome:** *M. truncatula* v5 (MtrunA17r5.0)

---

## Pipeline Overview

```
Raw reads (FASTQ)
    │
    ▼
FastP trimming (02.CleanData/)
    │
    ▼
STAR alignment ──► BAM/BAI (21 samples)
    │
    ├──► htseq-count ──► raw gene counts (49,421 genes)
    ├──► StringTie -A ──► FPKM / TPM per gene
    └──► StringTie -B ──► Ballgown ctab files
            │
            ▼
    Master expression tables + count/TPM/FPKM matrices
            │
            ▼
    Unified edgeR model (all 21 samples, ~0+group design)
    ──► glmTreat (lfc=1) ──► 6 genotype-vs-WT contrasts
            │
            ▼
    QC plots, per-comparison DE tables, volcano/MA/heatmaps
            │
            ▼
    Cross-genotype summaries: UpSet, Venn, enrichment, panel heatmaps
```

**Design choices:**
- **Unified fit** (not 6 independent pairwise fits): one model across all 21 samples gives 14 residual df for stable dispersion estimation
- **glmTreat** (lfc=1): tests significance *against* the fold-change threshold directly (more rigorous than post-hoc FC filtering)
- **filterByExpr**: 25,950 of 49,421 genes retained (52.5%) — the common gene universe for all comparisons
- **TMM normalization** for library-size correction

---

## 1. Quality Control

Global QC across all 21 samples confirms tight replicate clustering and clear genotype separation.

| | |
|:---:|:---:|
| ![MDS](figures/QC/01.MDS_all_samples.png) | ![PCA](figures/QC/02.PCA_logCPM.png) |
| **MDS plot** — multidimensional scaling of logCPM distances | **PCA** — principal component analysis of logCPM expression |
| ![Correlation](figures/QC/03.Sample_correlation.png) | ![Library sizes](figures/QC/04.Library_sizes.png) |
| **Sample correlation** (Spearman, logCPM) | **Library sizes** per sample |

**Key observations:**
- Replicates cluster tightly within each genotype
- WT and *dnf2* separate cleanly on PC1 (the largest source of variance)
- A86485 replicates show slightly more spread along PC1 (rep 3 sits apart from 1/2)
- A86256 clusters close to A256, and A86485 clusters close to *dnf2* — consistent with the Venn/heatmap results below

---

## 2. Differential Expression Results

Each genotype was compared to WT using the unified edgeR model with glmTreat (FDR < 0.05, |log₂FC| > 1).

### DEG Counts

| Genotype | UP | DOWN | Total |
|:---------|---:|-----:|------:|
| A86 | 438 | 197 | 635 |
| A256 | 1,162 | 959 | 2,121 |
| A485 | 1,292 | 1,015 | 2,307 |
| A86256 | 1,305 | 1,114 | 2,419 |
| A86485 | 279 | 1,122 | 1,401 |
| dnf2 | 1,017 | 1,685 | 2,702 |

**Internal positive control:** MtDNF2 (MtrunA17Chr4g0044681) is the #1 most significant gene in *dnf2*_vs_WT, DOWN (FDR 1.2×10⁻¹⁴) — the loss-of-function mutant's own gene is knocked down, validating sample labeling, contrast direction, and statistics.

### Volcano Plots

| | |
|:---:|:---:|
| ![A256 volcano](figures/DE_per_genotype/A256_vs_WT_Volcano.png) | ![A485 volcano](figures/DE_per_genotype/A485_vs_WT_Volcano.png) |
| **A256 vs WT** | **A485 vs WT** |
| ![A86 volcano](figures/DE_per_genotype/A86_vs_WT_Volcano.png) | ![A86256 volcano](figures/DE_per_genotype/A86256_vs_WT_Volcano.png) |
| **A86 vs WT** | **A86256 vs WT** |
| ![A86485 volcano](figures/DE_per_genotype/A86485_vs_WT_Volcano.png) | ![dnf2 volcano](figures/DE_per_genotype/dnf2_vs_WT_Volcano.png) |
| **A86485 vs WT** | **dnf2 vs WT** |

---

## 3. Cross-Genotype Comparisons

### 3.1 UpSet Plots

UpSet plots show every non-empty intersection across all 6 genotypes simultaneously (the right tool for >3 sets).

| UP DEGs | DOWN DEGs |
|:---:|:---:|
| ![UpSet UP](figures/Summary/UpSet_UP.png) | ![UpSet DOWN](figures/Summary/UpSet_DOWN.png) |

**Key pattern:** *dnf2*-only is the largest intersection (602 UP, 931 DOWN), consistent with its distinct Fix⁻ transcriptome. The A86485 + *dnf2* shared-DOWN intersection (954 genes) is striking — see [Section 3.4](#34-a86485-vs-dnf2--shared-symbiotic-failure-signature).

### 3.2 Combinatorial-Triple Venn Diagrams

Each double mutant is compared against its two single-mutant parents, revealing epistatic and additive effects.

**A86 / A256 / A86256 (double = A86 + A256):**

| UP | DOWN |
|:---:|:---:|
| ![Venn A86-A256-A86256 UP](figures/Summary/Venn_A86_A256_A86256_UP.png) | ![Venn A86-A256-A86256 DOWN](figures/Summary/Venn_A86_A256_A86256_DOWN.png) |

**A86 / A485 / A86485 (double = A86 + A485):**

| UP | DOWN |
|:---:|:---:|
| ![Venn A86-A485-A86485 UP](figures/Summary/Venn_A86_A485_A86485_UP.png) | ![Venn A86-A485-A86485 DOWN](figures/Summary/Venn_A86_A485_A86485_DOWN.png) |

**Interpretation:**
- A86256 ≈ A256 (redundancy): A86 alone adds few unique genes; A256 and A86256 share heavily (511 UP shared)
- A485 dominates A86485: 914 A485-only UP vs only 258 A86485-only, with A86485 showing a strong DOWN shift (1,001 DOWN-only)
- A86485 behaves more like the Fix⁻ *dnf2* control than like its A485 parent

### 3.3 Pairwise Venn Diagrams

Focused 2-way overlaps between key genotype pairs:

| A86 vs A485 | A86 vs dnf2 |
|:---:|:---:|
| ![Venn A86-A485 UP](figures/Summary/Venn_A86_vs_A485_UP.png) | ![Venn A86-dnf2 UP](figures/Summary/Venn_A86_vs_dnf2_UP.png) |
| ![Venn A86-A485 DOWN](figures/Summary/Venn_A86_vs_A485_DOWN.png) | ![Venn A86-dnf2 DOWN](figures/Summary/Venn_A86_vs_dnf2_DOWN.png) |

| A485 vs dnf2 |
|:---:|
| ![Venn A485-dnf2 UP](figures/Summary/Venn_A485_vs_dnf2_UP.png) |
| ![Venn A485-dnf2 DOWN](figures/Summary/Venn_A485_vs_dnf2_DOWN.png) |

### 3.4 A86485 vs dnf2 — Shared Symbiotic-Failure Signature

The PCA and data-driven heatmap showed A86485 clustering with *dnf2*. This Venn quantifies the overlap:

| UP | DOWN |
|:---:|:---:|
| ![Venn A86485-dnf2 UP](figures/Summary/Venn_A86485_vs_dnf2_UP.png) | ![Venn A86485-dnf2 DOWN](figures/Summary/Venn_A86485_vs_dnf2_DOWN.png) |

**DOWN overlap is the signature:** 954 shared genes = 85% of all A86485 DOWN-DEGs are also down in *dnf2* (vs only 168 A86485-specific). Discordant genes are negligible (2 + 3). The 954 shared-DOWN set was used for GO/KEGG enrichment (Section 6).

---

## 4. Summary Heatmaps

### Data-Driven Heatmap (Top DEGs)

Union of the top 50 DEGs per genotype (213 genes total), clustered by expression pattern:

![Summary heatmap top DEGs](figures/Summary/Summary_Heatmap_topDEGs.png)

**Pattern:** A86485 + *dnf2* cluster together (right), distinct from the A86/A256/A485/A86256 group (left). A large block of genes is regulated in opposite directions between these two clusters.

### DEG Count Summary

![DEG counts barplot](figures/Summary/DEG_counts_barplot.png)

---

## 5. Curated Nodulation-Gene Panel (Data11)

A curated panel of 41 nodulation/symbiosis genes (NCRs, leghemoglobins, Nod-factor receptors, defense/PR, regulatory genes), ported from the J002 analysis. Gene panel list: [`Data11_panel_genes.tsv`](scripts/Data11_panel_genes.tsv)

### Significance-Masked (faithful reproduction of original)

Only log₂FC values where the gene is a significant DEG in that genotype are shown; grey = not significant.

![Data11 panel sigmasked](figures/Summary/Data11_panel_heatmap_sigmasked.png)

### All log₂FC (shows underlying trends)

All log₂FC values colored regardless of significance:

![Data11 panel alllogFC](figures/Summary/Data11_panel_heatmap_alllogFC.png)

**Key patterns:**
- A86 is almost entirely grey (only 2 significant panel genes) — matching the original J002 rsp86 column
- Leghemoglobins (Lb1, Lb2, Lb5, Lb6, Lb7, Lb8, Lb10) go progressively down in higher-order mutants, strongest in *dnf2*
- PR/defense genes (PR10, MtNFS1) go strongly up — classic defense-up/symbiosis-down pattern
- Significance increases through the combinatorials to *dnf2* (16 DOWN panel genes)

---

## 6. GO/KEGG Enrichment Analysis

Over-representation analysis (ORA) of the 954 shared-DOWN genes (A86485 ∩ *dnf2*) using clusterProfiler. Background = 25,950 tested genes (correct universe).

### GO Enrichment

| Biological Process | Molecular Function |
|:---:|:---:|
| ![GO dotplot BP](figures/Enrichment/GO_dotplot_BP.png) | ![GO dotplot MF](figures/Enrichment/GO_dotplot_MF.png) |

| Combined GO Barplot | GO Network (BP) |
|:---:|:---:|
| ![GO barplot](figures/Enrichment/GO_barplot_combined.png) | ![GO cnet](figures/Enrichment/GO_cnet_BP.png) |

### KEGG Enrichment

| KEGG Dotplot | KEGG Network |
|:---:|:---:|
| ![KEGG dotplot](figures/Enrichment/KEGG_dotplot.png) | ![KEGG cnet](figures/Enrichment/KEGG_cnet.png) |

### Key Enriched Terms

| Category | Term | FDR | Key Genes |
|:---------|:-----|----:|:----------|
| GO-BP | Oxygen transport | 1.5×10⁻⁶ | MtLb1, MtLb2, MtLb5, MtLb6, MtLb7, MtLb8, MtLb10 |
| GO-MF | Oxygen carrier activity | 4.7×10⁻⁵ | (same leghemoglobins) |
| GO-BP | Gibberellin biosynthesis | 0.012 | MtGA20ox5 |
| GO-MF | Iron ion binding (trending) | 0.06 | — |
| KEGG | Nitrogen metabolism | 4.5×10⁻⁴ | MtGS1a, MtNRT2.1 |
| KEGG | Ala/Asp/Glu metabolism | 0.014 | — |

**Annotation coverage caveat:** Of the 954 shared-DOWN genes, 622 (65%) have GO and only 125 (13%) have a KEGG KO. Nodule-specific peptides (NCRs) are largely unannotated — ORA captures the metabolic/O₂/defense arms, while the symbiotic-peptide arm is annotation-invisible.

**Internal verification:** Top GO term p-value independently recomputed by hand (phyper) = 8.4257×10⁻⁹ = clusterProfiler's reported p-value (exact match). All member genes traced back to the master DE table.

---

## 7. Key Biological Findings

1. **A86485 phenocopies *dnf2***: 85% of A86485's DOWN-DEGs are shared with *dnf2*; the PCA and data-driven heatmap confirm they cluster together. The double mutant *rsp86*×*rsp485* triggers a symbiotic-failure response nearly identical to the Fix⁻ control.

2. **Nitrogen-fixation collapse signature**: The shared A86485∩*dnf2* DOWN set is dominated by 7 leghemoglobins (O₂-buffering for nitrogenase) + nitrogen metabolism genes (MtGS1a, MtNRT2.1) — the textbook molecular signature of Fix⁻ nodules.

3. **A485 dominates A86485**: In the A86/A485/A86485 Venn, A485 contributes the vast majority of unique DEGs (914 UP-only vs 258 A86485-only), with A86485 showing a strong DOWN shift — the "485" lesion is the primary effector.

4. **A86256 ≈ A256 (redundancy)**: A86 alone has the fewest DEGs (635); A256 and A86256 share heavily, suggesting A86 contributes minimally to the A86256 double-mutant phenotype.

5. **Defense-symbiosis trade-off**: PR/defense genes (PR10, MtNFS1) are strongly UP in the Fix⁻-like genotypes while leghemoglobins and NCRs are DOWN — consistent with the known antagonism between defense and symbiotic programs.

---

## Directory Structure

```
02.3wpi.RNAseq/
├── README.md                          # This file
├── scripts/                           # All analysis scripts
│   ├── config.sh                      # Environment configuration
│   ├── 00.Mapping.srun                # STAR alignment template
│   ├── 01-07.*.srun                   # Per-genotype alignment scripts
│   ├── 00.QC.srun                     # Alignment QC (qualimap)
│   ├── 05.build_expression_tables.R   # Master table + matrix builder
│   ├── 00.edgeR_unified.R             # Unified edgeR DE analysis
│   ├── 01.run_DE.sh                   # Runner script
│   ├── 02.summary_upset_heatmaps.R    # UpSet + summary heatmaps
│   ├── 03.venn_diagrams.py            # Combinatorial-triple Venns
│   ├── 04.venn_A86485_vs_dnf2.py      # A86485 vs dnf2 focused Venn
│   ├── 05.Data11_panel_heatmap.R      # Curated nodulation-gene panel
│   ├── 06.venn_pairwise.py            # Pairwise Venns
│   ├── 06.GO_KEGG_enrichment.R        # GO/KEGG ORA
│   ├── 07.verify_enrichment.R         # Independent verification
│   └── Data11_panel_genes.tsv         # Curated gene panel (41 genes)
├── figures/
│   ├── QC/                            # Quality control plots
│   ├── DE_per_genotype/               # Volcano + heatmap per genotype
│   ├── Summary/                       # Cross-genotype summaries
│   └── Enrichment/                    # GO/KEGG enrichment plots
└── tables/
    ├── MasterTable_with_DE.tsv        # 25,950 genes × all genotypes
    ├── All_comparisons_DE_summary.tsv # DEG counts per comparison
    ├── Data11_panel_heatmap_data.tsv  # Panel heatmap underlying data
    ├── GeneList_*.tsv                 # Venn region gene lists
    ├── GO_enrichment_all.tsv          # Full GO enrichment results
    ├── KEGG_enrichment.tsv            # Full KEGG enrichment results
    └── ...
```

---

## Dependencies

**R (v4.4.3):** edgeR, limma, statmod, data.table, ggplot2, pheatmap, RColorBrewer, ggrepel, UpSetR, clusterProfiler, GO.db, KEGGREST, AnnotationDbi, enrichplot, DOSE, grid

**Python 3:** matplotlib, matplotlib_venn, csv

**CLI tools:** STAR (≥2.7), samtools, htseq-count, stringtie (≥3.0), qualimap

---

## Usage

1. **Configure environment:** Edit `scripts/config.sh` to set paths for your system, then `source scripts/config.sh`

2. **Alignment:** Run per-genotype `.srun` scripts (STAR + htseq-count + StringTie):
   ```bash
   cd scripts
   bash 01.A256.srun 2>&1 | tee 01.A256.log
   # ... repeat for each genotype
   ```

3. **Build expression tables:**
   ```bash
   Rscript scripts/05.build_expression_tables.R
   ```

4. **DE analysis:**
   ```bash
   bash scripts/01.run_DE.sh    # runs 00.edgeR_unified.R for all 6 contrasts
   ```

5. **Summary visualizations:**
   ```bash
   Rscript scripts/02.summary_upset_heatmaps.R
   python3 scripts/03.venn_diagrams.py
   python3 scripts/04.venn_A86485_vs_dnf2.py
   python3 scripts/06.venn_pairwise.py
   Rscript scripts/05.Data11_panel_heatmap.R
   ```

6. **Enrichment:**
   ```bash
   Rscript scripts/06.GO_KEGG_enrichment.R
   Rscript scripts/07.verify_enrichment.R   # independent verification
   ```
