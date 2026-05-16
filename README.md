# HN00273522 — Bulk RNA-seq Analysis of BMDM × Intestinal Organoid Co-culture under TNFα Challenge

[![DOI](https://zenodo.org/badge/DOI/XXXXXXX.svg)](https://doi.org/XXXXXXX)
[![R](https://img.shields.io/badge/R-%E2%89%A54.3.0-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-green)]()

R-based bulk RNA-seq analysis pipeline for a bone marrow–derived macrophage (BMDM) × intestinal organoid co-culture system under TNFα challenge. The pipeline performs differential expression analysis, multi-platform validation, functional enrichment, co-expression network construction, and curated gene-panel visualization across 24 samples organized into eight biological groups (n = 3 per group).

---

## Overview

The experimental design comprises eight conditions:

| Group | TNFα | BMDM | Description |
| --- | :---: | :---: | --- |
| `Organoid` | − | − | Intestinal organoid monoculture |
| `Organoid_TNFa` | + | − | Organoid + TNFα |
| `Organoid_BMDM_1K` | − | 1K | Organoid + 1,000 BMDM |
| `Organoid_BMDM_5K` | − | 5K | Organoid + 5,000 BMDM |
| `Organoid_TNFa_BMDM_1K` | + | 1K | Organoid + TNFα + 1,000 BMDM |
| `Organoid_TNFa_BMDM_5K` | + | 5K | Organoid + TNFα + 5,000 BMDM |
| `Macrophage` | − | − | BMDM monoculture |
| `Macrophage_TNFa` | + | − | BMDM + TNFα |

Pairwise comparisons are organized into a **6-tier, 14-contrast framework** to dissect basal co-culture effects, TNFα-induced inflammation, and macrophage–organoid hybrid responses against distinct biological baselines.

| Tier | Baseline | # Contrasts | Biological question |
| --- | --- | :---: | --- |
| Tier 1 | Organoid | 2 | Basal effect of BMDM addition |
| Tier 2 | Macrophage | 2 | Basal effect of organoid co-culture on BMDM |
| Tier 3 | Organoid | 3 | Absolute TNFα-induced changes (Organoid reference) |
| Tier 4 | Organoid + TNFα | 2 | BMDM modulation of TNFα-stimulated organoid |
| Tier 5 | Macrophage | 3 | TNFα effect on BMDM ± organoid |
| Tier 6 | Macrophage + TNFα | 2 | Organoid modulation of TNFα-stimulated BMDM |

---

## Repository structure

```
.
├── README.md
├── LICENSE
├── renv.lock                         # Pinned package versions
├── sessionInfo.txt                   # Full R session record
│
├── scripts/
│   ├── 01_inspect_count_matrix.R     # Raw matrix structure & QC
│   ├── 02_DEG_master_pipeline.R      # DESeq2 + limma + edgeR + WGCNA + GSEA export
│   ├── 03_inflammation_bubble.R      # Inflammation gene-panel visualization
│   ├── 04_inflammation_bubble_noTNF.R
│   ├── 04_inflammation_bubble_TNF.R
│   ├── 05_stemness_bubble_TNF.R
│   ├── 06_stemness_bubble_noTNF.R
│   ├── 06_1_stemness_layoutA.R       # Curated GO BP bubble layout — basal
│   └── 06_2_stemness_layoutB.R       # Curated GO BP bubble layout — full
│
├── data/
│   └── README.md                     # Pointer to deposited raw FASTQ + count matrix
│
├── output/                           # Auto-generated, .gitignored
│   ├── QC/
│   ├── DEG_results/
│   ├── Volcano/
│   ├── GO_KEGG/
│   ├── Pathview/
│   ├── gProfiler2/
│   ├── limma/
│   ├── edgeR/
│   ├── WGCNA/
│   ├── GSEA_export/
│   └── tables/
│
└── docs/
    └── pipeline_overview.png         # Workflow diagram
```

---

## Requirements

### System

- R ≥ 4.3.0
- Operating system: Windows 10/11, macOS ≥ 12, or Linux (Ubuntu ≥ 20.04)
- RAM: ≥ 16 GB recommended (WGCNA step is memory-intensive)
- Internet connection required for KEGG / Reactome / g:Profiler API queries

### R packages

All package versions are pinned in `renv.lock` and can be restored in one step:

```r
install.packages("renv")
renv::restore()
```

**Core analysis (Bioconductor):** `DESeq2`, `limma`, `edgeR`, `ashr`, `apeglm`, `clusterProfiler`, `enrichplot`, `DOSE`, `pathview`, `gage`, `gageData`, `GSVA`, `EnhancedVolcano`, `WGCNA`, `topGO`, `KEGGREST`

**Annotation:** `org.Mm.eg.db`, `EnsDb.Mmusculus.v79`, `AnnotationDbi`

**Visualization & utilities:** `ggplot2`, `ggrepel`, `pheatmap`, `RColorBrewer`, `patchwork`, `corrplot`, `ggfortify`, `ggforce`, `gprofiler2`, `tidyverse`, `matrixStats`, `Hmisc`, `dynamicTreeCut`, `fastcluster`

Exact versions used in the published analysis are recorded in `sessionInfo.txt`.

---

## Upstream pipeline (sequencing → counts)

Library preparation and sequencing were performed by **Macrogen Inc. (Seoul, Republic of Korea)** using the **TruSeq Stranded Total RNA with Ribo-Zero H/M/R_Gold** kit, paired-end **101 bp × 2** on the Illumina platform.

Read processing was performed by the sequencing provider as follows:

```
FASTQ (raw)
   └── Trimmomatic (adapter & quality trimming)
        └── HISAT2 (alignment to GRCm38, strand-specific)
             └── StringTie (transcript assembly + quantification)
                  └── prepDE.py → gene_count_matrix.csv
                                  transcript_count_matrix.csv
```

This repository takes `gene_count_matrix.csv` as its primary input.

---

## Downstream analysis (this repository)

```
gene_count_matrix.csv
   │
   ├── 01: Matrix structure & integrity QC
   │
   ├── 02: Integrated master pipeline
   │   ├── DESeqDataSet construction (~ group, 8 levels)
   │   ├── Pre-filter (≥10 counts in ≥3 samples)
   │   ├── DESeq2 normalization + VST + PCA + sample-distance heatmap
   │   ├── 14 contrasts × ashr LFC shrinkage
   │   ├── EnhancedVolcano (×14)
   │   ├── clusterProfiler GO (BP/CC/MF) + KEGG (×14)
   │   ├── pathview KEGG overlays (selected contrasts)
   │   ├── gprofiler2 cross-validation (g:SCS)
   │   ├── limma + voom (TMM) cross-validation
   │   ├── edgeR + glmQLFit + glmTreat cross-validation
   │   ├── WGCNA co-expression network (signed, R² ≥ 0.80)
   │   └── GSEA-formatted exports (.gct + .cls)
   │
   ├── 03–04: Inflammation gene-panel visualization
   │   ├── Sub-panel LFC heatmaps per tier
   │   ├── Sub-panel module score (mean LFC + one-sample t-test)
   │   ├── KEGG ORA lollipop plots (Inflammation vs Apoptosis/Resolution)
   │   ├── Filtered KEGG dotplots (pathogen/disease blacklist)
   │   ├── GSEA enrichment plots for inflammation pathways
   │   └── GSVA sample-level scoring (supplementary)
   │
   └── 05–06: Stemness gene-panel visualization
       ├── ISC + lineage panel heatmaps (8 sub-panels)
       ├── Curated GO BP bubble plots (Layout A: basal 4 contrasts; Layout B: full 14)
       ├── KEGG / GO MF / GO CC bubble plots per baseline
       └── 6-way Venn intersection of stem-cell DEGs
```

---

## Usage

### 1. Restore environment

```r
renv::restore()
```

### 2. Place input data

Copy or symlink `gene_count_matrix.csv` into the project root, or modify the `file_path` variable at the top of `02_DEG_master_pipeline.R`.

### 3. Run the master pipeline

```r
source("scripts/02_DEG_master_pipeline.R")
```

This generates the full `output/` directory (~3–5 GB) and writes the following key RDS objects to the working directory:

- `dds_full_8groups.rds`
- `vsd_full_8groups.rds`
- `deg_results_full_8groups.rds`
- `gene_annotation.rds`
- `contrast_tiers.rds`

### 4. Run downstream visualization scripts

These depend on the RDS objects from step 3 and can be run in any order:

```r
source("scripts/03_inflammation_bubble.R")
source("scripts/05_stemness_bubble_TNF.R")
source("scripts/06_stemness_bubble_noTNF.R")
```

### Expected runtime

| Step | Approximate time (16 GB RAM, 8-core CPU) |
| --- | --- |
| 01 (matrix QC) | < 1 min |
| 02 (master pipeline, all 14 contrasts) | 45–90 min |
| 02 — WGCNA only | 5–15 min |
| 03–06 (visualization) | 5–15 min each |

---

## Key parameters

All thresholds are defined as variables at the top of `02_DEG_master_pipeline.R` and inherited by downstream scripts:

| Parameter | Value | Description |
| --- | :---: | --- |
| `padj_cutoff` | 0.05 | BH-adjusted p-value cutoff for DEG |
| `lfc_cutoff` | 1.0 | Minimum \|log₂FC\| for DEG |
| `min_dds_count` | 10 | Per-gene count threshold for pre-filter |
| `min_dds_count_min_rowSums` | 3 | Minimum samples meeting count threshold |
| `wgcna_seed` | 6689 | Random seed for `blockwiseModules()` |
| WGCNA scale-free R² | 0.80 | Minimum R² for soft-threshold selection |
| WGCNA `mergeCutHeight` | 0.25 | Module-merging threshold |
| GSEA `minGSSize` | 10 | Minimum gene-set size |
| GSEA `maxGSSize` | 500 | Maximum gene-set size |

---

## Data availability

- **Raw sequencing reads (FASTQ):** deposited at ArrayExpress under accession `E-MTAB-XXXXX`, with the corresponding ENA Study accession `PRJEB-XXXXX`.
- **Gene-level count matrix:** distributed with the ArrayExpress submission as processed data, and additionally archived at Zenodo (DOI: `10.5281/zenodo.XXXXXXX`).
- **Curated gene panels** (Inflammation: 9 sub-panels; Stemness / Lineage: 8 sub-panels) are defined in-script and exported as Supplementary Tables in the associated manuscript.

---

## Citation

If you use this pipeline or its outputs, please cite:

> Lee T.B., *et al.* (2026). [Manuscript title]. *[Journal]*, [volume]:[pages]. DOI: [xxxxxxx]

And the underlying tools — at minimum:

- Love M.I., Huber W., Anders S. (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* **Genome Biology** 15:550.
- Ritchie M.E., *et al.* (2015). *limma powers differential expression analyses for RNA-sequencing and microarray studies.* **Nucleic Acids Research** 43:e47.
- Robinson M.D., McCarthy D.J., Smyth G.K. (2010). *edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.* **Bioinformatics** 26:139–140.
- Wu T., *et al.* (2021). *clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.* **The Innovation** 2:100141.
- Langfelder P., Horvath S. (2008). *WGCNA: an R package for weighted correlation network analysis.* **BMC Bioinformatics** 9:559.

---

## Reproducibility checklist

- [x] All package versions pinned via `renv.lock`
- [x] Full `sessionInfo()` archived
- [x] Random seeds fixed for WGCNA and GSEA
- [x] Pre-filtering and DEG cutoffs declared as named variables
- [x] Raw counts and metadata publicly deposited
- [x] Pipeline runnable end-to-end from `gene_count_matrix.csv`

---

## License

Released under the MIT License — see `LICENSE`.

---

## Contact

**Tae Baek Lee (이태백)**
Laboratory of Veterinary Physiology, College of Veterinary Medicine, Jeju National University
- Email: [your.email@jejunu.ac.kr]
- ORCID: [0000-0000-0000-0000]
- GitHub: [@yourhandle]

**Principal Investigator:** Prof. Changhwan Ahn
