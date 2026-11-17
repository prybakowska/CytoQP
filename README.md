# CytoQP

**Cyto Quality Pipeline**

## Introduction

CytoQP is an R package designed for preprocessing, quality control, normalization, and analysis of high-dimensional cytometry data, with a particular focus on mass cytometry (CyTOF).

The package provides tools for reproducible data processing in large-scale, multi-batch, and multicenter studies. It is particularly suited for whole-blood samples acquired in multiple aliquots to minimize prolonged exposure to acquisition solutions.

## Main Features

- **Bead-based normalization:** Normalization of signal intensities using calibration beads and reference samples.
- **Flow-rate and signal cleaning:** Detection and removal of acquisition anomalies using flowAI-adapted methods, flowCut, and PeacoQC.
- **File quality assessment:** Identification of low-quality samples and acquisition outliers.
- **Sample debarcoding:** Separation of multiplexed samples using barcode information.
- **File deconvolution and aggregation:** Reconstruction and merging of data acquired in multiple aliquots.
- **Cell gating:** Identification of intact cells, singlets, CD45-positive cells, and viable populations, including bead exclusion.
- **Reference-based normalization:** Correction of batch-associated technical variation.
- **FlowSOM-based feature extraction:** Quantification of cluster and metacluster frequencies and marker signal intensities (MSI).
- **Data visualization:** Diagnostic plots for preprocessing, batch effects, marker distributions, and quality control.

## Required R Packages

CytoQP requires several dependencies available through CRAN, Bioconductor, and GitHub.

### CRAN Packages

- dplyr
- ggplot2
- ggpubr
- gridExtra
- magrittr
- pheatmap
- RColorBrewer
- stringr
- uwot

### Bioconductor Packages

- BiocParallel
- CATALYST
- flowCore
- flowCut
- flowDensity
- FlowSOM
- PeacoQC
- S4Vectors

### GitHub Packages

- CytoNorm — `saeyslab/CytoNorm`
- cytutils — `ismmshimc/cytutils`

The package also uses standard R packages, including `graphics`, `grDevices`, `grid`, `methods`, `stats`, and `utils`. These are distributed with R and do not require separate installation.

The flow-rate quality control implementation also includes functions adapted from **flowAI**. The `flowAI` package is acknowledged as a methodological source but is not required as a direct runtime dependency for the adapted implementation.

## Installation

### 1. Install CytoQP

```r
devtools::install_github("prybakowska/CytoQP")
```

### 2. Load CytoQP

```r
library(CytoQP)
```

Dependencies declared in the package `DESCRIPTION` file are also resolved automatically during installation where supported.

## Applications

CytoQP is particularly suitable for studies involving multiple acquisition batches, different cytometry instruments or centers, and large numbers of blood samples.

Its modular design enables users to perform preprocessing steps independently or combine them into comprehensive workflows while generating diagnostic outputs to evaluate data quality and processing performance.
