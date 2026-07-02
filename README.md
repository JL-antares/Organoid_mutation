# Supporting code for organoid similarity analysis

This repository contains the supporting R code used for organoid similarity analysis based on somatic mutation profiles. The current workflow focuses on mutation landscape visualization, sample-level mutation sharing, mutational spectrum/signature analysis, transition/transversion profiling, and optional copy-number heatmap visualization.

The code and project structure are organized for research reproducibility and GitHub deposition. Raw sequencing files and large intermediate results are not included in this repository.

## Introduction

### Background

Patient-derived organoids are commonly compared with matched tumor or tissue samples to evaluate whether the in vitro models preserve key genomic characteristics of the original lesion. Somatic mutation and copy-number profiles provide a direct way to assess this similarity across samples.

### Analysis scope

The workflow in this repository summarizes mutation profiles from MAF/VCF-style inputs and compares mutation overlap among samples. The current script includes the following analysis modules:

* recurrently mutated gene visualization with `maftools::oncoplot`
* variant allele frequency (VAF) extraction and plotting
* shared/private mutation site analysis with UpSet and Venn plots
* 96-channel mutational spectrum and NMF-based signature extraction
* COSMIC signature comparison and signature contribution visualization
* transition/transversion (Ti/Tv) summary
* optional copy-number heatmap for selected genes

### Expected outputs

The pipeline generates publication-oriented PDF/PNG figures and Excel tables, including oncoplots, mutation-site intersection plots, mutational spectrum profiles, mutational signature contributions, and Ti/Tv summaries.

## Contents of this repository

The repository is organized in a numbered layout similar to a research-supporting code repository:

```text
project
├── 00.Processed Data
├── 01.Organoid Similarity Analysis
├── RESULT
├── README.md
├── .gitattributes
└── .gitignore
```

* `00.Processed Data` stores example input notes and small metadata templates. Large raw/intermediate genomic files should not be committed.
* `01.Organoid Similarity Analysis` stores the main R analysis script for organoid mutation similarity analysis.
* `RESULT` stores generated figures and tables. The folder is kept in the repository, while generated result files are ignored by default.

## Input files

The main script expects user-provided mutation and metadata files. Based on the current code, the important inputs are:

| File or folder | Description |
| --- | --- |
| `data/allsamples.vep.maf` | MAF file containing annotated somatic variants for all samples. |
| `data/` | Folder containing per-sample mutation tables or VCF files used by the intersection and mutational spectrum modules. |
| `Groups.xlsx` | Sample grouping table. The script currently reads sheet 2 and expects columns such as `Samples` and `Group2`. |
| `GENE.xlsx` | Optional gene list used by the CNV heatmap module. |
| `.cns` files | Optional CNVkit segment files used by the CNV module. |

Before running the analysis, update the placeholder paths and filename patterns in `01.Organoid Similarity Analysis/Organoid_mutation.R` to match the local project data.

## R dependencies

The script uses the following major R packages:

```r
BSgenome
MutationalPatterns
maftools
ComplexHeatmap
tidyverse
data.table
reshape2
openxlsx
readxl
ggplot2
ggsci
patchwork
VennDiagram
venn
UpSetR
ggupset
ggplotify
ggimage
circlize
NMF
```

For the current human reference setting, the script uses:

```r
BSgenome.Hsapiens.UCSC.hg19
```

Install missing packages from CRAN or Bioconductor before running the workflow.

## How to run

From the repository root, place private input data in the expected local folders, then run:

```r
source("01.Organoid Similarity Analysis/Organoid_mutation.R")
```

or from a shell with R installed:

```bash
Rscript "01.Organoid Similarity Analysis/Organoid_mutation.R"
```

The default output directory used by the script is `RESULT/`.

## Notes for reuse

This repository currently preserves the analysis logic as a single research script. Some blocks contain project-specific placeholders, disabled exploratory sections, and local file assumptions. Before reuse on a new dataset, check:

* input paths and file patterns
* sample names and grouping metadata
* reference genome build
* selected genes for CNV or oncoplot visualization
* output filenames and figure sizes

## Citation and contact

If this repository accompanies a manuscript, project report, or public dataset, add the citation, data accession, and contact information here before publication.
