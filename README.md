# HoMiStats

Statistics tools for host-microbiome dual transcriptomics

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/sterrettJD/HoMiStats/actions/workflows/r.yaml/badge.svg)](https://github.com/sterrettJD/HoMiStats/actions/workflows/r.yaml)
  <!-- badges: end -->


# Installation

## Bioconductor (preferred)
**Pending Bioconductor approval** HoMiStats can be installed vi:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("HoMiStats")
```

## GitHub installation (not preferred)
To install HoMiStats with the latest updates from GitHub:
`devtools::install_github("sterrettJD/HoMiStats")`

# Package capabilities

## Statistical modules
HoMiStats consists of two main modules:

1.  Metatranscriptomics differential expression (mtxDE implemented in `run_mtxDE()`)

    -   Zero inflated beta regression, motivated by [this review](https://academic.oup.com/bib/article/24/5/bbad279/7239897) and using a model inspired by [Peng et al., 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6109378/)
        - Implemented via [GAMLSS BEZI](https://www.rdocumentation.org/packages/gamlss.dist/versions/6.1-1/topics/BEZI) and [ZIBR](https://github.com/PennChopMicrobiomeProgram/ZIBR) (supports random effects) 
    -   Linear regression and linear mixed-effects regression

2.  Construction of correlation network between host and microbial transcripts (HoMiCorr implemented in `run_HoMiCorr()`)
    - Also supports both zero inflated beta regression (gamlss or ZIBR for random effects) and linear regression or linear mixed-effects regression

## Utilities
HoMiStats contains tools for transforming metatranscriptomic data, as well as tools for filtering genes based on low expression, prevalence, or variance.
