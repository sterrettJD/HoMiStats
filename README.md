# HoMiStats

Statistics tools for host-microbiome dual transcriptomics

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/sterrettJD/CauDA/actions/workflows/r.yml/badge.svg)](https://github.com/sterrettJD/CauDA/actions/workflows/r.yml)
  <!-- badges: end -->


# Installation

`devtools::install_github("sterrettJD/HoMiStats")`

# Package capabilities

## Statistical modules
HoMiStats consists of two main modules:

1.  Metatranscriptomics differential expression (mtxDE)

    -   Zero inflated beta regression, motivated by [this review](https://academic.oup.com/bib/article/24/5/bbad279/7239897) and using a model inspired by [Peng et al., 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6109378/)
        - Implemented via [GAMLSS BEZI](https://www.rdocumentation.org/packages/gamlss.dist/versions/6.1-1/topics/BEZI) and [ZIBR](https://github.com/PennChopMicrobiomeProgram/ZIBR) (supports random effects) in the function `run_mtxDE()`
    -   Other models? Not yet implemented...

2.  Construction of correlation network between host and microbial transcripts - Not yet implemented...

## Utilities
HoMiStats additionally contains tools for transforming metatranscriptomic data
