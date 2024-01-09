# HoMiStats

Statstics tools for host-microbiome dual transcriptomics

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/sterrettJD/CauDA/actions/workflows/r.yml/badge.svg)](https://github.com/sterrettJD/CauDA/actions/workflows/r.yml)
  <!-- badges: end -->


# Installation

`devtools::install_github("sterrettJD/HoMiStats")`

# Package capabilities

HoMiStats consists of two main modules:

1.  Metatranscriptomics differential expression (mtxDE)

    -   Zero inflated beta regression, motivated by [this review](https://academic.oup.com/bib/article/24/5/bbad279/7239897) and using a model inspired by [Peng et al., 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6109378/)
    -   Other models? Not yet implemented...

2.  Construction of correlation network between host and microbial transcripts - Not yet implemented...
