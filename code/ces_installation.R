#!/usr/bin/env Rscript

## In debian it needs the packages r-base and r-base-dev (probably
## more than that)

.libPaths(c("./.Rlibs", .libPaths()))

options(timeout = 600)

install.packages("remotes")

remotes::install_github("Townsend-Lab-Yale/cancereffectsizeR@*release",
  dependencies = TRUE)

library(BiocManager)
install("BSgenome.Hsapiens.UCSC.hg38")

remotes::install_github("Townsend-Lab-Yale/ces.refset.hg38@*release")
