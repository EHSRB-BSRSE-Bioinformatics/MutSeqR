#!/usr/bin/env Rscript
# Install development tooling and all package dependencies (from
# DESCRIPTION) into the current R environment (e.g. the pixi env).
#
# Usage: Rscript scripts/setup-deps.R
#
# Set BIOC_VERSION to force a specific Bioconductor version
# (default: latest).

options(repos = c(CRAN = "https://cloud.r-project.org"))

## bootstrap dev tooling ---------------------------------------------
boot <- c("devtools", "pak", "remotes", "BiocManager", "rcmdcheck",
          "BiocCheck", "knitr", "rmarkdown", "pkgdown", "testthat")
missing <- setdiff(boot, rownames(installed.packages()))
if (length(missing)) install.packages(missing)

## Bioconductor --------------------------------------------------------
bioc_ver <- Sys.getenv("BIOC_VERSION", unset = "")
if (nzchar(bioc_ver)) {
  BiocManager::install(version = bioc_ver, ask = FALSE, update = FALSE)
} else if (!requireNamespace("BiocManager", quietly = TRUE)) {
  BiocManager::install(ask = FALSE, update = FALSE)
}

## package dependencies (hard + soft, from DESCRIPTION) ---------------
remotes::install_local(
  dependencies = TRUE,
  force_update = FALSE,
  upgrade = "never"
)

message("MutSeqR dependencies installed.")
