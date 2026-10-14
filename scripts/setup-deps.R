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
          "knitr", "rmarkdown", "pkgdown", "testthat", "yaml", "roxygen2")
missing <- setdiff(boot, rownames(installed.packages()))
if (length(missing)) install.packages(missing)

## Match the version used to generate committed Rd files ----------------
roxygen_version <- read.dcf("DESCRIPTION")[1, "Config/roxygen2/version"]
if (!identical(as.character(utils::packageVersion("roxygen2")), roxygen_version)) {
  remotes::install_version("roxygen2", version = roxygen_version,
                           repos = "https://cloud.r-project.org",
                           upgrade = "never")
}

## Bioconductor --------------------------------------------------------
bioc_ver <- Sys.getenv("BIOC_VERSION", unset = "")
if (nzchar(bioc_ver)) {
  BiocManager::install(version = bioc_ver, ask = FALSE, update = FALSE)
} else {
  BiocManager::install(ask = FALSE, update = FALSE)
}
BiocManager::install("BiocCheck", ask = FALSE, update = FALSE)

## package dependencies (hard + soft, from DESCRIPTION) ---------------
# Include the matching Bioconductor repositories so remotes can resolve
# Bioconductor Imports as well as CRAN dependencies.
repos <- BiocManager::repositories()

remotes::install_local(
  dependencies = TRUE,
  repos = repos,
  force_update = FALSE,
  upgrade = "never"
)

message("MutSeqR dependencies installed.")
