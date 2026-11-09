#!/usr/bin/env Rscript

mode <- commandArgs(trailingOnly = TRUE)
if (length(mode) != 1L || !mode %in% c("validate", "release-check")) {
  stop("Usage: Rscript scripts/run-validation.R {validate|release-check}")
}

scripts <- c(
  if (mode == "release-check") "check-release.R",
  "check-docs.R", "check-readme.R", "check-reference.R",
  "run-tests.R", "run-check.R", "run-bioccheck.R"
)
for (script in scripts) {
  message("==> ", script)
  status <- system2(file.path(R.home("bin"), "Rscript"), file.path("scripts", script))
  if (status != 0L) quit(save = "no", status = status)
}
