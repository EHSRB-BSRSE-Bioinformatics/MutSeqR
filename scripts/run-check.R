#!/usr/bin/env Rscript

# Include the PDF manual and vignettes. Keep this identical locally and in CI.
check_dir <- Sys.getenv("MUTSEQR_CHECK_DIR", "")
prefetch_status <- system2(
  file.path(R.home("bin"), "Rscript"),
  file.path("scripts", "prefetch-vignette-data.R")
)
if (prefetch_status != 0L) {
  stop("Failed to prefetch ExperimentHub data needed to build vignettes.")
}

rcmdcheck::rcmdcheck(
  args = character(), build_args = character(), error_on = "note",
  check_dir = if (nzchar(check_dir)) check_dir else NULL
)
