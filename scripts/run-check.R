#!/usr/bin/env Rscript

# Include the PDF manual and vignettes. Keep this identical locally and in CI.
check_dir <- Sys.getenv("MUTSEQR_CHECK_DIR", "")
rcmdcheck::rcmdcheck(
  args = character(), build_args = character(), error_on = "note",
  check_dir = if (nzchar(check_dir)) check_dir else NULL
)
