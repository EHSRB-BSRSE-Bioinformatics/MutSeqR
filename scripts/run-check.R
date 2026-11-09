#!/usr/bin/env Rscript

# Include the PDF manual and vignettes. Keep this identical locally and in CI.
check_dir <- Sys.getenv("MUTSEQR_CHECK_DIR", "")
skip_vignettes <- tolower(Sys.getenv("MUTSEQR_SKIP_VIGNETTES", "false")) %in%
  c("true", "1", "yes")

if (skip_vignettes) {
  message("Skipping vignette build/check as requested by MUTSEQR_SKIP_VIGNETTES.")
  build_args <- "--no-build-vignettes"
  check_args <- "--ignore-vignettes"
} else {
  prefetch_status <- system2(
    file.path(R.home("bin"), "Rscript"),
    file.path("scripts", "prefetch-vignette-data.R")
  )
  if (prefetch_status != 0L) {
    stop("Failed to prefetch ExperimentHub data needed to build vignettes.")
  }
  build_args <- character()
  check_args <- character()
}

rcmdcheck::rcmdcheck(
  args = check_args, build_args = build_args, error_on = "note",
  check_dir = if (nzchar(check_dir)) check_dir else NULL
)
