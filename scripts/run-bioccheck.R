#!/usr/bin/env Rscript

# BiocCheck installs the package; use a writable library for non-root containers.
check_library <- tempfile("mutseqr-bioccheck-library-")
dir.create(check_library)
.libPaths(c(check_library, .libPaths()))
BiocCheck::BiocCheck(
  ".", checkDir = Sys.getenv("MUTSEQR_RESULTS_DIR", tempdir()),
  `quit-with-status` = TRUE
)
