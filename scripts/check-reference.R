#!/usr/bin/env Rscript
# Use pkgdown's native reference coverage check, then build the reference
# without executing examples. The build goes to a temporary directory.

if (!requireNamespace("pkgdown", quietly = TRUE)) {
  stop("Install pkgdown to validate the reference index", call. = FALSE)
}

pkgdown::check_pkgdown(".")
output <- tempfile("mutseqr-reference-")
dir.create(output)
pkgdown::build_reference(pkg = ".", examples = FALSE, lazy = FALSE,
                          preview = FALSE,
                          override = list(destination = output))
message("Pkgdown reference coverage is current (examples disabled).")
