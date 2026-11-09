#!/usr/bin/env Rscript
# Validate the metadata of the version that will actually be shipped.

description <- read.dcf("DESCRIPTION")
version <- unname(description[1, "Version"])
if (!grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", version)) {
  stop("DESCRIPTION Version must have the form x.y.z")
}
news <- readLines("NEWS.md", warn = FALSE)
headings <- grep("^# ", news, value = TRUE)
if (!length(headings) || !identical(headings[1], paste("# MutSeqR", version))) {
  stop("The first NEWS.md heading must be '# MutSeqR ", version, "'")
}
message("Release metadata matches MutSeqR ", version, ".")
