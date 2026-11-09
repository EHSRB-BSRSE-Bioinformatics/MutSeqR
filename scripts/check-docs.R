#!/usr/bin/env Rscript
# Regenerate roxygen output and fail if doing so changes generated files.
# Only files changed by this invocation are considered, so unrelated existing
# worktree edits do not affect the result.

description <- read.dcf("DESCRIPTION")
expected <- unname(description[1, "Config/roxygen2/version"])
if (is.na(expected) || !nzchar(expected)) {
  stop("DESCRIPTION must set Config/roxygen2/version", call. = FALSE)
}
if (!requireNamespace("roxygen2", quietly = TRUE)) {
  stop("Install roxygen2 ", expected, " before checking documentation", call. = FALSE)
}
installed <- as.character(utils::packageVersion("roxygen2"))
if (!identical(installed, expected)) {
  stop("roxygen2 version mismatch: DESCRIPTION requires ", expected,
       ", but ", installed, " is installed", call. = FALSE)
}

generated_files <- function() {
  c("NAMESPACE", "DESCRIPTION", sort(list.files("man", pattern = "\\.Rd$",
                                                   recursive = TRUE, full.names = TRUE)))
}

file_state <- function(files) {
  exists <- file.exists(files)
  hashes <- rep(NA_character_, length(files))
  hashes[exists] <- unname(tools::md5sum(files[exists]))
  stats::setNames(hashes, files)
}

description_fields <- function() {
  current <- read.dcf("DESCRIPTION")
  fields <- c("Roxygen", "RoxygenNote", "Config/roxygen2/version")
  stats::setNames(vapply(fields, function(field) {
    if (field %in% colnames(current)) current[1, field] else NA_character_
  }, character(1)), fields)
}

before_files <- generated_files()
before_state <- file_state(before_files)
before_description <- description_fields()

devtools::document()

after_files <- generated_files()
after_state <- file_state(after_files)
all_files <- union(names(before_state), names(after_state))
before_values <- before_state[all_files]
after_values <- after_state[all_files]
changed <- all_files[is.na(before_values) != is.na(after_values) |
                       (!is.na(before_values) & !is.na(after_values) &
                          before_values != after_values)]
after_description <- description_fields()
changed_metadata <- names(before_description)[
  is.na(before_description) != is.na(after_description) |
    (!is.na(before_description) & !is.na(after_description) &
       before_description != after_description)
]

if (length(changed) || length(changed_metadata)) {
  if (length(changed)) {
    message("Documentation generation changed these files:")
    message(paste0("  ", changed, collapse = "\n"))
  }
  if (length(changed_metadata)) {
    message("Documentation generation changed DESCRIPTION metadata: ",
            paste(changed_metadata, collapse = ", "))
  }
  stop("Commit the generated roxygen output and rerun this check", call. = FALSE)
}

message("Roxygen documentation is current (roxygen2 ", expected, ").")
