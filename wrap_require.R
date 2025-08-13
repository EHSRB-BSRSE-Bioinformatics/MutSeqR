wrap_examples_with_require <- function(file) {
  lines <- readLines(file, warn = FALSE)
  example_starts <- grep("^#' @examples", lines)
  if (length(example_starts) == 0) return(FALSE)

  changed <- FALSE
  for (start in rev(example_starts)) {
    # Find the last line of this roxygen example block
    end <- start + 1
    while (end <= length(lines) && grepl("^#'", lines[end])) {
      end <- end + 1
    }
    end <- end - 1  # adjust because loop exits after first non-#' line

    # Skip if already wrapped
    if (start + 1 <= length(lines) &&
        grepl("^#' if \\(requireNamespace\\(\"MutSeqRData\"", lines[start + 1])) next

    # Insert opening brace line after @examples
    lines <- append(lines,
                    "#' if (requireNamespace(\"MutSeqRData\", quietly = TRUE)) {",
                    after = start)

    # Adjust end index due to inserted line
    end <- end + 1

    # Insert closing brace line after last roxygen example line
    lines <- append(lines,
                    "#' }",
                    after = end)

    changed <- TRUE
  }

  if (changed) writeLines(lines, file)
  return(changed)
}

# Apply to all R files in 'R' directory
files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
modified_files <- sapply(files, wrap_examples_with_require)
cat("Modified files:\n")
print(files[modified_files])