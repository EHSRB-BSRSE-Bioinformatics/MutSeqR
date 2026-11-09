#!/usr/bin/env Rscript

retry <- function(expr, label, attempts = 3L) {
  for (attempt in seq_len(attempts)) {
    result <- tryCatch(
      list(value = force(expr)),
      error = function(error) list(error = error)
    )

    if (!is.null(result$value) || is.null(result$error)) {
      return(result$value)
    }

    if (attempt == attempts) {
      stop(result$error)
    }

    message(
      label, " failed (attempt ", attempt, "/", attempts,
      "); retrying: ", conditionMessage(result$error)
    )
    Sys.sleep(5 * attempt)
  }
}

eh <- retry(ExperimentHub::ExperimentHub(), "ExperimentHub index")

for (accession in sprintf("EH%d", 9857:9865)) {
  retry({
    invisible(eh[[accession]])
    TRUE
  }, accession)
}

message("Vignette ExperimentHub resources are available in the local cache.")