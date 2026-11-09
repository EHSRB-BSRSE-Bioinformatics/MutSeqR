#!/usr/bin/env Rscript
# Render README.Rmd and fail only when this invocation changes README.md.

before <- if (file.exists("README.md")) tools::md5sum("README.md") else NA_character_
rmarkdown::render("README.Rmd", output_format = "github_document",
                  output_file = "README.md", clean = TRUE)
after <- if (file.exists("README.md")) tools::md5sum("README.md") else NA_character_
if (!identical(unname(before), unname(after))) {
  stop("README.md is stale relative to README.Rmd; commit the rendered output",
       call. = FALSE)
}
message("README.md is current.")
