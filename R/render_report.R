#' Read configuration file and render R Markdown document
#'
#' This function reads a configuration file in YAML format, extracts the parameters,
#' and renders an R Markdown document using the specified parameters.
#'
#' @param config_filepath The path to the configuration file.
#' @param output_file The path to the output file.
#'
#' @return None
#'
#' @examples
#' read_config_and_render("config.yaml", "output.html")
#'
#' @import yaml
#' @import rmarkdown
#' @importFrom utils install.packages
#' @export
render_report <- function(config_filepath, output_file) {
  config_filepath <- "./config.yml"
  output_file <- "./output.html"
  # Check if yaml package is available, install if not
  if (!requireNamespace("yaml", quietly = TRUE)) {
    install.packages("yaml")
  }

  # Check if rmarkdown package is available, install if not
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    install.packages("rmarkdown")
  }

  # Check if MutSeqR package is available, install if not
  if (!requireNamespace("MutSeqR", quietly = TRUE)) {
    install.packages("MutSeqR")
  }
  
  # take config_filepath and make sure it is a valid file, convert to path if needed
  config_filepath <- normalizePath(config_filepath)
  if (!file.exists(config_filepath)) {
    stop("Config file not found")
  }
  
  # Read the configuration file
  config <- yaml::yaml.load_file(config_filepath)
  
  # Use provided params or default to list if none provided
  #params <- ifelse(!is.null(config$params), config$params, list())
  # Construct the path to the .Rmd file within the installed package directory
  rmd_file <- "DS_summary_report.Rmd"
  rmd_path <- system.file("extdata", rmd_file, package = "MutSeqR", mustWork = TRUE)
  # Rendering the R Markdown document
  rmarkdown::render(
    input = rmd_path,
    output_file = output_file,
    params = config$params,
    envir = new.env()
  )
}