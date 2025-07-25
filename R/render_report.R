#' Read configuration file and render R Markdown document
#'
#' @description This function reads a configuration file in YAML format,
#' extracts the parameters, and renders an R Markdown document using the
#' specified parameters.
#'
#' @param config_filepath The path to the configuration file.
#' @param output_file The name of the output file. Will be saved to the
#' outputdir in config params.
#' @param output_format The format of the output file. Options are
#' "html_document" (default), "pdf_document", or "all".
#' @return None
#'
#' @importFrom utils install.packages
#' @importFrom here here
#' @export
render_report <- function(
  config_filepath,
  output_file = "./MutSeqR_Summary_Report.html",
  output_format = "html_document"
) {
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
    install.packages("MutSeqR") # Fix: Will depend on where it gets published.
  }

  # Validate config_filepath
  config_filepath <- normalizePath(config_filepath)
  if (!file.exists(config_filepath)) {
    stop("Config file not found")
  }

  # Read the configuration file
  config <- yaml::yaml.load_file(config_filepath)
  params <- config$params

  # Set directories
  if (is.null(params$projectdir)) {
    params$projectdir <- here::here()
  }
  if (!file.exists(normalizePath(params$projectdir))) {
    stop("Your project directory doesn't exist")
  }
  if (is.null(params$outputdir)) {
    params$outputdir <- params$projectdir
  }
  if (!file.exists(normalizePath(params$outputdir))) { # check for path from current wd
    output_path <- file.path(params$projectdir, params$outputdir)
    if (!file.exists(normalizePath(output_path))) { # check for path from projectdir
      dir.create(normalizePath(output_path)) # create as needed, within projectdir
    }
    params$outputdir <- normalizePath(output_path) # set the param
  }

  # Load the profile config
  profile_config <- yaml::yaml.load_file(system.file("extdata", "inputs", "profile_config.yaml", package = "MutSeqR"))
  if (!params$config_profile %in% c(
    "Duplex Sequencing Mouse Mutagenesis Panel",
    "Duplex Sequencing Human Mutagenesis Panel",
    "Duplex Sequencing Rat Mutagenesis Panel",
    "CODEC",
    "None"
  )) {
    stop("Invalid profile configuration in the YAML file.")
  }
  if (grepl("^Duplex Sequencing (Human|Mouse|Rat) Mutagenesis Panel$",
            params$config_profile)) {
    # Join the DS parameters with params list.
    params <- c(params, profile_config$Duplex_Sequencing)
    # Define the regions parameter for TS panels.
    params$regions <- paste0("TSpanel_",
      tolower(sub("Duplex Sequencing (Human|Mouse|Rat) Mutagenesis Panel",
      "\\1", params$config_profile, perl = TRUE))
    )
    message("Setting up parameters for Duplex Sequencing on ", params$regions)
    params$filtering_regions <- params$regions
  } else if (params$config_profile == "CODEC") {
    params <- c(params, profile_config$CODEC)
  } else if (params$config_profile == "None") {
    params <- c(params, config$Custom_Profile_Params)
  }

  # Construct the path to the .Rmd file within the installed package directory
  rmd_file <- "Summary_report.Rmd"
  rmd_path <- system.file("extdata", rmd_file,
                          package = "MutSeqR", mustWork = TRUE)
  message("project directory", params$projectdir)
  message("output directory:", params$outputdir)

  # Validate the Parameters
  if (!is.character(params$project_title)) stop("project_title must be a character string")
  if (!is.character(params$researcher_name)) stop("researcher_name must be a character string")
  if (!is.character(params$user_name)) stop("user_name must be a character string")
  if (!params$file_type %in% c("table", "vcf")) stop("file_type must be either 'table' or 'vcf'")

  # Rendering the R Markdown document
  rmarkdown::render(
    input = rmd_path,
    output_dir = params$outputdir,
    output_file = output_file,
    output_format = output_format,
    params = params,
    envir = new.env(),
    knit_root_dir = params$projectdir
  )
}
