# For this to work there are a few dependencies that must be met.
# You must first install python dependencies, OR use the venv approach that is the default method here.
# TO DO: SigProfiler packages are in suggests; write a feature that will ask users if they want to install them. 

#' Run COSMIC signatures comparison
#'
#' After cleaning the mutation data input, runs several Alexandrov Lab tools for COSMIC signature analysis (assigns signatures to best explain the input data).
#' @param mutation_data A data frame, imported from a .mut file
#' @param project_name The name of the project; used to get mutation data into the required .txt format for SigProfiler
#' @param project_genome A string describing the reference genome to use; 
#' e.g. GRCh37, GRCH38, mm10, mm9, rn6
#' output_path The directory where output results should be written. *not a parameter of the function
#' @param env_name The name of the virtual environment. This will be created on first use. 
#' @param group The column in the mutation data used to aggregate groups (e.g., sample ID, tissue, dose)
#' @param output_path The filepath to thedirectory in which the output folder will be created to store results.
#' @param python_version The version of python to be used. 
#' @param python_path The path to the version of python to be used with reticulate. It is important that this version of python meets the dependencies, including the SigProfiler python tools.
#' @param python_home The path to the conda virtual environment that contains the required python dependencies
#' @returns Creates a subfolder in the output directory with SigProfiler tools results.
#'  SigProfilerMatrixGeneratorR  SigProfilerMatrixGeneratorR install
#' @importFrom here here
#' @importFrom dplyr filter select rename mutate relocate
#' @importFrom utils write.table
#' @importFrom rlang .data
#' @import reticulate
#' @import stringr
#' @export
#' 
signature_fitting <- function(mutation_data,
                              project_name = "Default",
                              project_genome = "GRCh38",
                              env_name = "MutSeqR",
                              group = "sample",
                              output_path = NULL,
                              python_version = "3.11", # We should test with other versions
                              python_path = "~/../../AppData/Local/Programs/Python/Python310/python.exe", #"/usr/bin/python3.9"
                              python_home = "C:/Users/adodge/OneDrive - HC-SC PHAC-ASPC/Documents/.virtualenvs/r-reticulate") {
  if (!requireNamespace("reticulate")) {
    stop("Reticulate not installed: you need this to run SigProfiler tools in R.")
  }
  message("Note: This function requires python to be installed on the users 
          computer. If you do not have python installed, you can do so using: 
          reticulate::install_python().
          \n\nThis function will create a virtual environment using reticulate to 
          run python, as this is a requirement for the SigProfiler suite of tools.
          Note that it will also install several python dependencies using
          a conda virtual environment on first use. Please be aware of the 
          implications of this. For advanced use, it is suggested to
          use the SigProfiler python tools directly in python as described
          in their respective documentation.")

  # Only run this once, not every time that the function is called
  #python_version <- "3.11:latest

  installed_envs <- reticulate::virtualenv_list()
  # Check if MutSeqR virtualenv already exists
  if (env_name %in% installed_envs) {
    # Virtualenv exists, use it
    reticulate::use_virtualenv(env_name)
  } else {
    # Ask the user for confirmation
    user_input <- utils::menu("Do you want to create a virtual environment and 
                              install the required Python packages? This may 
                              take several minutes",
                              title = "Confirmation", choices = c("Yes", "No"))

    if (user_input == 1) {
      # User chose to install the packages
      # Virtualenv doesn't exist, set it up
      reticulate::virtualenv_create(env_name, python = reticulate::virtualenv_starter(python_version))
      # Install required packages
        # Patched version of pandas and scipy to avoid dependency errors: binom_test
      reticulate::virtualenv_install(env_name, c("SigProfilerMatrixGenerator", "SigProfilerAssignment", "SigProfilerExtractor", "pandas==1.5.3", "scipy==1.11.4"))
    } else {
      # User chose not to install the packages
      cat("Installation aborted by the user.\n")
      stop("Function terminated.")
    }
  }

  # reticulate::install_python(version = python_version)
  reticulate::use_virtualenv(env_name)
  
  # This will only install the genome if it's not found in the current env
  if (!requireNamespace("SigProfilerMatrixGeneratorR")) {
    stop("SigProfilerMatrixGeneratorR not installed: you need this to run SigProfiler tools in R. Install using devtools::install_github('AlexandrovLab/SigProfilerMatrixGeneratorR')")
  }
  SigProfilerMatrixGeneratorR::install(project_genome)
  signatures_python_code <- system.file('extdata', 'signatures.py',
                                        package = "MutSeqR")
  reticulate::source_python(signatures_python_code)

  message("Creating cleaned data for input into SigProfiler...")
  # Clean data into required format for Alexandrov Lab tools...
  #ID doesn't always exist. 
  signature_data <- as.data.frame(mutation_data)

  # Check if "id" column exists
  if (!"id" %in% colnames(signature_data)) {
    # If not, create "id" column and populate with "."
   signature_data$id <- "."
  }

  # Check if "seqnames" column exists (ie if it came from a GRanges)
  if ("seqnames" %in% colnames(signature_data)) {
    # If "seqnames" exists, rename it to "contig"
    signature_data <- dplyr::rename(signature_data, contig = seqnames)
  }

  signature_data <- signature_data %>%
    dplyr::filter(.data$variation_type %in% "snv") %>%
    dplyr::filter(.data$is_germline == FALSE) %>%   
    dplyr::select(all_of(group), "id", "variation_type", "contig", "start", "end", "ref", "alt") %>%
    dplyr::rename(
      "Samples" = all_of(group),
      "ID" = "id",
      "mut_type" = "variation_type",
      "chrom" = "contig",
      "pos_start" = "start",
      "pos_end" = "end"
    ) %>%
    dplyr::mutate(Samples = stringr::str_replace_all(.data$Samples, " ", "_")) %>%
    dplyr::mutate(chrom = stringr::str_replace(.data$chrom, "chr", "")) %>%
    dplyr::mutate(Project = project_name) %>%
    dplyr::mutate(Genome = project_genome) %>%
    dplyr::mutate(Type = "SOMATIC") %>%
    dplyr::relocate(.data$Project) %>%
    dplyr::relocate(.data$Genome, .after = .data$ID) %>%
    dplyr::mutate(mut_type = "SNP") # This should be fixed before using on other datasets.
  
# Make sure Samples column is NOT numeric
# Note that the values will be class character, but even if so,
  # number values will cause an issue
signature_data <- signature_data %>%
  dplyr::mutate(Samples = paste0(!!group, "_", Samples))

  message("Generating output path string...")
  if (is.null(output_path)) {
    output_path <- file.path(
      here::here(),
      "output",
      "SigProfiler",
      group
    )
  } else {
    output_path <- file.path(
      output_path,
      "SigProfiler",
      group
    )
  }

  message(paste0("Creating directory ", output_path))
  if (!dir.exists(file.path(output_path, "matrices"))) {
    dir.create(file.path(output_path, "matrices"), recursive = TRUE)
  }

  message("Writing mutation matrix to use as input to SigProfiler...")
  write.table(signature_data,
    file = file.path(output_path, "matrices", "mutations.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  message("Running SigProfilerMatrixGeneratorR...")
  signature_matrices <-
    SigProfilerMatrixGeneratorR::SigProfilerMatrixGeneratorR(
      project = project_name,
      genome = project_genome,
      matrix_path = file.path(output_path, "matrices"),
      plot = TRUE,
      exome = FALSE,
      bed_file = NULL,
      chrom_based = FALSE,
      tsb_stat = TRUE,
      seqInfo = TRUE,
      cushion = 100
    )

  message("Running COSMIC fitting...")
  cosmic_fit_MutSeqR(
    samples = file.path(output_path, "matrices", "output", "SBS",
                        paste0(project_name, ".SBS96.all")),
    output = file.path(output_path, "matrices", "output"),
    input_type = "matrix", # "vcf", "seg:TYPE", "matrix"
    context_type = "96", # Required for vcf input
    cosmic_version = 3.4,
    exome = FALSE,
    genome_build = project_genome,
    signature_database = NULL, #tab delimited file of signatures
    exclude_signature_subgroups = NULL,
    export_probabilities = TRUE,
    export_probabilities_per_mutation = FALSE, # Only for vcf input
    make_plots = TRUE,
    sample_reconstruction_plots = "png",
    verbose = TRUE
  )
}
