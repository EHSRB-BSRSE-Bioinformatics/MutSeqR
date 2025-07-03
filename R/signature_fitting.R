#' Run COSMIC signatures comparison using SigProfilerAssignment
#'
#' @details Assign COSMIC SBS signatures to mutation data using
#' SigProfilerAssignment. Data is cleaned and formatted for input into
#' SigProfiler tools. This function will create a virtual environment using
#' reticulate to run python, as this is a requirement for the SigProfiler suite
#' of tools. Note that it will also install several python dependencies using
#' a conda virtual environment on first use. Please be aware of the
#' implications of this. For advanced use, it is suggested to use the
#' SigProfiler python tools directly in python as described in their
#' respective documentation. Users must have python installed on their
#' computer to use this function.
#' @param mutation_data A data frame containing mutation data.
#' @param project_name The name of the project. This is used to format
#' the data into required .txt format for SigProfiler tools.
#' @param project_genome The reference genome to use. On first use, the
#' function will install the genome using SigProfilerMatrixGeneratorR::install.
#' e.x. GRCh37, GRCH38, mm10, mm9, rn6
#' @param env_name The name of the virtual environment. This will be created on
#' first use.
#' @param group The column in the mutation data used to aggregate groups.
#' Signature assignment will be performed on each group separately.
#' @param output_path The filepath to the directory in which the output folder
#' will be created to store results. Default is NULL. This will store results
#' in the current working directory.
#' @param python_version The version of python installed on the user's
#' computer.
#' @returns Creates a subfolder "SigProfiler" in the output directory with
#' SigProfiler tools results. For a complete breakdown of the results, see the
#' Readme file for MutSeqR. Most relevant results are stored in SigProfiler >
#' [group] > matrices > output > Assignment_Solution > Activities >
#' SampleReconstruction > WebPNGs.
#' These plots show a summary of the signature assignment results for each
#' group. In each plot, the top left panel represents the base_96 mutation
#' count for the group. The bottom left panel represents the reconstructed
#' profile. Below the reconstruction are the solution statistics that indicate
#' the goodness of fit of the reconstructed profile to the observed profile.
#' (Recommended cosine similarity > 0.9). The panels on the right represent the
#' SBS signatures that contribute to the reconstructed profile. The signature
#' name and its contribution % are shown in the panel. A high contribution
#' means a high association of the signature with the group's mutation
#' spectra.
#' @details Mutation data will be filtered to only include SNVs. Variants
#' flagged by the filter_mut column will be excluded.
#' @examples
#' \dontrun{
#' # Example data consists of 24 mouse bone marrow DNA samples imported
#' # using import_mut_data() and filtered with filter_mut as in Example 4.
#' # Sequenced on TS Mouse Mutagenesis Panel. Example data is
#' # retrieved from MutSeqRData, an ExperimentHub data package.
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' example_data <- eh[["EH9861"]]
#' 
#' signature_fitting(mutation_data = example_data,
#'                   project_name = "Example",
#'                   project_genome = "mm10",
#'                   env_name = "MutSeqR",
#'                   group = "dose",
#'                   python_version = "3.11")
#' }
#' @importFrom here here
#' @importFrom dplyr filter select rename mutate relocate
#' @importFrom utils write.table
#' @importFrom rlang .data
#' @import stringr
#' @export
#'
signature_fitting <- function(mutation_data,
                              project_name = "Default",
                              project_genome = "GRCh38",
                              env_name = "MutSeqR",
                              group = "sample",
                              output_path = NULL,
                              python_version) {
  if (!requireNamespace("reticulate")) {
    stop("Reticulate not installed: you need this to run SigProfiler tools in R.")
  }
  if (!requireNamespace("SigProfilerMatrixGeneratorR")) {
    stop("SigProfilerMatrixGeneratorR not installed: you need this to run SigProfiler tools in R. Install using devtools::install_github('AlexandrovLab/SigProfilerMatrixGeneratorR')")
  }
  message("Note: This function requires python to be installed on the users 
          computer. If you do not have python installed, you can do so using: 
          reticulate::install_python().
          \n\nThis function will create a virtual environment using reticulate
          to run python. Note that it will also install several python
          dependencies using a conda virtual environment on first use. Please
          be aware of the implications of this. For advanced use, it is
          suggested to use the SigProfiler python tools directly in python as
          described in their respective documentation.")
  setup_mutseqr_python(force = FALSE)
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
      # Create venv and install packages
      reticulate::virtualenv_create(env_name, python = reticulate::virtualenv_starter(python_version))
      # Install required packages
      # Patched version of pandas and scipy to avoid dependency errors: binom_test
      reticulate::virtualenv_install(env_name, c("SigProfilerMatrixGenerator",
                                                 "SigProfilerAssignment",
                                                 "SigProfilerExtractor",
                                                 "pandas==1.5.3",
                                                 "scipy==1.11.4",
                                                 "pypdf==4.3.1"))
    } else {
      # User chose not to install the packages
      cat("Installation aborted by the user.\n")
      stop("Function terminated.")
    }
  }

  # reticulate::install_python(version = python_version)
  reticulate::use_virtualenv(env_name)

  SigProfilerMatrixGeneratorR::install(project_genome)
  signatures_python_code <- system.file("extdata", "signatures.py",
                                        package = "MutSeqR")
  sig_py <- reticulate::source_python(signatures_python_code, envir = new.env())

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
    dplyr::filter(.data$filter_mut == FALSE) %>%
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
  sig_py$cosmic_fit_MutSeqR(
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
