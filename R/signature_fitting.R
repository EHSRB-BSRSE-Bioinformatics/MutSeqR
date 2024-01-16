# For this to work there are a few dependencies that must be met.
# You must first install python dependencies, OR use the venv approach that is the default method here.
# TO DO: SigProfiler packages are in suggests; write a feature that will ask users if they want to install them. 

#' Run COSMIC signatures comparison
#'
#' After cleaning the mutation data input, runs several Alexandrov Lab tools for COSMIC signature analysis (assigns signatures to best explain the input data).
#' @param mutations A data frame, imported from a .mut file
#' @param project_name The name of the project; used to get mutation data into the required .txt format for SigProfiler
#' @param project_genome A string describing the reference genome to use; e.g., GRCh38
#' output_path The directory where output results should be written. *not a parameter of the function
#' @param group The column in the mutation data used to aggregate groups (e.g., sample ID, tissue, dose)
#' @param python_path The path to the version of python to be used with reticulate. It is important that this version of python meets the dependencies, including the SigProfiler python tools.
#' @param python_home The path to the conda virtual environment that contains the required python dependencies
#' @returns Creates a subfolder in the output directory with SigProfiler tools results.
#'  SigProfilerAssignmentR cosmic_fit
#'  SigProfilerExtractorR sigprofilerextractor
#'  SigProfilerMatrixGeneratorR  SigProfilerMatrixGeneratorR install
#' @importFrom here here
#' @importFrom dplyr filter select rename mutate relocate 
#' @importFrom utils write.table
#' @importFrom rlang .data
#' @import reticulate
#' @import stringr
#' @export
signature_fitting <- function(mutations,
                              project_name = "Default",
                              project_genome = "GRCh38",
                              group = "sample",
                              python_version = "3.11:latest", # We should test with other versions
                              python_path = "~/../../AppData/Local/Programs/Python/Python310/python.exe", #"/usr/bin/python3.9"
                              python_home = "~/.virtualenvs/r-reticulate") {
  if (!requireNamespace("reticulate")) {
    stop("Reticulate not installed: you need this to run SigProfiler tools in R.")
  }
  message("This function will attempt to install python using reticulate, as
          this is a requirement for the SigProfiler suite of tools.
          Note that it will also install several python dependencies using
          a conda virtual environment on first use. Please be aware of the 
          implications of this. For advanced use, it is suggested to
          use the SigProfiler python tools directly in python as described
          in their respective documentation.")
  
  # Only run this once, not ever time function is called
  #python_version <- "3.11:latest"
  installed_envs <- reticulate::virtualenv_list()
  # Check if DupSeqR virtualenv already exists
  if ("DupSeqR" %in% installed_envs) {
    # Virtualenv exists, use it
    reticulate::use_virtualenv("DupSeqR")
  } else {
    # Virtualenv doesn't exist, set it up
    reticulate::virtualenv_create("DupSeqR", python = python_version)
    # Install required packages
    reticulate::virtualenv_install("DupSeqR", c("SigProfilerMatrixGenerator", "SigProfilerAssignment", "SigProfilerExtractor"))
  }
  
  # reticulate::install_python(version = python_version)
  reticulate::use_virtualenv("DupSeqR")
  
  # This will only install the genome if it's not found in the current env
  SigProfilerMatrixGeneratorR::install(project_genome)
  signatures_python_code <- system.file('extdata', 'signatures.py',
                                        package = "DupSeqR")
  reticulate::source_python(signatures_python_code)
  
  message("Creating cleaned data for input into SigProfiler...")
  # Clean data into required format for Alexandrov Lab tools...
  #ID doesn't always exist. 
  # New F(x) to create MAF/txt input for web sigprofiler
  signature_data <- as.data.frame(mutations) %>%
    dplyr::filter(.data$variation_type %in% "snv") %>%
    #dplyr::filter(.data$is_germline == FALSE) %>%   
    dplyr::select(all_of(group), .data$id, .data$variation_type, seqnames, .data$start, .data$end, .data$ref, .data$alt) %>%
    dplyr::rename(
      "Sample" = group,
      "ID" = "id",
      "mut_type" = "variation_type",
      "chrom" = "seqnames",
      "pos_start" = "start",
      "pos_end" = "end"
    ) %>%
    dplyr::mutate(Sample = stringr::str_replace_all(.data$Sample, " ", "_")) %>%
    dplyr::mutate(chrom = stringr::str_replace(.data$chrom, "chr", "")) %>%
    dplyr::mutate(Project = project) %>%
    dplyr::mutate(Genome = project_genome) %>%
    dplyr::mutate(Type = "SOMATIC") %>%
    dplyr::relocate(.data$Project) %>%
    dplyr::relocate(.data$Genome, .after = .data$ID) %>%
    dplyr::mutate(mut_type = "SNP") # This should be fixed before using on other datasets.
  
  message("Generating output path string...")
  output_path <- file.path(
    here::here(),
    "output",
    "SigProfiler",
    group
  )
  
  message(paste0("Creating directory ", output_path))
  if (!dir.exists(file.path(output_path, "matrices"))) {
    dir.create(file.path(output_path, "matrices"), recursive = T)
  }
  
  message("Writing mutation matrix to use as input to SigProfiler...")
  write.table(signature_data,
              file = file.path(output_path, "matrices", "mutations.txt"),
              sep = "\t", row.names = F, quote = F
  )
  
  
  message("Running SigProfilerMatrixGeneratorR...")
  signature_matrices <-
    SigProfilerMatrixGeneratorR::SigProfilerMatrixGeneratorR(
      project = project_name,
      genome = project_genome,
      matrix_path = file.path(output_path, "matrices"),
      plot = T, exome = F, bed_file = NULL,
      chrom_based = F, tsb_stat = T, seqInfo = T,
      cushion = 100
    )
  
  message("Running COSMIC fitting...")
  cosmic_fit_DupSeqR(
    samples = file.path(output_path, "matrices", "output", "SBS",
                        paste0(project_name, ".SBS96.all")),
    output = file.path(output_path, "matrices", "output"),
    genome_build=project_genome,
    cosmic_version=3.3,
    verbose=TRUE,
    exome=FALSE
  )
  
  # Errors here that I haven't been able to troubleshoot; maybe this shoudl be
  # a totally separate function anyway, putting on the backburner for now.
  # SigProfilerExtractorR::sigprofilerextractor(
  #   file.path(output_path, "matrices", "output", "SBS", paste0(project_name, ".SBS96.all")),
  #   file.path(output_path, "SigProfilerExtractor"),
  #   project_genome
  # )
}
