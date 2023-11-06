# For this to work there are a few dependenices that must be met.
# You must first install python dependencies.
# e.g., on the linux command line:
# pip install SigProfilerAssignment
# pip install SigProfilerMatrixGenerator
# pip install SigProfilerExtractor
# ...other python tools from Alexandrov lab
# You must also install a reference genome:
# SigProfilerMatrixGeneratorR::install('GRCh38') # Don't run this every time...
# TO DO: sigporfiler packages are in suggests; write a feature that will ask users if they want to install them. 

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
signature_decomposition <- function(mutations,
                                    project_name = "Default",
                                    project_genome = "GRCh38",
                                    group = "sample",
                                    python_path = "~/../../AppData/Local/Programs/Python/Python310/python.exe", #"/usr/bin/python3.9", #
                                    python_home = "~/.virtualenvs/r-reticulate") {
  if (!requireNamespace(reticulate)) {
    stop("reticulate not installed")
  }
  message("This function requires python to be installed on whichever
          platform you are using and to be available on the command line.
          Note that it will also install several python dependencies using
          a conda virtual environment on first use. Please be aware of the 
          implications of this. For advanced use, it is suggested to
          use the SigProfiler python tools directly in python as described
          in their respective documentation.")
  
  # TODO - only run this once, not ever time function is called
  # Use an if statement to determine if the dependencies are met already
  reticulate::install_python(version = "3.9:latest")
  reticulate::virtualenv_create("DupSeqR", python = "python3.9")
  reticulate::virtualenv_install("DupSeqR", "SigProfilerAssignment")
  reticulate::virtualenv_install("DupSeqR", "SigProfilerMatrixGenerator")
  reticulate::virtualenv_install("DupSeqR", "SigProfilerExtractor")
  reticulate::use_virtualenv("DupSeqR")
  
  # message(paste0("Creating a folder to store python dependencies at ",
  #                python_home, ". This can be avoided by setting it manually
  #                to a location of your choosing."))
  # 
  # if (!dir.exists(python_home)) { dir.create(python_home) }
  # 
  # options(reticulate.conda_binary = python_home)
  
  #Sys.setenv(RETICULATE_PYTHON = python_path)
  #Sys.setenv(RETICULATE_PYTHON =  py_config()$python)
  #Sys.setenv(RETICULATE_PYTHON_ENV =  py_config()$python)
  #reticulate::use_python(python_path)
  #reticulate::py_config()
  #cat(paste(py_config()))
  #use_virtualenv("~/.virtualenvs/r-reticulate/")
  #use_python("~/.virtualenvs/r-reticulate/Scripts/python.exe") 
  #cat(paste(py_config()))
  #have_SigProfilerAssignment <- py_module_available("SigProfilerAssignment")
  #have_SigProfilerExtractor <- py_module_available("SigProfilerExtractor")
  #have_SigProfilerMatrixGenerator <- py_module_available("SigProfilerMatrixGenerator")
  
  #if (!have_SigProfilerAssignment) {
    #reticulate::py_install("SigProfilerAssignment", pip = F) 
  #}
  #if (!have_SigProfilerExtractor) {
    #reticulate::py_install("SigProfilerExtractor", pip = F) 
  #}
  #if (!have_SigProfilerMatrixGenerator) {
    #reticulate::py_install("SigProfilerMatrixGenerator", pip = F) 
  #}

  # Clean data into required format for Alexandrov Lab tools...
  signature_data <- as.data.frame(mutations) %>%
    dplyr::filter(!.data$variation_type %in% "no_variant") %>%
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
    dplyr::mutate(Project = project_name) %>%
    dplyr::mutate(Genome = project_genome) %>%
    dplyr::mutate(Type = "SOMATIC") %>%
    dplyr::relocate(.data$Project) %>%
    dplyr::relocate(.data$Genome, .after = .data$ID) %>%
    dplyr::mutate(mut_type = "SNP") # This should be fixed before using on other datasets.

  output_path <- file.path(
    here::here(),
    "output",
    "SigProfiler",
    group
  )

  if (!dir.exists(output_path)) {
    dir.create(file.path(output_path, "matrices"), recursive = T)
  }

  write.table(signature_data,
    file = file.path(output_path, "matrices", "mutations.txt"),
    sep = "\t", row.names = F, quote = F
  )
  
  signatures_python_code <- system.file('extdata', 'signatures.py',
                                        package = "DupSeqR")
  #reticulate::use_python(Sys.which("python"))
  #use_python("C:/Users/MAMEIER/OneDrive - HC-SC PHAC-ASPC/Documents/.virtualenvs/r-reticulate/Scripts/python.exe")
  #reticulate::use_virtualenv(
  #  "C:/Users/MAMEIER/OneDrive - HC-SC PHAC-ASPC/Documents/.virtualenvs/r-reticulate")
  #reticulate::use_condaenv("myenv")
  
  reticulate::source_python(signatures_python_code)
  
  # only do this once as well?
  install_genome(project_genome)
  signature_matrices <-
    SigProfilerMatrixGeneratorR::SigProfilerMatrixGeneratorR(
      project = project_name,
      genome = project_genome,
      matrix_path = file.path(output_path, "matrices"),
      plot = T, exome = F, bed_file = NULL,
      chrom_based = F, tsb_stat = T, seqInfo = T,
      cushion = 100
    )

  

  cosmic_fit_DupSeqR(
    samples = file.path(output_path, "matrices", "output", "SBS",
                        paste0(project_name, ".SBS96.all")),
    genome_build="GRCh38",
    cosmic_version=3.3,
    verbose=True,
    exome=False
  )

  SigProfilerExtractorR::sigprofilerextractor(
    file.path(output_path, "matrices", "output", "SBS", paste0(project_name, ".SBS96.all")),
    file.path(output_path, "SigProfilerExtractor"),
    project_genome
  )
}
