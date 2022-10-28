# For this to work there are a few dependenices that must be met.
# You must first install python dependencies.
# e.g., on the linux command line:
# pip install SigProfilerAssignment
# pip install SigProfilerMatrixGenerator
# pip install SigProfilerExtractor
# ...other python tools from Alexandrov lab
# You must also install a reference genome:
# SigProfilerMatrixGeneratorR::install('GRCh38') # Don't run this every time...


#' Run COSMIC signatures comparison
#'
#' After cleaning the mutation data input, runs several Alexandrov Lab tools for COSMIC signature analysis (assigns signatures to best explain the input data).
#' @param mutation_data A data frame imported from a .mut file
#' @param project_name The name of the project; used to get mutation data into the required .txt format for SigProfiler
#' @param project_genome A string describing the reference genome to use; e.g., GRCh38
#' @param output_path The directory where output results should be written
#' @param group The column in the mutation data used to aggregate groups (e.g., sample ID, tissue, dose)
#' @param python_path The path to the version of python to be used with reticulate. It is important that this version of python meets the dependecies, including the SigProfiler python tools.
#' @returns Creates a subfolder in the output directory with SigProfiler tools results.
#' @export
signature_decomposition <- function(mutations = mutation_data,
                                    project_name = "Default",
                                    project_genome = "GRCh38",
                                    group = "sample",
                                    python_path = "/usr/bin/python3.9") {
  if (!require(reticulate)) {
    stop("reticulate not installed")
  }
  if (!require(SigProfilerMatrixGeneratorR)) {
    stop("SigProfilerMatrixGeneratorR not installed")
  }
  if (!require(tidyverse)) {
    stop("tidyverse not installed")
  }

  # reticulate::use_python(python_path)
  # reticulate::py_config()

  # Clean data into required format for Alexandrov Lab tools...
  signature_data <- as.data.frame(mutations) %>%
    dplyr::select(!!ensym(group), id, variation_type, seqnames, start, end, ref, alt) %>%
    dplyr::rename(
      "Sample" = group,
      "ID" = "id",
      "mut_type" = "variation_type",
      "chrom" = "seqnames",
      "pos_start" = "start",
      "pos_end" = "end"
    ) %>%
    dplyr::mutate(Sample = stringr::str_replace_all(Sample, " ", "_")) %>%
    dplyr::mutate(chrom = stringr::str_replace(chrom, "chr", "")) %>%
    dplyr::mutate(Project = project_name) %>%
    dplyr::mutate(Genome = project_genome) %>%
    dplyr::mutate(Type = "SOMATIC") %>%
    dplyr::relocate(Project) %>%
    dplyr::relocate(Genome, .after = ID) %>%
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

  signature_matrices <- SigProfilerMatrixGeneratorR(
    project = project_name,
    genome = project_genome,
    matrix_path = file.path(output_path, "matrices"),
    plot = T, exome = F, bed_file = NULL,
    chrom_based = F, tsb_stat = T, seqInfo = T,
    cushion = 100
  )

  source_python(file.path(here::here(), "inst", "signatures.py"))

  cosmic_fitR(
    file.path(output_path, "matrices", "output", "SBS", paste0(project_name, ".SBS96.all")),
    output_path
  )

  sigProfilerExtractorR(
    file.path(output_path, "matrices", "output", "SBS", paste0(project_name, ".SBS96.all")),
    file.path(output_path, "SigProfilerExtractor"),
    project_genome
  )
}
