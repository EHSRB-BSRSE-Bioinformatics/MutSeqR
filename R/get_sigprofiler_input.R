#' Write the mutation calling file to input into the SigProfiler Assignment
#' web application.
#'
#' Creates a .txt file from mutation data that can be used for mutational
#' signatures analysis using the SigProfiler Assignment web application.
#' Currently only supports SBS analysis i.e. snvs.
#'
#' @param mutation_data The object containing the mutation data.
#' The output of import_mut_data() or import_vcf_data().
#' @param project_name The name of the project. Default is "Example".
#' @param project_genome The reference genome to use.
#'  (e.g., Human: GRCh38, Mouse mm10: GRCm38)
#' @param output_path The path to save the output file. If NULL, files
#' will be saved in the current working directory. Default is NULL.
#' @returns a .txt file that can be uploaded to the SigProfiler Assignment web
#' application (https://cancer.sanger.ac.uk/signatures/assignment/) as a
#' "Mutational calling file".
#' @details Mutations will be be filtered for SNVs. Mutations flagged in
#' `filter_mut` will be excluded from the output.
#' @importFrom dplyr rename filter select mutate relocate
#' @examples
#' example_file <- system.file("extdata", "example_mutation_data_filtered.rds", package = "MutSeqR")
#' example_data <- readRDS(example_file)
#' temp_output <- tempdir()
#' write_mutation_calling_file(mutation_data = example_data,
#'                             project_name = "Example",
#'                             project_genome = "GRCm38",
#'                             output_path = temp_output)
#' list.files(temp_output)
#' # The file is saved in the temporary directory
#' # To view the file, use the following code:
#' ## output_file <- file.path(temp_output, "mutation_calling_file.txt")
#' ## file.show(output_file
#' @importFrom dplyr rename filter select mutate relocate
#' @importFrom here here
#' @importFrom utils write.table
#' @export

write_mutation_calling_file <- function(mutation_data,
                                        project_name = "Example",
                                        project_genome = "GRCm38",
                                        output_path = NULL) {

  # Check if data is provided as GRanges: if so, convert to data frame.
  if (inherits(mutation_data, "GRanges")) {
    mutation_data <- as.data.frame(mutation_data)
    mutation_data <- mutation_data %>%
      dplyr::rename(contig = "seqnames")
  }
  if (!inherits(mutation_data, "data.frame")) {
    warning("You should use a data frame as input here.")
  }

  signature_data <- mutation_data %>%
    dplyr::filter(.data$variation_type == "snv") %>%
    dplyr::filter(.data$filter_mut == FALSE)

  if ("id" %in% colnames(signature_data)) {
    signature_data <- signature_data %>%
      dplyr::select("sample", "id",
                    "variation_type",
                    "contig", "start",
                    "end", "ref", "alt") %>%
      dplyr::rename("Sample" = "sample",
        "ID" = "id",
        "mut_type" = "variation_type",
        "chrom" = "contig",
        "pos_start" = "start",
        "pos_end" = "end"
      ) %>%
      dplyr::mutate(Project = project_name) %>%
      dplyr::mutate(Sample = stringr::str_replace_all(.data$Sample, " ", "_")) %>%
      dplyr::mutate(chrom = stringr::str_replace(.data$chrom, "chr", "")) %>%
      dplyr::mutate(Genome = project_genome) %>%
      dplyr::relocate(.data$Project) %>%
      dplyr::relocate(.data$Genome, .after = .data$ID) %>%
      dplyr::mutate(mut_type = "SNP") # This should be fixed before using on other datasets.
  } else {
    signature_data <- signature_data %>%
      dplyr::select("sample", "variation_type", "contig",
                    "start", "end", "ref", "alt") %>%
      dplyr::rename(
        "Sample" = "sample",
        "mut_type" = "variation_type",
        "chrom" = "contig",
        "pos_start" = "start",
        "pos_end" = "end"
      ) %>%
      dplyr::mutate(Project = project_name) %>%
      dplyr::mutate(Sample = stringr::str_replace_all(.data$Sample, " ", "_")) %>%
      dplyr::mutate(ID = ".") %>%
      dplyr::mutate(chrom = stringr::str_replace(.data$chrom, "chr", "")) %>%
      dplyr::mutate(Genome = project_genome) %>%
      dplyr::relocate(.data$Project) %>%
      dplyr::relocate(.data$ID, .after = .data$Sample) %>%
      dplyr::relocate(.data$Genome, .after = .data$ID) %>%
      dplyr::mutate(mut_type = "SNP")
  }

  if (is.null(output_path)) {
    output_path <- file.path(
      here::here(),
      "output",
      "SigProfiler_input"
    )
  }

  if (!dir.exists(output_path)) {
    dir.create(file.path(output_path), recursive = TRUE)
  }
  utils::write.table(signature_data,
    file = file.path(output_path, "mutation_calling_file.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

#' Write a Mutational Matrix to input into the sigprofiler web application
#'
#' Creates a .txt file from mutation data that can be used for mutational
#' signatures analysis using the SigProfiler web application. Can handle
#' group analyses (ex dose, tissue, etc). Currently only supports SBS matrices
#' i.e. snvs.
#'
#' @param mutation_data The object containing the mutation data.
#' The output of import_mut_data() or import_vcf_data().
#' @param group The column in the mutation data used to aggregate groups
#' (e.g., sample, tissue, dose).
#' @param subtype_resolution The resolution of the mutation subtypes. Options
#' are "base_6" or "base_96". Default is "base_96".
#' @param mf_type The mutation counting method to use. Options are "min" or
#' "max". Default is "min".
#' @param output_path The path to save the output file. If not provided, the
#' file will be saved in the current working directory. Default is NULL.
#' @returns a .txt file that can be uploaded to the SigProfiler Assignment web
#' application (https://cancer.sanger.ac.uk/signatures/assignment/) as
#' a "Mutational Matrix".
#' @details Mutations will be be filtered for SNVs. Mutations flagged in
#' `filter_mut` will be excluded from the output. Mutations will be summed
#' across the groups specified in the `group` argument.
#' @examples
#' example_file <- system.file("extdata", "example_mutation_data_filtered.rds", package = "MutSeqR")
#' example_data <- readRDS(example_file)
#' temp_output <- tempdir()
#' write_mutational_matrix(mutation_data = example_data,
#'                         group = "dose_group",
#'                         subtype_resolution = "base_96",
#'                         mf_type = "min",
#'                         output_path = temp_output)
#' list.files(temp_output)
#' # The file is saved in the temporary directory
#' # To view the file, use the following code:
#' ## output_file <- file.path(temp_output, "dose_group_base_96_mutational_matrix.txt")
#' ## file.show(output_file)
#' @importFrom stats reshape
#' @importFrom dplyr rename filter group_by mutate ungroup select distinct
#' @importFrom here here
#' @export
#'
write_mutational_matrix <- function(mutation_data,
                                    group = "dose",
                                    subtype_resolution = "base_96",
                                    mf_type = "min",
                                    output_path = NULL) {

  if (!subtype_resolution %in% c("base_6", "base_96")) {
    stop("The subtype_resolution argument must be either 'base_6' or 'base_96'")
  }
  mut_matrix <- suppressWarnings(calculate_mf(mutation_data,
                                              cols_to_group = group,
                                              subtype_resolution = subtype_resolution,
                                              variant_types = "snv",
                                              calculate_depth = FALSE,
                                              precalc_depth_data = NULL))
  mut_matrix <- mut_matrix %>%
    dplyr::rename(mut_count = paste0("sum_", mf_type),
      MutationType = MutSeqR::subtype_dict[[subtype_resolution]]
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Group = stringr::str_c(!!!syms(group), sep = "_")) %>%
    dplyr::ungroup() %>%
    dplyr::select("MutationType", "Group", "mut_count")

  mut_matrix_wide <- tidyr::pivot_wider(mut_matrix,
                                        names_from = .data$Group,
                                        values_from = .data$mut_count)
  
  colnames(mut_matrix_wide) <- ifelse(colnames(mut_matrix_wide) == "MutationType",
                                      "MutationType",
                                      paste0("Group_", colnames(mut_matrix_wide)))
  if (is.null(output_path)) {
    output_path <- file.path(
      here::here(),
      "SigProfiler_input"
    )
  }

  if (!dir.exists(output_path)) {
    dir.create(file.path(output_path), recursive = TRUE)
  }
  filename <- paste0(group, "_", subtype_resolution, "_", "mutational_matrix.txt")
  utils::write.table(mut_matrix_wide,
    file = file.path(output_path, filename),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}