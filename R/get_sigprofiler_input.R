#' Write the mutation calling file to input into the sigprofiler web application
#' 
#' Creates a .txt file from mutation data that can be used for mutational signatures 
#' analysis using the SigProfiler web application.Cannot group higher than sample. 
#' 
#' @param mutation_data The GRanges object containing the mutation data.
#' The output of import_mut_data() or read_vcf().
#' @param project_name The name of the project; used to get mutation data
#' into the required .txt format for SigProfiler
#' @param project_genome A string describing the reference genome to use
#'  (e.g., GRCh38, GRCm38)
#' @param output_path The path to save the output file. Default is NULL.
#' See import_mut_data() or read_vcf() for more details.
#' @returns a .txt file that can be uploaded to the sigprofiler web application
#' as "Mutational calling file"
#' Filters out ostensibly germline mutations identified in mutation data.
#' @importFrom dplyr rename filter select mutate relocate
#' @importFrom here here
#' @importFrom utils write.table
#' @export

write_mutation_calling_file <- function(mutation_data,
                                  project_name = "test",
                                  project_genome = "GRCm38",
                                  output_path = NULL
                                  ) {
  
  # Check if data is provided as GRanges: if so, convert to data frame.
  if (inherits(mutation_data, "GRanges")) { 
    mutation_data <- as.data.frame(mutation_data)
    mutation_data <- mutation_data %>%
      dplyr::rename(contig = "seqnames")}
  if (!inherits(mutation_data, "data.frame")) { warning("You should use a 
                                             data frame as input here.")}
  
  signature_data <- mutation_data %>%
    dplyr::filter(.data$variation_type %in% "snv") %>%
    dplyr::filter(.data$is_germline == FALSE)

    
  if ("id" %in% colnames(signature_data)) {
  signature_data <- signature_data %>% 
    dplyr::select("sample", "id", "variation_type", "contig", "start", "end", "ref", "alt") %>%
    dplyr::rename(
      "Sample" = "sample",
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
      dplyr::select("sample", "variation_type", "contig", "start", "end", "ref", "alt") %>%
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
      dplyr::mutate(mut_type = "SNP") # This should be fixed before using on other datasets. 
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
#' Creates a .txt file from mutation data that can be used for mutational signatures 
#' analysis using the SigProfiler web application. Can handle group analyses
#' (ex dose, tissue, etc). Currently only supports snv. 
#' 
#' @param mutation_data The GRanges object containing the mutation data. 
#' The output of import_mut_data() or read_vcf().
#' @param group The column in the mutation data used to aggregate groups 
#' (e.g., sample, tissue, dose)
#' @param matrices At what resolution should the snv mutation counts be
#' calculated? Options are "base_96", 1536?, 384?, 6144?, DINUC? 6, 24, INDEL
#' TO DO: add more matrices options.
#' @param mf_type Mutation counting method. Options are "min" or "max".
#' Minimum Independent Frequency (min): All identical mutations within a
#' sample are assumed to be the result of clonal expansion and are thus
#' only counted once. Maximum Independent Frequency (max): All identical
#' mutations within a sample are assumed to be idenpendant mutational evens
#' and are included in the mutation frequency calculation. Note that this
#' does not apply for germline variants. Default is "min".
#' TO DO: clonality cut_off?
#' @param filter Parameter allows you to choose to filter for only somatic
#' or germline mutations.
#' @param output_path The path to save the output file. Default is NULL,
#' which will save the file in the current working directory. 
#' Values = c("somatic", "germline", "none). "none" will leave the data unfiltered. 
#' @returns a .txt file that can be uploaded to the sigprofiler web application
#' as "Mutational Matrix"
#' @importFrom stats reshape
#' @importFrom dplyr rename filter group_by mutate ungroup select distinct
#' @importFrom here here
#' @export
#' 
write_mutational_matrix <- function(mutation_data,
                                    group = "dose",
                                    matrices = "base_96",
                                    mf_type = "min",
                                    filter = "somatic",
                                    output_path = NULL) {
  
  
  if (!mf_type %in% c("min", "max")) {
    stop("Error: mf_type must be either 'min' or 'max'")
  }
  if (!filter %in% c("somatic", "germline", "none")) {
  stop("Error: filter must be either 'somatic', 'germline', 'none'")
}

  if (inherits(mutation_data, "GRanges")) { 
    mutation_data <- as.data.frame(mutation_data)
    mutation_data <- mutation_data %>%
      dplyr::rename(contig = "seqnames")}
  if (!inherits(mutation_data, "data.frame")) { warning("You should use a 
                                             data frame as input here.")}
  signature_data <- mutation_data %>%
    dplyr::filter(.data$variation_type %in% "snv")
    
  if (filter == "somatic") {
      signature_data <- dplyr::filter(signature_data, .data$is_germline == FALSE)
    } else if (filter == "germline") {
      signature_data <- dplyr::filter(signature_data, .data$is_germline == TRUE)
    } else if (filter == "none") {
       signature_data <- signature_data
      }
  
  # Numerator groups
  numerator_groups <- c(group, MutSeqR::subtype_dict[[matrices]])
  numerator_groups <- numerator_groups[!is.na(numerator_groups)]
  
  signature_data <- signature_data %>% 
    dplyr::group_by(dplyr::across(dplyr::all_of(c(numerator_groups)))) %>%
    dplyr::mutate(mut_count = ifelse(mf_type == "max",
                        sum(.data$alt_depth),
                        length(.data$alt_depth))) %>%
    dplyr::ungroup()
  
  # Create a list of vectors, one for each column in cols_to_group
  col_values <- lapply(group, function(col) unique(signature_data[[col]]))
  # Create a dataframe with every combination of subtype_resolution and values from cols_to_group
  mut_matrix <- do.call(expand.grid, c(subtype = list(MutSeqR::subtype_list[[matrices]]), col_values))
  
  col_names <- c(paste(MutSeqR::subtype_dict[[matrices]]), group)
  
  colnames(mut_matrix) <- col_names 
  
  summary_cols <- c(numerator_groups,"mut_count")
  
  summary_data <- signature_data %>%
    dplyr::select({{ summary_cols }})  %>%
    dplyr::distinct()
  mut_matrix <- merge(mut_matrix, summary_data, by = col_names, all = TRUE)
  mut_matrix$mut_count[is.na(mut_matrix$mut_count)] <- 0
  mut_matrix_wide <- stats::reshape(mut_matrix, idvar = paste(MutSeqR::subtype_dict[[matrices]]), timevar = paste(group), direction = "wide")
  colnames(mut_matrix_wide) <- gsub("mut_count.", "", colnames(mut_matrix_wide)) 
  colnames(mut_matrix_wide)[1] <- "MutationType"
   
   if (is.null(output_path)) {
    output_path <- file.path(
      here::here(),
      "output",
      "SigProfiler_input",
      group
    )
   }
  
  if (!dir.exists(output_path)) {
    dir.create(file.path(output_path), recursive = TRUE)
  }
  
  utils::write.table(mut_matrix_wide,
              file = file.path(output_path, "mutational_matrix.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE
  )  
}