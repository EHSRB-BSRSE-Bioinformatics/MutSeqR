#' Write the mutation calling file to input into the sigprofiler web application
#' 
#' Creates a .txt file from mutation data that can be used for mutational signatures 
#' analysis using the SigProfiler web application.Cannot group higher than sample. 
#' 
#' @param mutations The GRanges object containing the mutation data. 
#' The output of import_mut_data() or read_vcf().
#' @param project_name The name of the project; used to get mutation data 
#' into the required .txt format for SigProfiler
#' @param project_genome A string describing the reference genome to use
#'  (e.g., GRCh38, GRCm38)
#' See import_mut_data() or read_vcf() for more details. 
#' @returns a .txt file that can be uploaded to the sigprofiler web application
#' as "Mutational calling file"
#' Filters out ostensibly germline mutations identified in mutaiton data.. 
#' @export

write_mutation_calling_file <- function(mutations,
                                  project_name = "test",
                                  project_genome = "GRCm38"
                                  ) {
  
  # Check if data is provided as GRanges: if so, convert to data frame.
  if (inherits(mutations, "GRanges")) { 
    mutations <- as.data.frame(mutations)
    mutations <- mutations %>%
      dplyr::rename(contig = contig)}
  if (!inherits(mutations, "data.frame")) { warning("You should use a 
                                             data frame as input here.")}
  
  signature_data <- mutations %>%
    dplyr::filter(.data$variation_type %in% "snv") %>%
    dplyr::filter(is_germline == FALSE)

    
  if ("id" %in% colnames(signature_data)) {
  signature_data <- signature_data %>% 
    dplyr::select(.data$sample, .data$id, .data$variation_type, contig, .data$start, .data$end, .data$ref, .data$alt) %>%
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
      dplyr::select(.data$sample, .data$variation_type, contig, .data$start, .data$end, .data$ref, .data$alt) %>%
      dplyr::rename(
        "Sample" = sample,
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
  
  output_path <- file.path(
    here::here(),
    "output",
    "SigProfiler",
    "sample"
  )
  
  if (!dir.exists(output_path)) {
    dir.create(file.path(output_path, "matrices"), recursive = T)
  }
  
  write.table(signature_data,
              file = file.path(output_path, "matrices", "mutation_calling_file.txt"),
              sep = "\t", row.names = F, quote = F
  )
  
}

#' Write a Mutational Matrix to input into the sigprofiler web application
#' 
#' Creates a .txt file from mutation data that can be used for mutational signatures 
#' analysis using the SigProfiler web application. Can handle group analyses
#' (ex dose, tissue, etc). Currently only supports snv. 
#' 
#' @param mutations The GRanges object containing the mutation data. 
#' The output of import_mut_data() or read_vcf().
#' @param group The column in the mutation data used to aggregate groups 
#' (e.g., sample, tissue, dose)
#' @param matrices At what resolution should the snv mutation counts be
#' calculated? Options are "base_96", 1536?, 384?, 6144?, DINUC? 6, 24, INDEL
#' TO DO: add more matrices options. 
#' @param include_clonal Should the mutational matrix include possibly clonally expanded mutations. 
#' If FALSE, any mutation occuring more than once in a single sample will only be counted once.
#' TO DO: clonality cut_off? 
#' @param filter Parameter allows you to choose to filter for only somatic or germline mutations. 
#' Values = c("somatic", "germline", "none). "none" will leave the data unfiltered. 
#' #' @returns a .txt file that can be uploaded to the sigprofiler web application as "Mutational Matrix"
#' @export
#' 
write_mutational_matrix <- function(mutations,  
                                    group = "dose",
                                    matrices = "base_96",
                                    include_clonal = FALSE,
                                    filter = "somatic") {
  
  
  if (inherits(mutations, "GRanges")) { 
    mutations <- as.data.frame(mutations)
    mutations <- mutations %>%
      dplyr::rename(contig = contig)}
  if (!inherits(mutations, "data.frame")) { warning("You should use a 
                                             data frame as input here.")}
signature_data <- mutations %>%
  dplyr::filter(.data$variation_type %in% "snv")

# Remove prefix from metadata columns
# Identify and rename sample_data columns
columns_with_sample_data_prefix <- grep("^sample_data_", colnames(signature_data), value = TRUE)
columns_with_region_data_prefix <- grep("^region_data_", colnames(signature_data), value = TRUE)
# Remove the prefixes to get the clean column names
stripped_sample_data_columns <- sub("^sample_data_", "", columns_with_sample_data_prefix)
stripped_region_data_columns <- sub("^region_data_", "", columns_with_region_data_prefix)
# Apply new column names to signature_data
colnames(signature_data)[colnames(signature_data) %in% columns_with_sample_data_prefix] <- stripped_sample_data_columns
colnames(signature_data)[colnames(signature_data) %in% columns_with_region_data_prefix] <- stripped_region_data_columns
    
  if (filter == "somatic") {
      signature_data <- filter(signature_data, is_germline == FALSE)
    } else if (filter == "germline") {
      signature_data <- filter(signature_data, is_germline == TRUE)
    } else if (filter == "none") {
       signature_data <- signature_data
      }
  
  # Numerator groups
  numerator_groups <- c(group, DupSeqR::subtype_dict[[matrices]])
  numerator_groups <- numerator_groups[!is.na(numerator_groups)]
  
  signature_data <- signature_data %>% 
    dplyr::group_by(dplyr::across(dplyr::all_of(c(numerator_groups)))) %>%
    mutate(mut_count = ifelse(include_clonal, 
                        sum(.data$alt_depth), 
                        length(.data$alt_depth))) %>%
    dplyr::ungroup()
  
  # Create a list of vectors, one for each column in cols_to_group
  col_values <- lapply(group, function(col) unique(signature_data[[col]]))
  # Create a dataframe with every combination of subtype_resolution and values from cols_to_group
  mut_matrix <- do.call(expand.grid, c(subtype = list(subtype_list[[matrices]]), col_values))
  
  col_names <- c(paste(DupSeqR::subtype_dict[[matrices]]), group)
  
  colnames(mut_matrix) <- col_names 
  
  summary_cols <- c(numerator_groups,"mut_count")
  
  summary_data <- signature_data %>%
    dplyr::select({{ summary_cols }})  %>%
    dplyr::distinct()
  mut_matrix <- merge(mut_matrix, summary_data, by = col_names, all = TRUE)
  mut_matrix$mut_count[is.na(mut_matrix$mut_count)] <- 0
  mut_matrix_wide <- reshape(mut_matrix, idvar = paste(DupSeqR::subtype_dict[[matrices]]), timevar = paste(group), direction = "wide")
  colnames(mut_matrix_wide) <- gsub("mut_count.", "", colnames(mut_matrix_wide)) 
  colnames(mut_matrix_wide)[1] <- "MutationType"
   output_path <- file.path(
    here::here(),
    "output",
    "SigProfiler",
    group
  )
  
  if (!dir.exists(output_path)) {
    dir.create(file.path(output_path, "matrices"), recursive = T)
  }
  
  write.table(mut_matrix_wide,
              file = file.path(output_path, "matrices", "mutational_matrix.txt"),
              sep = "\t", row.names = F, quote = F
  )  
}