#' Calculate mutation frequency
#'
#' Calculates the mutation frequency for arbitrary groupings and creates a new
#' dataframe with the results. Mutation frequency is # mutations / total bases,
#' but this can be subset in different ways: e.g., by mutation context. In this
#' case, it is necessary to change the denominator of total bases to reflect
#' the sequencing depth at the proper reference bases under consideration.
#' 
#' Additionally, the operation is run by default using clonally expanded as well
#' as unique mutations.
#' @param data The data frame to be processed containing mutation data.
#' @param cols_to_group A vector of grouping variables: this should be the groups
#' of interest that you want to calculate a frequency for. For instance, getting
#' the frequency by `sample`. Other options might include `locus`, or,
#' `c("sample","locus")`. Must be a column in the mutation data table.
#' @param subtype_resolution At what resolution should the frequencies be
#' calculated? This is an important consideration because the reference bases
#' sequencing depth will differ depending on which mutation subtypes are being
#' considered. Options are "none", "type", base_6", "base_12", "base_96", and "base_192".
#' Could add 24 base option eventually. By selecting "none", this will use the
#' total depth across all the selected metadata columns.
#' @param vaf_cutoff Exclude rows from frequency calculations
#' with a variant allele fraction (VAF) greater than this value. Default is 0.1.
#' @param summary TRUE or FALSE, whether to return a summary table (i.e., where
#' only relevant columns for frequencies and groupings are returned). Setting 
#' this to false returns all columns in the original data, which might make
#' plotting more difficult, but may provide additional flexibility to power
#' users.
#' @param  clonality_cutoff NOT CURRENTLY IMPLEMENTED! Up for consideration.
#' This value determines the fraction of reads that
#' is considered a constitutional variant. If a mutation is present at a 
#' fraction higher than this value, the reference base will be swapped,
#' and the alt_depth recalculated. 0.3 (30%) would be a sane default?
#' @param variant_types Include these variant types. A vector of one or more
#'  "snv", "indel", "sv", "mnv", "no_variant"
#' @returns A data frame with the mutation frequency calculated.
#' @importFrom dplyr across all_of filter group_by mutate n row_number select distinct ungroup  
#' @importFrom magrittr %>%
#' @importFrom rlang := .data
#' @importFrom utils modifyList
#' @export
calculate_mut_freq <- function(data,
                               cols_to_group = c("sample", "description"),
                               subtype_resolution = "base_6",
                               vaf_cutoff = 0.1,
                              # clonality_cutoff = 0.3,
                               summary = TRUE,
                               variant_types = c("snv", "deletion", "insertion", "complex", "mnv","symbolic")) {
  
  if (subtype_resolution %in% c("base_6", "base_12", "base_96", "base_192") && !("snv" %in% variant_types)) {
    warning("Please include 'snv' in parameter 'variant_types' to calculate single-nucleotide variant subtype frequencies.")
  }
  
  freq_col_prefix <- paste(cols_to_group, collapse = "_")
  
  if (!subtype_resolution %in% names(DupSeqR::subtype_dict)) {
    stop(paste0("Error: you need to set subtype_resolution to one of: ",
                paste(names(DupSeqR::subtype_dict), collapse = ", ")))
  }

  numerator_groups <- c(cols_to_group, DupSeqR::subtype_dict[[subtype_resolution]])
  numerator_groups <- numerator_groups[!is.na(numerator_groups)]
  denominator_groups <- c(cols_to_group, DupSeqR::denominator_dict[[subtype_resolution]])
  denominator_groups <- denominator_groups[!is.na(denominator_groups)]
  # Check if data is provided as GRanges: if so, convert to data frame.
  if (inherits(data, "GRanges")) { data <- as.data.frame(data) }
  if (!inherits(data, "data.frame")) { warning("You should use a 
                                             data frame as input here.")}
 # Rename columns in data to default
   data <- rename_columns(data) 
# Un-duplicating the depth col
  # When there are +1 calls at the same position, modify total_depth such that
  # 1st call retains value, all other duplicates have total_depth = 0. 
  # Note that if users imported data using import_mut_file or read_vcf than total_depth
  # will be identical across all duplicates in a group (either take_del or take_mean)
  
 
  # TO DO : make into utlity function
  # Add parameter, default null
   mut_freq_table <- data %>%
    dplyr::group_by(.data$sample, .data$contig, .data$start) %>%
    dplyr::mutate(dup_id = row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(total_depth_unduplicated = ifelse(.data$dup_id > 1, 0, .data$total_depth)) %>%
    dplyr::select(-.data$total_depth, -.data$dup_id) %>%
    dplyr::rename(total_depth = total_depth_unduplicated)
  
    
 # Numerator groups
   mut_freq_table <- mut_freq_table %>% 
     dplyr::group_by(dplyr::across(dplyr::all_of(c(numerator_groups)))) %>%
    mutate(!!paste0(freq_col_prefix, "_sum_clonal") :=
        sum(.data$alt_depth[.data$variation_type %in% variant_types & .data$vaf < vaf_cutoff])) %>%
    mutate(!!paste0(freq_col_prefix, "_sum_unique") :=
             length(.data$alt_depth[.data$variation_type %in% variant_types & .data$vaf < vaf_cutoff])) %>%
    dplyr::ungroup()
    
  # Calculate denominator (same for clonal and unique mutations)
    # Calculate snv depth, which accounts for snv subtype resolution and 
   # group depth, which is the summed total_depth across the group for other variants. 
    mut_freq_table <- mut_freq_table %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(cols_to_group))) %>%
    mutate(!!paste0(freq_col_prefix, "_group_depth") := sum(.data$total_depth)) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(denominator_groups)))) %>%
    mutate(!!paste0(freq_col_prefix, "_snv_depth") := sum(.data$total_depth)) %>%
      dplyr::ungroup()
    
    # Calculate frequencies
    mut_freq_table <- mut_freq_table %>%
      mutate(!!paste0(freq_col_prefix, "_MF_clonal") := 
               ifelse(.data$variation_type == "snv", 
                      .data[[paste0(freq_col_prefix, "_sum_clonal")]] /
                      .data[[paste0(freq_col_prefix, "_snv_depth")]],
                        .data[[paste0(freq_col_prefix, "_sum_clonal")]] /
                        .data[[paste0(freq_col_prefix, "_group_depth")]])) %>%
    mutate(!!paste0(freq_col_prefix, "_MF_unique") :=
             ifelse(.data$variation_type == "snv", 
                    .data[[paste0(freq_col_prefix, "_sum_unique")]] /
                      .data[[paste0(freq_col_prefix, "_snv_depth")]],
                    .data[[paste0(freq_col_prefix, "_sum_unique")]] /
                      .data[[paste0(freq_col_prefix, "_group_depth")]])) %>%
    dplyr::ungroup()

  #######################################
# Create a summary table
    #We need to append the subtype list with the mnvs sv and indels as needed. 
  
if(subtype_resolution != "none"){    
  # Create a list of vectors, one for each column in cols_to_group
  col_values <- lapply(cols_to_group, function(col) unique(mut_freq_table[[col]]))
  # Grab and append list of subtypes

  if(subtype_resolution != "type"){
    subtype <- list(subtype_list[[subtype_resolution]])  
    subset_type <- subtype_list$type[subtype_list$type %in% variant_types]
    subset_type <- setdiff(subset_type, "snv")
    subtype[[1]] <- c(subtype[[1]], subset_type)
  } else {
    subtype <- list(subtype_list[[subtype_resolution]]) 
  }
  # Create a dataframe with every combination of subtype_resolution and values from cols_to_group
  summary_table <- do.call(expand.grid, c(muttype = subtype, col_values))
  # Set the column names
  col_names <- c(paste(DupSeqR::subtype_dict[[subtype_resolution]]), cols_to_group)
  colnames(summary_table) <- col_names
  } else {
  # Create a list of vectors, one for each column in cols_to_group
  col_values <- lapply(cols_to_group, function(col) unique(mut_freq_table[[col]]))
  # Create a dataframe with every combination of values from cols_to_group
  summary_table <- do.call(expand.grid, col_values)
  # Set the column names  
  col_names <- paste(cols_to_group)
  summary_table <- setNames(summary_table, col_names)
  }
    # Define columns of interest for summary table
    # Note if we want to include the snv depth, than the mnv and indels will be called twice. 
    # extra filtering afterwards?
  summary_cols <- c(
    numerator_groups,
    paste0(freq_col_prefix, "_sum_unique"),
    paste0(freq_col_prefix, "_sum_clonal"),
    paste0(freq_col_prefix, "_group_depth"),
    paste0(freq_col_prefix, "_snv_depth"),
    paste0(freq_col_prefix, "_MF_unique"),
    paste0(freq_col_prefix, "_MF_clonal")
  )

  # Make summary table of frequencies
  # This is also where subtype proportions are calculated
  summary_data <- mut_freq_table %>%
    dplyr::select({{ summary_cols }})  %>%
    dplyr::distinct(across(all_of(numerator_groups)), .keep_all = TRUE)
   
 # TO DO:
  # Fill in 0s with their snv depths (mut types that did not appear in the data)
summary_table <- dplyr::left_join(summary_table, summary_data, by = col_names)
summary_table[is.na(summary_table)] <- 0
summary_table <- summary_table %>% 
  dplyr::group_by(across(all_of(cols_to_group))) %>%
  dplyr::mutate(!!paste0(freq_col_prefix, "_group_depth") :=
                  ifelse(.data[[paste0(freq_col_prefix, "_group_depth")]] == 0, 
                         max(.data[[paste0(freq_col_prefix, "_group_depth")]]), 
                         .data[[paste0(freq_col_prefix, "_group_depth")]])) %>%
  dplyr::ungroup() %>%
  
     mutate(freq_clonal =
             .data[[paste0(freq_col_prefix, "_sum_clonal")]] /
             sum(.data[[paste0(freq_col_prefix, "_sum_clonal")]]) /
             .data[[paste0(freq_col_prefix, "_group_depth")]] ) %>%
    mutate(prop_clonal = .data$freq_clonal / sum(.data$freq_clonal)) %>%
    mutate(freq_unique =
             .data[[paste0(freq_col_prefix, "_sum_unique")]] /
             sum(.data[[paste0(freq_col_prefix, "_sum_unique")]]) /
             .data[[paste0(freq_col_prefix, "_group_depth")]] ) %>%
    mutate(prop_unique = .data$freq_unique / sum(.data$freq_unique))

  if (!summary) {
    return(mut_freq_table)
  } else if (summary) {
    return(summary_table)
  } else {
    stop("Error: summary must be TRUE or FALSE.")
  }
}

# For testing
#cols_to_group <- c("sample")
#cols_to_group <- c("sample","description")
#cols_to_group <- c("sample", "CpG_site")

# To do... locate and enumerate recurrent mutations?
