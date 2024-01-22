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
#' @param variant_types Include these variant types. A vector of one or more
#'  "snv", "complex", "deletion", "insertion", "mnv", "symbolic", "no_variant".
#'  Default includes all variants. 
#' @param vaf_cutoff Exclude rows from frequency calculations
#' with a variant allele fraction (VAF) greater than this value. No default provided.
#' @param summary TRUE or FALSE, whether to return a summary table (i.e., where
#' only relevant columns for frequencies and groupings are returned). Setting 
#' this to false returns all columns in the original data, which might make
#' plotting more difficult, but may provide additional flexibility to power
#' users.
#' @param retain_metadata_cols list the metadata columns that you would like to 
#' retain in the summary table. This may be useful for plotting your summary data.
#' @param  clonality_cutoff NOT CURRENTLY IMPLEMENTED! Up for consideration.
#' This value determines the fraction of reads that
#' is considered a constitutional variant. If a mutation is present at a 
#' fraction higher than this value, the reference base will be swapped,
#' and the alt_depth recalculated. 0.3 (30%) would be a sane default?

#' @returns A data frame with the mutation frequency calculated.
#' @importFrom dplyr across all_of filter group_by mutate n row_number select distinct ungroup  
#' @importFrom magrittr %>%
#' @importFrom rlang := .data
#' @importFrom utils modifyList
#' @export
calculate_mut_freq <- function(data,
                               cols_to_group = c("sample", "description"),
                               subtype_resolution = "base_96",
                               variant_types = c("snv", "deletion", "insertion", "complex", "mnv","symbolic"),
                               vaf_cutoff,
                              #clonality_cutoff = 0.3,
                               summary = TRUE,
                               retain_metadata_cols = NULL
                               ) {
  
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
   mut_freq_table <- data %>%
    dplyr::group_by(.data$sample, .data$contig, .data$start) %>%
    dplyr::mutate(dup_id = row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(total_depth_unduplicated = ifelse(.data$dup_id > 1, 0, .data$total_depth)) %>%
    dplyr::select(-.data$total_depth, -.data$dup_id) %>%
    dplyr::rename(total_depth = total_depth_unduplicated)
  
 # Calculate denominator Numerator groups
   mut_freq_table <- mut_freq_table %>% 
     dplyr::group_by(dplyr::across(dplyr::all_of(c(numerator_groups)))) %>%
    dplyr::mutate(!!paste0(freq_col_prefix, "_sum_clonal") :=
        sum(.data$alt_depth[.data$variation_type %in% variant_types & .data$vaf < vaf_cutoff])) %>%
    dplyr::mutate(!!paste0(freq_col_prefix, "_sum_unique") :=
             length(.data$alt_depth[.data$variation_type %in% variant_types & .data$vaf < vaf_cutoff])) %>%
    dplyr::ungroup()
    
  # Calculate denominator (same for clonal and unique mutations)
    # snv depth: depth across groups and snv subtype resolution
    # group depth: depth across groups. 
    mut_freq_table <- mut_freq_table %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(cols_to_group))) %>%
    dplyr::mutate(!!paste0(freq_col_prefix, "_group_depth") := sum(.data$total_depth)) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(denominator_groups)))) %>%
    dplyr::mutate(!!paste0(freq_col_prefix, "_snv_depth") := sum(.data$total_depth)) %>%
      dplyr::ungroup()
    
    # Calculate frequencies
    mut_freq_table <- mut_freq_table %>%
      dplyr::mutate(!!paste0(freq_col_prefix, "_MF_clonal") := 
               ifelse(.data$variation_type == "snv", 
                      .data[[paste0(freq_col_prefix, "_sum_clonal")]] /
                      .data[[paste0(freq_col_prefix, "_snv_depth")]],
                        .data[[paste0(freq_col_prefix, "_sum_clonal")]] /
                        .data[[paste0(freq_col_prefix, "_group_depth")]])) %>%
    dplyr::mutate(!!paste0(freq_col_prefix, "_MF_unique") :=
             ifelse(.data$variation_type == "snv", 
                    .data[[paste0(freq_col_prefix, "_sum_unique")]] /
                      .data[[paste0(freq_col_prefix, "_snv_depth")]],
                    .data[[paste0(freq_col_prefix, "_sum_unique")]] /
                      .data[[paste0(freq_col_prefix, "_group_depth")]])) %>%
    dplyr::ungroup()

#######################################
# Create a summary table
######################################
    # Grab all unique cols_to_group
  col_values <- lapply(cols_to_group, function(col) unique(mut_freq_table[[col]]))
    
    # Summary Data 
  summary_cols <- c(
    numerator_groups,
    paste0(freq_col_prefix, "_sum_unique"),
    paste0(freq_col_prefix, "_sum_clonal"),
    paste0(freq_col_prefix, "_group_depth"),
    paste0(freq_col_prefix, "_snv_depth"),
    paste0(freq_col_prefix, "_MF_unique"),
    paste0(freq_col_prefix, "_MF_clonal")
  )
summary_rows <- do.call(expand.grid, col_values)

subset_type <- subtype_list$type[subtype_list$type %in% variant_types] 

if(subtype_resolution != "none"){
summary_rows  <- merge(subset_type, summary_rows)
col_names <- c(paste(DupSeqR::subtype_dict[[subtype_resolution]]), cols_to_group)
} else {
col_names <- paste(cols_to_group) 
}

colnames(summary_rows) <- col_names 

# Grab the data and filter for mutations of interest including no_variants
summary_data <- mut_freq_table %>%
dplyr::filter(variation_type %in% subset_type | variation_type == "no_variant") %>%
dplyr::select({{ summary_cols }},
              if(!is.null(retain_metadata_cols))
               retain_metadata_cols
              ) %>%
dplyr::distinct(dplyr::across(dplyr::all_of(c(numerator_groups))), .keep_all = TRUE)

# Merge summary rows and data cols.    
summary_table <- merge(summary_rows, summary_data, by = c(col_names), all =  TRUE)

 # Make the subtype resolution list 
  # These are handled seperately in order to fill in the snv_depth for all rows
if(!is.na(DupSeqR::denominator_dict[[subtype_resolution]])){
  
  # Create a list of all snv subtypes at resolution for each group
  snv_subtype <- list(subtype_list[[subtype_resolution]]) 
  summary_rows_snv <- do.call(expand.grid, c(snv_subtype, col_values))
  colnames(summary_rows_snv) <- col_names
 # Add in the denominator reference column 
  summary_rows_snv <- summary_rows_snv %>%
    dplyr::rowwise() %>%
    dplyr::mutate(!!paste(denominator_dict[[subtype_resolution]]) := 
             get_ref_of_mut(get(subtype_dict[[subtype_resolution]])))

  # Grab the data and filter for snvs + no_variant  
  summary_data_snv <- mut_freq_table %>%
  dplyr::filter(variation_type == "snv" | variation_type == "no_variant") %>%
  dplyr::select(
    {{ summary_cols }},
    DupSeqR::denominator_dict[[subtype_resolution]],
    if(!is.null(retain_metadata_cols))
      retain_metadata_cols
  ) %>%
  dplyr::distinct(dplyr::across(dplyr::all_of(c(numerator_groups, 
                                  DupSeqR::denominator_dict[[subtype_resolution]]))), 
                  .keep_all = TRUE)
  
  # Merge summary rows and data cols.
  summary_table_snv <- merge(summary_rows_snv, summary_data_snv, all = TRUE, 
                            by = c(col_names, DupSeqR::denominator_dict[[subtype_resolution]]))
  
  # Fill in snv_depth for snv subtypes that were not in the data
  summary_table_snv <- summary_table_snv %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(denominator_groups))) %>%
    dplyr::mutate(
            dplyr::across(dplyr::all_of(!!paste0(freq_col_prefix, "_snv_depth")),
                          ~ dplyr::first(na.omit(.)))
    ) %>%
    dplyr::ungroup()

  # Bind together the types and snv data + summary_rows into one
   # In the types dataframes, create a column for the denominator reference to 
    # match snv dataframes & enable binding
   summary_rows <- summary_rows %>%
     dplyr::mutate(!!paste0(DupSeqR::denominator_dict[[subtype_resolution]]) := "N")
   
   summary_table <- summary_table %>%
     dplyr::mutate(!!paste0(DupSeqR::denominator_dict[[subtype_resolution]]) := "N")
   
   summary_rows <- rbind(summary_rows, summary_rows_snv)
   
   summary_table <- rbind(summary_table, summary_table_snv)
}  

  # Populate the group depth for rows that didn't exist in the data.
summary_table <- summary_table %>% 
  dplyr::group_by(dplyr::across(dplyr::all_of(cols_to_group))) %>%
  dplyr::mutate(
    if (!is.null(retain_metadata_cols)) { 
    dplyr::across(dplyr::all_of(retain_metadata_cols), ~ dplyr::first(na.omit(.)))},
    dplyr::across(dplyr::all_of(!!paste0(freq_col_prefix, "_group_depth")), ~ dplyr::first(na.omit(.)))
  ) %>%
  dplyr::ungroup() 

# left_join with summary_rows to keep only the (sub)types of interest
summary_table <- dplyr::left_join(summary_rows, summary_table)
# get rid of NA values
summary_table[is.na(summary_table)] <- 0

# For non-snv types, set snv_depth to NA since it is not relevant to them.
if(subtype_resolution %in% c("none", "type")) {
  summary_table <- summary_table %>%
    dplyr::select(-!!paste0(freq_col_prefix, "_snv_depth"))
} else {
  summary_table <- summary_table %>%
    dplyr::filter(.data[[paste0(DupSeqR::subtype_dict[[subtype_resolution]])]] != "snv") %>%
    dplyr::mutate(
      !!paste0(freq_col_prefix, "_snv_depth") :=
        ifelse(
          .data[[paste(subtype_dict[[subtype_resolution]])]] %in% variant_types,
          NA,
          .data[[paste0(freq_col_prefix, "_snv_depth")]]
        ))
# Create normalized proportion columns for the snv subtypes 

snv_proportions <- summary_table %>%
  filter(.data[[paste(DupSeqR::denominator_dict[[subtype_resolution]])]] != "N") %>% 
  group_by(across(all_of(c(cols_to_group)))) %>%
  mutate(total_group_mut_sum_unique = sum(.data[[paste0(freq_col_prefix, "_sum_unique")]]),
         total_group_mut_sum_clonal = sum(.data[[paste0(freq_col_prefix, "_sum_clonal")]])) %>%
  ungroup() %>%
  mutate(freq_unique = .data[[paste0(freq_col_prefix, "_sum_unique")]] / total_group_mut_sum_unique / .data[[paste0(freq_col_prefix, "_snv_depth")]],
         freq_clonal = .data[[paste0(freq_col_prefix, "_sum_clonal")]] / total_group_mut_sum_clonal / .data[[paste0(freq_col_prefix, "_snv_depth")]]) %>%
  group_by(across(all_of(c(cols_to_group)))) %>%
  mutate(total_freq_unique = sum(freq_unique),
         total_freq_clonal = sum(freq_clonal)) %>%
  ungroup() %>%
  mutate(norm_prop_unique = freq_unique/total_freq_unique,
         norm_prop_clonal = freq_clonal/total_freq_clonal) %>%
  select(-total_group_mut_sum_unique, -total_freq_unique, 
         -total_group_mut_sum_clonal, -total_freq_clonal, 
         )
  
summary_table <- dplyr::left_join(summary_table, snv_proportions)
}
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
