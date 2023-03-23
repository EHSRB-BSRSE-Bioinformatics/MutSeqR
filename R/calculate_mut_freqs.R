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
#' @param data The data frame to be processed
#' @param cols_to_group A vector of grouping variables: this should be the groups
#' of interest that you want to calculate a frequency for. For instance, getting
#' the frequency by `sample`. Other options might include `locus`, or,
#' `c("sample","locus")`. Must be a column in the mutation data table.
#' @param subtype_resolution At what resolution should the frequencies be
#' calculated? This is an important consideration because the reference bases
#' sequencing depth will differ depending on which mutation subtypes are being
#' considered. Options are "none", 6base", "12base", "96base", and "192base".
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
#' @import tidyverse
#' @importFrom rlang :=
#' @export
calculate_mut_freq <- function(data,
                               cols_to_group = "sample",
                               subtype_resolution = "6base",
                               vaf_cutoff = 0.1,
                               clonality_cutoff = 0.3,
                               summary = TRUE,
                               variant_types = "snv") {

  # Define internal objects
  
  freq_col_prefix <- paste0(cols_to_group, collapse = "_")
  
  subtype_dict <- c(
    "none" = NA,
    "6base" = "normalized_subtype",
    "12base" = "subtype",
    "96base" = "normalized_context_with_mutation",
    "192base" = "context_with_mutation"
    )
  
  denominator_dict <- c(
    "none" = NA,
    "6base" = "normalized_ref",
    "12base" = "short_ref",
    "96base" = "normalized_context",
    "192base" = "context"
    )
  
  if (!subtype_resolution %in% names(subtype_dict)) {
    stop("Error: you need to set subtype_resolution to one of \"none\",\"6base\",
         \"12base\", \"96base\" or \"192base\".")
  }

  numerator_groups <- c(cols_to_group, subtype_dict[[subtype_resolution]])
  numerator_groups <- numerator_groups[!is.na(numerator_groups)]
  denominator_groups <- c(cols_to_group, denominator_dict[[subtype_resolution]])
  denominator_groups <- denominator_groups[!is.na(denominator_groups)]
  
  # Calculate mutation frequencies
  mut_freq_table <- data %>%
    # Identify duplicate entries for depth calculation later on
    group_by(across(all_of(c(cols_to_group)))) %>%
    mutate(is_duplicate = duplicated(paste(!!sym(names(.)[1]), start))) %>% # NOTE THIS STILL CAN VARY BECAUSE THE DEPTH MAY BE DIFFERENT BETWEEN TWO DUPLICATES
    # Calculate numerators
    group_by(across(all_of(numerator_groups))) %>%
    mutate(
      !!paste0(freq_col_prefix,"_sum_clonal") :=
        sum(alt_depth[!alt =="." & VAF < vaf_cutoff])
      ) %>%
    mutate(
      !!paste0(freq_col_prefix, "_sum_unique") :=
        length(alt_depth[!alt =="." & VAF < vaf_cutoff])
      ) %>%
    # Calculate denominator (same for clonal and unique mutations)
    group_by(across(all_of(denominator_groups))) %>%
    mutate(
      !!paste0(freq_col_prefix, "_depth") := sum(total_depth[!is_duplicate])
      ) %>%
    # Calculate frequencies
    mutate(!!paste0(freq_col_prefix, "_MF_clonal") :=
             !!sym(paste0(freq_col_prefix, "_sum_clonal")) /
             !!sym(paste0(freq_col_prefix, "_depth"))) %>%
    mutate(!!paste0(freq_col_prefix, "_MF_unique") :=
             !!sym(paste0(freq_col_prefix, "_sum_unique")) /
             !!sym(paste0(freq_col_prefix, "_depth"))) %>%
    ungroup()

  # Define columns of interest for summary table
  cols <- c(
    numerator_groups,
    paste0(freq_col_prefix, "_sum_unique"),
    paste0(freq_col_prefix, "_sum_clonal"),
    paste0(freq_col_prefix, "_depth"),
    paste0(freq_col_prefix, "_MF_clonal"),
    paste0(freq_col_prefix, "_MF_unique")
  )

  # Make summary table of frequencies
  # This is also where subtype proportions are calculated
  summary_table <- mut_freq_table %>%
    filter(variation_type %in% variant_types) %>%
    dplyr::select({{ cols }}) %>%
    distinct() %>%
    mutate(freq_clonal =
             !!sym(paste0(freq_col_prefix, "_sum_clonal")) /
             sum(!!sym(paste0(freq_col_prefix, "_sum_clonal"))) /
             !!sym(paste0(freq_col_prefix, "_depth")) ) %>%
    mutate(prop_clonal = freq_clonal / sum(freq_clonal)) %>%
    mutate(freq_unique =
             !!sym(paste0(freq_col_prefix, "_sum_unique")) /
             sum(!!sym(paste0(freq_col_prefix, "_sum_unique"))) /
             !!sym(paste0(freq_col_prefix, "_depth")) ) %>%
    mutate(prop_unique = freq_unique / sum(freq_unique))

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
#cols_to_group <- c("sample","variation_type")
#cols_to_group <- c("sample", "CpG_site")

# To do... locate and enumerate recurrent mutations?
