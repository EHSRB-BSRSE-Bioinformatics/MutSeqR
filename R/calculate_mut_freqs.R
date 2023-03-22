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
#' `c("sample","locus")`
#' @param freq_col_prefix The prefix for column names to use for the new
#' calculated frequency. This is a prefix because we will append `_clonal` and
#' `_unique`.
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
#' @returns A data frame with the mutation frequency calculated.
#' @import tidyverse
#' @importFrom rlang :=
#' @export
calculate_mut_freq <- function(data,
                               cols_to_group = "sample",
                               freq_col_prefix = "sample",
                               subtype_resolution = "6base",
                               vaf_cutoff = 0.1,
                               summary = TRUE) {

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

  # If subtype_resolution is not in subtype_dict, throw an error... TODO
  
  numerator_groups <- c(cols_to_group, subtype_dict[[subtype_resolution]])
  numerator_groups <- numerator_groups[!is.na(numerator_groups)]
  denominator_groups <- c(cols_to_group, denominator_dict[[subtype_resolution]])
  denominator_groups <- denominator_groups[!is.na(denominator_groups)]
  mut_freq_table <- data %>%
    # Calculate numerators
    group_by(across(all_of(numerator_groups))) %>%
    mutate(
      !!paste0(freq_col_prefix,"_sum_clonal_included") :=
        sum(alt_depth[VAF < vaf_cutoff])
      ) %>%
    mutate(
      !!paste0(freq_col_prefix, "_sum_unique_only") :=
        length(alt_depth[!alt =="." & VAF < vaf_cutoff])
      ) %>%
    # Add grouping by contig and start to prevent double counts for denominator
    # Must be 1st and 2nd columns, i.e., the default scenario; names won't matter
    # this way.
    group_by(across(all_of(c(denominator_groups,  names(.)[1], names(.)[2])))) %>%
    # Remove duplicate sites for depth calculation
    distinct(!!sym(names(.)[1]), !!sym(names(.)[2]), .keep_all = TRUE) %>% # IS THIS DESTRUCTIVE?
    group_by(across(all_of(c(denominator_groups)))) %>%
    mutate(
      !!paste0(freq_col_prefix, "_depth") := sum(total_depth)
      ) %>%
    mutate(!!paste0(freq_col_prefix, "_MF_clonal") :=
             !!sym(paste0(freq_col_prefix, "_sum_clonal_included")) /
             !!sym(paste0(freq_col_prefix, "_depth"))) %>%
    mutate(!!paste0(freq_col_prefix, "_MF_unique_only") :=
             !!sym(paste0(freq_col_prefix, "_sum_unique_only")) /
             !!sym(paste0(freq_col_prefix, "_depth"))) %>%
    mutate(proportion_clonal_included =
             !!sym(paste0(freq_col_prefix, "_sum_clonal_included")) /
             sum(!!sym(paste0(freq_col_prefix, "_sum_clonal_included"))) /
             !!sym(paste0(freq_col_prefix, "_depth")) /
             sum(!!sym(paste0(freq_col_prefix, "_depth")))) %>%
    mutate(proportion_unique =
             !!sym(paste0(freq_col_prefix, "_sum_unique_only")) /
             sum(!!sym(paste0(freq_col_prefix, "_sum_unique_only")))) %>% # /
             #!!sym(paste0(freq_col_prefix, "_depth"))) %>%
    #mutate(proportion_unique = proportion_unique / sum(proportion_unique)) %>%
    ungroup()

  cols <- c(
    union(numerator_groups, denominator_groups),
    paste0(freq_col_prefix, "_sum_unique_only"),
    paste0(freq_col_prefix, "_sum_clonal_included"),
    paste0(freq_col_prefix, "_depth"),
    paste0(freq_col_prefix, "_MF_clonal"),
    paste0(freq_col_prefix, "_MF_unique_only"),
    "proportion_clonal_included",
    "proportion_unique"
  )
  summary_table <- mut_freq_table %>%
    dplyr::select({{ cols }}) %>%
    distinct()
  
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
