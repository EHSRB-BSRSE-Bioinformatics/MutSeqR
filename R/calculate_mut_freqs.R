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
#' @param custom_column_names A list of names to specify the meaning of column
#'  headers. Since column names can vary with data, this might be necessary to
#'  digest the mutation data table properly. Typical defaults are set, but can
#'  be substituted in the form of `list(total_depth = "my_custom_depth_name", 
#'  sample = "my_custom_sample_column_name")`. For a comprehensive list, see 
#'  examples. You can change one or more of these. The most likely column to
#'  need changing is the "chromosome name" (chr), which by default is "seqnames"
#'  but could be "contig", "chr", or others.  (TODO - MAKE EXAMPLES)
#' @returns A data frame with the mutation frequency calculated.
#' @import tidyverse
#' @import collapse
#' @importFrom rlang :=
#' @export
calculate_mut_freq <- function(data,
                               cols_to_group = "sample",
                               subtype_resolution = "6base",
                               vaf_cutoff = 0.1,
                               clonality_cutoff = 0.3,
                               summary = TRUE,
                               variant_types = c("snv","indel","mnv","sv"),
                               custom_column_names = list(chr = "seqnames")) {

  # Define internal objects
  default_columns <- list(
    start = "start",
    total_depth = "total_depth",
    n_depth = "no_calls",
    chr = "seqnames"
  )
  cols <- modifyList(default_columns, custom_column_names)
  
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
         \"12base\", \"96base\", \"192base\".")
  }

  numerator_groups <- c(cols_to_group, subtype_dict[[subtype_resolution]])
  numerator_groups <- numerator_groups[!is.na(numerator_groups)]
  denominator_groups <- c(cols_to_group, denominator_dict[[subtype_resolution]])
  denominator_groups <- denominator_groups[!is.na(denominator_groups)]
  
  # Calculate mutation frequencies
  mut_freq_table <- data %>%
    # Identify duplicate entries prior to depth calculation
    group_by(across(all_of(c(cols_to_group, cols$chr, cols$start)))) %>%
    mutate(num_dups = n(), 
           dup_id = row_number()) %>% 
    ungroup() %>%
    mutate(is_duplicated = num_dups > 1) %>%
    mutate(depth_undupes = ifelse(dup_id>1, 0, !!sym(cols$total_depth))) %>%
    #Identify values with the same start, sample, and depth
    group_by(across(all_of(c(cols_to_group, cols$chr, cols$start, cols$total_depth)))) %>%
    mutate(is_depth_duplicated = n() > 1) %>%
    ungroup() %>%
    # For sites that share start positions within a group, deduplicate depth.
    # For duplicated depths, set all instance but one to zero; for
    # indels/svs/mnvs, take the 'no_variant' depth instead.
    mutate(depth_final = ifelse(is_duplicated == TRUE & is_depth_duplicated == FALSE,
                                ifelse(variation_type == "no_variant", !!sym(cols$total_depth), 0),
                                ifelse(is_depth_duplicated == TRUE, depth_undupes, total_depth))) %>%
    group_by(across(all_of(numerator_groups))) %>%
    mutate(!!paste0(freq_col_prefix, "_sum_clonal") :=
        sum(alt_depth[!variation_type == "no_variant" & VAF < vaf_cutoff])) %>%
    mutate(!!paste0(freq_col_prefix, "_sum_unique") :=
             length(alt_depth[!variation_type == "no_variant" & VAF < vaf_cutoff])) %>%
    # Calculate denominator (same for clonal and unique mutations)
    group_by(across(all_of(denominator_groups))) %>%
    mutate(!!paste0(freq_col_prefix, "_depth") := sum(depth_final)) %>%
    #mutate(!!paste0(freq_col_prefix, "_depth") := sum(total_depth[!is_duplicate])) %>%
    # Calculate frequencies
    mutate(!!paste0(freq_col_prefix, "_MF_clonal") :=
             !!sym(paste0(freq_col_prefix, "_sum_clonal")) /
             !!sym(paste0(freq_col_prefix, "_depth"))) %>%
    mutate(!!paste0(freq_col_prefix, "_MF_unique") :=
             !!sym(paste0(freq_col_prefix, "_sum_unique")) /
             !!sym(paste0(freq_col_prefix, "_depth"))) %>%
    ungroup()

  # Define columns of interest for summary table
  summary_cols <- c(
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
    dplyr::select({{ summary_cols }}) %>%
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
