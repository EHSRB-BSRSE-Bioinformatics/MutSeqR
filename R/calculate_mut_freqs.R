#' Calculate mutation frequency
#'
#' Calculates the mutation frequency for arbitrary groupings and adds a
#' new named column to the data frame. 
#' @param data The data frame to be processed
#' @param group_vars A vector of grouping variables: this should be the groups
#' of interest that you want to calculate a frequency for. For instance, to get
#' the 
#' @param freq_col_name The column name to use for the new calculated frequency
#' @returns A data frame where an additional column with the new calculated 
#' frequency has been added
#' @import tidyverse
#' @importFrom rlang :=
#' @export
calculate_mut_freq <- function(data, group_vars, freq_col_name) {
  data %>%
    group_by(across(all_of(group_vars))) %>%
    mutate(!!freq_col_name := sum(alt_depth) / sum(total_depth)) %>%
    ungroup()
}


# To do... locate and enumerate recurrent mutations?
# Calculate depth for each of the 32 sequence contexts, by sample and by group
# Calculate frequency for each mouse within each 96 trinucleotide mutation

#   dplyr::group_by(normalized_context, !!sym(grouping_variable)) %>%
#   dplyr::mutate(group_depth = sum(total_depth)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(!!sym(grouping_variable)) %>%
#   dplyr::mutate(group_mut_count = sum(mut_depth)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(context_with_mutation, !!sym(grouping_variable)) %>%
#   dplyr::mutate(group_mut_count_by_type = sum(mut_depth)) %>%
#   dplyr::mutate(group_frequency = group_mut_count_by_type / group_mut_count / group_depth) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(normalized_context, sample) %>%
#   dplyr::mutate(sample_depth = sum(total_depth)) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(sample_frequency = (sum(mut_depth) / sample_depth)) %>%
#   dplyr::ungroup() %>%
#   dplyr::filter(variation_type == "snv") %>%
#   dplyr::filter(!mut_depth == 0) %>%

#group_depth = sum(total_depth),
#group_mut_count = sum(mut_depth),
#group_mut_count_by_type = sum(mut_depth),
#group_frequency = group_mut_count_by_type / group_mut_count / group_depth
