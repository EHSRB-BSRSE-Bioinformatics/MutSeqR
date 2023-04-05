#' Column names for mut tables
#'
#' A list of column specifications
#'
#' @format A list with potentially variable column names
#' @export
op <- list()
op$column <- list()
op$column$total_depth <- "total_depth"
op$column$n_depth <- "no_calls"
op$column$vaf <- "vaf"
op$column$depth <- "depth"
op$column$alt_depth <- "alt_depth"
op$column$is.snp <- "is_snp"
op$column$mut_depth <- "mut_depth"
op$column$mut_freq <- "mut_freq"
op$column$chr <- "contig"
op$column$start <- "start"
op$column$end <- "end"
op$column$sample <- "sample"
op$column$lower_ci <- "lower_ci"
op$column$upper_ci <- "upper_ci"
op$site$columns <- c("contig", "start")
op$mut_count_method <- "min"
op$processed_required_mut_cols <-
  c("mut_depth",
    "total_depth",
    "var_type",
    "subtype",
    "context",
    "vaf")
op$base_required_mut_cols <-
  c("contig",
    "start",
    "end",
    "sample",
    "ref",
    "alt",
    "alt_depth",
    "depth",
    "no_calls")
op$default_vaf_cutoffs <- c(0.3, 0.7, 0.9)

#' Values accepted for mutation subtypes
#'
#' These values are used to enable user input to translate to columns in a
#' mut file 
#'
#' @format A vector with corresponding values
#' @export
subtype_dict <- c(
  "none" = NA,
  "6base" = "normalized_subtype",
  "12base" = "subtype",
  "96base" = "normalized_context_with_mutation",
  "192base" = "context_with_mutation"
)

#' Values used for denominators in frequency calculations
#'
#' These values are used to cross reference base substitution types to their
#' appropriate denominators for calculations. That is, for example, the 6 base
#' substitution frequency should be subsetted based on the normalized_ref 
#' column which would contain only T or C (i.e., the pyrimidine context for
#' base substitutions).
#'
#' @format A vector with corresponding values
#' @export
denominator_dict <- c(
  "none" = NA,
  "6base" = "normalized_ref",
  "12base" = "short_ref",
  "96base" = "normalized_context",
  "192base" = "context"
)