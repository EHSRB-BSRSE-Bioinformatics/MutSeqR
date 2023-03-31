op <- list()
op$column.total_depth <- "total_depth"
op$column.n_depth <- "no_calls"
op$column.vaf <- "vaf"
op$column.depth <- "depth"
op$column.alt_depth <- "alt_depth"
op$column.is.snp <- "is_snp"
op$column.mut_depth <- "mut_depth"
op$column.mut_freq <- "mut_freq"
op$column.chr <- "contig"
op$column.start <- "start"
op$column.end <- "end"
op$column.sample <- "sample"
op$column.lower_ci <- "lower_ci"
op$column.upper_ci <- "upper_ci"
op$site.columns <- c("contig", "start")
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

#' Column names for mut tables
#'
#' A list of column specifications
#'
#' @format A list with potentially variable column names
#' @export
op