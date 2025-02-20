#' Add binomial confidence intervals to mutation frequencies.
#'
#' Uses the binomial distribution to create confidence intervals for
#' mutation frequencies calculated from a single point estimate.
#' Calculating binomial confidence intervals for mutation frequencies is
#' not part of MutSeqR's recommended workflow, but is provided here for
#' users who wish to use it.
#'
#' @param mf_data The data frame containing the mutation frequencies per sample.
#' Obtained as an output from `calculate_mf`.
#' @param sum_col Column name that specifies the mutation count (e.g., sum_min)
#' @param depth_col Column name that specifies the sequencing depth
#' (e.g., total_depth)
#' @param conf_level Confidence interval to calculate, default 95% (0.95)
#' @param method The method used by binom::binom.confint to calculate intervals.
#' Default is "wilson" (recommended).
#' @returns A mf data frame with added columns indicating the confidence
#' intervals.
#' @examples
#' example_file <- system.file("extdata", "example_mutation_data_filtered.rds", package = "MutSeqR")
#' example_data <- readRDS(example_file)
#' mf <- calculate_mf(example_data)
#' confint <- get_binom_ci(mf_data = mf,
#'                         sum_col = "sum_min",
#'                         depth_col = "group_depth")
#' @importFrom binom binom.confint
#' @importFrom dplyr bind_rows rename select
#' @export
get_binom_ci <- function(mf_data,
                         sum_col = "sum_min",
                         depth_col = "group_depth",
                         conf_level = 0.95,
                         method = "wilson") {
  if (length(method) != 1 || method == "all") {
    stop("Must select only one method.")
  }
  mf_data <- as.data.frame(mf_data)
  not_included <- setdiff(c(sum_col, depth_col), colnames(mf_data))
  if (length(not_included) > 0) {
    stop("Input dataframe does not include all required columns: ",
      paste(not_included, collapse = ", ")
    )
  }
  if (!is.numeric(mf_data[[sum_col]]) || !is.numeric(mf_data[[depth_col]])) {
    stop("sum_col (", sum_col, ", ", class(mf_data[[sum_col]]),
      ") and depth_col (", depth_col, ", ", class(mf_data[[depth_col]]),
      ") must be numeric."
    )
  }
  if (nrow(mf_data) == 0) {
    mf_data_ci <- data.frame(numeric(0), numeric(0), numeric(0))
    colnames(mf_data_ci) <- c("mean", "lower_ci", "upper_ci")

  } else {
    mf_data_ci <- dplyr::bind_rows(mapply(function(x_val, n_val) {
      if (is.na(x_val) || is.na(n_val)) {
        data.frame(
          method = NA_character_,
          sum_col = NA_integer_,
          depth_col = NA_integer_,
          mean = NA_integer_,
          lower = NA_integer_,
          upper = NA_integer_
        )
      } else {
        binom::binom.confint(x_val,
                             n_val,
                             conf.level = conf_level,
                             method = method)

      }
    }, mf_data[, sum_col], mf_data[, depth_col], SIMPLIFY = FALSE))

    mf_data_ci <- mf_data_ci %>%
      dplyr::select("mean", "lower", "upper") %>%
      dplyr::rename("lower_ci" = "lower", "upper_ci" = "upper")
  }
  cbind(mf_data, mf_data_ci)
}