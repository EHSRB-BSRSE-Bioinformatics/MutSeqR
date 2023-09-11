#' Add confidence intervals to mutation frequencies
#' 
#' Uses the binomial distribution to create confidence intervals for
#' mutation frequencies calculated from a single point estimate using DNA
#' sequencing
#' 
#' @param df The summary data frame containing the mutation frequencies
#' @param x Column name that specifies the mutation count (e.g., mut_depth) 
#' @param n Column name that specifies the sequencing depth (e.g., total_depth)
#' @param conf.level Confidence interval to calculate, default 95% (0.95)
#' @param method The method used by binom::binom.confint to calculate intervals.
#' Default is "wilson".
#' @returns A data frame with added columns indicating the confidence intervals.
#' @export
add_binom_conf_intervals <-
  function (df,
            x,
            n,
            conf.level = 0.95,
            method = "wilson") {
    if (length(method) != 1 || method == "all"){
      stop("Must select only one method.")
    }
    df <- as.data.frame(df)
    not_included <- setdiff(c(x, n), colnames(df))
    if (length(not_included) > 0) {
      stop(paste0(
        "Input dataframe does not include all required columns: ",
        paste(not_included, collapse = ", ")
      ))
    }
    if (!is.numeric(df[[x]]) | !is.numeric(df[[n]])) {
      stop(paste0(
        "x (",
        x,
        ", ",
        class(df[[x]]),
        ") and n (",
        n,
        ", ",
        class(df[[n]]),
        ") must be numeric."
      ))
    }
    if (nrow(df) == 0) {
      df_ci <- data.frame(numeric(0), numeric(0), numeric(0))
      colnames(df_ci) <- c("mean", op$column.lower_ci,
                           op$column.upper_ci)
    }
    else {
      df_ci <- bind_rows(mapply(function(x_val, n_val) {
        if (is.na(x_val) || is.na(n_val)) {
          data.frame(
            method = NA_character_,
            x = NA_integer_,
            n = NA_integer_,
            mean = NA_integer_,
            lower = NA_integer_,
            upper = NA_integer_
          )
        }
        else {
          binom.confint(x_val,
                               n_val,
                               conf.level = conf.level,
                               method = method)
          
        }
      }, df[, x], df[, n], SIMPLIFY = FALSE)) %>% select("mean",
                                                         "lower", "upper") %>% rename(`:=`(!!op$column.lower_ci,
                                                                                           "lower"),
                                                                                      `:=`(!!op$column.upper_ci,
                                                                                           "upper"))
    }
    cbind(df, df_ci)
  }