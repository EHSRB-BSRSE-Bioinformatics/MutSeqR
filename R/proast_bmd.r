#' proast_bmd
#' @description This function is an extension of the PROAST software that
#' calculates the BMD of continuous, individual-level data.
#' @param mf_data A data frame containing the data to be analyzed. This can
#' be created using calculate_mut_freq(cols_to_group = "sample",
#' retain_metadata_cols = "dose", summary = TRUE)
#' @param dose_col The name of the column in mf_data that contains the dose.
#' Must be a numeric value.
#' @param response_col The name of the column in mf_data that contains the
#' mutation frequency. Currently, only one response column at a time is
#' supported.
#' @param covariate_col The name of the column in mf_data that contains the
#' covariate. If no covariate is present, set to NULL.
#' @param CES The critical effect size to be used in the analysis. The Benchmark
#' Response is calculated as a relative increase of CES from the control.
#' @param adjust_CES_to_group_SD A logical value indicating whether the group
#' standard deviation should be used as the CES.
#' @param model_averaging A logical value indicating whether confidence
#' intervals should be calculated using model averaging.
#' @param num_bootstraps The number of bootstrap resamples to be used in the
#' model averaging. Default is 200.
#' @param summary A logical value indicating whether a summary of the results
#' should be returned. If FALSE, raw results from the PROAST analysis are
#' returned.
#' @return A list containing the results of the PROAST analysis.
#' @export
#' @details This function is a  modified vresion of the original interactive
#' PROAST software to allow for batch processing of data. The function is
#' designed to be used with the output of the calculate_mut_freq function
#' for the purpose of calculating the Benchmark Dose of mutation frequency
#' data. As such, some functionality of the original PROAST software has
#' been removed.
#' The function will fit model 3 or 5 from various families of models
#' (Exponential, Hill, Inverse Exponential, LogNormal). It will then compare
#' the fits of models 3 and 5 for each model family and select the model with
#' the lowest AIC. The BMD confidence intervals will be calculated for each
#' model family using its selected model (3 or 5). The BMD confidence
#' may also be calculated using the bootstrap method.
#'
#' var: represents the residual variance around the fitted curve on the natural log-scale
#' Plots: the fitted curves relate to the median at each dose.
#'
#' CED the calculated value of the CED is returned in the orginal dose units, while the 
#' legend to the plot is printined in the same dose units as used in the plot (thus they
#' may differ by the dose scalling factor)
#' The CI is calculated by the profile likelihood method (likelihood ratio method)
#' 
proast_bmd <- function(mf_data,
                       dose_col = "dose",
                       response_col = "sample_MF_min",
                       covariate_col = NULL,
                       CES = 0.5,
                       adjust_CES_to_group_SD = FALSE,
                       model_averaging = TRUE,
                       num_bootstraps = 200,
                       summary = TRUE) {


  CES_sd <- as.numeric(adjust_CES_to_group_SD) + 1

  if (is.null(covariate_col)) {
    covariate <- 0
  } else {
    covariate <- covariate_col
  }

  results <- f.proast(mf_data,
                      interactive_mode = FALSE,
                      datatype = "continuous, individual data",
                      model_choice = "select model 3 or 5 from various families of models",
                      setting_choice = "set of models",
                      indep_var_choice = dose_col,
                      Vyans_input = response_col,
                      covariates = covariate,
                      custom_CES = CES,
                      adjust_CES_to_group_SD = CES_sd,
                      model_selection = "previous option with lognormal DR model added",
                      lower_dd = NULL,
                      upper_dd = NULL,
                      selected_model = "exponential",
                      model_averaging = model_averaging,
                      num_bootstraps = num_bootstraps)

  if (summary == TRUE) {
    results_df <- results[[2]]
    results_df <- results_df %>%
      dplyr::mutate(across(c(CED, CEDL, CEDU, weights), as.numeric))
    lowest_AIC_model <- results_df[results_df$AIC == min(results_df$AIC), ]
    results_ls <- list(lowest_AIC_model = lowest_AIC_model,
                       summary = results_df)

    # Cleveland plot: all models w weights
    if (model_averaging) {
      model_order <- results_df %>%
        dplyr::filter(.data$Selected.Model != "Model averaging") %>%
        dplyr::arrange(weights) %>%
        dplyr::pull(Selected.Model)
      c.plot_df <- results_df %>%
        dplyr::mutate(Selected.Model = factor(Selected.Model,
                                              levels = c(model_order,
                                                         "Model averaging")))
      # assign dummy values to Model averaging,
      # making sure it is in range of the other CEDL and CEDU.
      c.plot_df$weights[c.plot_df$Selected.Model == "Model averaging"] <- NA
      c.plot_df$CED[c.plot_df$Selected.Model == "Model averaging"] <- with(subset(c.plot_df, Selected.Model == "Model averaging"), (CEDL + CEDU) / 2)

      c <- ggplot(c.plot_df, aes(x = CED, y = Selected.Model)) +
        geom_errorbar(aes(xmin = CEDL, xmax = CEDU),
                      color = "gray",
                      width = 0.1) +
        geom_point(aes(size = weights),
                   color = "red") +
        scale_size_continuous(guide = "none") +
        ggplot2::theme(panel.background = ggplot2::element_blank(),
                       axis.line = ggplot2::element_line(),
                       panel.grid = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_line(),
                       axis.ticks.y = ggplot2::element_line()) +
        ggplot2::xlab("BMD Estimate") +
        ggplot2::ylab("Model") +
        ggtitle("BMD by Selected Model (Sorted by Weights)")

      results_ls$Cleveland_plot <- c
    }

    #
    results_raw <- results[[1]]
    names(results_raw)
    expon <- results_raw$'Expon. m5-'
    names(expon)
    expon$regr.resid.raw

    return(results_ls)
  } else {
    return(results[[1]])
  }


}