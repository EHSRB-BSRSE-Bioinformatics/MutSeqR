#' proast_bmd
#' @description This function is an extension of the PROAST software that
#' calculates the BMD of continuous, individual-level data.
#' @param mf_data A data frame containing the data to be analyzed. This can
#' be created using calculate_mut_freq(cols_to_group = "sample",
#' retain_metadata_cols = "dose", summary = TRUE)
#' @param dose_col The name of the column in mf_data that contains the dose.
#' Must be a numeric value.
#' @param response_col The name of the column(s) in mf_data that contains the
#' mutation frequency. Currently, only one response column at a time is
#' supported. *FIX* this.
#' @param covariate_col The name of the column in mf_data that contains the
#' covariate. If no covariate is present, set to NULL.
#' @param CES The critical effect size to be used in the analysis. The Benchmark
#' Response is calculated as a relative increase of CES from the control.
#' @param adjust_CES_to_group_SD A logical value indicating whether the group
#' standard deviation should be used as the CES.
#' @param model_averaging A logical value indicating whether confidence
#' intervals should be calculated using model averaging.
#' @param num_bootstraps The number of bootstrap resamples to be used in the
#' model averaging.
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
#' CED the calculated value of the CED is retirn in the orginal dose units, while the 
#' legend to the plot is printined in the same dose units as used in the plot (thus they
#' may differ by the dose scalling factor)
#' The CI is calculated by the profile likelihood method (likelihood ratio method)
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

  if(is.null(covariate_col)) {
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
    if (!is.null(covariate_col)) {
      ced.list <- list()
      for (i in names(results[[1]])) {
        ced.list[[i]] <- results[[1]][[i]]$ced.table
      }
      ced.df.list <- lapply(names(ced.list), function(name) {
        df <- ced.list[[name]]
        df$model <- name
        return(df)
      })
      ced.df <- do.call(rbind, ced.df.list)
      results_df <- dplyr::left_join(results_df, ced.df)
    }
    BMD <- results_df[results_df$AIC == min(results_df$AIC), ]
    results <- list(BMD = BMD,
                    summary = results_df)
    # Plot Confidence Intervals
    results_bmd_df_plot <- BMD %>%
    dplyr::group_by(response) %>%
    dplyr::mutate(max = BMDU) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(across(where(is.numeric), round, 1)) %>%
    tidyr::pivot_longer(cols = c("BMD", "BMDL", "BMDU"))

  nudge_value <- 0.3

  g <- ggplot(results_bmd_df_plot,
              aes(x = value, y = response, color = name)) +
    geom_line(aes(group = response), color = "#b8b8b8", linewidth = 3.5) +
    geom_point(size = 3) +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.y = element_text(color = "black"),
          axis.text.x = element_text(color = "#000000"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          panel.grid = element_blank()) +
    scale_color_manual(values = c("black", "#BF2F24", "#436685")) +
    scale_x_continuous() +
    geom_text(aes(label = value, color = name),
              size = 3.25,
              nudge_x = dplyr::if_else(
                results_bmd_df_plot$value == results_bmd_df_plot$max, # if it's the larger value...
                nudge_value,   # move it to the right of the point
                -nudge_value), # otherwise, move it to the left of the point
              hjust = dplyr::if_else(
                results_bmd_df_plot$value==results_bmd_df_plot$max, #if it's the larger value
                0, # left justify
                1)# otherwise, right justify
    ) +
    ggplot2::labs(x = "BMD", y = "Response",
                  title = paste0("BMD with ",
                                  conf_int,
                                  "% Confidence Intervals"),
                  color = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0,
                                                       vjust = 0.5,
                                                       hjust = 0.5),
                    plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(axis.ticks = element_line(color = "black", size = 0.5))
  } else {
    return(results[[1]])
  }


}