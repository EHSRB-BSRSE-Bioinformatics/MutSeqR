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
#' @param plot_results A logical value indicating whether the results should
#' be plotted. Plots can either be saved as an svg to the output_path or
#' displayed in the R plot viewer.
#' @param output_path The filepath to save the plots .svg files. If NULL, the plots
#' will be saved in the current working directory.
#' @param output_type svg or none. If svg, the plots will be saved as .svg files
#' in the output_path. If none, the plots will NOT be saved, but be displayed in
#' the R plot viewer.
#' @return If summary is TRUE, a data frame of final results. If summary is
#' FALSE, a list of the raw results from the PROAST analysis.
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
                       summary = TRUE,
                       plot_results = FALSE,
                       output_path = NULL,
                       output_type = "none") {

  if (!dose_col %in% colnames(mf_data)) {
    stop("Dose column not found in mf_data")
  }
  if (!any(response_col %in% colnames(mf_data))) {
    stop("Response column not found in mf_data")
  }
  if (!is.null(covariate_col) && !covariate_col %in% colnames(mf_data)) {
    stop("Covariate column not found in mf_data")
  }
  if (model_averaging == TRUE) {
    message("Model averaging is set to TRUE. This may take some time to run.")
  }
  # ensure that dose is numeric #
  if (!is.numeric(mf_data[[dose_col]])) {
    stop("Dose column must be numeric")
  }
  if (plot_results == TRUE && model_averaging == TRUE && output_type == "svg") {
    if (!require("svglite", quietly = TRUE)) {
      stop("The 'svglite' package is required to save model averaging plots as svgs. Please install to use this functionality.")
}
  }

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
                      num_bootstraps = num_bootstraps,
                      display_plots = FALSE)

  if (plot_results == TRUE) {
    f.plot.result(results[[1]],
                  output_path = output_path,
                  output_type = output_type,
                  model_averaging = FALSE)

    if (model_averaging == TRUE) {
      f.plot.result(results[[1]],
                    output_path = output_path,
                    output_type = output_type,
                    model_averaging = TRUE)

      # Cleveland plot: all models w weights
      results_df <- results[[2]] %>%
        dplyr::mutate(CED = as.numeric(CED),
                      CEDL = as.numeric(CEDL),
                      CEDU = as.numeric(CEDU))

      for (i in unique(results_df$Response)) {
        c.plot.df <- results_df %>%
          dplyr::filter(Response == i)

        if (!is.null(covariate_col)) {
          c.plot.df <- c.plot.df %>%
            dplyr::rename(Model = "Selected.Model")
          c.plot.df$Selected.Model <- paste(c.plot.df$Model, c.plot.df$Covariate)
          model_order <- c.plot.df %>%
            dplyr::arrange(Covariates, weights) %>%
            dplyr::pull(Selected.Model)
          c.plot.df$Selected.Model <- factor(c.plot.df$Selected.Model,
                                             levels = model_order)
        } else {
          model_order <- c.plot.df %>%
            dplyr::filter(.data$Selected.Model != "Model averaging") %>%
            dplyr::arrange(weights) %>%
            dplyr::pull(Selected.Model)
          
          c.plot.df <- c.plot.df %>%
            dplyr::mutate(Selected.Model = factor(Selected.Model,
                                                  levels = c(unique(model_order),
                                                            "Model averaging")))
        }
        # assign dummy values to Model averaging,
        # making sure it is in range of the other CEDL and CEDU.
        c.plot.df$weights[c.plot.df$Model == "Model averaging"] <- NA
        c.plot.df$CED[c.plot.df$Model == "Model averaging"] <- with(subset(c.plot.df, Model == "Model averaging"), (CEDL + CEDU) / 2)

        c <- ggplot(c.plot.df, aes(x = CED, y = Selected.Model)) +
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
          ggplot2::xlab(paste("BMD Estimate for", i)) +
          ggplot2::ylab("Model") +
          ggtitle("BMD by Selected Model (Sorted by Weights)")
        if (output_type == "svg") {
          file_name <- file.path(output_path, paste0("PROAST_", i , "_cleveland.svg"))
          ggsave(filename = file_name, plot = c, device = "svg")
        } else if (output_type == "none") {
          # Determine the operating system to open the plot in the correct viewer
          os <- Sys.info()[["sysname"]]
          # Open the plot in the correct viewer
          if (os == "Windows") {
            windows(width = 8, height = 6)
          } else if (os == "Darwin") {  # Darwin is macOS
            quartz(width = 8, height = 6)
          } else if (os == "Linux") {
            X11(width = 8, height = 6)
          }
          print(c)
          dev.flush()
        }
      }
    }
  }
  if (summary == TRUE) {
    return(results[[2]])
  } else {
    return(results[[1]])
  }
}