#' BMD modeling using PROAST
#' @description Calculate the benchmark dose (BMD) of continuous,
#' individual-level data with optional model averaging. This function
#' is intended to model the dose-response of mutation frequency.
#' This function is an extension of the PROAST software (copyright RIVM
#' National Institute for Public Health and the Environment).
#' @param mf_data A data frame containing the data to be analyzed. Data should
#' be individual for each sample. Required columns are the column containing
#' the dose `dose_col` the column(s) containing the mutation frequency
#' `response_col`, and the column containing the covariate `covariate_col`,
#' if applicable.
#' @param dose_col The column in `mf_data` containing the dose data. Values
#' must be numeric. Default is "dose".
#' @param response_col The column(s) in `mf_data` containing the mutation
#' frequency. Multiple `response_col`s can be provided. Default is "mf_min".
#' @param covariate_col The column in `mf_data` containing the covariate.
#' If no covariate is present, set to \code{NULL} (default).
#' @param bmr The Benchmark Response value. The BMR is defined as a
#' bmr-percent change in mean response relative to the controls.
#' Default is 0.5 (50% change).
#' @param adjust_bmr_to_group_sd A logical value indicating whether the group
#' standard deviation should be used as the BMR. If TRUE, the BMR will be
#' bet set to one standard deviation above the control group mean. Default is
#' FALSE.
#' @param model_averaging A logical value indicating whether confidence
#' intervals should be calculated using model averaging. Default is TRUE
#' (recommended).
#' @param num_bootstraps The number of bootstrap resamples to be used in the
#' model averaging. Default is 200 (recommended).
#' @param raw_results A logical value indicating whether to return the raw
#' results from the PROAST analysis. If FALSE, data is returned as a
#' summary table.
#' @param plot_results A logical value indicating whether to plot the BMD models
#' and/or the Cleveland plots. Default is FALSE. Plots may be exported
#' directly to an output_path, or returned within a list to the user.
#' @param output_path The file path indicating where to save the plots.
#' If NULL, the plots will automatically be displayed to the graphics
#' window and then returned as a list alongside the bmd results.
#' @return A summmary data frame of final results. If plots or raw results
#' are selected, all data will be returned within a list.
#'
#' The summary will include the following for each response variable and
#' covariate subgroup (if applicable):
#'  \itemize{
#'    \item `Model`: The m3 or m5 model selected for each model family
#' (Exponential, Hill, Inverse Exponential, LogNormal).
#'    \item `Response`: The response variable.
#'    \item `Covariate`: The covariate subgroup, if applicable.
#'    \item `bmr`: The specified Benchmark Response.
#'    \item `BMD`: The Benchmark Dose, in original dose units, estimated for
#' the given model.
#'    \item `BMDL`: The lower bound of the 90% confidence interval for the BMD,
#' calculated by the profile likelihood method.
#'    \item `BMDU`: The upper bound of the 90% confidence interval for the BMD,
#' calculated by the profile likelihood method.
#'   \item `AIC`: The Akaike Information Criterion for the selected model.
#' Lower values indicate a better fit. It is advised to choose the BMD value
#' from the model with the lowest AIC.
#'  \item `weights`: The weight of the model in the model averaging process,
#' if applicable.
#' \item `Model averaging`: The BMDL and BMDU calculated by the bootstrap
#' method if \code{model_averaging = TRUE}.
#' }
#' If there is no significant response in the data, the function will return an
#' empty data frame.
#'
#' If \code{plot_results = TRUE} the function will create the following plots
#' for each response variable. The plots will be saved to the output_path.
#' If no output_path is provided, then they will be returned within a list
#' alongside the summary data frame.
#'  \itemize{
#'   \item Model Plots. The following plot will be created for each model
#' family (Exponential, Hill, Inverse Exponential, LogNormal): The fitted curve
#' of the selected (3 or 5) model. Data is log-transformed. Individual data
#' points are plotted using small triangles. The geometric mean (median) at
#' each dose is plotted as a large triangle. The BMD is indicated by the
#' dotted line. If applicable, the covariate subgroup is indicated by color.
#' \item bootstrap_curves If \code{model_averaging = TRUE}, the bootstrap
#' curves based on model averaging. The geometric mean (median) at each dose
#' is plotted as a large triangle. Data is log-transformed.
#' \item cleveland plot if \code{model_averaging = TRUE} The BMD estimate
#' for each model is plotted as a red point alongside the 90% confidence
#' intervals. The size of the BMD point represents the model weight assigned
#' during model averaging, based on the AIC.
#' }
#'
#' If \code{raw_results = TRUE}, the function will return the raw results of
#' the PROAST analysis alongside the summary data frame. PROAST raw_results
#' is a list of variables and data that is continuously modified as it is
#' passed through the proast functions. It can be given to f.proast() to
#' resume analysis.
#' @details This function is a  modified version of the original interactive
#' PROAST software (\url{https://www.rivm.nl/en/proast} that allows for batch
#' processing of data. The function is designed to be used with the output of
#' \code{calculate_mf} for the purpose of calculating the
#' Benchmark Dose of mutation frequency data. As such, some functionality of
#' the original PROAST software has been removed.
#'
#' This function will accept continuous data, with an observation for each
#' individual subject. It is assumed that data are lognormally distributed.
#' The response data is log-transformed, then back-transformed after the
#' statistical analysis. The function will fit model 3 or 5 from various
#' families of models (Exponential, Hill, Inverse Exponential, LogNormal).
#' It will then compare the fits of models 3 and 5 for each model family and
#' select the model with the lowest AIC. The BMD 90% confidence intervals will
#' be calculated based on the selected model (3 or 5) for each model family
#' using the profile likelihood method. The BMD 90% confidence interval may
#' also be calculated using the bootstrap method if
#' \code{model_averaging = TRUE}. It is recommended to use 200 bootstraps for
#' model averaging.
#'
#' To replicate these results in the PROAST interactive software,
#' select the following menu options:
#' \enumerate{
#'    \item f.proast(mf_data)
#'    \item What type of response data do you want to consider?
#' \emph{1: continuous, individual data}
#'    \item Do you want to fit a single model or fit various nested families
#' of models? \emph{3: select model 3 or 5 from various families of models}
#'    \item Q1: Which variable do you want to consider as the independent
#' variable? \emph{# : dose_col}
#'    \item Give number(s) of the response(s) you want to analyse.
#' \emph{# : response_col}
#'    \item Give number of factor serving as potential covariate (e.g.sex)
#' type 0 if none. \emph{# : covariate_col}
#'    \item Do you want to adjust CES to within group SD?
#' \emph{1: no, 2: yes | adjust_bmr_to_group_sd: FALSE/TRUE}
#'    \item Give value for CES (always positive) type 0 to avoid calculation
#' of CIs. \emph{bmr}
#'    \item Do you want to calculate the BMD confidence interval by model
#' averaging? \emph{1: no 2: yes | model_averaging: FALSE/TRUE}
#'    \item give number of bootstrap runs for calculating BMD confidence
#' interval based on MA (e.g. 200) \emph{num_bootstraps}
#'    \item Which models do you want to be fitted?
#' \emph{4 : previous option with lognormal DR model added}
#' }
#' @examples
#' # Calculate the BMD for a 50% increase in mutation frequency from control
#' # With Model averaging.
#' # For the purpose of this example, num_bootstraps is set to 5 to reduce
#' # run time. 200 bootstraps is recommended.
#' example_file <- system.file("extdata", "Example_files",
#'                             "example_mutation_data_filtered.rds",
#'                             package = "MutSeqR")
#' example_data <- readRDS(example_file)
#' mf <- calculate_mf(example_data, retain_metadata_cols = "dose")
#' bmd <- bmd_proast(mf_data = mf,
#'                   dose_col = "dose",
#'                   response_col = c("mf_min", "mf_max"),
#'                   bmr = 0.5,
#'                   model_averaging = TRUE,
#'                   num_bootstraps = 5)
#' # Plot the Model Averaging 90% CI using plot_ci()
#' plot_df <- bmd %>%
#'  dplyr::filter(Model == "Model averaging") %>%
#'  dplyr::select(Response, BMD, BMDL, BMDU)
#' plot <- plot_ci(plot_df, order = "asc", log_scale = FALSE)
#' @export
#' @importFrom dplyr arrange filter mutate pull rename
bmd_proast <- function(
  mf_data,
  dose_col = "dose",
  response_col = "mf_min",
  covariate_col = NULL,
  bmr = 0.5,
  adjust_bmr_to_group_sd = FALSE,
  model_averaging = TRUE,
  num_bootstraps = 200,
  plot_results = FALSE,
  output_path = NULL,
  raw_results = FALSE
) {

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
  # ensure that dose is numeric
  if (!is.numeric(mf_data[[dose_col]])) {
    stop("Dose column must be numeric")
  }

  bmr_sd <- as.numeric(adjust_bmr_to_group_sd) + 1

  if (is.null(covariate_col)) {
    covariate <- 0
  } else {
    covariate <- covariate_col
  }

  results <- f.proast(
    mf_data,
    interactive_mode = FALSE,
    datatype = "continuous, individual data",
    model_choice = "select model 3 or 5 from various families of models",
    setting_choice = "set of models",
    indep_var_choice = dose_col,
    Vyans_input = response_col,
    covariates = covariate,
    custom_CES = bmr,
    adjust_CES_to_group_SD = bmr_sd,
    model_selection = "previous option with lognormal DR model added",
    lower_dd = NULL,
    upper_dd = NULL,
    selected_model = "exponential",
    model_averaging = model_averaging,
    num_bootstraps = num_bootstraps,
    display_plots = FALSE
  )

  if (plot_results == TRUE) {
    if (!is.null(output_path)) { # Export the plots
      output_path <- file.path(output_path)
      if (!exists(output_path)) {
        dir.create(output_path)
      }
      f.plot.result(
        results[[1]],
        output_path = output_path,
        output_type = "svg",
        model_averaging = FALSE
      )
      if (model_averaging == TRUE) {
        f.plot.result(
          results[[1]],
          output_path = output_path,
          output_type = "svg",
          model_averaging = TRUE
        )
        cleveland_plot(
          results,
          covariate_col = covariate_col,
          output_path = output_path
        )
      }
    } else { # return the plots
      plots <- f.plot.result(
        results[[1]],
        output_path = output_path,
        output_type = "none",
        model_averaging = FALSE
      )
      if (model_averaging == TRUE) {
        ma_plots <- f.plot.result(
          results[[1]],
          output_path = output_path,
          output_type = "none",
          model_averaging = TRUE
        )
        c_plot <- cleveland_plot(
          results,
          covariate_col = covariate_col,
          output_path = output_path
        )
        for (i in seq_along(c_plot)) {
          grDevices::dev.new()
          print(c_plot[[i]])
        }
        plots <- c(plots, ma_plots, c_plot)
      }
    }
  }
  # Select the BMD with the lowest AIC for each response
  # If they have the same AIC, take the mean.
  summary <- results[[2]] %>%
    dplyr::rename(BMD = "CED",
                  BMDL = "CEDL",
                  BMDU = "CEDU",
                  Model = "Selected.Model",
                  BMR = "CES") %>%
    dplyr::select(-"Log.Likelihood", -"Var", -"a", -"d")
  results_list <- list(summary = summary)
  if (plot_results && is.null(output_path)) {
    results_list <- c(results_list, plots)
  }
  if (raw_results) {
    results_list$raw_results <- results[[1]]
  }
  if (length(results_list) > 1) {
    return(results_list)
  } else {
    return(summary)
  }
}