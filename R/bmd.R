#' Fit a BMD model for mutation frequency data
#' 
#' @description This function fits a continuous BMD model to
#' dose-response data for mutation frequency.
#' @param mf_data A data frame with columns "dose" and "MFmin"
#' @param model_type A string specifying the model type.
#' Options are "hill", "exp-3","exp-5", "power", "polynomial"
#' @param BMR A numeric value specifying the benchmark response.
#' Default is 0.5.
#' @param model_avg A logical value specifying whether to
#' average the model fits. Default is TRUE.
#' @return A list with the following components:
#' @importFrom ToxicR single_continuous_fit
#' @export
mf_bmd <- function(mf_data,
                   dose_col = "dose",
                   response_col = c("sample_MF_unique", "sample_MF_clonal"),
                   model_types = c("hill",
                                   "exp-3",
                                   "exp-5",
                                   "power",
                                   "polynomial"),
                   BMR = 0.5) {

  if (!requireNamespace("ToxicR")) {
    stop("ToxicR is not installed. Please install from https://github.com/NIEHS/ToxicR")
  }
  result_list <- list()
  # Run the function for each response column
  for (response in response_col) {
    fit_results <- list() # Initialize a list to store the fit results
    plot_list <- list()  # Initialize a list to store the plots

    # Run the function for each model type
    for (model_type in model_types) {
      fit_results[[model_type]] <-
        ToxicR::single_continuous_fit(mf_data[, dose_col],
                                      mf_data[, response],
                                      model_type = model_type,
                                      BMR = BMR,
                                      BMR_TYPE = "rel",
                                      fit_type = "laplace")

      # Plot the model
      plot_list[[model_type]] <- plot(fit_results[[model_type]])
    }

    # Extract the bmd values
    bmd_values <- lapply(fit_results, function(x) x$bmd)
    # Create a dataframe from the list and add a column for the names
    result_df <- dplyr::bind_rows(bmd_values, .id = "dataset")

    # Store the result_df and plot_list in the result_list
    result_list[[response]] <- list(data = result_df, plots = plot_list)
  }
# Calculate the AIC for each model and add it to the fit_results_df
  return(result_list)

# AIC <- model max + 2 x DF
# AIC<- -model_indiv$maximum + 2*summary(model_indiv)$GOF[1,2]

}

#' Fit a model averaged continuous BMD model
#' @description This function fits a model averaged continuous BMD model
#' to dose-response data for mutation frequency.
#' @param mf_data A data frame containing the dose-response data. Data may
#' be individual for each sample or averaged over dose groups.
#' Required columns for individual data are "dose" and "response".
#' Multiple response columns are allowed. Required columns for summarised
#' data are "dose", "mean response", "sample size", and "standard deviation".
#' Only one response column is allowed for summarised data.
#' @param data_type A string specifying the type of response data.
#' Data may be response per individual or summarised across dose groups.
#' ("individual", "summary").
#' @param dose_col A character string specifying the column in mf_data
#' containing the dose data.
#' @param response_cols A character vector specifying the columns in mf_data
#' containing the response data. For summarised data types, this should be
#' the mean response for each dose group.
#' @param sd_col A character string specifying the column in mf_data containing
#' the standard deviation of the response data. This is only required for
#' summarised data types.
#' @param n_col A character string specifying the column in mf_data containing
#' the sample size of each dose group. This is only required for summarised
#' data types.
#' @param bmr_type A string specifying the type of benchmark response.
#' For continuous models, there are four types of BMD definitions that
#' are commonly used:
#'    \itemize{
#'      \item Relative deviation (default; 'BMR_TYPE = "rel"'). This defines the
#' BMD as the dose that changes the control mean/median a certain
#' percentage from the background dose, i.e. it is the dose, BMD that
#' solves \eqn{\mid f(dose) - f(0) \mid = (1 \pm BMR) f(0)}
#'      \item Standard deviation ('BMR_TYPE = "sd"'). This defines the BMD as
#' the dose associated with the mean/median changing a specified number of
#' standard deviations from the mean at the control dose., i.e., it is the
#' dose, BMD, that solves \eqn{\mid f(dose)-f(0) \mid = BMR \times \sigma}
#'      \item Hybrid deviation ('BMR_TYPE = "hybrid"'). This defines the
#' BMD that changes the probability of an adverse event by a stated amount
#' relative to no exposure (i.e 0).  That is, it is the dose, BMD, that solves
#' \eqn{\frac{Pr(X > x| dose) - Pr(X >x|0)}{Pr(X < x|0)} = BMR}.
#' For this definition,
#' \eqn{Pr(X < x|0) = 1 - Pr(X > X|0) = \pi_0}, where \eqn{0 \leq \pi_0 < 1}
#' is defined by the user as "point_p," and it defaults to 0.01.  Note: this
#' discussion assumed increasing data.  The fitter determines the direction
#' of the data and inverts the probability statements for decreasing data.
#'      \item Absolute deviation ('BMR_TYPE="abs"'). This defines the BMD
#' as an absolute change from the control dose of zero by a specified amount.
#' That is the BMD is the dose that solves the equation
#' \eqn{\mid f(dose) - f(0) \mid = BMR}.
#'    }
#' @param bmr A numeric value specifying the benchmark response. The BMR is
#' defined in relation to the calculation requested in bmr_type. Default is 0.5.
#' @param fit A string specifying the method used to fit the model.
#' Options are ("laplace", "mle", or "mcmc"). Default is "laplace".
#' @param a The specified nominal coverage rate for computation of
#' the lower bound on the BMDL and BMDU, i.e., one computes a
#' \eqn{100\times(1-\alpha)\%} confidence interval.  For the
#' interval (BMDL,BMDU) this is a \eqn{100\times(1-2\alpha)\% }.
#' By default, it is set to 0.025 for a CI of 95%.
#' @param ... Additional arguments to be passed to the model fitting function.
#' See \link[ToxicR]{ma_continuous_fit} for details.

#' @details Model averaging is done over the default models described
#' in The European Food Safety Authority's (2022) Guidance on the use
#' of the benchmark dose approach in risk assessment. These models are:
#'  \itemize{
#'      \item \code{"exp-aerts"}: \eqn{f(x) = a(1 + (c-1)(1-\exp(-bx^{d}))) }
#'      \item \code{"invexp-aerts"}: \eqn{f(x) = a(1 + (c-1)(\exp(-bx^{-d})))}
#'      \item \code{"hill-aerts"}:
#' \eqn{f(x) = a(1 + (c-1)(1-\frac{b^d}{b^d + x^d}))}
#'      \item \code{"lognormal-aerts"}:
#' \eqn{f(x) = a\left\{1 + (c-1)\left(\Phi( \ln(b) + d\times \ln(x))\right) \right\}}
#'      \item \code{"gamma-efsa"}: \eqn{f(x) = a(1 + (c-1)(\Gamma(bx; d))) }
#'      \item \code{"LMS"}: \eqn{f(x) = a(1 + (c-1)(1 - \exp(-bx - dx^2))) }
#'      \item \code{"probit-aerts"}:
#' \eqn{f(x) = c\left(\Phi(a + b\times x^d)\right) }
#'      \item \code{"logistic-aerts"}:
#' \eqn{f(x) = \frac{c}{1 + \exp(-a - b\times x^d)} }
#'    }
#'   Here: \eqn{\Phi(\cdot)} is the standard normal distribution and
#'         \eqn{\Phi_{SN}(\cdot;\cdot)} is the skew-normal distribution
#'
#' @return A list with the following components:
#' \itemize{
#' \item BMD: A data frame containing the BMD values and
#' the \eqn{100\times(1-2\alpha)\% } confidence intervals
#' for each response column.
#' \item summary: A list containing the summary statistics for each
#' repsonse column.
#' \item model_plots: A list containing the model plots for each
#' response column.
#' \item cleveland_plots: A list containing the Cleveland
#' plots for each response column.
#' \item models: A list containing the model fit values for
#' each response column. Values include:
#' \itemize{
#' \item Individual_Model_X: Individual model fits for each model,
#' X, used for model averaging. See \link[ToxicR]{single_continuous_fit}
#' for details on the model object class structure.
#' \item ma_bmd: The CDF of the model averaged BMD distribution.
#' \item posterior_probs: the posterior probabilities used in the
#' model averaging.
#' }
#' }
#' @importFrom dplyr select rename
#' @import ggplot2
#' @export
bmd_ma <- function(mf_data,
                   data_type = c("individual", "summary"),
                   dose_col,
                   response_cols = c("sample_MF_unique", "sample_MF_clonal"),
                   sd_col = NULL,
                   n_col = NULL,
                   bmr_type = "rel",
                   bmr = 0.5,
                   fit = "laplace",
                   a = 0.025,
                   ...) {

  if (!requireNamespace("ToxicR", quietly = TRUE)) {
    stop("ToxicR is not installed. Please install from https://github.com/NIEHS/ToxicR")
  }


  if (data_type == "individual") {
    # Initialize empty lists to store results
    results_bmd <- list()
    results_summary <- list()
    results_model <- list()
    results_model_plots <- list()
    results_cleveland_plots <- list()

    for (i in seq_along(response_cols)) {
      model <- ToxicR::ma_continuous_fit(D = mf_data[, dose_col],
                                         Y = mf_data[, response_cols[i]],
                                         fit_type = fit,
                                         BMR_TYPE = bmr_type,
                                         BMR = bmr,
                                         alpha = a,
                                         ...
                                         )

      results_bmd[[response_cols[i]]] <- summary(model)$BMD
      results_summary[[response_cols[i]]] <- summary(model)
      results_model[[response_cols[i]]] <- model
      results_model_plots[[response_cols[i]]] <- plot(model)
      results_cleveland_plots[[response_cols[i]]] <- ToxicR::cleveland_plot(model)
    }
    results_bmd_df <- do.call(rbind, results_bmd)
  } else if (data_type == "summary") {
    response_mat <- mf_data %>%
      dplyr::select({{response_cols}}, {{n_col}}, {{sd_col}})
    response_mat <- response_mat %>% rename(Mean = {{response_cols}},
                                            N = {{n_col}},
                                            SD = {{sd_col}})
    response_mat <- as.matrix(response_mat)

    model <- ToxicR::ma_continuous_fit(D = mf_data[, dose_col],
                                       response_mat = response_mat,
                                       fit_type = fit,
                                       BMR_TYPE = bmr_type,
                                       BMR = bmr,
                                       alpha = a,
                                       ...)
    results_bmd <- summary(model)$BMD
    results_bmd_df <- data.frame(BMDL = results_bmd[1],
                                 BMD = results_bmd[2],
                                 BMDU = results_bmd[3])
    row.names(results_bmd_df) <- response_cols
    results_summary <- summary(model)
    results_model <- model
    results_model_plots <- plot(model)
    results_cleveland_plots <- ToxicR::cleveland_plot(model)
  }

  results_bmd_df <- as.data.frame(results_bmd_df) %>%
    dplyr::mutate(response = row.names(results_bmd_df))

  conf_int <- 100 * (1 - 2 * a)
  g <- ggplot2::ggplot(results_bmd_df,
                        ggplot2::aes(x = results_bmd_df$response,
                                     y = results_bmd_df$BMD)) +
    ggplot2::geom_point(size = 5) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = results_bmd_df$BMDL,
                                        ymax = results_bmd_df$BMDU),
                           width = 0.2) +
    ggplot2::labs(x = "Response", y = "BMD",
                  title = paste0("BMD with ",
                                 conf_int,
                                 "% Confidence Intervals")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                       vjust = 0.5,
                                                       hjust = 1),
                   axis.line = ggplot2::element_line(colour = "black"),
                   plot.title = ggplot2::element_text(hjust = 0.5))

  results_list <- list(BMD = results_bmd_df,
                       BMD_plot = g,
                       summary = results_summary,
                       model_plots = results_model_plots,
                       cleveland_plots = results_cleveland_plots,
                       models = results_model)

# TO DO:
 # calculate some kind of goodness of fit for models
return(results_list)
}