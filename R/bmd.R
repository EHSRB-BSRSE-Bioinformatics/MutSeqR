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
#' 
#' @return A list with the following components:
#' @importFrom ToxicR single_continuous_fit
mf_bmd <- function(mf_data,
                   dose_col = "dose",
                   response_col = c("sample_MF_unique", "sample_MF_clonal"),
                   model_types = c("hill", "exp-3", "exp-5", "power", "polynomial"),
                   BMR = 0.5) {

  result_list <- list() # Initialize a list to store the fit results for each response

  # Run the function for each response column
  for (response in response_col) {
    fit_results <- list() # Initialize a list to store the fit results
    plot_list <- list()  # Initialize a list to store the plots

    # Run the function for each model type
    for (model_type in model_types) {
      fit_results[[model_type]] <- ToxicR::single_continuous_fit(mf_data[, dose_col],
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
# AIC
}

#' Fit a model averaged continuous BMD model
#' @description This function fits a model averaged continuous BMD model
#' to dose-response data for mutation frequency.
#' @param mf_data A data frame containing the dose-response data. Data may
#' be individual for each sample or averaged over dose groups.
#' @param response_cols A character vector specifying the columns in mf_data
#' containing the response data.
#' @param dose_col A character string specifying the column in mf_data
#' containing the dose data.
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
#' @param BMR A numeric value specifying the benchmark response. The BMR is
#' defined in relation to the calculation requested in bmr_type. Default is 0.5.
#' @param fit A string specifying the method used to fit the model.
#' Options are ("laplace", "mle", or "mcmc"). Default is "laplace".
#' @param alpha The specified nominal coverage rate for computation of
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
#'      \item \code{"exp-aerts"}:         \eqn{f(x) = a(1 + (c-1)(1-\exp(-bx^{d}))) }
#'      \item \code{"invexp-aerts"}:      \eqn{f(x) = a(1 + (c-1)(\exp(-bx^{-d})))}
#'      \item \code{"hill-aerts"}:        \eqn{f(x) = a(1 + (c-1)(1-\frac{b^d}{b^d + x^d}))}
#'      \item \code{"lognormal-aerts"}:   \eqn{f(x) = a\left\{1 + (c-1)\left(\Phi( \ln(b) + d\times \ln(x))\right) \right\}}
#'      \item \code{"gamma-efsa"}:        \eqn{f(x) = a(1 + (c-1)(\Gamma(bx; d))) }
#'      \item \code{"LMS"}:               \eqn{f(x) = a(1 + (c-1)(1 - \exp(-bx - dx^2))) }
#'      \item \code{"probit-aerts"}:      \eqn{f(x) = c\left(\Phi(a + b\times x^d)\right) }
#'      \item \code{"logistic-aerts"}:    \eqn{f(x) = \frac{c}{1 + \exp(-a - b\times x^d)} }
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
#' \itemize {
#' \item Individual_Model_X: Individual model fits for each model,
#' X, used for model averaging. See \link[ToxicR]{single_continuous_fit}
#' for details on the model object class structure.
#' \item ma_bmd: The CDF of the model averaged BMD distribution.
#' \item posterior_probs: the posterior probabilities used in the
#' model averaging.
#' }}
#' @importFrom ToxicR single_continuous_fit
bmd_ma <- function(mf_data,
                  response_cols = c("sample_MF_unique", "sample_MF_clonal"),
                  dose_col = "dose",
                  bmr_type = "rel",
                  bmr = 0.5,
                  fit = "laplace",
                  alpha = 0.025,
                  ...) {
# Initialize empty lists to store results
results_bmd <- list()
results_summary <- list()
results_model <- list()
results_model_plots <- list()
results_cleveland_plots <- list()
# Loop over response columns
for (i in seq_along(response_cols)) {
  # Fit model
  model <- ToxicR::ma_continuous_fit(D = mf_data[, dose_col],
                                     Y = mf_data[, response_cols[i]],
                                     fit_type = fit,
                                     BMR_TYPE = bmr_type,
                                     BMR = bmr,
                                     alpha = alpha,
                                     ...
                                     )
  # Grab results
  results_bmd[[response_cols[i]]] <- summary(model)$BMD
  results_summary[[response_cols[i]]] <- summary(model)
  results_model[[response_cols[i]]] <- model
  # Create plot
  results_model_plots[[response_cols[i]]] <- plot(model)
  results_cleveland_plots[[response_cols[i]]] <- ToxicR::cleveland_plot(model)
}
# Combine all summaries into one dataframe
results_bmd_df <- do.call(rbind, results_bmd)
results_list <- list(BMD = results_bmd_df, 
                     summary = results_summary,
                     model_plots = results_model_plots,
                     cleveland_plots = results_cleveland_plots,
                     models = results_model)

# TO DO:
 # plot all BMDs with ggplot
 # calculate some kind of goodness of fit for models
                     
#g <- ggplot2::ggplot(mf_dat, ggplot2::aes(x = Response, y = BMD)) +
#  ggplot2::geom_point() +  # BMD points
#  ggplot2::geom_errorbar(ggplot2::aes(ggplot2::ymin = BMDL, ggplot2::ymax = BMDU), ggplot2::width = 0.2) +  # Error bars for confidence intervals
#  ggplot2::labs(x = "Response", y = "BMD", title = "BMD with Confidence Intervals") +  # Labels
#  ggplot2::theme_minimal()   +  # Minimal theme
#  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  # Rotate x-axis labels


return(results_list)
}
#results_list$BMD
#results_list$plots$sample_MF_unique
#results_list$plots$sample_MF_clonal
#results_list$summary$sample_MF_unique
#results_list$models$sample_MF_unique$`Indiv_exp-aerts_normal`$maximum



                          

