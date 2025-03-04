#' BMD modeling using ToxicR
#' @description Calculate the benchmark dose (BMD) for continuous dose-response
#' data with optional model averaging. This function is intended to model the
#' dose-response of mutation frequency using the ToxicR software.
#' @param mf_data A data frame containing the dose-response data. Data may
#' be individual for each sample or averaged over dose groups.
#' Required columns for individual data are the column containing the dose
#' `dose_col` and the column(s) containing the mutation frequency data
#' `response_col`(s). Summary data must include the `dose_col`, the
#' `response_col`(s) containing the mean response for each dose group,
#' the `sd_col` containing the standard deviation of the response data,
#' and the `n_col` containing the sample size for each dose group.
#' @param data_type A string specifying the type of response data.
#' Data may be response per individual or summarised across dose groups.
#' Options are ("individual", "summary"). Default is "individual".
#' @param dose_col The column in `mf_data` containing the dose data. Values
#' must be numeric. Default is "dose".
#' @param response_col The column(s) in mf_data containing the mutation
#' frequency data. For summarised data types, this should be the mean response
#' for each dose group. Multiple `response_col`s can be provided.
#' @param sd_col The column in mf_data containing the standard deviation of
#' the summarised response data. This is only required for
#' `data_type = "summary"`. If multiple response columns are provided,
#' multiple `sd_col`s should be provided in the same order. Default is NULL.
#' @param n_col The column in mf_data containing the sample size of each dose
#' group. This is only required for `data_type = "summary"`. If multiple
#' response columns are provided, multiple `n_col`s should be provided in the
#' same order. Default is NULL.
#' @param bmr_type The type of benchmark response. Options are: "rel", "sd",
#' "hybrid", "abs". Default is "rel". See details for more information.
#' @param bmr A numeric value specifying the benchmark response. The bmr is
#' defined in relation to the calculation requested in bmr_type. Default is 0.5.
#' @param model The model type to use. Options are "all" or
#' a vector of model types. Default is "exp-aerts", the Exponential model.
#' See details for available models. Note that model averaging will use
#' a pre-defined model set. See details for more information.
#' @param model_averaging A logical value indicating whether to use model
#' averaging. Default is TRUE (recommended).
#' @param alpha The specified nominal coverage rate for computation of
#' the lower and upper confidence intervals for the benchmark dose
#' (BMDL, BMDU). The confidence level is calculated as
#' \eqn{100\times(1-2\alpha)\% }. The default is 0.05 (90% CI).
#' @param plot_results A logical value indicating whether to plot the BMD models
#' and/or the Cleveland plots. Default is FALSE. If TRUE, the function will
#' save plots to the `output_path`.
#' @param ma_summary A logical value indicating whether to return the summary
#' of the model averaging results. Default is FALSE.
#' @param output_path The file path indicating where to save the plots.
#' If NULL, the plots will be saved to the working directory. Default is NULL.
#' @param ... Additional arguments to be passed to the model fitting function.
#' For more information, see \link[ToxicR]{single_continuous_fit} or
#' \link[ToxicR]{ma_continuous_fit} if model averaging.
#' @details Available model types for single model fitting are:
#' \itemize{
#'    \item \code{"exp-aerts"}: \eqn{f(x) = a(1 + (c-1)(1-\exp(-bx^{d}))) }
#'    \item \code{"invexp-aerts"}: \eqn{f(x) = a(1 + (c-1)(\exp(-bx^{-d})))}
#'    \item \code{"gamma-aerts"}: \eqn{f(x) = a(1 + (c-1)(Gamma(bx^d;xi)))}
#'    \item \code{"invgamma-aerts"}:
#' \eqn{f(x) = a(1 + (c-1)(1-Gamma(bx^{-d};xi)))}
#'    \item \code{"hill-aerts"}:
#' \eqn{f(x) = a(1 + (c-1)(1-\frac{b^d}{b^d + x^d}))}
#'    \item \code{"lomax-aerts"}:
#' \eqn{f(x) = a\left\{1 + (c-1)(1-\left(\frac{b}{b+x^d} \right)^\xi) \right\}}
#'    \item \code{"invlomax-aerts"}:
#'\eqn{f(x) = a\left\{1 + (c-1)(\left(\frac{b}{b+x^{-d}} \right))^\xi \right\}}
#'    \item \code{"lognormal-aerts"}:
#' \eqn{f(x) = a\left\{1 + (c-1)\left(\Phi( \ln(b) + d\times \ln(x))\right) \right\}}
#'    \item \code{"logskew-aerts"}:
#' \eqn{f(x) = a\left\{1 + (c-1)\left(\Phi_{SN}( \ln(b) + d\times \ln(x); \xi )\right) \right\}}
#'    \item \code{"invlogskew-aerts"}:
#' \eqn{f(x) = a\left\{1 + (c-1)\left(1 - \Phi_{SN}( \ln(b) - d\times \ln(x); \xi )\right) \right\}}
#'    \item \code{"logistic-aerts"}:
#' \eqn{f(x) = \frac{c}{1 + \exp(-a - b\times x^d)} }
#'    \item \code{"probit-aerts"}:
#' \eqn{f(x) = c\left(\Phi(a + b\times x^d)\right) }
#'    \item \code{"LMS"}: \eqn{f(x) = a(1 + (c-1)(1 - \exp(-bx - dx^2)))}
#'    \item \code{"gamma-efsa"}: \eqn{f(x) = a(1 + (c-1)(Gamma(bx; d)))}
#' }
#' Here: \eqn{\Phi(\cdot)} is the standard normal distribution and
#'        \eqn{\Phi_{SN}(\cdot;\cdot)} is the skew-normal distribution.
#' See \link[ToxicR]{single_continuous_fit} for more details.
#'
#'
#' Model averaging is done over the the model set described in The European
#' Food Safety Authority's (2022) Guidance on the use of the benchmark dose
#' approach in risk assessment. These models are (normal then lognormal for
#' each model): exp-aerts, invexp-aerts, hill-aerts, lognormal-aerts,
#' gamma-efsa, LMS, probit-aerts, and logistic-aerts. See
#' \link[ToxicR]{ma_continuous_fit} for more details.
#'
#'
#'  BMR types for continuous models:
#'  \itemize{
#'    \item Relative deviation (default; `bmr_type = "rel"`). This defines the
#' BMD as the dose that changes the control mean/median a certain
#' percentage from the background dose. It is the dose that
#' solves \eqn{\mid f(dose) - f(0) \mid = (1 \pm BMR) f(0)}
#'    \item Standard deviation (`bmr_type = "sd"`). This defines the BMD as
#' the dose associated with the mean/median changing a specified number of
#' standard deviations from the mean at the control dose. It is the
#' dose that solves \eqn{\mid f(dose)-f(0) \mid = BMR \times \sigma}
#'    \item Absolute deviation (`bmr_type="abs"`). This defines the BMD
#' as an absolute change from the control dose of zero by a specified amount.
#' That is the BMD is the dose that solves the equation
#' \eqn{\mid f(dose) - f(0) \mid = BMR}.
#'    \item Hybrid deviation (`bmr_type = "hybrid"`). This defines the
#' BMD that changes the probability of an adverse event by a stated amount
#' relative to no exposure (i.e 0).  That is, it is the dose that solves
#' \eqn{\frac{Pr(X > x| dose) - Pr(X >x|0)}{Pr(X < x|0)} = BMR}.
#' For this definition,
#' \eqn{Pr(X < x|0) = 1 - Pr(X > X|0) = \pi_0}, where \eqn{0 \leq \pi_0 < 1}
#' is defined by the user as "point_p," and it defaults to 0.01.  Note: this
#' discussion assumed increasing data. The fitter determines the direction
#' of the data and inverts the probability statements for decreasing data.
#' }
#' @return If `model_averaging = FALSE`, the function returns a data frame
#' with the BMD values and the \eqn{100\times(1-2\alpha)\% } confidence
#' intervals (BMDL, BMDU)for each response column and each model listed. The
#' AIC value is calculated for each model to compare fits. The AIC is
#' calculated as maximum likelihood + 2 * degrees of freedom. If
#' `plot_results = TRUE`, the function will plot all fitted models to the
#' data and save them to the `output_path`.
#'
#'
#'  If `model_averaging = TRUE`, the function returns a data frame with the
#' BMD values and the \eqn{100\times(1-2\alpha)\% } confidence intervals
#' (BMDL, BMDU) for each response column calculated using model averaging.
#' If `ma_summary = TRUE`, the function will return the posterior probabilities
#' used in the model averaging. If `plot_results = TRUE`, the function will plot
#' the model averaged model to the data and save it to the `output_path`.
#' The function will also make a Cleveland plot, saved to the  `output_path`.
#' Here, the BMDs are plotted for each model in the set alongside the model
#' averaged BMD. The BMD is represented by a red dot. The size of the dot is
#' scaled on the model probability with the Model Average having a value of
#' 100%. The BMDL and BMDU are expressed as interval bars.
#' @examples
#' # Calculate the BMD for a 50% increase in mutation frequency from control
#' # Individual data with Model averaging.
#' example_file <- system.file("extdata", "example_mutation_data_filtered.rds", package = "MutSeqR")
#' example_data <- readRDS(example_file)
#' mf <- calculate_mf(example_data, retain_metadata_cols = "dose")
#' bmd <- bmd_toxicr(mf_data = mf,
#'                   dose_col = "dose",
#'                   response_col = c("mf_min", "mf_max"))
#' # Plot the results using plot_ci()
#' plot <- plot_ci(bmd, order = "asc", log_scale = FALSE)
#' 
#' # Summary data with Model averaging.
#' mf_sum <- mf %>%
#'  dplyr::group_by(dose) %>%
#'  dplyr::summarise(mean_mf_min = mean(mf_min), sd_min = sd(mf_min), n_min = n(),
#'                   mean_mf_max = mean(mf_max), sd_max = sd(mf_max), n_max = n())
#' bmd <- bmd_toxicr(mf_data = mf_sum,
#'                   dose_col = "dose",
#'                   response_col = c("mean_mf_min", "mean_mf_max"),
#'                   sd_col = c("sd_min", "sd_max"),
#'                   n_col = c("n_min", "n_max"))
#' @importFrom dplyr select rename if_else
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @export
bmd_toxicr <- function(mf_data,
                       data_type = "individual",
                       dose_col = "dose",
                       response_col = c("mf_min", "mf_max"),
                       sd_col = NULL,
                       n_col = NULL,
                       bmr_type = "rel",
                       bmr = 0.5,
                       model = "exp-aerts",
                       alpha = 0.05,
                       model_averaging = TRUE,
                       plot_results = FALSE,
                       ma_summary = FALSE,
                       output_path = NULL,
                       ...) {
  if (!requireNamespace("ToxicR", quietly = TRUE)) {
    stop("ToxicR is not installed. Please install from https://github.com/NIEHS/ToxicR")
  }
  # Output directory for plots
  if (plot_results) {
    if (is.null(output_path)) {
      output_dir <- file.path(here::here())
    } else {
      output_dir <- file.path(output_path)
    }
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
  }
  # Data: Individual or Summary
  if (data_type == "summary") {
    create_matrices_list <- function(df, response_col, sd_col, n_col) {
      # Ensure the parameters are provided as character vectors
      response_cols <- as.character(response_col)
      sd_cols <- as.character(sd_col)
      n_cols <- as.character(n_col)

      matrices_list <- setNames(
        lapply(seq_along(response_cols), function(i) {
          matrix <- df %>%
            dplyr::select(dplyr::all_of(response_cols[i]),
                          dplyr::all_of(n_cols[i]),
                          dplyr::all_of(sd_cols[i])) %>%
            dplyr::rename(Mean = dplyr::all_of(response_cols[i]),
                          N = dplyr::all_of(n_cols[i]),
                          SD = dplyr::all_of(sd_cols[i])) %>%
            as.matrix()

          return(matrix)
        }),
        response_cols
      )

      return(matrices_list)
    }
    matrices_list <- create_matrices_list(mf_data, response_col, sd_col, n_col)
  }

  if (!model_averaging) {
    if (model == "all") {
      model_type <- c("exp-aerts", "invexp-aerts", "gamma-aerts",
                      "invgamma-aerts", "hill-aerts", "lomax-aerts",
                      "invlomax-aerts", "lognormal-aerts", "logskew-aerts",
                      "invlogskew-aerts", "logistic-aerts", "probit-aerts",
                      "LMS", "gamma-efsa")
    } else {
      model_type <- model
    }

    result_list <- lapply(response_col, function(response) {
      fit_results <- list() # Initialize a list to store the fit results

      # Select the appropriate y values based on data_type
      y <- if (data_type == "summary") matrices_list[[response]] else mf_data[, response]

      # Fit models using lapply
      fit_results <- lapply(model_type, function(modeltype) {
        fit <- ToxicR::single_continuous_fit(mf_data[, dose_col],
                                             y,
                                             model_type = modeltype,
                                             BMR = bmr,
                                             BMR_TYPE = bmr_type,
                                             alpha = alpha,
                                             ...)
        # Plot results if needed
        if (plot_results) {
          plot_file <- file.path(output_dir,
                                 paste0(response, "_", modeltype, ".svg"))
          plot <- plot(fit)
          ggsave(plot_file, plot, width = 7, height = 5, dpi = 300)
        }

        return(fit)
      })

      names(fit_results) <- model_type

      # Extract BMD values and create a data frame
      # Using lapply and dplyr::bind_rows to include model_type and AIC
      bmd_values <- lapply(model_type, function(modeltype) {
        fit <- fit_results[[modeltype]]
        bmd_df <- data.frame(Model = modeltype,
                             BMD = fit$bmd[1],
                             BMDL = fit$bmd[2],
                             BMDU = fit$bmd[3],
                             AIC = fit$maximum + 2 * summary(fit)$GOF[1, 2])
        # AIC <- model_indiv$maximum + 2*summary(model_indiv)$GOF[1,2] ## I don't remember where I found this....
        rownames(bmd_df) <- NULL
        return(bmd_df)
      })
      result_df <- dplyr::bind_rows(bmd_values)
      result_df$Response <- response
      result_df$confidence_level <- paste0((1 - 2 * alpha) * 100, "%")
      return(result_df)
    })
    results_single <- do.call(rbind, result_list)
    return(results_single)
  } else {
    results_list <- lapply(response_col, function(response) {
     y <- if (data_type == "summary") matrices_list[[response]] else mf_data[, response]
      fit <- ToxicR::ma_continuous_fit(mf_data[, dose_col],
                                       y,
                                       BMR = bmr,
                                       BMR_TYPE = bmr_type,
                                       alpha = alpha,
                                       ...)
      bmd <- data.frame(model_names = "Model Averaging",
                        BMD = fit$bmd[1],
                        BMDL = fit$bmd[2],
                        BMDU = fit$bmd[3])
      if (ma_summary) {
        bmd$post_p <- NA
        summary <- summary(fit)$fit_table
        results <- rbind(summary, bmd)
      } else {
        results <- bmd
      }
      results$Response <- response
      results$confidence_level <- paste0((1 - 2 * alpha) * 100, "%")
      results <- dplyr::rename(results, "Model" = "model_names")
      rownames(results) <- NULL
      # Plot results: MA model + cleveland plot
      if (plot_results) {
        plot_file <- file.path(output_dir, paste0(response, "_ma.svg"))
        plot <- plot(fit) +
          ggplot2::theme(panel.background = ggplot2::element_blank(),
                         axis.line = ggplot2::element_line(),
                         panel.grid = ggplot2::element_blank(),
                         axis.ticks.x = ggplot2::element_line())
        ggsave(plot_file, plot, width = 7, height = 5, dpi = 300)
        cplot_file <- file.path(output_dir, paste0(response, "_cleveland.svg"))
        cplot <- ToxicR::cleveland_plot(fit) +
          ggplot2::theme(panel.background = ggplot2::element_blank(),
                         axis.line = ggplot2::element_line(),
                         panel.grid = ggplot2::element_blank(),
                         axis.ticks.x = ggplot2::element_line(),
                         axis.ticks.y = ggplot2::element_line()) +
          ggplot2::xlab("BMD Estimate") +
          ggplot2::ylab("Model") +
          ggplot2::ggtitle("BMD per Model")
        ggsave(cplot_file, cplot, width = 7, height = 5, dpi = 300)
      }
      return(results)
    })
    results_ma <- do.call(rbind, results_list)
    return(results_ma)
  }
}