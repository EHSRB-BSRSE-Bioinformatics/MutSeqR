#' Perform linear modelling on mutation frequency for given
#' fixed and random effects
#'
#' `model_mf` will fit a linear model to analyse the effect(s) of given
#' factor(s)  on mutation frequency and perform specified pairwise comparisons.
#' This function will fit either a generalized linear model (\link[stats]{glm})
#' or, if supplied random effects, a generalized linear mixed-effects model
#' (\link[lme4]{glmer}). Pairwise comparisons are conducted using the doBy
#' library (\link[doBy]{esticon}) and estimates are then back-transformed. The
#' delta method is  employed to approximate the  back-transformed
#' standard-errors. A Sidak correction is applied to adjust p-values for
#' multiple comparisons.
#' @param mf_data The data frame containing the mutation frequency data.
#' Mutation counts and total sequencing depth should be summarized per sample
#' alongside columns for your fixed effects.
#' This data can be obtained using `calculate_mf(summary=TRUE)`.
#' @param fixed_effects The name(s) of the column(s) that will act as the
#' fixed_effects (factor/independent variable) for modelling mutation frequency.
#' @param test_interaction a logical value. Whether or not your model should
#' include the interaction between the `fixed_effects`.
#' @param random_effects The name of the column(s) to be analysed as a
#' random effect in the model. Providing this effect will cause the function to
#' fit a generalized linear mixed-effects model.
#' @param reference_level Refers to one of the levels within each of your
#' fixed_effects. The coefficient for the reference level will represent the
#' baseline effect. The coefficients of the other levels will be interpreted in
#' relation to the reference_level as deviations from the baseline effect.
#' @param muts The column containing the mutation count per sample.
#' @param total_count The column containing the sequencing depth per sample.
#' @param contrasts a data frame or a  filepath to a file that will
#' provide the information necessary to make pairwise comparisons between
#' groups. The table must consist of two columns. The first column will be a
#' group within your fixed_effects and the second column must be the group that
#' it will be compared to.  The values must correspond to entries in your
#' mf_data column for each fixed effect. Put the group that you expect to have
#' the higher mutation frequency in the 1st column and the group that you expect
#' to have a lower mutation frequency in the second column. For multiple fixed
#' effects, separate the levels of each `fixed_effect` of a group with a colon.
#' Ensure that all `fixed_effects` are represented in each entry for the table.
#' See `details` for examples.
#' @param cont_sep The delimiter for importing the contrast table file.
#' Default is tab-delimited.
#' @param ... Extra arguments for \link[stats]{glm}  or \link[lme4]{glmer}. The
#' `glmer` function is used when a `random_effect` is supplied, otherwise, the
#' model uses the `glm` function.
#' @details
#' `fixed_effects` are variables that have a direct and constant effect on the
#' dependent variable (ie mutation frequency).They are typically the
#' experimental factors or covariates of interest for their impact on the
#' dependent variable. One or more fixed_effect may be provided. If you are
#' providing more than one fixed effect, avoid using correlated variables;
#' each fixed effect must independently predict the dependent variable.
#' Ex. `fixed_effects = c("dose", "genomic_target", "tissue", "age", etc)`.
#'
#' Interaction terms enable you to examine whether the relationship between the
#' dependent and independent variable changes based on the value of another
#' independent variable. In other words, if an interaction is significant, then
#' the relationship between the fixed effects is not constant across all levels
#' of each variable. Ex. Consider investigating the effect of dose group and
#' tissue on mutation frequency. An interaction between dose and tissue would
#' capture whether the dose response differs between tissues.
#'
#' `random_effects` account for the unmeasured sources of statistical variance that
#' affect certain groups in the data. They help account for unobserved
#' heterogeneity or correlation within groups. Ex. If your model uses repeated
#' measures within a sample, `random_effects = "sample"`.
#'
#' Setting a `reference_level` for your fixed effects enhances the interpretability
#' of the model. Ex. Consider a `fixed_effect` "dose" with levels 0, 25, 50, and 100 mg/kg.
#' Intuitively, the reference_level would refer to the  negative control dose, "0"
#' since we are interested in testing how the treatment might change mutation
#' frequency relative to the control.
#'
#' Examples of `contrasts`:
#'
#' If you have a `fixed_effect` "dose" with dose groups 0, 25, 50, 100,
#' then the first column would contain the treated groups (25, 50, 100), while
#' the second column would be 0, thus comparing each treated group to the control group.
#'
#' 25 0
#'
#' 50 0
#'
#' 100 0
#'
#'
#' Alternatively, if you would like to compare mutation frequency between treated
#' dose groups, then the contrast table would look as follows, with the lower
#' dose always in the second column, as we expect it to have a lower mutation
#' frequency. Keeping this format aids in interpretability of the estimates for
#' the pairwise comparisons. Should the columns be reversed, with the higher
#' group in the second column, then the model will compute the fold-decrease
#' instead of the fold-increase.
#'
#' 100 25
#'
#' 100 50
#'
#' 50 25
#'
#'
#' Ex. Consider the scenario where the `fixed_effects `
#' are "dose" (0, 25, 50, 100) and "genomic_target" ("chr1", "chr2"). To compare
#' the three treated dose groups to the control for each genomic target, the
#' contrast table would look like:
#'
#' 25:chr1	0:chr1
#'
#' 50:chr1	0:chr1
#'
#' 100:chr1	0:chr1
#'
#' 25:chr2	0:chr2
#'
#' 50:chr2	0:chr2
#'
#' 100:chr2	0:chr2
#'
#' Troubleshooting: If you are having issues with convergence for your
#' generalized linear mixed-effects model, it may be advisable to increase
#' the tolerance level for convergence checking during model fitting. This
#' is done through the `control` argument for the `lme4::glmer` function. The
#' default tolerance is tol = 0.002. Add this argument as an extra argument
#' in the `model_mf` function.
#' Ex. `control = lme4::glmerControl(check.conv.grad = lme4::.makeCC("warning",
#'                                                              tol = 3e-3,
#'                                                              relTol = NULL))`
#' Alternate approach:
#' `control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))`
#' Similar approaches may be taken for glm models.
#' @returns Model results are output as a list. Included are:
#' - model_data: the supplied mf_data with added column for the Pearson's
#' residuals of the model.
#' - summary: the summary of the model.
#' - anova: the analysis of variance for models with two or more effects. \link[car]{Anova}`(model) `
#' - residuals_histogram: the Pearson's residuals plotted as a histogram. This is
#' used to check whether the variance is normally distributed. A symmetric
#' bell-shaped histogram, evenly distributed around zero indicates that the
#' normality assumption is likely to be true.
#' - residuals_qq_plot: the Pearson's residuals plotted in a quantile-quantile plot.
#'  For a normal distribution, we expect points to roughly follow the y=x line.
#' - point_estimates_matrix: the contrast matrix used to generate point-estimates for the fixed effects.
#' - point_estimates: the point estimates for the fixed effects.
#' - pairwise_comparisons_matrix: the contrast matrix used to conduct the pairwise comparisons specified in the `contrasts`.
#' - pairwise_comparisons: the results of pairwise comparisons specified in the `contrasts`.
#' @examples
#' # Example data  consists of 24 mouse bone marrow
#' # samples exposed to three doses of BaP alongside vehicle controls.
#' # Libraries were sequenced with Duplex Sequencing using
#' # the TwinStrand Mouse Mutagenesis Panel which consists of 20 2.4kb
#' # targets = 48kb of sequence. Example data can be retrieved from
#' # MutSeqRData, an ExperimentHub data package:
#' ## library(ExperimentHub)
#' ## eh <- ExperimentHub()
#' ## query(eh, "MutSeqRData")
#' # Mutation frequency data was precalculated using
#' ## mf_data_global <- calculate_mf(mutation_data = eh[["EH9861"]],
#' ##   cols_to_group = "sample",
#' ##   retain_metadata_cols = c("dose_group", "dose"))
#'
#' # We will model the effect of dose on mutation frequency min.
#' mf_example <- readRDS(system.file("extdata/Example_files/mf_data_global.rds",
#'   package = "MutSeqR"
#' ))
#' # We will compare all BaP dose groups to the control group
#' # using pairwise comparisons.
#' contrasts <- data.frame(
#'   col1 = c("12.5", "25", "50"),
#'   col2 = c("0", "0", "0")
#' )
#' # Fit the model
#' model <- model_mf(
#'   mf_data = mf_example,
#'   fixed_effects = "dose",
#'   reference_level = "0",
#'   muts = "sum_min",
#'   total_count = "group_depth",
#'   contrasts = contrasts
#' )
#' @importFrom magrittr %>%
#' @importFrom graphics abline boxplot hist par
#' @importFrom stats as.formula model.matrix qqnorm relevel residuals qqline glm
#' @importFrom dplyr mutate select case_when
#' @importFrom utils read.delim
#' @export
model_mf <- function(
    mf_data,
    fixed_effects,
    test_interaction = TRUE,
    random_effects = NULL,
    reference_level,
    muts = "sum_min",
    total_count = "group_depth",
    contrasts = NULL,
    cont_sep = "\t",
    ...) {

  # --- 1. Input Validation & Setup ---
  muts <- match.arg(muts, choices = colnames(mf_data))
  total_count <- match.arg(total_count, choices = colnames(mf_data))

  if (!requireNamespace("doBy", quietly = TRUE)) stop("Package doBy is required.")
  if (length(fixed_effects) > 1 && !requireNamespace("car", quietly = TRUE)) stop("Package car is required.")
  if (!is.null(random_effects) && !requireNamespace("lme4", quietly = TRUE)) stop("Package lme4 is required.")

  # numeric conversion to prevent integer overflow
  mf_data[c(muts, total_count)] <- lapply(mf_data[c(muts, total_count)], as.numeric)
  
  # factor conversion
  mf_data[fixed_effects] <- lapply(mf_data[fixed_effects], as.factor)

  # --- 2. Reference Level Validation & Assignment ---
  # Check validity of reference levels
  for (i in seq_along(fixed_effects)) {
    f_name <- fixed_effects[i]
    ref <- reference_level[i]
    if (!ref %in% levels(mf_data[[f_name]])) {
      stop("Invalid reference level '", ref, "' for factor '", f_name, "'.")
    }
    message("Reference level for factor ", f_name, ": ", ref)
  }

  # Releveling
  # Map applies relevel to each column in the list using the corresponding ref
  # We force as.character() to ensure relevel matches by Name, not by Index
  mf_data[fixed_effects] <- Map(stats::relevel,
                                x = mf_data[fixed_effects],
                                ref = as.list(as.character(reference_level)))

  # --- 3. Model Fitting ---
  # Construct formula
  op <- if (test_interaction) "*" else "+"
  rhs <- paste(fixed_effects, collapse = op)

  if (!is.null(random_effects)) {
    # GLMM
    rhs <- paste(rhs, paste0("(1|", random_effects, ")"), sep = " + ")
    formula_str <- paste0("cbind(", muts, ", ", total_count, ") ~ ", rhs)
    message("Fitting GLMM: lme4::glmer(", formula_str, ", family = binomial)")
    model <- lme4::glmer(stats::as.formula(formula_str), family = "binomial", data = mf_data, ...)
  } else {
    # GLM
    formula_str <- paste0("cbind(", muts, ", ", total_count, ") ~ ", rhs)
    message("Fitting GLM: glm(", formula_str, ", family = quasibinomial)")

    model_formula <- stats::as.formula(formula_str)
    model <- stats::glm(model_formula, family = "quasibinomial", data = mf_data, ...)
    # Check dispersion
    if (summary(model)$dispersion < 1) {
      warning("Dispersion < 1. Switching to binomial distribution.")
      model <- stats::glm(model_formula, family = "binomial", data = mf_data, ...)
    }
  }

  model_summary <- summary(model)
  model_anova <- if (length(fixed_effects) > 1) car::Anova(model) else NULL

  # --- 4. Residuals ---
  mf_data$residuals <- stats::residuals(model, type = "pearson")

  max_resid_idx <- which.max(abs(mf_data$residuals))
  message("Max absolute residual: ", mf_data$residuals[max_resid_idx], " (Row ", max_resid_idx, ")")

  # Plots
  graphics::par(las = 1, xaxs = "i", yaxs = "i")
  hist_obj <- graphics::hist(mf_data$residuals, plot = FALSE) # Capture data
  graphics::hist(mf_data$residuals, main = "Residuals", col = "yellow", ylim = c(0, max(hist_obj$counts) + 1))

  qq <- stats::qqnorm(mf_data$residuals, main = "QQ Plot of Residuals")
  stats::qqline(mf_data$residuals, col = "red")

  # --- 5. Point Estimates (Vectorized) ---
  fixed_effects_levels <- lapply(mf_data[fixed_effects], levels)
  design_matrix <- do.call(expand.grid, fixed_effects_levels)

  # Create simplified model formula for model.matrix
  fixed_formula <- stats::as.formula(paste("~", paste(fixed_effects, collapse = op)))
  model_matrix_est <- stats::model.matrix(fixed_formula, data = design_matrix)

  # Create readable rownames (e.g., "LevelA:LevelB")
  row_names <- do.call(paste, c(design_matrix, sep = ":"))
  rownames(model_matrix_est) <- row_names

  # Compute estimates
  est_obj <- doBy::esticon(obj = model, L = model_matrix_est)

  # Process results
  model_estimates <- as.data.frame(est_obj)
  model_estimates <- model_estimates %>%
    dplyr::mutate(
      estimate_raw = .data$estimate, # keep raw log-odds if needed
      Estimate = exp(.data$estimate),
      Lower = exp(.data$lwr),
      Upper = exp(.data$upr),
      # Delta method approx for SE back-transform: SE_p = p * SE_logit
      Std.Err = .data$Estimate * .data$std.error 
    ) %>%
    dplyr::select("Estimate", "Std.Err", "Lower", "Upper")

  # Add fixed effect levels as indiv columns to the p estimate table
  if (length(fixed_effects) > 1) {
    meta_cols <- do.call(rbind, strsplit(rownames(model_estimates), ":"))
  } else {
    meta_cols <- matrix(rownames(model_estimates), ncol = 1)
  }
  colnames(meta_cols) <- fixed_effects
  model_estimates <- cbind(model_estimates, as.data.frame(meta_cols))


  # --- 6. Pairwise Comparisons (Vectorized) ---
  pairwise_results <- NULL
  result_matrix <- NULL

  if (!is.null(contrasts)) {
    # Load contrasts
    if (is.data.frame(contrasts)) {
      contrast_table <- contrasts
    } else {
      contrast_table <- read.delim(contrasts, sep = cont_sep, header = FALSE)
    }
    
    # Check dimensions
    if (ncol(contrast_table) != 2) stop("Contrast table must have exactly 2 columns.")
    
    # Validate levels
    all_contrast_levels <- unique(unlist(contrast_table))
    # Note: row_names variable from above contains all valid combinations (e.g. "25:chr1")
    if (!all(all_contrast_levels %in% row_names)) {
      stop("Contrast table contains values not found in model fixed effects combinations.")
    }

    # MATRIX SUBTRACTION (Vectorized)
    V1 <- as.character(contrast_table[, 1])
    V2 <- as.character(contrast_table[, 2])

    # Direct matrix algebra instead of loop
    result_matrix <- model_matrix_est[V1, , drop = FALSE] - model_matrix_est[V2, , drop = FALSE]
    rownames(result_matrix) <- paste(V1, "vs", V2)

    # Perform comparisons
    pc_obj <- doBy::esticon(obj = model, L = result_matrix)

    # Clean results
    pairwise_results <- as.data.frame(pc_obj) %>%
      dplyr::mutate(
        Fold.Change = exp(.data$estimate),
        FC.Lower = exp(.data$lwr),
        FC.Upper = exp(.data$upr),
        FC.Std.Err = .data$Fold.Change * .data$std.error, # Delta method
        adj_p.value = MutSeqR::sidak(.data$p.value)$SidakP,
        Significance = dplyr::case_when(
          adj_p.value <= 0.001 ~ "***",
          adj_p.value <= 0.01 ~ "**",
          adj_p.value <= 0.05 ~ "*",
          TRUE ~ ""
        )
      ) %>%
      dplyr::select("Fold.Change", "FC.Std.Err", Obs.T = "statistic", "p.value", "df", "FC.Lower", "FC.Upper", "adj_p.value", "Significance")

    # Parse Contrast Names back into indiv columns
    # Split "A vs B"
    groups <- do.call(rbind, strsplit(rownames(pairwise_results), " vs "))

    # Split Group 1 parts (e.g., "25:chr1" -> "25", "chr1")
    if (length(fixed_effects) > 1) {
      g1_parts <- do.call(rbind, strsplit(groups[, 1], ":"))
      g2_parts <- do.call(rbind, strsplit(groups[, 2], ":"))
    } else {
      g1_parts <- matrix(groups[, 1], ncol = 1)
      g2_parts <- matrix(groups[, 2], ncol = 1)
    }

    colnames(g1_parts) <- paste0(fixed_effects, "_1")
    colnames(g2_parts) <- paste0(fixed_effects, "_2")

    pairwise_results <- cbind(pairwise_results, as.data.frame(g1_parts), as.data.frame(g2_parts))
  }

  # --- 7. Return ---
  out <- list(
    model = model,
    model_data = mf_data,
    model_formula = formula_str,
    summary = model_summary,
    residuals_histogram = hist_obj,
    residuals_qq_plot = qq,
    point_estimates_matrix = model_matrix_est,
    point_estimates = model_estimates
  )

  if (!is.null(model_anova)) out$anova <- model_anova
  if (!is.null(pairwise_results)) {
    out$pairwise_comparisons_matrix <- result_matrix
    out$pairwise_comparisons <- pairwise_results
  }

  return(out)
}