#' Perform linear modelling on mutation frequency for given
#' fixed and random effects
#'
#' `model_mf` will fit a linear model to analyse the effect(s) of given factor(s) 
#' on mutation frequency and perform specified pairwise comparisons. This function
#' will fit either a generalized linear model (\link[stats]{glm}) or, if supplied 
#' random effects, a generalized linear mixed-effects model (\link[lme4]{glmer}) . 
#' Pairwise comparisons are conducted using the doBy library (\link[doBy]{esticon}) 
#' and estimates are then back-transformed. The delta method is  employed to
#' approximate the  back-transformed standard-errors. A Sidak  correction
#' is applied to adjust p-values for multiple comparisons.
#' @param mf_data The data frame containing the mutation frequency data.
#' Mutation counts and total sequencing depth should be summarized per sample
#' alongside columns for your fixed effects.
#' This data can be obtained using `calculate_mut_freq(summary=TRUE)`.
#' @param fixed_effects The name(s) of the column(s) that will act as the
#' fixed_effects (factor/independent variable) for modelling mutation frequency.
#' @param test_interaction a logical value. Whether or not your model should include 
#' the interaction between the `fixed_effects`.
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
#' provide the information necessary to make pairwise comparisons between groups. 
#' The table must consist of two columns. The first column will be a group within 
#' your fixed_effects and the second column must be the group that it will be 
#' compared to.  The values must correspond to entries in your 
#' mf_data column for each fixed effect. Put the group that you expect to have the higher 
#' mutation frequency in the 1st column and the group that you expect to have a 
#' lower mutation frequency in the second column. For multiple fixed effects, 
#' separate the levels of each `fixed_effect` of a group with a colon. Ensure 
#' that all `fixed_effects` are represented in each entry for the table. 
#' See `details` for examples.
#' @param cont_sep The delimiter for importing the contrast table file. Default is tab-delimited.
#' @param ... Extra arguments for \link[stats]{glm}  or \link[lme4]{glmer}. The 
#' `glmer` function is used when a `random_effect` is supplied, otherwise, the 
#' model uses the `glm` function.
#' @export
#' 
#' @details 
#' 
#' `fixed_effects` are variables that have a direct and constant effect on the 
#' dependent variable (ie mutation frequency).They are typically the experimental 
#' factors or covariates of interest for their impact on the dependent variable.
#' One or more fixed_effect may be provided. If you are providing more than one 
#' fixed effect, avoid using correlated variables; each fixed effect must 
#' independently predict the dependent variable. 
#' Ex. `fixed_effects = c("dose", "genomic_target", "tissue", "age", etc)`.
#' 
#' Interaction terms enable you to examine whether the relationship between the
#' dependent and independent variable changes based on the value of another 
#' independent variable. In other words, if an interaction is significant, then 
#' the relationship between the fixed effects is not constant across all levels 
#' of each variable. Ex. Consider investigating the effect of dose group and tissue
#' on mutation frequency. An interaction between dose and tissue would capture 
#' whether the dose response differs between tissues. 
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
#' Troubleshooting: If you are having issues with convergence for your generalized linear mixed-
#' effects model, it may be advisable to increase the tolerance level for convergence
#' checking during model fitting. This is done through the `control` argument for 
#' the `lme4::glmer` function. The default tolerance is tol = 0.002. Add this
#' argument as an extra argument in the `model_mf` function.
#' Ex. `control = lme4::glmerControl(check.conv.grad = lme4::.makeCC("warning",
#'                                                              tol = 3e-3,
#'                                                              relTol = NULL))`
#' @returns Model results are output as a list. Included are:
#' - model_data: the supplied mf_data with added column for model residuals.
#' - summary: the summary of the model.
#' - anova: the analysis of variance for models with two or more effects. \link[car]{Anova}`(model) `
#' - residuals_histogram: the model residuals plotted as a histogram. This is
#' used to check whether the variance is normally distributed. A symmetric
#' bell-shaped histogram, evenly distributed around zero indicates that the
#' normality assumption is likely to be true.
#' - residuals_qq_plot: the model residuals plotted in a quantile-quantile plot.
#'  For a normal distribution, we expect points to roughly follow the y=x line.  
#' - point_estimates_matrix: the contrast matrix used to generate point-estimates for the fixed effects. 
#' - point_estimates: the point estimates for the fixed effects.
#' - pairwise_comparisons_matrix: the contrast matrix used to conduct the pairwise comparisons specified in the `contrasts`.
#' - pairwise_comparisons: the results of pairwise comparisons specified in the `contrasts`.
#' @importFrom magrittr %>%
#' @importFrom doBy esticon
#' @importFrom lme4 glmer
#' @importFrom car Anova
#' @importFrom graphics abline boxplot hist par
#' @importFrom stats as.formula model.matrix qqnorm relevel residuals
#' @export
model_mf <- function(mf_data,
  fixed_effects,
  test_interaction = TRUE,
  random_effects = NULL,
  reference_level,
  muts = "sample_sum_min",
  total_count = "sample_group_depth",
 # family = NULL, # Note: back-transforming the estimates requires the family to be binomial
  contrasts = NULL,
  cont_sep = "\t",
  ...
) {
  # Convert muts and total_count to numeric to avoid integer overflow
  mf_data[[muts]] <- as.numeric(mf_data[[muts]])
  mf_data[[total_count]] <- as.numeric(mf_data[[total_count]])

  # Convert specified columns to factors
  mf_data[, fixed_effects] <- lapply(mf_data[, fixed_effects, drop = FALSE], as.factor)

    # Check that the reference levels are valid levels of the factors
  for (factor_name in fixed_effects) {
    factor_levels <- levels(mf_data[[factor_name]])
    reference_level_char <- as.character(reference_level[fixed_effects == factor_name])
    invalid_levels <- reference_level_char[!reference_level_char %in% factor_levels]
    if (length(invalid_levels) > 0) {
      stop(paste("Invalid reference level(s) for factor", factor_name, ":", paste(invalid_levels, collapse = ", ")))
    } else {
    message(paste("Reference level for factor", factor_name, ":", reference_level_char))
  }
}
  # Set the reference level for each factor in fixed_effects
  if (length(fixed_effects) == 1) {
    mf_data[[fixed_effects]] <- stats::relevel(mf_data[[fixed_effects]], ref = as.character(reference_level))
  } else {
    for (factor_name in as.list(fixed_effects)) {
    reference_level_for_factor <- reference_level[match(factor_name, fixed_effects)]
    mf_data[[factor_name]] <- relevel(mf_data[[factor_name]], ref = reference_level_for_factor)
  }
}
  # Construct the model formula
  if (test_interaction) {
    formula_str <- paste("cbind(", muts, ",", total_count, ") ~ ", paste(fixed_effects, collapse = "*"))
  } else {
    formula_str <- paste("cbind(", muts, ",", total_count, ") ~ ", paste(fixed_effects, collapse = "+"))
  }
  
  # Add random effects to the formula
  if (!is.null(random_effects)) {
   # Add random effect to model formula
    random_formula <- paste(paste("(1|", random_effects, ")", collapse = "+"))
    formula_str <- paste(formula_str, "+", random_formula)
    model_formula <- stats::as.formula(formula_str)
    
  # GLMM
  message(paste0("Fitting generalized linear mixed-effects model. lme4::glmer(", formula_str, ", family = ", family_param, ")"))

    model <- lme4::glmer(model_formula,
        family = "binomial",
        data = mf_data,
        ...
      )

  } else {
    model_formula <- stats::as.formula(formula_str)

    #GLM
    message(paste0("Fitting generalized linear model. glm(", formula_str, ", family = quasibinomial"))
    model <- stats::glm(model_formula,
                        family = "quasibinomial",
                        data = mf_data,
                        weights = get(total_count),
                        ...
                        )
  }
  
  model_summary <- summary(model)
  
  if(length(fixed_effects) > 1){
    model_anova <-car::Anova(model)}
 #################################
  # Check residuals
  if(!is.null(random_effects)){
   mf_data$residuals <- stats::residuals(model)
  } else {
    mf_data$residuals <- model$residuals
  }
 
   # Print the row with the maximum residual
  max_residual_index <- which.max(abs(mf_data$residuals))
  max_residual_row <- mf_data[max_residual_index, ]
  message("The row with the maximum residual in absolute value is:\n")
  print(max_residual_row)
  
  ## TO DO: Let's walk users through what the residuals should look like and give advice on them
  
  par(las = 1)
  hist <- hist(mf_data$residuals, main = "Residuals", col = "yellow")

  qqplot <- stats::qqnorm(mf_data$residuals, main = "QQ Plot of Residuals")
  stats::qqline(mf_data$residuals, col = "red")

## TO DO: Have a message about what these should look like.

  #########################################
  # Point Estimates By Fixed Effects
  ###########################################################
  # Define your levels for each factor
  fixed_effects_levels <- lapply(fixed_effects, function(factor) levels(mf_data[[factor]]))
  
  # Define all combinations of each factor
  design_matrix <- do.call(expand.grid, fixed_effects_levels)
  colnames(design_matrix) <- fixed_effects
  
   # Create simplified model formula
  if(test_interaction) {
    fixed_effect_formula <- stats::as.formula(paste(" ~ ", paste(fixed_effects, collapse = "*")))
  } else {
    fixed_effect_formula <- stats::as.formula(paste(" ~ ", paste(fixed_effects, collapse = "+")))
  }

  # Create the model matrix using model.matrix
  model_matrix <- stats::model.matrix(fixed_effect_formula, data = design_matrix)
  row_names <- apply(design_matrix, 1, paste, collapse = ":")
  rownames(model_matrix) <- row_names

  # Computed estimates
  model_estimates <- doBy::esticon(obj = model, L = model_matrix)

  model_estimates <- as.data.frame(model_estimates)
  model_estimates$estimate <- exp(model_estimates$estimate)
  delta <- model_estimates$estimate^2
  model_estimates$lwr <- exp(model_estimates$lwr)
  model_estimates$upr <- exp(model_estimates$upr)
  model_estimates$std.error <- sqrt(delta * model_estimates$std.error^2)
  # Clean data
  model_estimates <- model_estimates[, -c(3, 4, 5, 6)]
  colnames(model_estimates) <- c("Estimate", "Std.Err", "Lower", "Upper")

# Split the string into individual characters
chars <- strsplit(fixed_effects, " ")
# Count the characters
count <- length(unlist(chars))
for (i in seq_along(chars)) {
  # Assign each element to a separate variable in the global environment
  assign(paste("var", i, sep = ""), chars[[i]])
}
    # Extract rownames into a column for each contrast variable
for (i in 1:count) {
    var_i <- get(paste0("var", i))
    model_estimates[[var_i]] <- sapply(strsplit(rownames(model_estimates), ":"), "[", i)
}


  ##############################################################
 # Pairwise Comparisons  
 ##################################################################
  # load contrast table file and do checks
  if (!is.null(contrasts)) {
    if (is.data.frame(contrasts)) {
        contrast_table <- contrasts
    } else {
        contrast_table <- read.delim(file.path(contrasts), sep = cont_sep, header = F)
    if (ncol(contrast_table) <= 1) {
      stop("Your contrast_table only has one column. Make sure to set the proper delimiter with cont_sep.")
    }
  }
  model_matrix <- as.data.frame(model_matrix)
  contrast_table <- as.data.frame(contrast_table)  # Convert to data frame if needed
  
  # Create an empty list to store the result of matrix subtractions
  result_list <- list()
  
  # Loop through each row in contrast_table
  for (i in seq_len(nrow(contrast_table))) {
    
    # Get the model_row values
    V1 <- as.character(contrast_table[i, 1])
    V2 <- as.character(contrast_table[i, 2])
    
    # Perform matrix subtraction and store the result in the list
    result_list[[i]] <- model_matrix[V1, ] - model_matrix[V2, ]
  }
  
  # Convert the list of matrices to a single matrix
  result_matrix <- as.matrix(do.call(rbind, result_list))
  # Set row names for result_matrix
  rownames(result_matrix) <- paste(contrast_table[, 1], "vs", contrast_table[, 2])
  
 # Perform comparisons
  pairwise_comparisons <- doBy::esticon(obj = model, L = result_matrix)

  # Clean results  
  pairwise_comparisons <- as.data.frame(pairwise_comparisons)
  pairwise_comparisons$estimate <- exp(pairwise_comparisons$estimate)
  delta <- pairwise_comparisons$estimate^2
  pairwise_comparisons$lwr <- exp(pairwise_comparisons$lwr)
  pairwise_comparisons$upr <- exp(pairwise_comparisons$upr)
  pairwise_comparisons$std.error <- sqrt(delta*pairwise_comparisons$std.error^2)
  pairwise_comparisons <- pairwise_comparisons[,-5]
  colnames(pairwise_comparisons) <- c("Fold.Change", "FC.Std.Err", "Obs.T", "p.value", "df", "FC.Lower", "FC.Upper")
  
  pairwise_comparisons$adj_p.value <- MutSeqR::sidak(pairwise_comparisons$p.value)$SidakP
  pairwise_comparisons <- pairwise_comparisons %>%
  dplyr::mutate(
    Significance = case_when(
      adj_p.value <= 0.001 ~ "***",
      adj_p.value <= 0.01 ~ "**",
      adj_p.value <= 0.05 ~ "*",
      TRUE ~ ""
    )
  )

pairwise_comparisons$contrast_group1 <- sapply(strsplit(rownames(pairwise_comparisons), " vs "), "[", 1)
pairwise_comparisons$contrast_group2 <- sapply(strsplit(rownames(pairwise_comparisons), " vs "), "[", 2)
for (i in 1:count) {
    var_i <- get(paste0("var", i))
    pairwise_comparisons[[paste0(var_i, "_1")]] <- sapply(strsplit(pairwise_comparisons$contrast_group1, ":"), "[", i)
    pairwise_comparisons[[paste0(var_i, "_2")]] <- sapply(strsplit(pairwise_comparisons$contrast_group2, ":"), "[", i)
}
pairwise_comparisons <- pairwise_comparisons %>%
  dplyr::select(-contrast_group1, -contrast_group2)
  }
  
  model_results <- list(model = model,
                        model_data = mf_data,
                        model_formula = model_formula,
                        summary = model_summary,
                        residuals_histogram = hist,
residuals_qq_plot = qqplot,
                        point_estimates_matrix = model_matrix,
                        point_estimates = model_estimates)
  if (length(fixed_effects) > 1) {
    model_results$anova <- model_anova
  }
  if(!is.null(contrasts)){
    model_results$pairwise_comparisons_matrix <- result_matrix
    model_results$pairwise_comparisons <- pairwise_comparisons
  }
  
  return(model_results)
  
}
  