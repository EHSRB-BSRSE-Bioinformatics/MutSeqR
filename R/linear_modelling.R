#' Perform General Linear Modelling on DupSeq Mutation Frequency
#'
#' Some notes: 
#' GLM fixed effect ex. dose, tissue
#' GLMM: fixed effect(s) + random effects ex. sample
#' If you have repeated measures for the same sample, use glmm w sample as a random effect. 
#' Ex. MF ~ dose*target + (1|sample) (multiple targets per sample = random effect)
#' Ex. MF ~ dose*tissue = glm since each sample is only in one of each group
#' 
#' 
#' @param mf_data The data frame containing the mutation frequency data.
#' Mutation counts and total_depth should be summarized per sample alongside a
#' column for your factor of interest. 
#' This data can be obtained using the calculate_mut_freq(summary = TRUE). 
#' @param factor This should be the name of the column that will act as the 
#' factor/fixed effect for modelling mutation frequency. Only one factor can be analysed within 
#' the general linear model. Ex. "dose", or "time_point".
#' @param reference_level Refers to one of the levels within your factor. 
#' Remaining levels will be contrasted against the reference_level.
#' Ex. if my factor was a dose column, then my reference_level would refer to my
#' negative control dose. 
#' @param muts The column containing the mutation count per sample.
#' @param total_count The column containing the sequencing depth per sample.
#' @param contrast_table_file a filepath to a tab-delimited .txt file that will 
#' provide the information necessary to make comparisons between samples. 
#' The first column must be a level within your factor and the second column 
#' must be the level within your factor that it will be compared to. The names 
#' must correspond to entries in your mf_data factor column. Put the factor that
#' you expect to have the higher mutation frequency in the 1st column and the 
#' factor level that you expect to have a lower mutation frequency in the second 
#' column. For example, if your factor is "dose" with dose groups 0, 25, 50, 100, 
#' then the first column would contain the treated groups (25, 50, 100), while 
#' the second column would be 0, thus comparing each treated group to the 
#' control group. The contrast table for this example would be: 
#' 
#' 25 0 
#' 
#' 50 0 
#' 
#' 100 0 
#' 
#' Alternatively, if you would like to compare mutation frequency between treated
#' dose groups, than the contrast table would look as follows, with the lower
#' dose always in the second column, as we expect it to have a lower mutation 
#' frequency.
#' 
#' 100 25
#' 
#' 100 50
#' 
#' 50 25
#' 
#' 
#' 
#' 
#' @returns A data frame with the model estimates and specified pairwise comparisons
#' @importFrom magrittr %>%
#' @importFrom doBy esticon
#' @importFrom dplyr bind_rows

glm_mf_by_factor <- function(mf_data,
                             factor = "dose",
                             reference_level = 0,
                             muts = "sample_sum_unique",
                             total_count = "sample_group_depth",
                             contrast_table_file =  NULL) {
  
  mf_data[[factor]] <- as.factor(mf_data[[factor]])
  # Extract unique levels of the factor
  factor_levels <- unique(mf_data[[factor]])
  # Check that the reference level is a level of the factor
  reference_valid <- reference_level %in% factor_levels
  if (!reference_valid) {
    stop(paste(reference_level, "is not a valid reference_level. Please ensure that your 
                  reference_level corresponds to a level of your factor column"))
  }
  
  # Set the reference_level; reorders the factor levels with reference_level first
  mf_data[[factor]] <- relevel(mf_data[[factor]], ref = as.character(reference_level))
  
  # Construct glm formula 
  glm_formula <- as.formula(sprintf("cbind(%s, %s) ~ %s", muts, total_count, factor))
  
  # quasibinomial in order to account for over dispersion
  model <- glm(glm_formula,
           family = "quasibinomial", data = mf_data,
           weights = get(total_count))
  
  mf_data$glm_residuals <- model$residuals
  
  # Create a subset based on residuals condition
  subset_df <- subset(mf_data, abs(glm_residuals) > 4)
    
  # Check if the subset is not empty
  if (nrow(subset_df) > 0) {
    # Print message and subset
    message("It is advised to remove observations with residuals greater than 4 in absolute value.",
        "The following observations have abs(glm_residuals) > 4:\n")
    print(subset_df)
  } else {
    max_residual_index <- which.max(abs(mf_data$glm_residuals))
    
    # Extract the row with the maximum residual
    max_residual_row <- mf_data[max_residual_index, ]
    
    # Print the entire row with the maximum residual
    message("The row with the maximum residual is:\n")
    print(max_residual_row)
    
  }
  
 # Plot the residuals
  hist_res <- hist(mf_data$glm_residuals, main = "Residuals", col = "yellow")
  message(" \n Printing histogram of model residuals: \n")
  print(hist_res)

  qqnorm_res <- qqnorm(mf_data$glm_residuals, main = "QQ Plot of Residuals")
  qqline_res <- qqline(mf_data$glm_residuals, col = "red")
  message(" \n Printing qq plot for model residuals  \n")
  print(qqnorm_res)
  print(qqline_res)    
  
  ##################################################################
  # Contrast matrix for the point estimates
  
  # Create Contrast matrix
    # Initialize an empty contrast matrix based on number of factor levels
  contrast_matrix <- matrix(0, nrow = length(factor_levels), ncol = length(factor_levels))
  # Set the contrast matrix based on the factor levels
  for (i in seq_along(factor_levels)) {
    contrast_matrix[i, i] <- 1
  }
  # Set the first column to all 1s
  contrast_matrix[, 1] <- 1
  
  #esticon Contrasts for glm. Computes weighted sums of the estimated regression parameters
  model_estimates <- doBy::esticon(obj = model,L = contrast_matrix)
  
  ###Wald Statistics
  model_estimates <- as.data.frame(model_estimates)
  
  #take the exponential function of the rows?
  model_estimates$estimate <- exp(model_estimates$estimate)
  delta <- model_estimates$estimate^2
  model_estimates$lwr <- exp(model_estimates$lwr)
  model_estimates$upr <- exp(model_estimates$upr)
  model_estimates$std.error <- sqrt(delta*model_estimates$std.error^2)
  #remove extra columns (statistic, p value, beta0, df)
  model_estimates <- model_estimates[,-c(3,4,5,6)]
  # Rename table
  colnames(model_estimates) <- c("Estimate", "Std.Err", "Lower", "Upper")
  rownames(model_estimates) <- levels(mf_data[[factor]])
  
##################################################################
  # Contrast matrix for the pairwise comparisons
  
  # load  and do checks
  if (!is.null(contrast_table_file)) {
    contrast_table <- read.delim(file.path(contrast_table_file), sep = "\t",
                             header = F)
  } else {
    cat("Please supply the file path to your contrast table")
  }
  if (ncol(contrast_table) <= 1) {
    stop("Your contrast_table only has one column. Make sure your file is tab-delimited.")
  }
  contrast_table <- as.matrix(contrast_table)
  all_valid <- all(contrast_table %in% factor_levels)
  if (!all_valid) {
    invalid_values <- contrast_table[!(contrast_table %in% factor_levels)]
    stop(paste("Invalid values in contrast_table:", paste(invalid_values, collapse = ","), 
                  ". Please ensure that values in the contrast_table correspond to values in your factor column."))
  }
  
  # Identify the index of the reference level in factor_levels
  ref_level_index <- match(reference_level, factor_levels)
  
  # Initialize an empty contrast matrix
  contrast_matrix <- matrix(0, nrow = nrow(contrast_table), ncol = length(factor_levels))
  
  # Identify the indices of the levels to be compared
  comp_indices <- matrix(0, nrow = nrow(contrast_table), ncol = 2)
  for (i in 1:nrow(contrast_table)) {
    comp_indices[i, ] <- match(contrast_table[i, ], factor_levels)
  }
  
  # Set the appropriate values in the contrast matrix
  for (i in 1:nrow(contrast_table)) {
    # Check if the first column index is not the reference level index, if yes, put 1, otherwise, put 0
    contrast_matrix[i, comp_indices[i, 1]] <- ifelse(comp_indices[i, 1] != ref_level_index, 1, 0)
    
    # Check if the second column index is not the reference level index, if yes, put -1, otherwise, put 0
    contrast_matrix[i, comp_indices[i, 2]] <- ifelse(comp_indices[i, 2] != ref_level_index, -1, 0)
  }
  
  #esticon Contrasts for glm.
  pairwise_comparisons <- doBy::esticon(obj = model, L = contrast_matrix )
  pairwise_comparisons <- as.data.frame(pairwise_comparisons)
  
  pairwise_comparisons$estimate <- exp(pairwise_comparisons$estimate)
  delta <- pairwise_comparisons$estimate^2
  pairwise_comparisons$lwr <- exp(pairwise_comparisons$lwr)
  pairwise_comparisons$upr <- exp(pairwise_comparisons$upr)
  pairwise_comparisons$std.error <- sqrt(delta*pairwise_comparisons$std.error^2)
  # Remove necessary column
  pairwise_comparisons <- pairwise_comparisons[,-5]
  # Rename the estimates table
  colnames(pairwise_comparisons) <- c("Estimate", "Std.Err", "Obs.T", "p.value", "df", "Lower", "Upper")
  new_row_names <- paste(contrast_table[, 1], "vs.", contrast_table[ , 2])
  # Rename the rows in the estimates table
  rownames(pairwise_comparisons) <- new_row_names
  
  ###############################################
  # Holm-Sidak Correction for multiple comparison
  my.holm.sidak <- function(P){
    m <- length(P)
    if(m > 1){
      Psort <- matrix(c(P, 1:m), 2, m, byrow = TRUE)
      Psort <- Psort[, order(Psort[1, ])]
      for (i in 1:m) {
        adjust <- m + 1 - i
        Psort[1, i] <- pmin(1, (1 - ((1 - Psort[1, i])^adjust)))
      }
      Psort <- Psort[, order(Psort[2, ])]
      P.adjust <- Psort[1, ]
      P.adjust
    } else{
      #No adjustment as there is only one comparison
      P.adjust <- P
      P.adjust
    }}
    
    my.holm.sidak(pairwise_comparisons$p.value)

 ############################################### 
  # Return output file
  glm_output <- dplyr::bind_rows(model_estimates, pairwise_comparisons)
  return(glm_output)
    
## To DO:  Results as list
    
}

########################################################################
########################################################################

#' Perform General Linear Mixed Modelling on DupSeq Mutation Frequency based on
#' fixed and random effects
#' 
#' @param mf_data The data frame containing the mutation frequency data.
#' Mutation counts and total_depth should be summarized per sample alongside a
#' columns for your fixed effects. 
#' This data can be obtained using the calculate_mut_freq(summary = TRUE). 
#' @param fixed_effects This should be the name of the columns that will act as the 
#' factor/fixed effect for modelling mutation frequency. Fixed effects encapsulate 
#' the tendencies/trends that are consistent at the levels of primary interest. 
#' These effects are considered fixed because they are non-random and assumed to 
#' be constant for the samples being studied. Ex. 'fixed_effects = c("dose", "target")'.
#' @param test_interaction TRUE or FALSE. Whether or not your model should test 
#' the interaction between fixed_effects on top of the effect of each factor.
#' @param random_effects This should be the name of the column to be analysed as 
#' a random effect in the model. Random effects introduce statistical variability 
#' at different levels of the data hierarchy. These account for the unmeasured 
#' sources of variance that affect certain groups in the data. Ex. If your model
#' uses repeated measures of sample, 'random_effects = "sample"'. 
#' @param reference_level Refers to one of the levels within your fixed_effects. 
#' Remaining levels will be contrasted against the reference_level.
#' Ex. if my factor was a dose column, then my reference_level would refer to my
#' negative control dose. 
#' @param muts The column containing the mutation count per sample.
#' @param total_count The column containing the sequencing depth per sample.
#' @param contrast_table_file a filepath to a tab-delimited .txt file that will 
#' provide the information necessary to make comparisons between samples. 
#' The first column must be a level within your factor and the second column 
#' must be the level within your factor that it will be compared to. The names 
#' must correspond to entries in your mf_data factor column. Put the factor that
#' you expect to have the higher mutation frequency in the 1st column and the 
#' factor level that you expect to have a lower mutation frequency in the second 
#' column. For example, if your factor is "dose" with dose groups 0, 25, 50, 100, 
#' then the first column would contain the treated groups (25, 50, 100), while 
#' the second column would be 0, thus comparing each treated group to the 
#' control group. The contrast table for this example would be: 
#' 
#' 25 0 
#' 
#' 50 0 
#' 
#' 100 0 
#' 
#' Alternatively, if you would like to compare mutation frequency between treated
#' dose groups, than the contrast table would look as follows, with the lower
#' dose always in the second column, as we expect it to have a lower mutation 
#' frequency.
#' 
#' 100 25
#' 
#' 100 50
#' 
#' 50 25
#' 
#' 
#' 
#' 
#' @returns A data frame with the model estimates and specified pairwise comparisons
#' @importFrom magrittr %>%
#' @importFrom doBy esticon
#' @importFrom dplyr bind_rows
#' @import lme4 glmer fixef
#' @import HLMdiag
#' @import Cairo
#' @import robustlmm
#' @import DHARMA
#' @import car Anova
#' @import multcomp


glmm_mf <- function(mf_data,
                    fixed_effects = c("dose", "label"),
                    test_interaction = TRUE,
                    random_effects = "sample",
                    reference_level = 0,
                    muts = "sample_label_sum_unique",
                    total_count = "sample_label_group_depth",
                    contrast_table_file =  NULL) {
  
  # Convert specified columns to factors
  mf_data[, fixed_effects] <- lapply(mf_data[, fixed_effects], as.factor) 
  
  # Extract unique levels from each factor and combine them
  factor_levels <- unique(unlist(lapply(fixed_effects, function(factor_name) {
    factor_variable <- mf_data[[factor_name]]
    levels(factor_variable)
  })))
  
  # Check that the reference level is a level of the factor
  reference_valid <- reference_level %in% factor_levels
  if (!reference_valid) {
    stop(paste(reference_level, "is not a valid reference_level. Please ensure that your 
                  reference_level corresponds to a level of your factor column"))
  }
  
  # factor_level <- lapply(fixed_effects, function(factor_name) {
  #  if(reference_level %in% levels(mf_data[[factor_name]])){
  #    mf_data[[factor_name]] <- relevel(mf_data[[factor_name]], ref = as.character(reference_level))
  #  }
  #  levels(mf_data[[factor_name]])
  # })
  
  # Set the reference_level; reorders the factor levels with reference_level first
  # mf_data[[factor]] <- relevel(mf_data[[factor]], ref = as.character(reference_level))
  
  # Construct the formula
  if(test_interaction){
    formula_str <- paste("cbind(", muts, ",", total_count, ") ~ ", paste(fixed_effects, collapse = "*"))
  } else {
    formula_str <- paste("cbind(", muts, ",", total_count, ") ~ ", paste(fixed_effects, collapse = "+"))
  }
  
  # Add random effects to the formula
  if (length(random_effects) > 0) {
    random_formula <- paste(paste("(1|", random_effects, ")", collapse = "+"))
    formula_str <- paste(formula_str, "+", random_formula)
  }
  
  # Convert the string formula to an actual formula object
  glmm_formula <- as.formula(formula_str)
  
  # We should have a think about the control and how it can be made flexible. 
  model <- lme4::glmer(glmm_formula, 
                 family = "binomial", 
                 data = mf_data, 
                 control = lme4::glmerControl(check.conv.grad = lme4::.makeCC("warning", 
                                                                  tol = 3e-3, 
                                                                  relTol = NULL)))
  
  model_summary <- summary(model)
  
  model_anova <-car::Anova(model)
 #################################
  # Check residuals
  mf_data$Resid <- residuals(model)
  # Print the entire row with the maximum residual
  cat("The row with the maximum residual is:\n")
  mf_data[abs(residuals(model)) == max(abs(residuals(model))),]
  
  ## TO DO: Let's walk users through what the residuals should look like and give advice on them
  
  # Create a data frame with the variables needed for the boxplot
  plot_data <- mf_data %>%
    dplyr::mutate(mutation_frequency = mf_data[[muts]] / mf_data[[total_count]])
  
  # Create a boxplot to visualise data
  # Constructing formula for boxplot
  plot_formula <- paste("mutation_frequency ~", paste(fixed_effects, collapse = "*"))
  plot_formula <- as.formula(plot_formula)
  par(las = 2)
  boxplot(plot_formula,
          data = plot_data,
          main = "Mutation Frequency by Fixed Effects",
          xlab = paste(fixed_effects[1],"x", fixed_effects[2]),
          ylab = "Mutation Frequency")
  
  
  par(las = 1)
  hist <- hist(residuals(model), main = "Residuals")
  cat(" \n Printing histogram of model residuals: \n")
  print(hist)
  
  qqnorm(residuals(model))
  abline(h = -3, col = "red", lwd = 2)
  abline(h = 3, col = "red", lwd = 2)
  
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
  if(test_interaction){
    fixed_effect_formula <- as.formula(paste(" ~ ", paste(fixed_effects, collapse = "*")))
  } else {
    fixed_effect_formula <- as.formula(paste(" ~ ", paste(fixed_effects, collapse = "+")))
  }
  
  # Create the model matrix using model.matrix  
  model_matrix <- model.matrix(fixed_effect_formula, data = design_matrix)
  row_names <- apply(design_matrix, 1, paste, collapse = ":")
  rownames(model_matrix) <- row_names 
  
  # Computed estimates
  model_estimates <-esticon(model,lambda)

  model_estimates <- as.data.frame(model_estimates)
  model_estimates$estimate <- exp(model_estimates$estimate)
  delta <- model_estimates$estimate^2
  model_estimates$lwr <- exp(model_estimates$lwr)
  model_estimates$upr <- exp(model_estimates$upr)
  model_estimates$std.error <- sqrt(delta*model_estimates$std.error^2)
  # Clean data
  model_estimates <- model_estimates[,-c(3,4,5,6)]
  colnames(model_estimates) <- c("Estimate", "Std.Err", "Lower", "Upper")

}
  