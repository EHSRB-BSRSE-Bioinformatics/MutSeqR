#' Perform General Linear Modelling on Mutation Frequency
#' 
#' @param mf_data The data frame containing the mutation data summary. 
#' This can be obtained using the calculate_mut_freq(summary = TRUE). 
#' @param factor The factor used for the model. This should be the name of 
#' the column that will act as the factor for the general linear model.
#' @param reference_level Refers to one of the levels within your factor. 
#' Remaining levels will be contrasted against the reference_level.
#' Ex. if my factor was a dose column, then my reference_level would refer to my
#' negative control dose. 
#' @param muts The column containing the mutation count per individual.
#' @param total_count The column containing the depth per individual.
#' 
#' @importFrom dplyr  
#' @importFrom magrittr %>%
#' @importFrom doby esticon

glm_mf_by_factor <- function(mf_data,
                             factor = "dose",
                             reference_level = 0,
                             muts = "sample_sum_unique",
                             total_count = "sample_group_depth") {
  
  mf_data[[factor]] <- as.factor(mf_data[[factor]])
  # Extract unique levels of the factor
  factor_levels <- unique(mf_data[[factor]])
  
  # Reorder levels with reference_level first
  factor_levels <- c(reference_level, setdiff(factor_levels, reference_level))
  mf_data[[factor]] <- factor(mf_data[[factor]], levels = factor_levels)
  
  # Construct glm formula 
  glm_formula <- as.formula(sprintf("cbind(%s, %s) ~ %s", muts, total_count, factor))
  
  
  model <- glm(glm_formula,
           family = "quasibinomial", data = mf_data,
           weights = get(total_count))
  
  mf_data$glm_residuals <- model$residuals
  
  # Create a subset based on residuals condition
  subset_df <- subset(mf_data, abs(glm_residuals) > 4)
    
  # Check if the subset is not empty
  if (nrow(subset_df) > 0) {
    # Print message and subset
    cat("It is advised to remove observations with residuals greater than 4 in absolute value.",
        "The following observations have abs(glm_residuals) > 4:\n")
    print(subset_df)
  } else {
    max_residual_index <- which.max(abs(mf_data$glm_residuals))
    
    # Extract the row with the maximum residual
    max_residual_row <- mf_data[max_residual_index, ]
    
    # Print the entire row with the maximum residual
    cat("The row with the maximum residual is:\n")
    print(max_residual_row)
    
  }
  
 # Plot the residuals
  hist_res <- hist(mf_data$glm_residuals, main = "Residuals", col = "yellow")
  hist_res

  qqnorm_res <- qqnorm(mf_data$glm_residuals, main = "QQ Plot of Residuals")
  qqline_res <- qqline(mf_data$glm_residuals, col = "red")
  qqnorm_res
  qqline_res    
  
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
  rownames(model_estimates) <- levels(mf_data$dose)
  
  message(paste("MF estimates per", factor, ":\n"))
  return(model_estimates)
}
  