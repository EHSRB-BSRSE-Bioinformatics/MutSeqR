#' Perform General Linear Modelling on DupSeq Mutation Frequency
#' 
#' @param mf_data The data frame containing the mutation frequency data.
#' Mutation counts and total_depth should be summarized per sample alongside a
#' column for your factor of interest. 
#' This data can be obtained using the calculate_mut_freq(summary = TRUE). 
#' @param factor This should be the name of the column that will act as the 
#' factor for modelling mutation frequency. Only one factor can be analysed within 
#' the general linear model. Ex. "dose", or "time_point".
#' @param reference_level Refers to one of the levels within your factor. 
#' Remaining levels will be contrasted against the reference_level.
#' Ex. if my factor was a dose column, then my reference_level would refer to my
#' negative control dose. 
#' @param muts The column containing the mutation count per sample.
#' @param total_count The column containing the sequencing depth per sample.
#' @param contrast_table a .txt file that will provide the information necessary 
#' to make comparisons between samples. The first column must be a level within
#' your factor and the second column must be the level within your factor that
#' it will be compared to. The names must correspond to entries in your mf_data
#' factor column. For example, if your factor is "dose" with dose groups 
#' 0, 25, 50, 100, then the first column would contain the treated groups 
#' (25, 50, 100), while the second column would be 0, thus comparing each treated
#' group to the control group. The contrast table for this example would be: 
#' 
#' 25 0 
#' 
#' 50 0 
#' 
#' 100 0 
#' @returns A data frame with the model estimates and pairwise comparisons
#' @importFrom magrittr %>%
#' @importFrom doBy esticon
#' @importFrom dplyr bind_rows

glm_mf_by_factor <- function(mf_data,
                             factor = "dose",
                             reference_level = 0,
                             muts = "sample_sum_unique",
                             total_count = "sample_group_depth",
                             contrast_table =  rbind(c(6.25, 0),
                                                    c(12.5, 0),
                                                    c(25, 0),
                                                    c(12.5, 6.25),
                                                    c(25, 6.25),
                                                    c(25, 12.5))) {
  
  mf_data[[factor]] <- as.factor(mf_data[[factor]])
  # Extract unique levels of the factor
  factor_levels <- unique(mf_data[[factor]])
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
  print(hist_res)

  qqnorm_res <- qqnorm(mf_data$glm_residuals, main = "QQ Plot of Residuals")
  qqline_res <- qqline(mf_data$glm_residuals, col = "red")
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
    
}
  