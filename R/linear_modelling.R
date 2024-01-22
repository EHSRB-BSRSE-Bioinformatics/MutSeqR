#' Perform General Linear Modelling on Mutation Frequency
#' 
#' @param mf_data The data frame containing the mutation data summary. 
#' This can be obtained using the calculate_mut_freq(summary = TRUE). 
#' @param factor The factor used for the model. This should be the name of 
#' the column that will act as the factor for the general linear model. 
#' @param muts The column containing the mutation count per individual.
#' @param total_count The column containing the depth per individual.
#' 
#' @importFrom dplyr  
#' @importFrom magrittr %>%
#' @importFrom doby

glm_mf_by_factor <- function(mf_data,
                             factor = "dose",
                             muts = "sample_sum_unique",
                             total_count = "sample_group_depth") {
  
  mf_data[[factor]] <- as.factor(mf_data[[factor]])
  
  glm_formula <- as.formula(paste0(muts,"/(",muts,"+",total_count,") ~ ", factor))
  
  model <- glm(glm_formula,
           family = "quasibinomial", data = mf_data,
           weights = get(muts) + get(total_count))
  
  
}