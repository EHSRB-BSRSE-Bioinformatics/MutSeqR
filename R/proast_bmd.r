#' proast_bmd
#' @description This function is an extension of the PROAST package that
#' calculates the BMD of continuous, individual-level data. The function
#' 
#' 
proast_bmd <- function(mf_data,
                       dose_col = "dose",
                       response_col = "sample_MF_min",
                       covariate_col = NULL,
                       CES = 0.5,
                       adjust_CES_to_group_SD = FALSE,
                       model_averaging = TRUE,
                       num_bootstraps = 200,
                       summary = TRUE) {

  CES_sd <- as.numeric(adjust_CES_to_group_SD) + 1
  
  if(is.null(covariate_col)) {
    covariate <- 0
  } else {
    covariate <- covariate_col
  }

  results <- f.proast(mf_data,
                      interactive_mode = FALSE,
                      datatype = "continuous, individual data",
                      model_choice = "select model 3 or 5 from various families of models",
                      setting_choice = "set of models",
                      indep_var_choice = dose_col,
                      Vyans_input = response_col,
                      covariates = covariate,
                      custom_CES = CES,
                      adjust_CES_to_group_SD = CES_sd,
                      model_selection = "previous option with lognormal DR model added",
                      lower_dd = NULL,
                      upper_dd = NULL,
                      selected_model = "exponential",
                      model_averaging = model_averaging,
                      num_bootstraps = num_bootstraps)

  if (summary == TRUE) {
    results_df <- results[[2]]
    if (model_averaging) {
      ma_results <- results[[1]]$model_averaging$conf.int.ma

      f.plot.con(bmd[[1]]$"Expon. m3-")
    }
  } else {
    return(results[[1]])
  }


}