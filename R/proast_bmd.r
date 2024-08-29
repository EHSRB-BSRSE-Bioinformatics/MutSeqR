#' proast_bmd
#' @description This function is an extension of the PROAST package that
#' calculates the BMD of continuous, individual-level data. The function
#' 
#' 
proast_bmd <- function(mf_data,
                       dose_col = "dose",
                       response_col = "sample_MF_min",
                       CES = 0.5,
                       adjust_CES_to_group_SD = FALSE,
                       model_options = "all",
                       model_selection = "m3 or m5",
                       model_averaging = TRUE,
                       num_bootstraps = 200,
                       summary = TRUE) {

  CES_sd <- as.numeric(adjust_CES_to_group_SD) + 1

  # Model Options
  if (model_options == "expon") {
    model.options <- "Exponential model only"
  } else if (model_options == "expon_hill") {
    model.options <- "Exponential and Hill model"
  } else if (model_options == "expon_hill_invexpon") {
    model.options <- "previous option with inverse exponential model added"
  } else if (model_options == "all") {
    model.options <- "previous option with lognormal DR model added"
  } else {
    stop("Invalid model_options argument. Please choose from 'expon', 'expon_hill', 'expon_hill_invexpon', or 'all'.")
  }

  # Model Selection
  if (model_selection == "m3") {
    model.selection <- "select model 3 from various nested families of models"
  } else if (model_selection == "m5") {
    model.selection <- "select model 5 from various nested families of models"
  } else if (model_selection == "m3_m5") {
    model.selection <-  "select model 3 or 5 from various families of models"
  } else {
     stop("Invalid model_selection argument. Please choose from 'm3', 'm5', or 'm3_m5'.")
  }

  results <- f.proast(mf_data, interactive_mode = FALSE,
                     datatype = "continuous, individual data",
                     model_choice = model.selection,
                     setting_choice = "set of models",
                     indep_var_choice = dose_col,
                     Vyans_input = response_col,
                     covariates = 0,
                     CES = CES,
                     adjust_CES_to_group_SD = CES_sd,
                     model_selection = model.options,
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