#' Cleveland Plot
#' @description Make a Cleveland plot for the PROAST results. Matches ToxicR.
#' @param results PROAST results object.
#' @param covariate_col Covariate column name.
#' @param output_path Output path for the plot.
#' @return ggplot object.
#' @importFrom dplyr arrange filter mutate pull rename
#' @import ggplot2
#' 
cleveland_plot <- function(results,
                           covariate_col = NULL,
                           output_path = NULL) {
  # Cleveland plot: all models w weights
  results_df <- results[[2]] %>%
    dplyr::mutate(CED = as.numeric(CED),
                  CEDL = as.numeric(CEDL),
                  CEDU = as.numeric(CEDU))

  for (i in unique(results_df$Response)) {
    c.plot.df <- results_df %>%
      dplyr::filter(.data$Response == i)

    if (!is.null(covariate_col)) {
      c.plot.df <- c.plot.df %>%
        dplyr::rename(Model = "Selected.Model")
      c.plot.df$Selected.Model <- paste(c.plot.df$Model, c.plot.df$Covariate)
      model_order <- c.plot.df %>%
        dplyr::arrange(Covariates, weights) %>%
        dplyr::pull(Selected.Model)
      c.plot.df$Selected.Model <- factor(c.plot.df$Selected.Model,
                                         levels = model_order)
    } else {
      model_order <- c.plot.df %>%
        dplyr::filter(.data$Selected.Model != "Model averaging") %>%
        dplyr::arrange(weights) %>%
        dplyr::pull(Selected.Model)

      c.plot.df <- c.plot.df %>%
        dplyr::mutate(Selected.Model = factor(Selected.Model,
                                              levels = c(unique(model_order),
                                                         "Model averaging")))
      c.plot.df$Model <- c.plot.df$Selected.Model
    }
    # assign dummy values to Model averaging,
    # making sure it is in range of the other CEDL and CEDU.
    c.plot.df$weights[c.plot.df$Model == "Model averaging"] <- NA
    c.plot.df$CED[c.plot.df$Model == "Model averaging"] <- with(subset(c.plot.df, Model == "Model averaging"), (CEDL + CEDU) / 2)

    c <- ggplot(c.plot.df, aes(x = CED, y = Selected.Model)) +
      geom_errorbar(aes(xmin = CEDL, xmax = CEDU),
                    color = "gray",
                    width = 0.1) +
      geom_point(aes(size = weights),
                 color = "red") +
      ggplot2::scale_size_continuous(guide = "none") +
      ggplot2::theme(panel.background = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(),
                     panel.grid = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_line(),
                     axis.ticks.y = ggplot2::element_line()) +
      ggplot2::xlab(paste("BMD Estimate for", i)) +
      ggplot2::ylab("Model") +
      ggtitle("BMD by Selected Model (Sorted by Weights)")
      
      file_name <- file.path(output_path, paste0("PROAST_", i , "_cleveland.svg"))
      ggsave(filename = file_name, plot = c, device = "svg")
  }
  return(c)
}