#' Cleveland Plot
#' @description Make a Cleveland plot for the PROAST results. Matches ToxicR.
#' @param results PROAST results object.
#' @param covariate_col Covariate column name.
#' @param output_path Output path for the plot. If the output_path doesn't
#' exist, it will be created. If NULL, the plots will not be exported.
#' @return A list of ggplot objects for each response in results.
#' @importFrom dplyr arrange filter mutate pull rename
#' @import ggplot2
#'
cleveland_plot <- function(results,
                           covariate_col = NULL,
                           output_path = NULL) {
  # Cleveland plot: all models w weights
  results_df <- results[[2]] %>%
    dplyr::mutate(CED = as.numeric(.data$CED),
                  CEDL = as.numeric(.data$CEDL),
                  CEDU = as.numeric(.data$CEDU))
  plots <- list()
  for (i in unique(results_df$Response)) {
    c.plot.df <- results_df %>%
      dplyr::filter(.data$Response == i)

    if (!is.null(covariate_col)) {
      c.plot.df$weights[c.plot.df$Selected.Model == "Model averaging"] <- 1
      c.plot.df <- c.plot.df %>%
        dplyr::rename(Model = "Selected.Model")
      c.plot.df$Selected.Model <- paste(c.plot.df$Model, c.plot.df$Covariate)
      model_order <- c.plot.df %>%
        dplyr::arrange(.data$Covariates, .data$weights) %>%
        dplyr::pull(.data$Selected.Model)
      c.plot.df$Selected.Model <- factor(c.plot.df$Selected.Model,
                                         levels = model_order)
    } else {
      c.plot.df$weights[c.plot.df$Selected.Model == "Model averaging"] <- 1
      model_order <- c.plot.df %>%
        dplyr::arrange(.data$weights) %>%
        dplyr::pull(.data$Selected.Model)

      c.plot.df <- c.plot.df %>%
        dplyr::mutate(Selected.Model = factor(.data$Selected.Model,
                                              levels = model_order))
      c.plot.df$Model <- c.plot.df$Selected.Model
    }

    c <- ggplot(c.plot.df, aes(x = CED, y = Selected.Model)) +
      geom_errorbar(aes(xmin = CEDL, xmax = CEDU),
                    color = "gray",
                    width = 0.1) +
      geom_point(aes(size = weights),
                 color = "red") +
      ggplot2::scale_size_continuous(guide = "none", range = c(0.5, 3)) +
      ggplot2::theme(panel.background = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(),
                     panel.grid = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_line(),
                     axis.ticks.y = ggplot2::element_line()) +
      ggplot2::xlab(paste("BMD Estimate for", i)) +
      ggplot2::ylab("Model") +
      ggtitle("BMD by Selected Model (Sorted by Weights)")
    plot_name <- paste0(i, "_cleveland")
    plots[[plot_name]] <- c
    if (!is.null(output_path)) {
      file_name <- file.path(paste0("PROAST_", i, "_cleveland.svg"))
      ggsave(
        filename = file_name,
        plot = c,
        device = "svg",
        path = output_path,
        create.dir = TRUE
      )
    }
  }
  return(plots)
}