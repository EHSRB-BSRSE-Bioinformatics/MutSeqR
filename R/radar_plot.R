#' Create a radar plot
#' @description Create a radar plot
#' @param mf_data A data frame with the data to plot
#' @param response_col The column with the response values
#' @param label_col The column with the labels for the radar plot.
#' @param facet_col The column with the group to facet the radar plots.
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr select filter
#' 
#' @return A radar plot
#' @export


# radar chart
radar_plot <- function(mf_data,
                      response_col,
                      label_col,
                      facet_col) {

  if (!requireNamespace("fmsb")) {
    stop("You need the package fmsb to run this function.")
  }
  plot_data <- mf_data %>%
    dplyr::select({{response_col}}, {{label_col}}, {{facet_col}}) %>%
    tidyr::pivot_wider(names_from = {{label_col}}, values_from = {{response_col}})

    # Convert dose column to a factor to handle arbitrary values
    plot_data[[facet_col]] <- as.factor(plot_data[[facet_col]])
    # Get levels
    facet_levels <- levels(plot_data[[facet_col]])
  
  # Set up the layout for the plots
  n_plots <- length(facet_levels)
  n_cols <- 2  # Number of columns in the grid
  n_rows <- ceiling(n_plots / n_cols)  # Number of rows needed
  graphics::layout(matrix(1:n_plots, nrow = n_rows, ncol = n_cols, byrow = TRUE))

  for (i in seq_along(facet_levels)) {
    facet <- facet_levels[i]
    df_i <- plot_data %>%
      dplyr::filter(.data[[paste0(facet_col)]] == facet) %>%
      select(-{{facet_col}})

    count <- ncol(df_i)
    # Add rows for max and min values
    max <- max(df_i, na.rm = TRUE) * 1.1
    max_row <- rep(max, count)
    min_row <- rep(0, count)
    df <- rbind(max_row, min_row, df_i)
    axis_labels <- seq(from = 0, to = max, length.out = 5)
    axis_labels <- sprintf("%.1e", axis_labels)
    title <- paste(facet_col, facet)
    plot <- fmsb::radarchart(df = df,
                axistype = 1,
                caxislabels = axis_labels,
                sep = 5, # number of axis ticks
                caxiscol = "grey",
                vlabels = NULL, # variable labels
                axislabcol = "grey",
                vlcex = 0.7,  # variable label font size
                title = title,
                pcol = "black", # color of the polygon
                pfcol = NULL, # color of the polygon fill
                plwd = 2, # width of the polygon line
                plty = 1, # type of the polygon line
                cglcol = "grey",
                cglty = 1,
                cglwd = 0.7)
    }
  graphics::layout(1)
  grDevices::dev.off()
}
