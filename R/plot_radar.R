#' Create a radar plot
#' @description Create a radar plot
#' @param mf_data A data frame with the data to plot
#' @param response_col The column with the response values
#' @param label_col The column with the labels for the radar plot.
#' @param facet_col The column with the group to facet the radar plots.
#' @param indiv_y A logical indicating whether to use individual y-axis scales
#' for each plot.
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr select filter pull
#' @examples
#' # Plot the mean MFmin of each genomic target per dose group
#' # Order the genomic targets by their genic context.
#' #Load the example data and calculate MF
#' example_file <- system.file("extdata", "Example_files",
#'                             "example_mutation_data_filtered.rds",
#'                             package = "MutSeqR")
#' example_data <- readRDS(example_file)
#' mf <- calculate_mf(mutation_data = example_data,
#'                    cols_to_group = c("sample", "label"),
#'                    retain_metadata_cols = c("dose_group", "genic_context"))
#' # Define the order of the genomic targets
#' label_order <- mf %>% dplyr::arrange(genic_context) %>%
#'   dplyr::pull(label) %>%
#'   unique()
#' # Calculate the mean MF per dose_group for each target.
#' mean <- mf %>%
#'   dplyr::group_by(dose_group, label) %>%
#'   dplyr::summarise(mean = mean(mf_min))
#' # Set the order of each column
#' mean$dose_group <- factor(mean$dose_group,
#'                           levels = c("Control",
#'                                      "Low",
#'                                      "Medium",
#'                                      "High"))
#' mean$label <- factor(mean$label,
#'                      levels = label_order)
#' # Plot
#' plot <- plot_radar(mf_data = mean,
#'                    response_col = "mean",
#'                    label_col = "label",
#'                    facet_col = "dose_group",
#'                    indiv_y = FALSE)
#' @return A radar plot
#' @export


# radar chart
plot_radar <- function(mf_data,
                       response_col,
                       label_col,
                       facet_col,
                       indiv_y = TRUE) {

  if (!requireNamespace("fmsb", quietly = TRUE)) {
    stop("You need the package fmsb to run this function.")
  }
  plot_data <- mf_data %>%
    dplyr::select({{response_col}}, {{label_col}}, {{facet_col}}) %>%
    tidyr::pivot_wider(names_from = {{label_col}}, values_from = {{response_col}})

  if (is.factor(mf_data[[rlang::as_name(enquo(label_col))]])) {
    label_levels <- levels(mf_data[[rlang::as_name(enquo(label_col))]])
    plot_data <- plot_data %>%
      dplyr::select({{facet_col}}, all_of(label_levels))
  }
    # Convert dose column to a factor to handle arbitrary values
    plot_data[[facet_col]] <- as.factor(plot_data[[facet_col]])
    # Get levels
    facet_levels <- levels(plot_data[[facet_col]])

    global_max <- if (!indiv_y) {
      max(plot_data %>% dplyr::ungroup() %>% dplyr::select(-{{facet_col}}), na.rm = TRUE) * 1.1
    } else {
      NULL
    }
  
  # Set up the layout for the plots
  n_plots <- length(facet_levels)
  n_cols <- 2  # Number of columns in the grid
  n_rows <- ceiling(n_plots / n_cols)  # Number of rows needed
  graphics::layout(matrix(1:n_plots, nrow = n_rows, ncol = n_cols, byrow = TRUE))

  for (i in seq_along(facet_levels)) {
    facet <- facet_levels[i]
    df_i <- plot_data %>%
      dplyr::filter(.data[[paste0(facet_col)]] == facet) %>%
      dplyr::ungroup() %>%
      dplyr::select(-{{facet_col}})

    count <- ncol(df_i)
    # Add rows for max and min values
    max_value <- if (is.null(global_max)) max(df_i, na.rm = TRUE) * 1.1 else global_max
    max_row <- rep(max_value, count)
    min_row <- rep(0, count)
    df <- rbind(max_row, min_row, df_i)
    axis_labels <- seq(from = 0, to = max_value, length.out = 5)
    axis_labels <- sprintf("%.1e", axis_labels)
    title <- paste(facet_col, facet)
    plot <- fmsb::radarchart(df = df,
                axistype = 1,
                caxislabels = axis_labels,
                sep = 5, # number of axis ticks
                caxiscol = "grey",
                vlabels = NULL, # variable labels
                axislabcol = "#797777",
                vlcex = 1.2,  # variable label font size
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
}
