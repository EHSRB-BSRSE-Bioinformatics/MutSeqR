#' Plot the Mutation Frequency
#' @description This function creates a plot of the mutation frequency.
#' @param mf_data A data frame containing the mutation frequency data. This is
#' obtained from the calculate_mut_freq function with SUMMARY = TRUE.
#' @param group_col The name of the column containing the sample/group names
#' for the x-axis.
#' @param plot_type The type of plot to create. Options are "bar" or "point".
#' @param mf_type The type of mutation frequency to plot. Options are "min",
#' "max", "both", or "stacked". If "both", the min and max mutation
#' frequencies are plotted side by side. "stacked" can be chosen for bar
#' plot_type only. If "stacked", the difference between the min and max
#' MF is stacked on top of the min MF such that the total height of both
#' bars represent the max MF.
#' @param fill_col The name of the column containing the fill variable.
#' @param custom_palette A character vector of colour codes to use for the plot.
#' If NULL, a default palette is used
#' @param group_order The order of the samples/groups along the x-axis.
#' ' Options include:
#' \itemize{
#'   \item `none`: No ordering is performed. Default.
#'   \item `smart`: Samples are ordered based on the sample names.
#'   \item `arranged`: Samples are ordered based on one or more factor column(s)
#' in mf_data. Factor column names are passed to the function using the
#' `group_order_input`.
#'  \item `custom`: Samples are ordered based on a custom vector of sample
#' names. The custom vector is passed to the function using the
#' `group_order_input`.
#' }
#' @param group_order_input The order of the samples/groups if group_order is
#' "custom". The column name by which to arrange samples/groups if group_order
#' is "arranged"
#' @param add_labels The data labels to display on the plot. Either "count", "MF", or
#' "none". Count labels display the number of mutations, MF labels display
#' the mutation frequency.
#' @param scale_y_axis The scale of the y axis. Either "linear" or "log".
#' @param x_lab The label for the x axis.
#' @param y_lab The label for the y axis.
#' @param title The title of the plot.
#' @return A ggplot object
#' @import ggplot2
#' @importFrom dplyr arrange across all_of rename
#' @export
plot_mean_mf <- function(mf_data,
                         group_col,
                         mf_type = "both",
                         plot_type = c("bar", "line"),
                         plot_error_bars = TRUE,
                         plot_indiv_vals = TRUE,
                         fill_col = NULL, # may want to color bars
                         group_order = "none",
                         group_order_input = NULL,
                         add_labels = c("indiv_count", "indiv_MF", "mean_count", "mean_MF", "none"),
                         scale_y_axis = "linear",
                         x_lab = NULL,
                         y_lab = NULL,
                         title = NULL,
                         custom_palette = NULL) {

  if (group_order == "smart" && !requireNamespace("gtools", quietly = TRUE)) {
      stop("Package gtools is required when using the 'smart' group_order option. Please install the package using 'install.packages('gtools')'")
  }
  if (!plot_indiv_vals && add_labels %in% c("indiv_count", "indiv_MF")) {
    stop("plot_indiv_vals must be TRUE when add_labels is set to 'indiv_count' or 'indiv_MF'")
  }

   # axis_labels
  if (!is.null(x_lab)) {
    x_lab <- x_lab
  } else {
    x_lab <- group_col
  }
  if (!is.null(y_lab)) {
     y_lab <- y_lab
  } else {
    y_lab <- "Mutation Frequency (mutations/bp)"
  }

 # Title
  if (is.null(title)) {
    if (mf_type %in% c("stacked", "both")) {
      title <- paste0("min and max Mutation Frequency per ", group_col)
    } else if (mf_type %in% c("min", "max")) {
      title <- paste0(mf_type, " mutation frequency per ", group_col)
    }
  } else {
    title <- title
  }

    # x-axis order
  if (group_order == "none") {
    order <- as.vector(unique(mf_data[[group_col]]))
    mf_data[[group_col]] <- factor(mf_data[[group_col]])
  } else if (group_order == "smart") {
    order <- gtools::mixedsort(as.vector(unique(mf_data[[group_col]])))
    mf_data[[group_col]] <- factor(mf_data[[group_col]], levels = order)
  } else if (group_order == "arranged") {
    mf_data <- mf_data %>%
      dplyr::arrange(dplyr::across(dplyr::all_of({{group_order_input}})))
    order <- as.vector(unique(mf_data[[group_col]]))
    mf_data[[group_col]] <- factor(mf_data[[group_col]], levels = order)
  } else if (group_order == "custom") {
    mf_data[[group_col]] <- factor(mf_data[[group_col]],
                                   levels = group_order_input)
  }

  # Plot data
    # identify the MF columns
    MF_min_column_pattern <- paste0(".*(_MF_min", ")$")
    min_mf_col <- names(mf_data)[grepl(MF_min_column_pattern, names(mf_data))]
    MF_max_column_pattern <- paste0(".*(_MF_max", ")$")
    max_mf_col <- names(mf_data)[grepl(MF_max_column_pattern, names(mf_data))]

    # identify the sum columns for label making
    sum_min_column_pattern <- paste0(".*(_sum_min", ")$")
    min_count_col <- names(mf_data)[grepl(sum_min_column_pattern, names(mf_data))]
    sum_max_column_pattern <- paste0(".*(_sum_max", ")$")
    max_count_col <- names(mf_data)[grepl(sum_max_column_pattern, names(mf_data))]

  # Make indiv plot data
  plot_data <- mf_data %>%
    dplyr::rename(group_col = dplyr::all_of(group_col)) %>%
    dplyr::rename(mf_min = dplyr::all_of(min_mf_col)) %>%
    dplyr::rename(mf_max = dplyr::all_of(max_mf_col)) %>%
    dplyr::rename(sum_min = dplyr::all_of(min_count_col)) %>%
    dplyr::rename(sum_max = dplyr::all_of(max_count_col))

  # Make Group Mean data
  group_mean <- plot_data %>%
    dplyr::group_by(.data$group_col) %>%
    summarize(min_Mean = mean(.data$mf_min, na.rm = TRUE),
              min_SE = sd(.data$mf_min, na.rm = TRUE) / sqrt(n()),
              max_Mean = mean(.data$mf_max, na.rm = TRUE),
              max_SE = sd(.data$mf_max, na.rm = TRUE) / sqrt(n()))
  group_mean <- as.data.frame(group_mean)

  if (mf_type == "min") {
    plot_data <- rename(plot_data, mf_col = mf_min, sum_col = sum_min)
    group_mean <- rename(group_mean, Mean = min_Mean, SE = min_SE)
  }
  if (mf_type == "max") {
    plot_data <- rename(plot_data, mf_col = mf_max, sum_col = sum_max)
    group_mean <- rename(group_mean, Mean = max_Mean, SE = max_SE)
  }

  # Change MFmax for stacked
  if (mf_type == "stacked") {
    plot_data <- transform(plot_data,
                           mf_max = plot_data$mf_max - plot_data$mf_min)
    # set the max y for stacked
    max_y <- max(plot_data$mf_min + plot_data$mf_max) * 1.1
    # Change the MFmax for group mean data
    group_mean <- transform(group_mean,
                          max_Mean = group_mean$max_Mean - group_mean$min_Mean)
  }

  if (mf_type %in% c("both", "stacked")) {
    # mf_type = both or stacked, pivot the data to long format
    plot_data <- reshape(plot_data,
                         varying = list(c("sum_min", "sum_max"),
                                        c("mf_min", "mf_max")),
                         v.names = c("sum_col", "mf_col"),
                         times = c("min", "max"),
                         timevar = "mf_type",
                         direction = "long")
    # pivot the mean data to long format
    group_mean <- reshape(group_mean,
                          varying = list(c("min_Mean", "max_Mean"),
                                         c("min_SE", "max_SE")),
                          v.names = c("Mean", "SE"),
                          times = c("min", "max"),
                          timevar = "mf_type",
                          direction = "long")
  }
  # set the max_y
  if (mf_type != "stacked") {
    max_y <- max(plot_data$mf_col) * 1.1
  }
    # set the mf_type order for both/stacked
    if (mf_type == "both") {
      plot_data$mf_type <- factor(plot_data$mf_type, levels = c("min", "max"))
      group_mean$mf_type <- factor(group_mean$mf_type, levels = c("min", "max"))
    }
    if (mf_type == "stacked") {
      plot_data$mf_type <- factor(plot_data$mf_type, levels = c("max", "min"))
      group_mean$mf_type <- factor(group_mean$mf_type, levels = c("max", "min"))
    }
# end plot data

  # the Fill Column
  if (!is.null(fill_col)) { # if fill col exists
    if (fill_col == group_col) { # if fill col is the same as group col
      plot_data$fill_col <- plot_data$group_col
    } else { # if fill col is different from group col
      plot_data <- dplyr::rename(plot_data, fill_col = dplyr::all_of(fill_col))
    }

    if (mf_type %in% c("both", "stacked")) { # if mf_type is both or stacked
      fill <- interaction(plot_data$mf_type, plot_data$fill_col)
      fill_label <- paste("MF Type and", fill_col)
    } else { # if mf_type is min or max
      fill <- plot_data$fill_col
      fill_label <- paste(fill_col)
    }
  } else { # if fill col is NULL
    if (mf_type %in% c("both", "stacked")) {
      fill <- plot_data$mf_type
      fill_label <- "MF Type"
      plot_data$fill_col <- ""
    } else { # if mf_type is min or max
      fill <- plot_data$fill_col <- ""
      fill_label <- NULL
    }
  } # end fill column

  # scale y axis
  if (scale_y_axis == "log") {
    yscale <- ggplot2::scale_y_log10()
  } else {
    yscale <- ggplot2::scale_y_continuous(limits = c(0, max_y))
  }

  # palette
  if (is.null(custom_palette)) {
    if (mf_type %in% c("both", "stacked")) {
      n_colors <- length(unique(plot_data$fill_col)) * 2
    } else if (mf_type %in% c("min", "max")) {
      n_colors <- length(unique(plot_data$fill_col))
    }
    gradient <- colorRampPalette(colors = c("#c5e5fc",
                                            "#5ab2ee",
                                            "#12587b",
                                            "#263247",
                                            "#ffedef",
                                            "#ffb9c1",
                                            "#ff5264",
                                            "#b23946"))
    palette <- gradient(n_colors)
  } else {
    palette <- custom_palette
  }
  
  # Position
  if (mf_type == "both") {
    position <- "dodge"
  } else if (mf_type == "stacked") {
    position <- "stack"
  } else {
    position <- "identity"
  }
  
  # Plot type: mean value
  if (plot_type == "bar") {
    mean_value <- ggplot2::geom_bar(data = group_mean,
                                    aes(x = group_col,
                                        y = Mean),
                                    stat = "identity",
                                    position = position,
                                    color = "black")
  } else if (plot_type == "line") {
    mean_value <- ggplot2::geom_point(data = group_mean,
                                      aes(x = group_col,
                                          y = Mean),
                                      shape = "\U2014",
                                      size = 9,
                                      position = position,
                                      color = "black")
  }
  # Error bars
  if (plot_error_bars) {
    error_bars <- ggplot2::geom_errorbar(data = group_mean,
                                         aes(x = group_col,
                                             ymin = Mean - SE,
                                             ymax = Mean + SE),
                                         position = position,
                                         width = 0.2)
  } else {
    error_bars <- NULL
  }
  # Individual values
  if (plot_indiv_vals) {
    indiv_vals <- ggplot2::geom_point(data = plot_data,
                                      aes(x = group_col,
                                          y = mf_col,
                                          fill = fill),
                                      shape = 21,
                                      size = 3,
                                      color = "black",
                                      position = ggplot2::position_jitter(width = 0.1,
                                                                          height = 0,
                                                                          seed = 123),
                                      inherit.aes = FALSE)
  } else {
    indiv_vals <- NULL
  }
  # Labels
  if (add_labels == "indiv_count") {
    label <- plot_data$sum_col
  } else if (add_labels == "indiv_MF") {
    label <- sprintf("%.2e", plot_data$mf_col)
  } else if (add_labels == "mean_count") {
    label <- group_mean$Mean # doesn't make sense.
  } else if (add_labels == "mean_MF") {
    label <- sprintf("%.2e", group_mean$Mean)
  } else if (add_labels == "none") {
    label <- ""
  }

  if (add_labels %in% c("indiv_count", "indiv_MF")) {

     labels <- ggrepel::geom_text_repel(data = plot_data,
                                        ggplot2::aes(x = group_col,
                                            y = mf_col,
                                            label = label),
                                       size = 3,
                                       color = "black",
                                       position = ggplot2::position_jitter(width = 0.1,
                                               height = 0,
                                               seed = 123),
                                       max.overlaps = Inf,
                                       inherit.aes = FALSE)
  } else if (add_labels %in% c("mean_count", "mean_MF")) {
    if (mf_type %in% c("min", "max")) {
      label_position <- "identity"
    } else if (mf_type == "both") {
      label_position <- ggplot2::position_dodge(width = 0.9)
    } else {
      label_position <- ggplot2::position_stack(vjust = 0.5)
    }
    labels <- ggplot2::geom_text(data = group_mean,
                               ggplot2::aes(x = group_col,
                                            y = Mean,
                                            label = label),
                               position = label_position,
                               vjust = -0.5,
                               size = 4,
                               color = "black")
  } else {
    labels <- NULL
  }

  # plot
  p <- ggplot2::ggplot() +
    mean_value +
    error_bars +
   # indiv_vals +
   # labels +
    yscale +
    ggplot2::xlab(x_lab) +
    ggplot2::ylab(y_lab) +
    ggplot2::ggtitle(title) +
  #  ggplot2::scale_fill_manual(values = palette) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line())
p
return(p)
}