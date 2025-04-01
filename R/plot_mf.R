#' Plot the Mutation Frequency
#' @description This function creates a plot of the mutation frequency.
#' @param mf_data A data frame containing the mutation frequency data. This is
#' obtained from the calculate_mf function with SUMMARY = TRUE.
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
#' @param labels The data labels to display on the plot. Either "count", "MF",
#' or "none". Count labels display the number of mutations, MF labels display
#' the mutation frequency.
#' @param scale_y_axis The scale of the y axis. Either "linear" or "log".
#' @param x_lab The label for the x axis.
#' @param y_lab The label for the y axis.
#' @param title The title of the plot.
#' @return A ggplot object
#' @examples
#' example_file <- system.file("extdata", "Example_files",
#'                             "example_mutation_data_filtered.rds",
#'                             package = "MutSeqR")
#' example_data <- readRDS(example_file)
#' example_data$dose_group <- factor(example_data$dose_group,
#'                                   levels = c("Control", "Low",
#'                                              "Medium", "High"))
#' mf <- calculate_mf(mutation_data = example_data,
#'                    cols_to_group = "sample",
#'                    subtype_resolution = "none",
#'                    retain_metadata_cols = "dose_group")
#' plot <- plot_mf(mf_data = mf,
#'                 group_col = "sample",
#'                 plot_type = "bar",
#'                 mf_type = "min",
#'                 fill_col = "dose_group",
#'                 group_order = "arranged",
#'                 group_order_input = "dose_group",
#'                 labels = "count",
#'                 title = "Mutation Frequency per Sample")
#' @import ggplot2
#' @importFrom dplyr arrange across all_of rename
#' @export
plot_mf <- function(mf_data,
                    group_col,
                    plot_type = "bar",
                    mf_type = "min",
                    fill_col = NULL,
                    custom_palette = NULL,
                    group_order = "none",
                    group_order_input = NULL,
                    labels = "count",
                    scale_y_axis = "linear",
                    x_lab = NULL,
                    y_lab = NULL,
                    title = NULL) {
  
  if (group_order == "smart" && !requireNamespace("gtools", quietly = TRUE)) {
      stop("Package gtools is required when using the 'smart' group_order option. Please install the package using 'install.packages('gtools')'")
  }
  if (mf_type == "stacked" && plot_type == "point") {
   stop("The 'stacked' mutation frequency type is not compatible with the 'point' plot type.")
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

  # Sample order
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

  if (mf_type %in% c("min", "max")) {
    # response column
    response_col <- paste0("mf_", mf_type)

    # sum column
    found_count_col <- paste0("sum_", mf_type)

    plot_data <- mf_data %>%
      dplyr::rename(group_col = dplyr::all_of(group_col)) %>%
      dplyr::rename(mf_col = dplyr::all_of(response_col)) %>%
      dplyr::rename(sum_col = dplyr::all_of(found_count_col))
    max_y <- max(plot_data$mf_col) * 1.1
  } else {
    plot_data <- mf_data %>%
      dplyr::rename(group_col = dplyr::all_of(group_col))

    if (mf_type == "stacked") {
      plot_data <- transform(plot_data, mf_max = plot_data$mf_max - plot_data$mf_min)
      max_y <- max(plot_data$mf_min + plot_data$mf_max) * 1.1
    } else {
      max_y <- max(plot_data$mf_max) * 1.1
    }
    # pivot long
    plot_data <- reshape(plot_data,
                         varying = list(c("sum_min", "sum_max"),
                                        c("mf_min", "mf_max")),
                         v.names = c("sum_col", "mf_col"),
                         times = c("min", "max"),
                         timevar = "mf_type",
                         direction = "long")
    if (mf_type == "both") {
      plot_data$mf_type <- factor(plot_data$mf_type, levels = c("min", "max"))
    }
    if (mf_type == "stacked") {
      plot_data$mf_type <- factor(plot_data$mf_type, levels = c("max", "min"))
    }
  }

  # fill column
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
  }

  # labels
  if (labels == "count") {
    label <- plot_data$sum_col
  } else if (labels == "MF") {
    label <- sprintf("%.2e", plot_data$mf_col)
  } else if (labels == "none") {
    label <- ""
  }

  # scale y axis
  if (scale_y_axis == "log") {
    yscale <- ggplot2::scale_y_log10()
  } else {
    yscale <- ggplot2::scale_y_continuous(limits = c(0, max_y))
  }

  # Position
  if (mf_type == "both" && plot_type == "bar") {
    position <- "dodge"
    label_position <- ggplot2::position_dodge(width = 0.9)
  } else if (mf_type == "stacked" && plot_type == "bar") {
    position <- "stack"
    label_position <- ggplot2::position_stack(vjust = 0.5)
  } else {
    position <- "identity"
    label_position <- "identity"
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

 # define the plot type
  if (plot_type == "bar") {
    type <- ggplot2::geom_bar(stat = "identity",
                              position = position,
                              color = "black")
  labels <- ggplot2::geom_text(ggplot2::aes(label = label),
                               position = label_position,
                               vjust = -0.5,
                               size = 3,
                               color = "black")
  } else if (plot_type == "point") {
    pos <- ggplot2::position_jitter(width = 0.1,
                                    height = 0,
                                    seed = 123)
    type <- ggplot2::geom_point(shape = 21,
                                size = 3,
                                color = "black",
                                position = pos)
    labels <- ggrepel::geom_text_repel(aes(label = label),
                                       size = 3,
                                       color = "black",
                                       position = pos,
                                       max.overlaps = Inf)
  }

  # Create the plot
  plot <- ggplot2::ggplot(plot_data,
                          ggplot2::aes(x = plot_data$group_col,
                                       y = plot_data$mf_col,
                                       fill = factor(fill))) +
    type +
    labels +
    yscale +
    ggplot2::labs(title = title,
                  fill = fill_label) +
    ggplot2::ylab(y_lab) +
    ggplot2::xlab(x_lab) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                       hjust = 1,
                                                       vjust = 0.5),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line()) +
    ggplot2::scale_fill_manual(values = palette)

  return(plot)
}