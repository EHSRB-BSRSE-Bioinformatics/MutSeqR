#' Plot the Mutation Frequency
#' @description This function creates a bar plot of the mutation frequency
#' @param mf_data A data frame containing the mutation frequency data. This is
#' obtained from the calculate_mut_freq function with SUMMARY = TRUE.
#' @param group_col The name of the column containing the sample/group names.
#' @param mf_type The type of mutation frequency to plot. Options are "min",
#' "max", "both", or "stacked". If "both", the min and max mutation
#' frequencies are plotted side by side. If "stacked", the difference between
#' the min and max MF is stacked on top of the min MF such that the
#' total height of both bars represent the max MF.
#' @param fill_col The name of the column containing the fill variable.
#' @param custom_palette A character vector of colour codes to use for the plot.
#' If NULL, a default palette is used
#' @param group_order The order of the samples.
#' ' Options include:
#' \itemize{
#'   \item `none`: No ordering is performed. Default.
#'   \item `smart`: Samples are ordered based on the sample names.
#'   \item `arranged`: Samples are ordered based on one or more factor column(s)
#' in mf_data. Column names are passed to the function using the
#' `group_order_input`.
#'  \item `custom`: Samples are ordered based on a custom vector of sample
#' names. The custom vector is passed to the function using the
#' `group_order_input`.
#' }
#' @param group_order_input The order of the samples if group_order is
#' "custom". The column name by which to arrange samples if group_order
#' is "arranged"
#' @param labels The labels to use for the bars. Either "count", "MF", or
#' "none". Count labels display the number of mutations, MF labels display
#' the mutation frequency.
#' @param scale_y_axis The scale of the y axis. Either "linear" or "log".
#' @param x_lab The label for the x axis.
#' @param y_lab The label for the y axis.
#' @return A ggplot object
#' @import ggplot2
#' @importFrom dplyr arrange across all_of rename
#' @export
plot_mf <- function(mf_data,
                    group_col,
                    mf_type = c("min", "max", "both", "stacked"),
                    fill_col = NULL,
                    custom_palette = NULL,
                    group_order = c("none", "smart", "arranged", "custom"),
                    group_order_input = NULL,
                    labels = c("count", "MF", "none"),
                    scale_y_axis = "linear",
                    x_lab = NULL,
                    y_lab = NULL) {
  
  if (group_order == "smart") {
    if (!requireNamespace("gtools", quietly = TRUE)) {
      stop("Package gtools is required when using the 'smart' group_order option. Please install the package using 'install.packages('gtools')'")
    }
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
    order <- as.vector(unique(mf_data[[group_col]]))
    order <- gtools::mixedsort(order)
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
    MF_column_pattern <- paste0(".*(_MF_", mf_type, ")$")
    response_col <- names(mf_data)[grepl(MF_column_pattern, names(mf_data))]

    # sum column
    sum_column_pattern <- paste0(".*(_sum_", mf_type, ")$")
    found_count_col <- names(mf_data)[grepl(sum_column_pattern, names(mf_data))]

    plot_data <- mf_data %>%
      dplyr::rename(group_col = dplyr::all_of(group_col)) %>%
      dplyr::rename(mf_col = dplyr::all_of(response_col)) %>%
      dplyr::rename(sum_col = dplyr::all_of(found_count_col))
    max_y <- max(plot_data$mf_col) * 1.1
  } else {
    # response columns
    MF_min_column_pattern <- paste0(".*(_MF_min", ")$")
    min_mf_col <- names(mf_data)[grepl(MF_min_column_pattern, names(mf_data))]
    MF_max_column_pattern <- paste0(".*(_MF_max", ")$")
    max_mf_col <- names(mf_data)[grepl(MF_max_column_pattern, names(mf_data))]

    # sum columns
    sum_min_column_pattern <- paste0(".*(_sum_min", ")$")
    min_count_col <- names(mf_data)[grepl(sum_min_column_pattern, names(mf_data))]
    sum_max_column_pattern <- paste0(".*(_sum_max", ")$")
    max_count_col <- names(mf_data)[grepl(sum_max_column_pattern, names(mf_data))]

    plot_data <- mf_data %>%
      dplyr::rename(group_col = dplyr::all_of(group_col)) %>%
      dplyr::rename(mf_min = dplyr::all_of(min_mf_col)) %>%
      dplyr::rename(mf_max = dplyr::all_of(max_mf_col)) %>%
      dplyr::rename(sum_min = dplyr::all_of(min_count_col)) %>%
      dplyr::rename(sum_max = dplyr::all_of(max_count_col))

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
  if (!is.null(fill_col)) {
        plot_data <- dplyr::rename(plot_data, fill_col = dplyr::all_of(fill_col))
        if (mf_type %in% c("both", "stacked")) {
          fill <- interaction(plot_data$mf_type, plot_data$fill_col)
          fill_label <- paste("MF Type and", fill_col)
    } else {
      fill <- plot_data$fill_col
      fill_label <- paste(fill_col)
    }
  } else {
    if (mf_type %in% c("both", "stacked")) {
      fill <- plot_data$mf_type
      fill_label <- "MF Type"
      plot_data$fill_col <- ""
    } else {
      fill <- plot_data$fill_col <- ""
      fill_label <- NULL
    }
  }

  # bar labels
  if (labels == "count") {
    label <- plot_data$sum_col
  } else if (labels == "MF") {
    label <- sprintf("%.2e", plot_data$mf_col)
  } else if (labels == "none") {
    label <- ""
  }

  # scale y axis
  if (scale_y_axis == "log") {
    yscale <- scale_y_log10()
  } else {
    yscale <- scale_y_continuous(limits = c(0, max_y))
  }

  # Position
  if (mf_type == "both") {
    position <- "dodge"
    label_position <- position_dodge(width = 0.9)
  } else if (mf_type == "stacked") {
    position <- "stack"
    label_position <- position_stack(vjust = 0.5)
  } else {
    position <- "identity"
    label_position <- "identity"
  }

  # Title
  if (mf_type %in% c("stacked", "both")) {
    title <- paste0("Mininimum and Maximum Mutation frequency per ", group_col)
  } else if (mf_type %in% c("min", "max")){
    title <- paste0(mf_type, " mutation frequency per ", group_col)
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
                                            "#b23946",
                                            "#ff5264",
                                            "#ffb9c1",
                                            "#ffedef"))
    palette <- gradient(n_colors)
  } else {
    palette <- custom_palette
  }

  plot <- ggplot(plot_data, aes(x = plot_data$group_col,
                                y = plot_data$mf_col,
                                fill = factor(fill))) +
    geom_bar(stat = "identity", position = position, color = "darkgrey") +
    geom_text(aes(y = plot_data$mf_col, label = label),
              position = label_position,
              vjust = -0.5,
              size = 4,
              color = "black") +
    yscale +
    labs(title = title,
         fill = fill_label) +
    ylab(y_lab) +
    xlab(x_lab) +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          panel.background = element_blank(),
          axis.line = element_line()) +
    scale_fill_manual(values = palette)

  return(plot)
}