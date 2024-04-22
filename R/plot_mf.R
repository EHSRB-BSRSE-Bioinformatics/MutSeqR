#' Plot the Mutation Frequency
#' @description This function creates a bar plot of the mutation frequency
#' @param mf_data A data frame containing the mutation frequency data. This is
#' obtained from the calculate_mut_freq function with SUMMARY = TRUE.
#' @param sample_col The name of the column containing the sample names.
#' @param mf_type The type of mutation frequency to plot. Options are "unique",
#' "clonal", "both", or "stacked". If "both", the unique and clonal mutation
#' frequencies are plotted side by side. If "stacked", the difference between
#' the unique and clonal MF is stacked on top of the unique MF such that the
#' total height of both bars represent the clonal MF.
#' @param fill_col The name of the column containing the fill variable.
#' @param custom_palette A character vector of colour codes to use for the plot.
#' If NULL, a default palette is used
#' @param sample_order The order of the samples.
#' ' Options include:
#' \itemize{
#'   \item `none`: No ordering is performed. Default.
#'   \item `smart`: Samples are ordered based on the sample names.
#'   \item `arranged`: Samples are ordered based on one or more factor column(s)
#' in mf_data. Column names are passed to the function using the
#' `sample_order_input`.
#'  \item `custom`: Samples are ordered based on a custom vector of sample
#' names. The custom vector is passed to the function using the
#' `sample_order_input`.
#' }
#' @param sample_order_input The order of the samples if sample_order is
#' "custom". The column name by which to arrange samples if sample_order
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
#' @importFrom gtools mixedsort
#' 
#' @export
plot_mf <- function(mf_data,
                    sample_col,
                    mf_type = c("unique", "clonal", "both", "stacked"),
                    fill_col = NULL,
                    custom_palette = NULL,
                    sample_order = c("none", "smart", "arranged", "custom"),
                    sample_order_input = NULL,
                    labels = c("count", "MF", "none"),
                    scale_y_axis = "linear",
                    x_lab = NULL,
                    y_lab = NULL) {
  # axis_labels
  if (!is.null(x_lab)) {
    x_lab <- x_lab
  } else {
    x_lab <- sample_col
  }
  if (!is.null(y_lab)) {
     y_lab <- y_lab
  } else {
    y_lab <- "Mutation Frequency (mutations/bp)"
  }

  # Sample order
  if (sample_order == "none") {
    order <- as.vector(unique(mf_data[[sample_col]]))
    mf_data[[sample_col]] <- factor(mf_data[[sample_col]])
  } else if (sample_order == "smart") {
    order <- as.vector(unique(mf_data[[sample_col]]))
    order <- gtools::mixedsort(order)
    mf_data[[sample_col]] <- factor(mf_data[[sample_col]], levels = order)
  } else if (sample_order == "arranged") {
    mf_data <- mf_data %>%
      dplyr::arrange(dplyr::across(dplyr::all_of({{sample_order_input}})))
    order <- as.vector(unique(mf_data[[sample_col]]))
    mf_data[[sample_col]] <- factor(mf_data[[sample_col]], levels = order)
  } else if (sample_order == "custom") {
    mf_data[[sample_col]] <- factor(mf_data[[sample_col]],
                                    levels = sample_order_input)
  }


  if (mf_type %in% c("unique", "clonal")) {
    # response column
    MF_column_pattern <- paste0(".*(_MF_", mf_type, ")$")
    response_col <- names(mf_data)[grepl(MF_column_pattern, names(mf_data))]

    # sum column
    sum_column_pattern <- paste0(".*(_sum_", mf_type, ")$")
    found_count_col <- names(mf_data)[grepl(sum_column_pattern, names(mf_data))]

    plot_data <- mf_data %>%
      dplyr::rename(sample_col = dplyr::all_of(sample_col)) %>%
      dplyr::rename(mf_col = dplyr::all_of(response_col)) %>%
      dplyr::rename(sum_col = dplyr::all_of(found_count_col))
    max_y <- max(plot_data$mf_col) * 1.1
  } else {
    # response columns
    MF_unique_column_pattern <- paste0(".*(_MF_unique", ")$")
    unique_mf_col <- names(mf_data)[grepl(MF_unique_column_pattern, names(mf_data))]
    MF_clonal_column_pattern <- paste0(".*(_MF_clonal", ")$")
    clonal_mf_col <- names(mf_data)[grepl(MF_clonal_column_pattern, names(mf_data))]

    # sum columns
    sum_unique_column_pattern <- paste0(".*(_sum_unique", ")$")
    unique_count_col <- names(mf_data)[grepl(sum_unique_column_pattern, names(mf_data))]
    sum_clonal_column_pattern <- paste0(".*(_sum_clonal", ")$")
    clonal_count_col <- names(mf_data)[grepl(sum_clonal_column_pattern, names(mf_data))]

    plot_data <- mf_data %>%
      dplyr::rename(sample_col = dplyr::all_of(sample_col)) %>%
      dplyr::rename(mf_unique = dplyr::all_of(unique_mf_col)) %>%
      dplyr::rename(mf_clonal = dplyr::all_of(clonal_mf_col)) %>%
      dplyr::rename(sum_unique = dplyr::all_of(unique_count_col)) %>%
      dplyr::rename(sum_clonal = dplyr::all_of(clonal_count_col))

    if (mf_type == "stacked") {
      plot_data <- transform(plot_data, mf_clonal = plot_data$mf_clonal - plot_data$mf_unique)
      max_y <- max(plot_data$mf_unique + plot_data$mf_clonal) * 1.1
    } else {
      max_y <- max(plot_data$mf_clonal) * 1.1
    }
    # pivot long
    plot_data <- reshape(plot_data,
                         varying = list(c("sum_unique", "sum_clonal"),
                                        c("mf_unique", "mf_clonal")),
                         v.names = c("sum_col", "mf_col"),
                         times = c("unique", "clonal"),
                         timevar = "mf_type",
                         direction = "long")
    if (mf_type == "both") {
      plot_data$mf_type <- factor(plot_data$mf_type, levels = c("unique", "clonal"))
    }
    if (mf_type == "stacked") {
      plot_data$mf_type <- factor(plot_data$mf_type, levels = c("clonal", "unique"))
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
    title <- paste0("Unique and Clonal Mutation frequency per ", sample_col)
  } else if (mf_type %in% c("unique", "clonal")){
    title <- paste0(mf_type, " mutation frequency per ", sample_col)
  }

  # palette
  if (is.null(custom_palette)) {
    if (mf_type %in% c("both", "stacked")) {
      n_colors <- length(unique(plot_data$fill_col)) * 2
    } else if (mf_type %in% c("unique", "clonal")) {
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

  plot <- ggplot(plot_data, aes(x = plot_data$sample_col,
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