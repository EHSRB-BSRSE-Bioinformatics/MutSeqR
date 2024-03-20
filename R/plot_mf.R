#' Plot the Mutation Frequency
#' @description This function creates a bar plot of the mutation frequency
#' @param mf_data A data frame containing the mutation frequency data
#' @param sample_col The name of the column containing the sample names
#' @param mf_type The type of mutation frequency to plot. Either "unique" or "clonal"
#' @param fill_col The name of the column containing the fill variable.
#' @param sample_order The order of the samples. Either "none", "smart", "arranged", or "custom"
#' @param sample_order_input The order of the samples if sample_order is "custom". The column name
#' by which to arrange samples if sample_order is "arranged"
#' @param labels The labels to use for the bars. Either "count", "MF", or "none"
#' @param scale_y_axis The scale of the y axis. Either "linear" or "log"
#' @param x_lab The label for the x axis
#' @param y_lab The label for the y axis
#' @return A ggplot object
#' @import ggplot2
#' @importFrom dplyr arrange across all_of rename
#' @export
#' 
plot_mf <- function(mf_data,
                    sample_col,
                    mf_type = c("unique", "clonal"),
                    fill_col = NULL,
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
    y_lab <- paste0(mf_type, "Mutation Frequency (mutations/bp)")
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

# fill column
if (!is.null(fill_col)) {
    plot_data <- dplyr::rename(plot_data, fill_col = dplyr::all_of(fill_col))
} else {
    plot_data$fill_col <- "all"
}

# bar labels
if (labels == "count") {
    label <- plot_data$sum_col
} else if (labels == "MF") {
    label <- plot_data$mf_col
} else if (labels == "none") {
    label <- NULL
}

# scale y axis
if (scale_y_axis == "log") {
    yscale <- scale_y_log10()
} else {
    yscale <- scale_y_continuous(limits = c(0, max(plot_data$mf_col) * 1.1))
}

plot <- ggplot(plot_data, aes(x = sample_col, y = mf_col, fill = factor(fill_col))) +
      geom_bar(stat = "identity", alpha = 0.7) +
      geom_text(aes(y = mf_col, label = label), vjust = -0.5) +
      yscale +
      labs(title = paste0(mf_type, " mutation frequency per ", sample_col),
           fill = paste(fill_col)) +
      ylab(y_lab) +
      xlab(x_lab) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
return(plot)
}