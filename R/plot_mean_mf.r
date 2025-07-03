#' Plot the Mean Mutatation Frequency
#' @description This function calculates the mean mutation frequency across
#' samples for given groups and plots the results.
#' @param mf_data A data frame containing the mutation frequency data. This is
#' obtained from the calculate_mf function with SUMMARY = TRUE.
#' @param group_col The column(s) in mf_data by which to calculate the mean.
#' When supplying more than one column, the values of all group columns will
#' be concatenated into a single value by which to calculate the mean. Values
#' will be displayed along the x-axis. Ex. "dose" or c("dose", "tissue").
#' @param fill_col An optional column name in the data used to define the fill
#' aesthetic in the plot. If fill_col has multiple levels within each group_col
#' level, the mean will be calculated for each level of fill_col
#' (recommend plot_type = "line" for this use case). Default is NULL.
#' @param mf_type The type of mutation frequency to plot. Options are "min",
#' "max", "both", or "stacked". If "both", the min and max mutation
#' frequencies are plotted side by side. "stacked" can be chosen for bar
#' plot_type only. If "stacked", the difference between the min and max
#' MF is stacked on top of the min MF such that the total height of both
#' bars represent the max MF. Default is "min".
#' @param plot_type The type of plot to create. Options are "bar" or "line".
#' Default is "bar".
#' @param plot_error_bars Whether to plot the error bars. Default is TRUE.
#' Error bars are standard error of the mean.
#' @param plot_indiv_vals Whether to plot the individual values as data points.
#' Default is FALSE.
#' @param group_order The order of the groups along the x-axis.
#' ' Options include:
#' \itemize{
#'   \item `none`: No ordering is performed. Default.
#'   \item `smart`: Groups are ordered based on the sample names.
#'   \item `arranged`: Groups are ordered based on one or more factor column(s)
#' in mf_data. Factor column names are passed to the function using the
#' `group_order_input`.
#'  \item `custom`: Groups are ordered based on a custom vector of group
#' names. The custom vector is passed to the function using the
#' `group_order_input`.
#' }
#' @param group_order_input The order of the groups if group_order is
#' "custom". The column name by which to arrange groups if group_order
#' is "arranged". If "custom", and using more than one group_col, values
#' are concatenated in the order listed, separated by a "_".
#' @param add_labels The data labels to display on the plot. Either
#' "indiv_count", "indiv_MF", "mean_count", "mean_MF", or "none".
#' Count labels display the number of mutations, MF labels display the mutation
#' frequency. Mean plots the mean value. Indiv plots the labels for individual
#' data points (only if plot_indiv_vals = TRUE). Default is "none".
#' @param plot_title The title of the plot. Default is
#' "Mean Mutation Frequency".
#' @param x_lab The x-axis label. Default is the value of group_col.
#' @param y_lab The y-axis label. Default is "Mutation Frequency
#' (mutations/bp)".
#' @param scale_y_axis The scale of the y axis. Either "linear" or "log".
#' Default is "linear".
#' @param custom_palette A custom color palette to use for the plot. Input a
#' character vector of colours. Input a named character vector to specify
#' olours to specific groups. Fill labels will be constructed by the following
#' components
#' 1. "Mean/Individual" if plot_indiv_vals = TRUE, fill labels will specify
#' Mean/Individual values.
#' 2. "min/max" if mf_type = "both" or "stacked", fill labels will specify
#' min/max values.
#' 3. fill_col value. Name colours to match the fill labels. Default is NULL.
#' If no custom_palette, a rainbow palette is generated. Min/Max values and
#' Mean/Individual values will be the same colour, different shades.
#' @param plot_legend Logical. Whether to show the fill (and color) legend.
#' Default is TRUE.
#' @param rotate_labels A logical value indicating whether data labels should
#' be rotated 90 degrees. Default is FALSE.
#' @param label_size A numeric value that controls the size of the data labels.
#' @return a ggplot object
#' @examples
#' # Example data consists of 24 mouse bone marrow DNA samples imported
#' # using import_mut_data() and filtered with filter_mut as in Example 4.
#' # Sequenced on TS Mouse Mutagenesis Panel. Example data is
#' # retrieved from MutSeqRData, an ExperimentHub data package.
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' example_data <- eh[["EH9861"]]
#' 
#' example_data$dose_group <- factor(example_data$dose_group,
#'                                   levels = c("Control", "Low",
#'                                              "Medium", "High"))
#' mf <- calculate_mf(mutation_data = example_data,
#'                    cols_to_group = "sample",
#'                    subtype_resolution = "none",
#'                    retain_metadata_cols = "dose_group")
#' plot <- plot_mean_mf(mf_data = mf,
#'                      group_col = "dose_group",
#'                      mf_type = "min",
#'                      plot_type = "line",
#'                      fill_col = "dose_group",
#'                      plot_error_bars = TRUE,
#'                      plot_indiv_vals = TRUE,
#'                      add_labels = "none")
#' @import ggplot2
#' @importFrom dplyr across all_of arrange rename group_by summarize
#' @importFrom stats sd setNames
#' @export

plot_mean_mf <- function(mf_data,
                         group_col = "dose",
                         fill_col = NULL,
                         mf_type = "both",
                         plot_type = "line",
                         plot_error_bars = TRUE,
                         plot_indiv_vals = TRUE,
                         group_order = "none",
                         group_order_input = NULL,
                         add_labels = "mean_count",
                         scale_y_axis = "linear",
                         x_lab = NULL,
                         y_lab = NULL,
                         plot_title = NULL,
                         custom_palette = NULL,
                         plot_legend = TRUE,
                         rotate_labels = FALSE,
                         label_size = 3) {
  # load required packages
  if (group_order == "smart" && !requireNamespace("gtools", quietly = TRUE)) {
      stop("Package gtools is required when using the 'smart' group_order option. Please install the package using 'install.packages('gtools')'")
  }
  if (add_labels %in% c("indiv_count", "indiv_MF") && !requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package ggrepel is required when using the 'indiv_count' or 'indiv_MF' add_labels options. Please install the package using 'install.packages('ggrepel')'")
  }

  if (is.null(custom_palette) && !requireNamespace("colorspace", quietly = TRUE)) {
    stop("Package colorspace is required when using default color palette. Please install the package using 'install.packages('colorspace')'")
  }

  if (!plot_indiv_vals && add_labels %in% c("indiv_count", "indiv_MF")) {
  stop("plot_indiv_vals must be TRUE when add_labels is set to 'indiv_count' or 'indiv_MF'")
  }

  if (mf_type == "stacked" && plot_error_bars && plot_type == "bar") {
    stop("Error bars are currently not supported with mf_type 'stacked', plot_type 'bar'. Sorry!")
  }
  # Concat group_col (if applicable)
  mf_data$group_col <- do.call(paste, c(mf_data[group_col], sep = "_"))

  # axis_labels
  if (!is.null(x_lab)) {
    x_lab <- x_lab
  } else {
    x_lab <- paste(group_col, collapse = " & ")
  }
  if (!is.null(y_lab)) {
    y_lab <- y_lab
  } else {
    y_lab <- "Mutation Frequency (mutations/bp)"
  }
  # title
  if (is.null(plot_title)) {
    if (mf_type == "min") {
      title <- paste("MFmin per", paste(group_col, collapse = " & "))
    } else if (mf_type == "max") {
      title <- paste0("MFmax per ", paste(group_col, collapse = " & "))
    } else {
      title <- paste0("MF per ", paste(group_col, collapse = " & "))
    }
  } else {
    title <- plot_title
  }

  # x-axis order
  if (group_order == "none") {
    mf_data$group_col <- factor(mf_data$group_col)
  } else if (group_order == "smart") {
    order <- gtools::mixedsort(as.vector(unique(mf_data$group_col)))
    mf_data$group_col <- factor(mf_data$group_col, levels = order)
  } else if (group_order == "arranged") {
    mf_data <- mf_data %>%
      dplyr::arrange(dplyr::across(dplyr::all_of({{group_order_input}})))
    order <- as.vector(unique(mf_data$group_col))
    mf_data$group_col <- factor(mf_data$group_col, levels = order)
  } else if (group_order == "custom") {
    mf_data$group_col <- factor(mf_data$group_col,
                                levels = group_order_input)
  }

  # Plot Data
  ## Make indiv plot data
  indiv_data <- mf_data

  if (!is.null(fill_col)) {
    if (fill_col %in% group_col) {
      indiv_data$fill_col <- indiv_data$group_col
    } else {
      indiv_data <- dplyr::rename(indiv_data,
                                  fill_col = dplyr::all_of(fill_col))
    }
  } else {
    indiv_data$fill_col <- "f1ll_c0l"
  }
  ## Make Group Mean data
  mean_data <- indiv_data %>%
    dplyr::group_by(.data$group_col, .data$fill_col) %>%
    dplyr::summarize(
      min_Mean = mean(.data$mf_min, na.rm = TRUE),
      min_SE = sd(.data$mf_min, na.rm = TRUE) / sqrt(dplyr::n()),
      min_sum_mean = mean(.data$sum_min, na.rm = TRUE),
      max_Mean = mean(.data$mf_max, na.rm = TRUE),
      max_SE = sd(.data$mf_max, na.rm = TRUE) / sqrt(dplyr::n()),
      max_sum_mean = mean(.data$sum_max, na.rm = TRUE),
      .groups = "drop"
    )
  mean_data <- as.data.frame(mean_data)

  ## Set the max y value
  if (mf_type == "min") {
    y_max <- max(indiv_data$mf_min) * 1.05
  } else {
    y_max <- max(indiv_data$mf_max) * 1.05
  }

  ## Rename the columns we want to plot to generic names
  if (mf_type == "min") {
    indiv_data <- dplyr::rename(indiv_data,
                                mf_col = "mf_min",
                                sum_col = "sum_min")
    mean_data <- dplyr::rename(mean_data,
                               Mean = "min_Mean",
                               SE = "min_SE",
                               mean_sum = "min_sum_mean")
  }
  if (mf_type == "max") {
    indiv_data <- dplyr::rename(indiv_data,
                                mf_col = "mf_max",
                                sum_col = "sum_max")
    mean_data <- dplyr::rename(mean_data,
                               Mean = "max_Mean",
                               SE = "max_SE",
                               mean_sum = "max_sum_mean")
  }
  ### Change MFmax value for stacked bar
  if (mf_type == "stacked" && plot_type == "bar") {
    mean_data$max_og <- mean_data$max_Mean
    mean_data$min_og <- mean_data$min_Mean
    mean_data <- transform(
      mean_data,
      max_Mean = mean_data$max_Mean - mean_data$min_Mean
    )
  }
  ### Pivot the data to long format for both and stacked options
  if (mf_type %in% c("both", "stacked")) {
    indiv_data <- reshape(indiv_data,
                          varying = list(c("sum_min", "sum_max"),
                                         c("mf_min", "mf_max")),
                          v.names = c("sum_col", "mf_col"),
                          times = c("min", "max"),
                          timevar = "mf_type",
                          direction = "long")
    if (mf_type == "stacked" && plot_type == "bar") {
      mean_data <- reshape(mean_data,
                           varying = list(c("min_Mean", "max_Mean"),
                                          c("min_SE", "max_SE"),
                                          c("min_sum_mean", "max_sum_mean"),
                                          c("min_og", "max_og")),
                           v.names = c("Mean", "SE", "mean_sum", "og_mf"),
                           times = c("min", "max"),
                           timevar = "mf_type",
                           direction = "long")
    } else {
      mean_data <- reshape(mean_data,
                           varying = list(c("min_Mean", "max_Mean"),
                                          c("min_SE", "max_SE"),
                                          c("min_sum_mean", "max_sum_mean")),
                           v.names = c("Mean", "SE", "mean_sum"),
                           times = c("min", "max"),
                           timevar = "mf_type",
                           direction = "long")
    }
  } # end plot data

  # scale y axis
  if (scale_y_axis == "log") {
    yscale <- ggplot2::scale_y_log10()
  } else {
    yscale <- ggplot2::scale_y_continuous(limits = c(0, y_max))
  }

  # Set the fill column: min/max
  if (mf_type %in% c("both", "stacked")) {
    mean_data$mean_fill_col <- paste(mean_data$mf_type, mean_data$fill_col)
    # Remove trailing spaces
    indiv_data$indiv_fill_col <- paste(indiv_data$mf_type, indiv_data$fill_col)
  } else { # if mf_type is min or max
    mean_data$mean_fill_col <- paste(mean_data$fill_col)
    indiv_data$indiv_fill_col <- paste(indiv_data$fill_col)
  }
  if (plot_indiv_vals) {
    indiv_data$indiv_fill_col <- paste("Individual", indiv_data$indiv_fill_col)
    mean_data$mean_fill_col <- paste("Mean", mean_data$mean_fill_col)
  }
  # Remove the f1ll_c0l placeholder
  mean_data$mean_fill_col <- sub(" f1ll_c0l$", "", mean_data$mean_fill_col)
  indiv_data$indiv_fill_col <- sub(" f1ll_c0l$", "", indiv_data$indiv_fill_col)
  mean_data$mean_fill_col <- stringr::str_trim(mean_data$mean_fill_col,
                                               side = "both")
  indiv_data$indiv_fill_col <- stringr::str_trim(indiv_data$indiv_fill_col,
                                                 side = "both")
  # end fill column

  # Palette
  if (is.null(custom_palette)) {
    fill_values <- unique(indiv_data$fill_col)
    palette <- grDevices::rainbow(length(fill_values))
    palette <- setNames(palette, fill_values)
    # Function to generate lighter/darker shades
    generate_shades <- function(color, steps = 2) {
      # Blend the input color with white to generate lighter shades
      shades <- colorspace::lighten(
        color,
        amount = seq(0, 0.5, length.out = steps)
      )
      return(shades)
    }
    # Generate shades for min/max
    if (mf_type == "both" || mf_type == "stacked") {
      palette <- lapply(palette, function(color) {
        generate_shades(color)  # Lighter for min, darker for max
      })
      palette <- unlist(lapply(names(palette), function(name) {
        setNames(palette[[name]], paste(c("max", "min"), name))
      }), recursive = FALSE)
    }
    # Generate shades for mean/indiv
    if (plot_indiv_vals) {
      palette <- lapply(palette, function(color) {
        generate_shades(color)  # Lighter for min, darker for max
      })
      palette <- unlist(lapply(names(palette), function(name) {
        setNames(palette[[name]], paste(c("Individual", "Mean"), name))
      }), recursive = FALSE)
    }
    if (is.null(fill_col)) {
      names(palette) <- sub(" ?f1ll_c0l$", "", names(palette))
    }
    if (plot_type == "line") {
      palette[grepl("^Mean", names(palette))] <- "#000000"
    }
  } else {
    palette <- custom_palette
  }
  mean_data$mean_fill_col <- sub(" f1ll_c0l$", "", mean_data$mean_fill_col)
  indiv_data$indiv_fill_col <- sub(" f1ll_c0l$", "", indiv_data$indiv_fill_col)

  # set the fill_col order for both/stacked.
  if (mf_type == "both") { # 1. min, 2. max so max is to the right of min
    sorted_levels <- unique(indiv_data$indiv_fill_col)[order(grepl("max", unique(indiv_data$indiv_fill_col)), unique(indiv_data$indiv_fill_col))]
    indiv_data$indiv_fill_col <- factor(indiv_data$indiv_fill_col, levels = sorted_levels)
    sorted_levels <- unique(mean_data$mean_fill_col)[order(grepl("max", unique(mean_data$mean_fill_col)), unique(mean_data$mean_fill_col))]
    mean_data$mean_fill_col <- factor(mean_data$mean_fill_col,
                                      levels = sorted_levels)
  }
  if (mf_type == "stacked") { # 1. max, 2. min so that max is stacked on top of min
    sorted_levels <- unique(indiv_data$indiv_fill_col)[order(grepl("min", unique(indiv_data$indiv_fill_col)), unique(indiv_data$indiv_fill_col))]
    indiv_data$indiv_fill_col <- factor(indiv_data$indiv_fill_col, levels = sorted_levels)
    sorted_levels <- unique(mean_data$mean_fill_col)[order(grepl("min", unique(mean_data$mean_fill_col)), unique(mean_data$mean_fill_col))]
    mean_data$mean_fill_col <- factor(mean_data$mean_fill_col,
                                      levels = sorted_levels)
  }
  # Position
  if (mf_type == "both") {
    bar_position <- ggplot2::position_dodge(width = 0.9)
    line_position <- ggplot2::position_dodge(width = 0.5)
    if (plot_type == "line") {
      error_position <- ggplot2::position_dodge(width = 0.5)
      indiv_position <- ggplot2::position_jitterdodge(dodge.width = 0.5,
                                                      jitter.width = 0.1,
                                                      jitter.height = 0,
                                                      seed = 123)
      indiv_label_position <- ggplot2::position_jitterdodge(dodge.width = 0.5,
                                                            jitter.width = 0.1,
                                                            jitter.height = 0,
                                                            seed = 123)
      mean_label_position <- ggplot2::position_dodge(width = 0.5)
    } else {
      error_position <- ggplot2::position_dodge(width = 0.9)
      indiv_position <- ggplot2::position_jitterdodge(dodge.width = 0.9,
                                                      jitter.width = 0.1,
                                                      jitter.height = 0,
                                                      seed = 123)
      indiv_label_position <- ggplot2::position_jitterdodge(dodge.width = 0.9,
                                                            jitter.width = 0.1,
                                                            jitter.height = 0,
                                                            seed = 123)
      mean_label_position <- ggplot2::position_dodge(width = 0.9)
    }
  } else if (mf_type == "stacked") {
    bar_position <- "stack"
    line_position <- "identity"
    error_position <- "identity"
    indiv_position <- ggplot2::position_jitter(width = 0.1,
                                               height = 0,
                                               seed = 123)
    indiv_label_position <- ggplot2::position_jitter(width = 0.1,
                                                     height = 0,
                                                     seed = 123)
    mean_label_position <- ggplot2::position_stack(vjust = 0.5)
  } else {
    bar_position <- "identity"
    line_position <- "identity"
    error_position <- "identity"
    indiv_position <- ggplot2::position_jitter(width = 0.1,
                                               height = 0,
                                               seed = 123)
    indiv_label_position <- ggplot2::position_jitter(width = 0.1,
                                                     height = 0,
                                                     seed = 123)
    mean_label_position <- "identity"
  }

  # Plot type: mean value
  if (plot_type == "bar") {
    mean_bar <- ggplot2::geom_bar(data = mean_data,
                                  aes(x = mean_data$group_col,
                                      y = mean_data$Mean,
                                      fill = mean_data$mean_fill_col),
                                  stat = "identity",
                                  position = bar_position,
                                  color = "black")
    error_bar_size <- 0.2
  } else {
    mean_bar <- NULL
  }
  if (plot_type == "line") {
    mean_line <- ggplot2::geom_point(
      data = mean_data,
      aes(x = mean_data$group_col,
          y = mean_data$Mean,
          fill = "black",
          color = mean_data$mean_fill_col),
      shape = "\U2014",
      size = 11,
      position = line_position
    )
    if (mf_type == "both") {
      error_bar_size <- 0.2
    } else {
      error_bar_size <- 0.1
    }

  } else {
    mean_line <- NULL
  }

  # Error bars
  if (plot_error_bars) {
    error_bars <- ggplot2::geom_errorbar(
      data = mean_data,
      aes(x = mean_data$group_col,
          ymin = mean_data$Mean - mean_data$SE,
          ymax = mean_data$Mean + mean_data$SE,
          group = interaction(mean_data$group_col,
                              mean_data$mean_fill_col)),
      position = error_position,
      width = error_bar_size,
      inherit.aes = FALSE,
      show.legend = FALSE
    )
  } else {
    error_bars <- NULL
  }

  # Individual values
  if (plot_indiv_vals) {
    indiv_vals <- ggplot2::geom_point(data = indiv_data,
                                      aes(x = indiv_data$group_col,
                                          y = indiv_data$mf_col,
                                          fill = indiv_data$indiv_fill_col),
                                      shape = 21,
                                      size = 2,
                                      color = "black",
                                      position = indiv_position,
                                      inherit.aes = FALSE)
  } else {
    indiv_vals <- NULL
  }

  # Data Labels
  if (add_labels == "indiv_count") {
    label <- indiv_data$sum_col
  } else if (add_labels == "indiv_MF") {
    label <- sprintf("%.2e", indiv_data$mf_col)
  } else if (add_labels == "mean_count") {
    label <- round(mean_data$mean_sum)
  } else if (add_labels == "mean_MF") {
    if (mf_type == "stacked" && plot_type == "bar") {
      label <- sprintf("%.2e", mean_data$og_mf)
    } else {
      label <- sprintf("%.2e", mean_data$Mean)
    }
  } else if (add_labels == "none") {
    label <- ""
  }

  if (add_labels %in% c("indiv_count", "indiv_MF")) {
    # label parameters
    if (rotate_labels) {
      label_angle <- 90
    } else {
      label_angle <- 0
    }
    labels <- ggrepel::geom_text_repel(
      data = indiv_data,
      ggplot2::aes(x = indiv_data$group_col,
                   y = indiv_data$mf_col,
                   label = label,
                   color = indiv_data$indiv_fill_col),
      size = label_size,
      angle = label_angle,
      position = indiv_label_position,
      max.overlaps = Inf,
      inherit.aes = FALSE
    )
  } else if (add_labels %in% c("mean_count", "mean_MF")) {

    if (mf_type == "stacked" && plot_type == "bar") {
      mean_data$label_position <- mean_data$og_mf
    } else {
      mean_data$label_position <- mean_data$Mean
    }

    if (plot_error_bars) {
      mean_data <- transform(
        mean_data,
        label_position = mean_data$label_position + mean_data$SE
      )
    }
    # Set label params
    if (rotate_labels) {
      label_angle <- 90
      vjust <- 0.5
      hjust <- -0.5
    } else {
      label_angle <- 0
      vjust <- -0.5
      hjust <- 0.5
    }
    labels <- ggplot2::geom_text(
      data = mean_data,
      ggplot2::aes(x = mean_data$group_col,
                   y = mean_data$label_position,
                   label = label,
                   group = interaction(mean_data$group_col, mean_data$mean_fill_col)),
      position = mean_label_position,
      vjust = vjust,
      hjust = hjust,
      angle = label_angle,
      size = label_size,
      color = "black"
    )
  } else {
    labels <- NULL
  }
  # plot
  p <- ggplot2::ggplot() +
    mean_bar +
    indiv_vals +
    error_bars +
    labels +
    mean_line +
    yscale +
    ggplot2::xlab(x_lab) +
    ggplot2::ylab(y_lab) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line()) +
    ggplot2::scale_fill_manual(values = palette) +
    labs(fill = NULL,
         color = NULL)

  if (plot_type == "line") {
    p <- p + ggplot2::scale_color_manual(values = palette)
  }
  if (!plot_legend) {
    p <- p + guides(fill = "none", color = "none")
  }

  return(p)
}
