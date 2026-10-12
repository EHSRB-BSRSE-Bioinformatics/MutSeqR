#' Create a radar plot
#' @description Create a radar plot of mutation frequencies or other metrics
#' across categories (labels), optionally grouped or faceted by a categorical variable.
#' @param mf_data A data frame with the data to plot.
#' @param response_col The column with the response values (numeric).
#' @param label_col The column with the labels for the radar plot axes.
#' @param group_col The column with the group variable to plot or facet across. Default is NULL.
#' @param facet A logical indicating whether to create separate faceted plots for each
#' group in \code{group_col} (TRUE, default) or plot all groups on a single radar plot
#' (FALSE).
#' @param custom_palette A custom color palette to distinguish groups. Can be a named character
#' vector of colors (with names corresponding to levels/values in \code{group_col}),
#' an unnamed vector of colors, or the name of an RColorBrewer palette. Default is NULL.
#' @param label_colors A custom color specification for the radar axis labels. Can be
#' a named vector of colors (with names corresponding to values in \code{label_col}),
#' an unnamed character vector of colors matching the order of labels, or a single color string.
#' Any labels not specified in a named vector default to "black". Default is NULL ("black").
#' @param indiv_y A logical indicating whether to use individual y-axis scales
#' for each plot when \code{facet = TRUE}. When \code{facet = FALSE}, all groups share the
#' same scale. Default is TRUE.
#' @param plot_legend A logical indicating whether to show a legend when multiple groups
#' are plotted on the same radar chart (\code{facet = FALSE}). Default is TRUE.
#' @param legend_pos The position of the legend when \code{plot_legend = TRUE}. Can be a
#' keyword such as "topright", "topleft", "bottomright", "bottomleft", etc. Default is "topright".
#' @param vlcex A numeric value specifying character expansion for variable labels. Default is 1.2.
#' @param plot_title An optional title for the plot. If NULL, a default title is generated.
#' @param facet_col Deprecated. Use \code{group_col} instead.
#' @param ... Additional arguments passed to \code{fmsb::radarchart}.
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr select all_of
#' @importFrom grDevices colorRampPalette rainbow
#' @importFrom graphics layout legend par text
#' @examples
#' # Plot the MFmin for each variation type, including the 6 SNV subtypes
#' # Facet the plots by dose group
#' mf_ex <- readRDS(system.file("extdata", "Example_files", "mf_data_6.rds",
#'      package = "MutSeqR")
#' )
#' plot_radar(
#'      mf_data = mf_ex,
#'      response_col = "mf_min",
#'      label_col = "normalized_subtype",
#'      group_col = "dose_group",
#'      facet = TRUE,
#'      indiv_y = TRUE
#' )
#'
#' # Plot all groups on the same plot with custom palette and label colors
#' plot_radar(
#'      mf_data = mf_ex,
#'      response_col = "mf_min",
#'      label_col = "normalized_subtype",
#'      group_col = "dose_group",
#'      facet = FALSE,
#'      custom_palette = c("Control" = "#45567f", "Low" = "#6a74a7",
#'                         "Medium" = "#c3c9e9", "High" = "#e2aab9"),
#'      label_colors = c("C>A" = "red", "C>T" = "blue")
#' )
#' @return Invisible NULL (generates a base R graphics plot)
#' @export
plot_radar <- function(mf_data,
                       response_col,
                       label_col,
                       group_col = NULL,
                       facet = TRUE,
                       custom_palette = NULL,
                       label_colors = NULL,
                       indiv_y = TRUE,
                       plot_legend = TRUE,
                       legend_pos = "topright",
                       vlcex = 1.2,
                       plot_title = NULL,
                       facet_col = NULL,
                       ...) {
  # Backward compatibility: support facet_col if group_col is not provided
  if (is.null(group_col) && !is.null(facet_col)) {
    group_col <- facet_col
  }

  # --- 1. Validation ---
  stopifnot(
    !missing(mf_data) && is.data.frame(mf_data),
    !missing(response_col), !missing(label_col)
  )
  if (!requireNamespace("fmsb", quietly = TRUE)) {
    stop("Package 'fmsb' required. Please install it using install.packages('fmsb').")
  }

  # Validate columns exist
  cols_to_check <- c(response_col, label_col)
  if (!is.null(group_col)) cols_to_check <- c(cols_to_check, group_col)

  if (!all(cols_to_check %in% colnames(mf_data))) {
    stop("Columns not found: ", paste(setdiff(cols_to_check, colnames(mf_data)), collapse = ", "))
  }

  # --- 2. Data Preparation ---
  # Handle the optional group column
  if (is.null(group_col)) {
    dummy_group_name <- "All_Data"
    mf_data[[dummy_group_name]] <- "All Data"
    group_col <- dummy_group_name
    facet <- FALSE
    indiv_y <- FALSE
  }

  # Ensure consistent ordering if label_col is a factor
  if (is.factor(mf_data[[label_col]])) {
    label_levels <- levels(mf_data[[label_col]])
  } else {
    label_levels <- sort(unique(mf_data[[label_col]]))
  }

  # Pivot to Wide Format (Rows=Groups, Cols=Labels)
  plot_data <- mf_data %>%
    dplyr::select(dplyr::all_of(c(group_col, label_col, response_col))) %>%
    tidyr::pivot_wider(
      names_from = dplyr::all_of(label_col),
      values_from = dplyr::all_of(response_col),
      values_fill = 0
    )

  # Reorder columns to match label levels
  valid_levels <- intersect(label_levels, colnames(plot_data))
  plot_data <- plot_data[, c(group_col, valid_levels)]

  # Ensure group is a factor and preserve ordering
  if (!is.factor(plot_data[[group_col]])) {
    group_levels <- unique(plot_data[[group_col]])
    plot_data[[group_col]] <- factor(plot_data[[group_col]], levels = group_levels)
  } else {
    group_levels <- levels(plot_data[[group_col]])
    group_levels <- intersect(group_levels, unique(plot_data[[group_col]]))
  }
  # Arrange plot_data rows to match group_levels
  plot_data <- plot_data[match(group_levels, plot_data[[group_col]]), , drop = FALSE]

  data_matrix <- as.matrix(plot_data[, -1, drop = FALSE])
  mode(data_matrix) <- "numeric"
  n_groups <- length(group_levels)

  # --- 3. Resolve Colors ---
  # Resolve group colors (custom_palette)
  if (!is.null(custom_palette)) {
    is_brewer_name <- is.character(custom_palette) &&
      length(custom_palette) == 1 &&
      requireNamespace("RColorBrewer", quietly = TRUE) &&
      (custom_palette %in% rownames(RColorBrewer::brewer.pal.info))

    if (is_brewer_name) {
      max_colors <- RColorBrewer::brewer.pal.info[custom_palette, "maxcolors"]
      if (n_groups > max_colors) {
        group_colors <- grDevices::colorRampPalette(
          RColorBrewer::brewer.pal(n = max_colors, name = custom_palette)
        )(n_groups)
      } else {
        group_colors <- RColorBrewer::brewer.pal(n = max(3, n_groups), name = custom_palette)[seq_len(n_groups)]
      }
      names(group_colors) <- group_levels
    } else if (!is.null(names(custom_palette)) && any(group_levels %in% names(custom_palette))) {
      # Named vector
      if (requireNamespace("RColorBrewer", quietly = TRUE) && n_groups <= 8) {
        fallback_cols <- RColorBrewer::brewer.pal(max(3, n_groups), "Dark2")[seq_len(n_groups)]
      } else {
        fallback_cols <- grDevices::rainbow(n_groups)
      }
      names(fallback_cols) <- group_levels
      group_colors <- character(n_groups)
      names(group_colors) <- group_levels
      for (lvl in group_levels) {
        if (lvl %in% names(custom_palette) && !is.na(custom_palette[lvl])) {
          group_colors[lvl] <- custom_palette[lvl]
        } else {
          group_colors[lvl] <- fallback_cols[lvl]
        }
      }
    } else {
      # Unnamed vector
      group_colors <- rep_len(custom_palette, n_groups)
      names(group_colors) <- group_levels
    }
  } else {
    # Default group colors
    if (!facet && n_groups > 1) {
      if (requireNamespace("RColorBrewer", quietly = TRUE) && n_groups <= 8) {
        group_colors <- RColorBrewer::brewer.pal(max(3, n_groups), "Dark2")[seq_len(n_groups)]
      } else if (requireNamespace("RColorBrewer", quietly = TRUE) && n_groups <= 12) {
        group_colors <- RColorBrewer::brewer.pal(n_groups, "Set3")
      } else {
        group_colors <- grDevices::rainbow(n_groups)
      }
      names(group_colors) <- group_levels
    } else {
      group_colors <- rep("black", n_groups)
      names(group_colors) <- group_levels
    }
  }

  # Resolve label colors
  if (is.null(label_colors)) {
    resolved_label_colors <- rep("black", length(valid_levels))
  } else if (!is.null(names(label_colors))) {
    resolved_label_colors <- rep("black", length(valid_levels))
    names(resolved_label_colors) <- valid_levels
    matched <- intersect(names(label_colors), valid_levels)
    resolved_label_colors[matched] <- label_colors[matched]
    resolved_label_colors <- unname(resolved_label_colors)
  } else if (length(label_colors) == 1) {
    resolved_label_colors <- rep(label_colors, length(valid_levels))
  } else {
    resolved_label_colors <- rep_len(label_colors, length(valid_levels))
  }

  # Helper to draw colored variable labels around the perimeter
  draw_radar_labels <- function(label_names, colors, cex) {
    n <- length(label_names)
    theta <- seq(90, 450, length = n + 1) * pi / 180
    theta <- theta[1:n]
    xx <- cos(theta)
    yy <- sin(theta)
    graphics::text(
      x = xx * 1.2,
      y = yy * 1.2,
      labels = label_names,
      col = colors,
      cex = cex,
      xpd = TRUE
    )
  }

  # Capture user ... arguments
  user_args <- list(...)

  # --- 4. Plotting ---
  if (facet && n_groups > 1) {
    # Faceted plots: one plot per group
    if (!indiv_y) {
      global_max <- max(data_matrix, na.rm = TRUE) * 1.1
      if (is.na(global_max) || global_max == 0) global_max <- 1
    } else {
      row_maxes <- apply(data_matrix, 1, max, na.rm = TRUE) * 1.1
      row_maxes[is.na(row_maxes) | row_maxes == 0] <- 1
      global_max <- NULL
    }

    n_cols <- 2
    n_rows <- ceiling(n_groups / n_cols)
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit({
      graphics::layout(1)
      graphics::par(old_par)
    }, add = TRUE)
    graphics::layout(matrix(seq_len(n_rows * n_cols), nrow = n_rows, ncol = n_cols, byrow = TRUE))

    lapply(seq_along(group_levels), function(i) {
      group_name <- group_levels[i]
      current_vals <- data_matrix[i, ]
      current_max <- if (indiv_y) row_maxes[i] else global_max

      df_radar <- rbind(
        rep(current_max, length(current_vals)),
        rep(0, length(current_vals)),
        current_vals
      )
      df_radar <- as.data.frame(df_radar)
      colnames(df_radar) <- colnames(data_matrix)

      axis_labels <- sprintf("%.1e", seq(0, current_max, length.out = 5))

      if (is.null(plot_title)) {
        current_title <- if (group_col == "All_Data") "Radar Plot" else paste(group_col, group_name)
      } else {
        current_title <- paste(plot_title, "-", group_name)
      }

      current_pcol <- group_colors[i]

      radar_args <- list(
        df = df_radar,
        axistype = 1,
        caxislabels = axis_labels,
        seg = 4,
        vlabels = rep("", ncol(df_radar)),
        axislabcol = "#797777",
        vlcex = vlcex,
        title = current_title,
        pcol = current_pcol,
        plwd = 2,
        plty = 1,
        cglcol = "grey",
        cglty = 1,
        cglwd = 0.7
      )
      radar_args[names(user_args)] <- user_args
      radar_args$vlabels <- rep("", ncol(df_radar))

      do.call(fmsb::radarchart, radar_args)
      draw_radar_labels(colnames(df_radar), resolved_label_colors, vlcex)
    })

  } else {
    # Single plot: all groups on the same radar plot
    global_max <- max(data_matrix, na.rm = TRUE) * 1.1
    if (is.na(global_max) || global_max == 0) global_max <- 1

    df_radar <- rbind(
      rep(global_max, ncol(data_matrix)),
      rep(0, ncol(data_matrix)),
      data_matrix
    )
    df_radar <- as.data.frame(df_radar)
    colnames(df_radar) <- colnames(data_matrix)
    rownames(df_radar) <- c("max", "min", group_levels)

    axis_labels <- sprintf("%.1e", seq(0, global_max, length.out = 5))

    if (is.null(plot_title)) {
      current_title <- if (group_col == "All_Data") "Radar Plot" else paste("Radar Plot by", group_col)
    } else {
      current_title <- plot_title
    }

    radar_args <- list(
      df = df_radar,
      axistype = 1,
      caxislabels = axis_labels,
      seg = 4,
      vlabels = rep("", ncol(df_radar)),
      axislabcol = "#797777",
      vlcex = vlcex,
      title = current_title,
      pcol = unname(group_colors),
      plwd = 2,
      plty = 1,
      cglcol = "grey",
      cglty = 1,
      cglwd = 0.7
    )
    radar_args[names(user_args)] <- user_args
    radar_args$vlabels <- rep("", ncol(df_radar))

    do.call(fmsb::radarchart, radar_args)
    draw_radar_labels(colnames(df_radar), resolved_label_colors, vlcex)

    if (plot_legend && n_groups > 1) {
      legend_lty <- if (!is.null(radar_args$plty)) radar_args$plty else 1
      legend_lwd <- if (!is.null(radar_args$plwd)) radar_args$plwd else 2
      graphics::legend(
        x = legend_pos,
        legend = group_levels,
        col = unname(group_colors),
        lty = legend_lty,
        lwd = legend_lwd,
        bty = "n",
        cex = 0.9,
        xpd = TRUE
      )
    }
  }

  invisible(NULL)
}