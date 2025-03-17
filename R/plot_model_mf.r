#' Plot your mf model
#' @description Provide a visualization of the point estimates derived using
#' model_mf()
#' @param model A model object created using model_mf()
#' @param plot_type The type of plot to create. Options are "bar" or "point".
#' @param x_effect If there are multiple fixed effects in the model, specify
#' the fixed effect to plot on the x-axis. The other will be used in the fill
#' aesthetic. Currently, only 2 fixed effects are supported.
#' @param ref_effect The fixed effect to use as the reference level when adding
#' significance labels. Only applicable if using two fixed effects.
#' @param plot_error_bars Logical. If TRUE, the estimated standard error will
#' be added to the plot.
#' @param plot_signif Logical. If TRUE, will add significance labels based on
#' the pairwise_comparisons data frame in the model object. This is only valid
#' if you supplied a contrasts table to model_mf(). Symbols will be applied
#' to plotted values that are significantly different from the reference. Your
#' contrasts table is structured as a data frame with two columns, each
#' containing levels of the fixed effects to be contrasted. When adding
#' significance labels, symbols will be added to the values defined in the
#' first column, while the second column will represent the reference. A
#' different symbol will be used for each unique reference level. If a single
#' plotted value has been contrasted against multiple references, then it will
#' gain multiple symbols for each significance difference.
#' @param x_order A character vector indicating the order of the levels for the
#' x_effect.
#' @param fill_order A character vector indicating the order of the levels for
#' the fill aesthetic, if applicable.
#' @param x_label The label for the x-axis.
#' @param y_label The label for the y-axis.
#' @param plot_title The title of the plot.
#' @param fill_label The label for the fill aesthetic, if applicable.
#' @param custom_palette A vector of colors to use for the fill and color
#' aesthetics. If not provided, a default palette will be used. When plotting
#' a model that has a single fixed effect, you can specify colors for "fill" and
#' "color" using a named vector. Likewise, when plotting a model with two fixed
#' effects, you can specify colors for the levels within your fill variable.
#' @details See model_mf() for examples.
#' @return A ggplot object.
#' @export
#' @importFrom dplyr select ends_with filter mutate group_by summarize left_join
#' rename across all_of
#' @importFrom ggplot2 aes geom_bar geom_errorbar geom_point
#' geom_text ggplot guides labs theme position_dodge scale_fill_manual
#' element_blank element_line
#' @importFrom grDevices colorRampPalette

plot_model_mf <- function(model,
                          plot_type = "point",
                          x_effect = NULL,
                          plot_error_bars = TRUE,
                          plot_signif = TRUE,
                          ref_effect = NULL,
                          x_order = NULL,
                          fill_order = NULL,
                          x_label = NULL,
                          y_label = NULL,
                          plot_title = NULL,
                          fill_label = NULL,
                          custom_palette = NULL) {

  # Check if point_estimates exist in the model_object
  if (!"point_estimates" %in% names(model)) {
    stop("The model object does not contain 'point_estimates'")
  }

  # Validate plot_type
  if (plot_type != "bar" && plot_type != "point") {
    stop("Invalid plot type. Choose either 'bar' or 'point'.")
  }

  plot_data <- model$point_estimates

  # Determine the x-axis variable
  nfixef <- ncol(plot_data) - 4
  if (nfixef == 1) {
    x_var <- names(plot_data)[5]
    other_effect <- NULL
  } else if (nfixef > 1) {
    if (is.null(x_effect)) {
      stop("Multiple fixed effects found. Please specify the fixed effect to plot on the x-axis using x_effect.")
    }
    if (!x_effect %in% names(plot_data)) {
    stop(paste("The fixed effect", x_effect, "is not found in point_estimates dataframe."))
    }
    x_var <- x_effect
    other_effect <- setdiff(names(plot_data)[5:(4 + nfixef)], x_effect)
  }

  # Order the data
  # x-axis order
  if (!is.null(x_order)) {
    plot_data[[x_var]] <- factor(plot_data[[x_var]], levels = x_order)
  }
  # Fill order
  if (!is.null(other_effect) && !is.null(fill_order)) {
    plot_data[[other_effect]] <- factor(plot_data[[other_effect]],
                                        levels = fill_order)
  }

  # Initialize the plot
  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = !!sym(x_var),
                                    y = plot_data$Estimate))

  # Add fill for other_effect
  if (!is.null(other_effect)) {
    p <- p + ggplot2::aes(color = !!sym(other_effect), fill = !!sym(other_effect))
    n_colors <- length(unique(plot_data[[other_effect]]))
  } else {
    plot_data$fill <- "fill"
    plot_data$color <- "color"
    p <- p + ggplot2::aes(fill = plot_data$fill, color = plot_data$color) +
      ggplot2::guides(fill = "none", color = "none")
    n_colors <- 1
  }
  # Palette
  if (is.null(custom_palette)) {
    gradient <- grDevices::colorRampPalette(colors = c("#c5e5fc",
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
    if (length(palette) < n_colors) {
      stop(paste0("The number of colors in the custom_palette is less than the number of levels in the fill aesthetic. ", n_colors, " needed but ", length(custom_palette), " provided."))
    }
  }
  # Estimate: bars or points
  if (plot_type == "bar") {
    pos <- ggplot2::position_dodge(0.9)
    bar <- ggplot2::geom_bar(stat = "identity",
                             position = pos)
    point <- NULL
  } else {
    pos <- ggplot2::position_dodge(width = 0.5)
    point <- ggplot2::geom_point(position = pos)
    bar <- NULL
  }
  # Error bars
  if (plot_error_bars) {
    plot_data$y_position_label <- plot_data$Estimate + plot_data$Std.Err
    error <- ggplot2::geom_errorbar(ggplot2::aes(ymin = plot_data$Estimate - plot_data$Std.Err,
                                                 ymax = plot_data$Estimate + plot_data$Std.Err),
                                    color = "black",
                                    width = 0.2,
                                    position = pos)
  } else {
    plot_data$y_position_label <- plot_data$Estimate
    error <- NULL
  }

  if (plot_signif) {
    if (!"pairwise_comparisons" %in% names(model)) {
      stop("The model object does not contain 'pairwise_comparisons'.
           Cannot plot significance labels. Please supply a 'contrast' table to
           model_mf() when creating your model object.")
    }
    signif_data <- model$pairwise_comparisons %>%
      dplyr::select(dplyr::ends_with("1"),
                    dplyr::ends_with("2"),
                    "Significance") %>%
      dplyr::filter(Significance != "")
    names(signif_data) <- sub("_1$", "", names(signif_data))

    # Define the reference effect
    if (is.null(ref_effect)) {
      split_names <- strsplit(row.names(signif_data), " vs ")
      signif_data$ref_level <- sapply(split_names, function(x) x[2])
    } else {
      signif_data <- dplyr::rename(signif_data,
                                   ref_level = paste0(ref_effect, "_2"))
    }
    # Assign a unique symbol based on the reference level
    unique_ref_levels <- unique(signif_data$ref_level)
    nref <- length(unique_ref_levels)
    if (nref == 1) {
      symbols_pool <- c("*")
    } else if (nref > 1 && nref <= 13) {
      symbols_pool <- c(strsplit("*†‡¥¢¤●Ø@#$%&", "")[[1]])
    } else if (nref > 13 && nref <= 26) {
      symbols_pool <- c(letters)
    } else {
      warning("Too many reference levels. Cannot add distinct significance labels.")
      symbols_pool <- rep("*", nref)
    }
    random_symbols <- sample(symbols_pool, nref, replace = FALSE)
    random_symbols <- setNames(random_symbols, unique_ref_levels)

    signif_data <- signif_data %>%
      dplyr::mutate(Symbol = random_symbols[ref_level]) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(x_var, other_effect)))) %>%
      dplyr::summarize(Significance = paste(unique(Symbol),
                                            collapse = ""),
                       .groups = 'drop')

    plot_data <- dplyr::left_join(plot_data, signif_data)
    plot_data$Significance <- ifelse(is.na(plot_data$Significance),
                                     "",
                                     plot_data$Significance)

    # plot the significance labels
    signif <- ggplot2::geom_text(ggplot2::aes(y = plot_data$y_position_label,
                                              label = plot_data$Significance),
                                 position = pos,
                                 vjust = -0.5,
                                 size = 3,
                                 color = "black")

    # Create a legend for the significance labels using dummy points
    significance_labels <- paste(random_symbols,
                                 "Significant differences from",
                                 names(random_symbols))
    plot_data$significance_labels <- rep(significance_labels,
                                         length.out = nrow(plot_data))
    plot_data$dummy_labels_y <- 0
    # Plot an invisible shape
    dummy <- ggplot2::geom_point(aes(y = plot_data$dummy_labels_y,
                                     shape = plot_data$significance_labels),
                                 size = 0,
                                 color = "white",
                                 show.legend = TRUE)
  } else {
    signif <- NULL
    dummy <- NULL
  }

  # Labels
  if (!is.null(x_label)) {
    xlab <- x_label
  } else {
    xlab <- x_var
  }
  if (!is.null(y_label)) {
    ylab <- y_label
  } else {
    ylab <- "Mean MF Model Estimate"
  }

  if (!is.null(plot_title)) {
    title <- plot_title
  } else {
    title <- "Model Estimates"
  }

  if (!is.null(fill_label)) {
    flabel <- fill_label
  } else {
    if (is.null(other_effect)) {
      flabel <- ""
    } else {
      flabel <- other_effect
    }
  }

  p <- p +
    dummy +
    bar +
    error +
    point +
    signif +
    ggplot2::labs(x = xlab,
                  y = ylab,
                  fill = flabel,
                  color = flabel,
                  title = title,
                  shape = "Significance") +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line()) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_color_manual(values = palette)

  return(p)
}
