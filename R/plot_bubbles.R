#' Generate Bubble Plots
#'
#' Produces a ggplot object of bubble plots from given mutation data.
#' Optionally, bubble plots can be facetted and
#' coloured by a specified column.
#'
#' @param mutation_data Data frame containing the mutation data.
#' @param size_by The column name by which to size the circles.
#' Recommended values are "alt_depth" or "vaf".
#' @param facet_col The column name by which to facet .
#' If NULL, no facetting will be done. Default is NULL.
#' @param color_by The column name by which to colour the mutations. Default is
#' "normalized_subtype".
#' @param circle_spacing Numerical value to adjust the spacing between circles.
#' Default is 1.
#' @param circle_outline Colour for the circle outline. Default is "none",
#' resulting in no outline colour. Other accepted values are colours in the R
#' language.
#' @param circle_resolution Number of points to use for the circle resolution.
#' Default is 50.
#' @param circle_palette Circle colours are generated using
#' /code{RColorBrewer::brewer.pal(n=8, name = "Set1")}. If you want to use a
#' different `brewer.pal` palette, you can specify it here. Default is NULL;
#' "Set1". `color_by =` "normalized_subtype" and "subtype" have a default color
#' palette to match the MutSeqR plot_spectra() function but these can also be
#' changed to the specified brewer.pal palette. See
#' \code{\link[RColorBrewer]{brewer.pal}} for palette options. You may
#' visualize the palettes at the ColorBrewer website:
#' \url{https://colorbrewer2.org/}.
#' @details The function will plot a circle for each mutation in
#' `mutation_data`. Mutations flagged by the `filter_mut` column will be
#' excluded from the plot. The size of the circle is determined by the
#' `size_by` parameter. Sizing by the "alt_depth" or the "vaf" will give users
#' the ability to visualize the the distribution of recurrent mutations within
#' their data with large multiplets having a large circle.
#' @return A ggplot object with the bubble plot, facetted if specified.
#' @examples
#' example_file <- system.file("extdata", "example_mutation_data_filtered.rds", package = "MutSeqR")
#' example_data <- readRDS(example_file)
#' plot <- plot_bubbles(mutation_data = example_data,
#'                      facet_col = "dose_group")
#' @importFrom dplyr arrange filter left_join
#' @import ggplot2
#' @export
## fix awkward legend facetting
plot_bubbles <- function(mutation_data,
                         size_by = "alt_depth",
                         facet_col = NULL,
                         color_by = "normalized_subtype",
                         circle_spacing = 1,
                         circle_outline = "none",
                         circle_resolution = 50,
                         circle_palette = NULL) {

  if (!requireNamespace("RColorBrewer")) {
    stop("You need the package RColorBrewer to run this function.")
  }
  if (!requireNamespace("packcircles")) {
    stop("You need the package packcircles to run this function.")
  }

  if (color_by == "normalized_subtype" && is.null(circle_palette)) {
    plotcolors <- c("C>A" = "#3288BD",
                    "C>G" = "#99D594",
                    "C>T" = "#E6F598",
                    "T>A" = "#FEE08B",
                    "T>C" = "#FC8D59",
                    "T>G" = "#D53E4F",
                    "mnv" = "pink",
                    "deletion" = "black",
                    "insertion" = "grey",
                    "sv" = "purple",
                    "ambiguous" = "darkgrey",
                    "uncategorized" = "white")
  } else if (color_by == "subtype" && is.null(circle_palette)) {
    plotcolors <- c("A>C" = "limegreen",
                    "A>G" = "forestgreen",
                    "A>T" = "darkgreen",
                    "C>A" = "skyblue1",
                    "C>G" = "dodgerblue2",
                    "C>T" = "darkblue",
                    "G>A" = "grey28",
                    "G>C" = "grey25",
                    "G>T" = "grey0",
                    "T>A" = "red2",
                    "T>C" = "red3",
                    "T>G" = "red4",
                    "mnv" = "hotpink",
                    "deletion" = "yellow",
                    "insertion" = "purple",
                    "sv" = "azure2",
                    "ambiguous" = "darkgrey",
                    "uncategorized" = "white")
  } else {
    plotcolors <- NULL
    if (is.null(circle_palette)) {
      palette <- "Set1"
    } else {
      palette <- circle_palette
    }
  }

  if (!is.null(facet_col) && !is.character(facet_col)) {
    stop("facet_col must be a character vector")
  }

  # Scale the size_by column if needed:
  # scale so smallest value is just above 1
  sizeby <- mutation_data[[size_by]]
  if (min(sizeby) < 1) {
    min_non_zero_value <- min(sizeby[sizeby > 0])
    scale_factor <- 10^(-floor(log10(min_non_zero_value)))
    scale_factor <- ifelse(min_non_zero_value * scale_factor < 1,
                           scale_factor * 10,
                           scale_factor)
  } else {
    scale_factor <- 1
  }
  mutation_data[[size_by]] <- mutation_data[[size_by]] * scale_factor
  x <- mutation_data %>%
    dplyr::filter(.data$variation_type != "no_variant" &
                    .data$filter_mut == FALSE) %>%
    dplyr::arrange(!!rlang::sym(color_by))
  data <- data.frame(group = paste(x$sample,
                                   x$contig,
                                   x$start,
                                   x$end,
                                   x$ref,
                                   x$alt,
                                   sep = "_"),
                     response = x[[size_by]],
                     color_column = x[[color_by]])

  data <- data %>% dplyr::arrange("color_column")
  data$color_column <- factor(data$color_column)

  if (!is.null(facet_col)) {
    if (is.factor(x[[facet_col]])) {
      original_facet_levels <- levels(x[[facet_col]])
    } else {
      original_facet_levels <- unique(x[[facet_col]])
    }
    all_facet_levels <- c(original_facet_levels, "Legend")
    data$facet <- factor(x[[facet_col]], levels = all_facet_levels)

    # Pack circles for each facet level
    circles <- lapply(seq_along(original_facet_levels), function(i) {
      facet_level <- original_facet_levels[[i]]
      filtered_data <- data %>% dplyr::filter(.data$facet == facet_level)
      circle_layout <- packcircles::circleProgressiveLayout(filtered_data, sizecol = "response", sizetype = 'area')
      circle_layout$radius <- circle_spacing * circle_layout$radius
      data2 <- cbind(filtered_data, circle_layout)
      return(data2)
    })
    data2 <- do.call(rbind, circles)
  } else {
    circles <- packcircles::circleProgressiveLayout(data,
                                                    sizecol = "response",
                                                    sizetype = "area")
    circles$radius <- circle_spacing * circles$radius
    data2 <- cbind(data, circles)
  }

  vertices <- packcircles::circleLayoutVertices(data2,
                                                npoints = circle_resolution,
                                                idcol = "group",
                                                xysizecols = c("x", "y", "radius"))

  plot_data <- dplyr::left_join(vertices,
                                data2,
                                by = c("id" = "group"),
                                suffix = c(".circle_coords", ".circle_centres"))

  # Create the main plot with bubbles
  p <- ggplot() +
    geom_polygon(data = plot_data,
                 aes(x = x.circle_coords,
                     y = y.circle_coords,
                     group = id,
                     fill = color_column),
                 colour = if(is.null(circle_outline) || circle_outline == "none") NA else circle_outline) +
    theme_void() +
    labs(fill = color_by) +
    theme(legend.position = "right") +
    coord_equal()

  if (is.null(plotcolors)) {
    num_colors <- length(unique(plot_data$color_column))
    color_palette <- colorRampPalette(RColorBrewer::brewer.pal(8, name = palette))(num_colors)
    p <- p + scale_fill_manual(values = color_palette)
  } else {
    p <- p + scale_fill_manual(values = plotcolors, breaks = names(plotcolors))
  }

  if (!is.null(facet_col)) {
    p <- p + facet_wrap(~facet)
  }

  # Generate a size legend within the same plot
  unique_responses <- sort(unique(data$response))
  num_legend_circles <- min(5, length(unique_responses))  # Up to 5 elements
  selected_indices <- seq(1,
                          length(unique_responses),
                          length.out = num_legend_circles) %>%
    round()
  selected_sizes <- unique_responses[selected_indices]

  # Calculate the range of x-values in the data and determine the center point of the x-axis
  x_range <- range(data2$x, na.rm = TRUE)
  x_center <- mean(x_range)

  # Calculate the total length of the legend
  legend_gap <- max(data2$radius) * 2  # Adjust gap based on maximum circle radius
  legend_width <- legend_gap * (num_legend_circles - 1)

  # Adjust legend positions to be centered about x_center
  legend_start <- x_center - (legend_width / 2)
  legend_positions <- seq(legend_start,
                          by = legend_gap,
                          length.out = num_legend_circles)

  size_legend_df <- data.frame(
    x = legend_positions,
    y = rep(min(data2$y) - legend_gap, num_legend_circles),
    radius = circle_spacing * sqrt(selected_sizes / pi),
    label = selected_sizes / scale_factor
  )
  size_legend_df <- size_legend_df %>%
    dplyr::mutate(label = ifelse(.data$label < 1,
                                 sprintf("%#.2e", label),
                                 as.character(label)))
  legend_data <- packcircles::circleLayoutVertices(size_legend_df, idcol = NULL, xysizecols = 1:3)
  if (!is.null(facet_col)) {
    legend_data$facet <- factor("Legend", levels = all_facet_levels)
    size_legend_df$facet <- factor("Legend", levels = all_facet_levels)
  }
  # Overlay the size legend on the main plot
  q <- p +
    geom_polygon(data = legend_data,
                 aes(x = x, y = y, group = "label"),
                 fill = "grey80") +
    geom_text(data = size_legend_df,
              aes(x = x, y = y + max(radius) * 1.2, label = label),
              vjust = 0.5)

  return(q)
}