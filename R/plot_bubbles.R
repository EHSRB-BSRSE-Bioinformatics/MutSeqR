#' Generate Bubble Plots
#'
#' Produces a ggplot object of bubble plots from given mutation data.
#' Optionally, bubble plots can be facetted by a specified column.
#'
#' @param mutation_data Data frame containing the mutation data.
#' @param size_by A string with the column name by which to size the circles.
#' Recommended values are "alt_depth" or "vaf".
#' @param facet_col A string with the column name by which to facet .
#' If NULL, no facetting will be done. Default is NULL.
#' @param color_by Character vector specifying how to color the mutations.
#' Accepted values are "normalized_subtype", "subtype", and
#' "trinucleotide_subtype". NOT FULLY IMPLEMENTED - open to any column name. Ex. contig
#' @param circle_spacing Numerical value to adjust the spacing between circles.
#' @param circle_outline Color for the circle outline. Default is "none", resulting in no outline color. Other accepted values are colors in the R language.
#' @param circle_resolution Number of points to use for the circle resolution. Default is 50.
#'
#' @return A ggplot object with the bubble plot, facetted if specified.
#'
#' @importFrom dplyr arrange filter left_join
#' @import ggplot2
#' @export
plot_bubbles <- function(mutation_data,
                                  size_by = "alt_depth",
                                  facet_col = NULL,
                                  color_by = "normalized_subtype",
                                  circle_spacing = 1,
                                  circle_outline = "none",
                                  circle_resolution = 50) {

  if (!requireNamespace("RColorBrewer")) {
    stop("You need the package RColorBrewer to run this function.")
  }
  
  if (color_by == "normalized_subtype") {
    plotcolors <- c("C>A" = "#3288BD",
                    "C>G" = "#99D594",
                    "C>T" = "#E6F598",
                    "T>A" = "#FEE08B",
                    "T>C" = "#FC8D59",
                    "T>G" = "#D53E4F",
                    "mnv" = "pink",
                    "deletion" = "black",
                    "insertion" = "grey",
                    "symbolic" = "purple")
  } else if (color_by == "subtype") {
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
                    "symbolic" = "azure2")
  } else if (color_by == "trinucleotide_subtype") {
    plotcolors <- NULL
  } else {
    plotcolors <- NULL
  }

  if (!is.null(facet_col) && !is.character(facet_col)) {
    stop("facet_col must be a character vector")
  }

  x <- mutation_data %>% 
    dplyr::filter(!.data$variation_type %in% "no_variant" &
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

  if (!is.null(facet_col)) {
    data$facet <- x[[facet_col]]
  }

  data <- data %>% dplyr::arrange("color_column")

  if (!is.null(facet_col)) {
    data$facet <- as.factor(data$facet)
    facet_levels <- levels(data$facet)
    # Pack circles for each facet level
    circles <- lapply(seq_along(facet_levels), function(i) {
      facet_level <- facet_levels[[i]]
      filtered_data <- data %>% filter(facet == facet_level)
      circle_layout <- packcircles::circleProgressiveLayout(filtered_data, sizecol = "response", sizetype = 'area')
      circle_layout$radius <- circle_spacing * circle_layout$radius
      data2 <- cbind(filtered_data, circle_layout)
      return(data2)
    })
    data2 <- do.call(rbind, circles)
  } else {
    circles <- packcircles::circleProgressiveLayout(data, sizecol = "response", sizetype = 'area')
    circles$radius <- circle_spacing * circles$radius
    data2 <- cbind(data, circles)
  }
  
  vertices <- packcircles::circleLayoutVertices(data2, npoints = circle_resolution, idcol = "group", xysizecols = c("x", "y", "radius"))
  
  plot_data <- left_join(vertices, data2, by = c("id" = "group"), suffix = c(".circle_coords", ".circle_centres"))
  
  # Create the main plot with bubbles
  p <- ggplot() +
    geom_polygon(data = plot_data,
                 aes(x = x.circle_coords, y = y.circle_coords, group = id, fill = color_column),
                 colour = if(is.null(circle_outline) || circle_outline == "none") NA else circle_outline) +
    theme_void() +
    labs(fill = "Mutation type") +
    theme(legend.position = "right") +
    coord_equal()
  
  if (is.null(plotcolors)) {
    num_colors <- length(unique(plot_data$color_column))
    color_palette <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(num_colors)
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
  selected_indices <- seq(1, length(unique_responses), length.out = num_legend_circles) %>% round()
  selected_sizes <- unique_responses[selected_indices]
  
  # Calculate the range of x-values in the data and determine the center point of the x-axis
  x_range <- range(data2$x, na.rm = TRUE)
  x_center <- mean(x_range)
  
  # Calculate the total length of the legend
  legend_width <- legend_gap * (num_legend_circles - 1)
  legend_gap <- max(data2$radius) * 2  # Adjust gap based on maximum circle radius
  
  # Adjust legend positions to be centered about x_center
  legend_start <- x_center - (legend_width / 2)
  legend_positions <- seq(legend_start, by = legend_gap, length.out = num_legend_circles)
  
  size_legend_df <- data.frame(
    x = legend_positions,
    y = rep(min(data2$y) - legend_gap, num_legend_circles),
    radius = circle_spacing * sqrt(selected_sizes / pi),
    label = as.character(selected_sizes)
  )
  legend_data <- packcircles::circleLayoutVertices(size_legend_df, idcol = NULL, xysizecols = 1:3)
  if (!is.null(facet_col)) {
    legend_data$facet <- levels(plot_data$facet)[length(levels(plot_data$facet))]
    size_legend_df$facet <- levels(plot_data$facet)[length(levels(plot_data$facet))]
  }
  # Overlay the size legend on the main plot
  q <- p +
    geom_polygon(data = legend_data,
                 aes(x = x, y = y, group = "label"),
                 fill = "grey80") +
    geom_text(data = size_legend_df, aes(x = x, y = y, label = label), vjust = 0.5)
  q
  return(q)
}