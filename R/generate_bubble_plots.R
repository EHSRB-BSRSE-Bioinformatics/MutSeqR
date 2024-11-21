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
generate_bubble_plots <- function(mutation_data,
                                  size_by = "alt_depth",
                                  facet_col = NULL,
                                  color_by = "normalized_subtype",
                                  circle_spacing = 1,
                                  circle_outline = "none",
                                  circle_resolution = 50) {
  if (!requireNamespace("fmsb")) {
    stop("You need the package fmsb to run this function.")
  }
  if (!requireNamespace("RColorBrewer")) {
    stop("You need the package RColorBrewer to run this function.")
  }
  # TO DO: I haven't fully implemented the color_by option
  # ideally it would conditionally build the color pallet depending on the value of this variable
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
    plotcolors <- NULL  # Use NULL to indicate automatic color generation
  } else {
    plotcolors <- NULL
  }

  # TODO - allow some text labels on the circles? Perhaps for clonally expanded mutants above a certain size (e.g., >1), we could show a number?

  # Ensure additional_columns is a character vector
  if (!is.null(facet_col) && !is.character(facet_col)) {
    stop("facet_col must be a character vector")
  }

  x <- mutation_data %>% 
    dplyr::filter(!.data$variation_type %in% "no_variant" &
                    .data$is_germline == FALSE) %>% # TO DO: add a parameter to filter out germline mutations
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
    # Convert facet column to a factor to handle arbitrary values
    data$facet <- as.factor(data$facet)
    # Get levels
    facet_levels <- levels(data$facet)
    # Loop through facets
    circles <- lapply(
                      seq_along(facet_levels),
                      function(i) {
                        facet_level <- facet_levels[[i]]
                        # Filter data for the current facet_level
                        filtered_data <- data %>% filter(.data$facet == facet_level)
                        
                        # Generate the layout for the filtered data with scaled response
                        circles <- packcircles::circleProgressiveLayout(filtered_data,
                                                          sizecol = "response", sizetype = 'area')
                        circles$radius <- circle_spacing * circles$radius
                        data2 <- cbind(filtered_data, circles)
                        return(list(data2 = data2))
                      })
    data2 <- do.call(rbind, lapply(circles, `[[`, "data2"))

  } else {

    # Generate the layout for the filtered data with scaled response
    circles <- packcircles::circleProgressiveLayout(data,
                                                    sizecol = "response",
                                                    sizetype = 'area')
    circles$radius <- circle_spacing * circles$radius
    data2 <- cbind(data, circles)
  }

  # Select the specified columns from data2 for vertices
  vertices <- packcircles::circleLayoutVertices(data2,
                                                npoints = circle_resolution,
                                                idcol = "group",
                                                xysizecols = (ncol(data2) - 2):(ncol(data2)))

  # Perform the left join with the suffixes to handle different x and y columns
  plot_data <- left_join(vertices, data2, by = c("id" = "group"),
                         suffix = c(".circle_coords", ".circle_centres"))

  # Create the plot
  p <- ggplot() +
    geom_polygon(data = plot_data,
                 aes(x = x.circle_coords, y = y.circle_coords,
                     group = id,
                     fill = color_column),
                 colour = if(is.null(circle_outline) || circle_outline == "none") NA else circle_outline) +
    theme_void() +
    labs(fill = "Mutation type") +
    theme(legend.position = "right") +
    coord_equal()

  if (is.null(plotcolors)) {
    num_colors <- length(unique(plot_data$color_column))
    color_palette <- colorRampPalette(RColorBrewer::Rcolorbrewer.pal(8, "Set1"))(num_colors)
    p <- p + scale_fill_manual(values = color_palette)
  } else {
    p <- p + scale_fill_manual(values = plotcolors,
                               breaks = names(plotcolors))
  }

  if (!is.null(facet_col)) {
    p <- p + facet_wrap(~facet)
    return(p)
  } else {
    return(p)
  }
}