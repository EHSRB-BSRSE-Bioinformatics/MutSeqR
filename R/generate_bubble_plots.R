#' Generate Bubble Plots
#'
#' Produces a ggplot object of bubble plots from given mutation data.
#' Optionally, bubble plots can be facetted by a specified column.
#'
#' @param mutation_data Data frame containing the mutation data.
#' @param facet_col A string with the column name to facet by. If NULL or not provided, no faceting is performed.
#' @param circle_spacing Numerical value to adjust the spacing between circles.
#' @param color_by Character vector specifying how to color the mutations. Accepted values are "normalized_subtype", "subtype", and "trinucleotide_subtype". NOT FULLY IMPLEMENTED.
#' @param circle_outline Color for the circle outline. Default is "none", resulting in no outline color. Other accepted values are colors in the R language.
#' @param circle_resolution Number of points to use for the circle resolution. Default is 50.
#' 
#' @return A ggplot object with the bubble plot, facetted if specified.
#' 
#' @importFrom dplyr arrange filter left_join
#' @import ggplot2
#' @export
generate_bubble_plots <- function(mutation_data,
                                  facet_col = "dose",
                                  circle_spacing = 1,
                                  color_by = "normalized_subtype",
                                  circle_outline = "none",
                                  circle_resolution = 50) {
  if (!requireNamespace("fmsb")) {
    stop("You need the package fmsb to run this function.")
  }
  # I haven't fully implemented the color_by option
  # ideally it would conditionally build the color pallet depending on the value of this variable
  if (color_by == "normalized_subtype") {
    plotcolors <- c( "C>A" = "#05C3EF",
                     "C>G" = "#000000",
                     "C>T" = "#E62F29",
                     "T>A" = "#D0CFCF",
                     "T>C" = "#A9D46C",
                     "T>G" = "#EECDCC",
                     "mnv" = "#66C2A5",
                     "deletion" = "#FC8D62",
                     "insertion" = "#8DA0CB",
                     "symbolic" = "#FFD92F")
  } else if (color_by == "subtype") {
      plotcolors <- NA
  } else if (color_by == "trinucleotide_subtype") {
      plotcolors <- NA
  } else {
    stop("Invalid color_by value. Use one of: 'normalized_subtype', 'subtype', or 'trinucleotide_subtype'.")
  }

  # TODO - allow some text labels on the circles? Perhaps for clonally expanded mutants above a certain size (e.g., >1), we could show a number?
  
  # Ensure additional_columns is a character vector
  if (!is.null(facet_col) && !is.character(facet_col)) {
    stop("facet_col must be a character vector")
  }
  
  #mutation_data <- readRDS("./Bbf_BM_mut_dat.RDS")
  
  x <- mutation_data %>%
    dplyr::filter(!.data$variation_type %in% "no_variant" &
                    .data$is_germline == FALSE) %>%
    dplyr::arrange(!!rlang::sym(color_by))
  
  data <- data.frame(group = paste0(x$sample, "_",
                                    x$contig, "_",
                                    x$start, "_",
                                    x$end, "_",
                                    x$ref, "_",
                                    x$alt),
                     alt_depth = x$alt_depth,
                     normalized_subtype = x$normalized_subtype,
                     subtype = x$subtype,
                     trinucleotide_subtype = x$normalized_context_with_mutation,
                     color_column = x[[color_by]],
                     dose = x$dose) %>%
    dplyr::arrange(color_column)
  
  if (!is.null(facet_col)) {
    # Convert dose column to a factor to handle arbitrary values
    data[[facet_col]] <- as.factor(data[[facet_col]])
    # Get levels
    facet_levels <- levels(data[[facet_col]])
    # Loop through facets
    circles <- lapply(
                      seq_along(facet_levels),
                      function(i) {
                        facet_level <- facet_levels[[i]]
                        # Filter data for the current facet_level
                        filtered_data <- data %>% filter(.[[facet_col]] == facet_level)
                        
                        # Generate the layout for the filtered data with scaled alt_depth
                        circles <- packcircles::circleProgressiveLayout(filtered_data,
                                                          sizecol = "alt_depth", sizetype = 'area')
                        circles$radius <- circle_spacing * circles$radius
                        data2 <- cbind(filtered_data, circles)
                        return(list(data2 = data2))
                      })
    data2 <- do.call(rbind, lapply(circles, `[[`, "data2"))
    
  } else {

    # Generate the layout for the filtered data with scaled alt_depth
    circles <- packcircles::circleProgressiveLayout(data,
                                       sizecol = "alt_depth", sizetype = 'area')
    circles$radius <- circle_spacing * circles$radius
    data2 <- cbind(data, circles)
  }

  # Select the specified columns from data2 for vertices
  vertices <- packcircles::circleLayoutVertices(data2, npoints = circle_resolution, idcol = "group",
                                   xysizecols = (ncol(data2) - 2):(ncol(data2)))
  
  # Perform the left join with the suffixes to handle different x and y columns
  plot_data <- left_join(vertices, data2, by = c("id" = "group"),
                         suffix = c(".circle_coords", ".circle_centres"))
  
  # Create the plot
  p <- ggplot() +
    geom_polygon(data = plot_data,
                 aes(x = x.circle_coords, y = y.circle_coords,
                     group = id,
                     fill = !!rlang::sym(color_by)),
                 colour = if(is.null(circle_outline) || circle_outline == "none") NA else circle_outline) +
    scale_fill_manual(values = plotcolors,
                      breaks = names(plotcolors)) +
    theme_void() +
    labs(fill = "Mutation type") +
    theme(legend.position = "right") +
    coord_equal()
  
  if (!is.null(facet_col)) {
    p <- p + facet_wrap(as.formula(paste0("~", facet_col)))
    return(p)
  } else {
    return(p)
  }
}