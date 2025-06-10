 #' Plot the trinucleotide spectrum
#' @description Creates barplots of the trinucleotide spectrum for all levels of
#' a given group.
#' @param mf_96 A data frame containing the mutation frequency data at the
#' 96-base resolution. This should be obtained using the 'calculate_mf' with
#' subtype_resolution set to 'base_96'. Generally, cols_to_group should be the
#' same as 'group_col'.
#' @param response A character string specifying the type of response to plot.
#' Must be one of 'frequency', 'proportion', or 'sum'.
#' @param mf_type A character string specifying the mutation count method to
#' plot. Must be one of 'min' or 'max'. Default is 'min'.
#' @param group_col A character string specifying the column(s) in 'mf_96'
#' to group the data by. Default is 'sample'. The sum, proportion, or frequency
#' will be plotted for all unique levels of this
#' group. You can specify more than one column to group by. Generally the same
#' as the 'cols_to_group' parameter in 'calculate_mf' when generating mf_96.
#' @param indiv_y A logical value specifying whether the the max response value
#' for the y-axis should be scaled independently for each group (TRUE) or scaled
#' the same for all groups (FALSE). Default is FALSE.
#' @param output_path An optional file path to an output directory. If provided,
#' the plots will be automatically exported using the graphics device
#' specified in output_type. The function will create the output directory if it
#' doesn't already exist. If NULL, plots will not be exported. Default is NULL.
#' @param output_type A character string specifying the type of output file.
#' Options are  'eps', 'ps', 'tex', 'pdf', or 'jpeg', 'tiff', 'png', 'bmp',
#' 'svg', or 'wmf' (windows only). Default is 'svg'.
#' @param sum_totals A logical value specifying whether to display the total
#' sum of mutations in the mutation labels. Default is TRUE.
#' @return A named list containing ggplots.
#' @importFrom dplyr arrange group_by mutate summarise
#' @import ggplot2
#' @importFrom stringr str_extract str_c
#' @details The function plots the trinucleotide spectrum for all levels of a
#' given group from the provided mf_96 data; the output of calculate_mf with
#' subtype_resolution = "base_96".
#' @examples
#' # Load example data
#' example_file <- system.file(
#'  "extdata", "Example_files",
#'  "example_mutation_data_filtered.rds",
#'  package = "MutSeqR"
#' )
#' example_data <- readRDS(example_file)
#'
#' # Calculate the mutation frequency data at the 96-base resolution
#' mf_96 <- calculate_mf(
#'  mutation_data = example_data,
#'  cols_to_group = "dose_group",
#'  subtype_resolution = "base_96",
#'  variant_types = "snv"
#' )
#' # Plot the trinucleotide proportions for each dose group
#' # Scale y-axis the same for all groups
#' plots <- plot_trinucleotide(
#'  mf_96 = mf_96,
#'  response = "proportion",
#'  mf_type = "min",
#'  group_col = "dose_group",
#'  indiv_y = FALSE,
#'  output_path = NULL
#' )
#' @export
plot_trinucleotide <- function(
  mf_96,
  response = "proportion",
  mf_type = "min",
  group_col = "dose",
  indiv_y = FALSE,
  sum_totals = TRUE,
  output_path = NULL,
  output_type = "svg"
) {
  mf_96 <- dplyr::filter(mf_96,
    !.data$normalized_context_with_mutation %in% setdiff(MutSeqR::subtype_list$type, "snv"))

  if (response == "proportion") {
    response_col <- paste0("proportion_", mf_type)
  } else if (response == "frequency") {
    response_col <- paste0("mf_", mf_type)
  } else if (response == "sum") {
    response_col <- paste0("sum_", mf_type)
  } else {
    stop("response must be one of 'frequency', 'proportion', or 'sum'")
  }

  if (length(group_col) > 1) {
    mf_96$group <- apply(mf_96[group_col], 1, paste, collapse = "_")
  } else {
    mf_96$group <- mf_96[[group_col]]
  }

  data <- mf_96 %>%
    dplyr::select(
      "group",
      "normalized_context_with_mutation",
      "normalized_context",
      all_of(paste0("sum_", mf_type)),
      all_of(response_col)
    ) %>%
    dplyr::rename(
      context = "normalized_context",
      subtype = "normalized_context_with_mutation",
      response = !!response_col,
      sum = paste0("sum_", mf_type)
    ) %>%
    dplyr::mutate(
      mutation = str_extract(.data$subtype, "(?<=\\[)[^\\]]+(?=\\])")
    ) %>%
    dplyr::arrange(.data$mutation, .data$context) %>%
    dplyr::mutate(
      subtype = factor(.data$subtype, levels = unique(.data$subtype)),
      mutation = factor(.data$mutation,
                        levels = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
    )

  if (response == "sum") {
    data$response <- data$sum
  }
  group_levels <- unique(data$group)
  plot_list <- setNames(vector("list", length(group_levels)), group_levels)

  # Predefine mutation colors
  plotcolours <- c(
    "C>A" = "#4DB6E9",
    "C>G" = "#000000",
    "C>T" = "#E74C43",
    "T>A" = "#CCCCCC",
    "T>C" = "#AAC96F",
    "T>G" = "#F7B6B5"
  )

  n_mut <- 6
  block_len <- 16

  for (i in seq_along(group_levels)) {
    # Subset the data into the group to be plotted
    plot_data <- dplyr::filter(data, .data$group == group_levels[i])

    # Sum totals and modify the labels as needed.
    mut_levels <- levels(plot_data$mutation)
    if (sum_totals) {
      mut_counts <- plot_data %>%
        dplyr::group_by(.data$mutation) %>%
        dplyr::summarise(nrmuts = sum(sum), .groups = "drop_last")
      # Create a character vector with the label for each mutation block
      labels_vec <- stringr::str_c(mut_counts$mutation, "\n(n = ", mut_counts$nrmuts, ")")
      labels <- setNames(labels_vec, mut_counts$mutation)
    } else {
      labels <- setNames(mut_levels, mut_levels)
    }

    # y-axis scaling
    if (indiv_y) {
      y_max <- max(plot_data$response, na.rm = TRUE)
    } else {
      y_max <- max(data$response, na.rm = TRUE)
    }
    if (response == "proportion") {
      y_lab <- "Proportion of Mutations"
      y_max <- ceiling(y_max * 10) / 10
    } else if (response == "frequency") {
      y_max_string <- format(y_max, scientific = TRUE)
      split_frequency <- strsplit(y_max_string, "e", fixed = TRUE)[[1]]
      coefficient <- as.numeric(split_frequency[1])
      rounded_coefficient <- ceiling(coefficient)
      round_str <- paste(rounded_coefficient, "e", split_frequency[2], sep = "")
      y_max <- as.numeric(round_str)
      y_lab <- "Frequency of Mutations"
    } else {
      y_lab <- "Sum of Mutations"
      y_max <- ceiling(y_max / 5) * 5
    }

    # assign x position to each subtype (for annotation)
    plot_data <- plot_data %>%
      dplyr::arrange(.data$mutation, .data$subtype) %>%
      dplyr::mutate(x_pos = as.numeric(subtype))
    # Change the x-labels to the Context for legibility.
    subtype_levels <- levels(plot_data$subtype)
    subtype_to_context <- setNames(plot_data$context, plot_data$subtype)[subtype_levels]

    # build rectangle/label dataframe
    gap <- 0.5 # gap between rectangles
    box_gap <- 0.01 * y_max
    box_height <- 0.02 * y_max
    rect_ymin <- y_max + box_gap
    rect_ymax <- rect_ymin + box_height
    text_y <- rect_ymax + 2 * box_gap
    rects <- data.frame(
      xmin = seq(0.5 + gap / 2, by = block_len, length.out = n_mut),
      xmax = seq(block_len + 0.5 - gap / 2, by = block_len, length.out = n_mut),
      mutation = mut_levels,
      label = labels[mut_levels],
      ymin = rect_ymin,
      ymax = rect_ymax
    )
    rects$fill <- plotcolours[rects$mutation]
    rects$xcenter <- (rects$xmin + rects$xmax) / 2

    n_bars <- length(levels(plot_data$subtype))

    # Make the ggplot
    p <- ggplot(plot_data, aes(x = subtype, y = response, fill = mutation)) +
      annotate("segment",
        x = 0.5, xend = 0.5,
        y = 0, yend = y_max,
        color = "gray80", linewidth = 0.6) +
      geom_rect(
        data = rects,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = mutation),
        inherit.aes = FALSE,
        show.legend = FALSE
      ) +
      annotate(
        "text",
        x = rects$xcenter, y = text_y,
        label = rects$label,
        color = "black",
        size = 4.5, fontface = 2,
        vjust = 0
      ) +
      coord_cartesian(ylim = c(0, y_max), clip = "off") +
      geom_col(width = 0.5, color = NA, show.legend = FALSE) +
      annotate("segment",
        x = 0.5, xend = n_bars + 0.5,
        y = 0, yend = 0,
        color = "gray80", linewidth = 0.6
      ) +
      scale_fill_manual(values = plotcolours) +
      scale_x_discrete(
        breaks = subtype_levels,
        labels = subtype_to_context,
        drop = FALSE,
        expand = c(0.002, 0.002)
      ) +
      labs(
        x = "Trinucleotide Context",
        y = y_lab,
        title = as.character(group_levels[i])
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 90,
          vjust = 1,
          hjust = 1,
          family = "mono", size = rel(0.75),
          margin = margin(t = -14) # decreases the gab between x-axis and the labels.
        ),
        axis.title.x = element_text(margin = margin(t = 5, b = 0)),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line("gray80"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(72, 24, 24, 48),
        legend.position = "none",
        panel.border = element_blank(),
        plot.background = element_blank(),
        plot.caption = element_text(hjust = 0)
      ) + ggtitle(as.character(group_levels[i])) +
      theme(
        plot.title.position = "plot",  # ← moves it above the plot panel
        plot.title = element_text(hjust = 0.5, margin = margin(b = 30))  # add bottom margin
      )
    plot_list[[i]] <- p

    # save the plots
    if (!is.null(output_path)) {
      output_dir <- file.path(output_path)
      output_filename <- paste0(
        "trinucleotide_plot_",
        group_col, "_", names(plot_list)[i],
        ".", output_type)
      ggsave(
        filename = output_filename,
        plot = p,
        device = output_type,
        path = output_path,
        create.dir = TRUE,
        width = 12,
        height = 6)
    }
  }
  return(plot_list)
}
