#' Plot the trinucleotide spectrum
#'
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
#' @importFrom stats setNames
#' @import ggplot2
#' @importFrom stringr str_extract str_c
#' @details The function plots the trinucleotide spectrum for all levels of a
#' given group from the provided mf_96 data; the output of calculate_mf with
#' subtype_resolution = "base_96".
#' @examples
#'  # Calculate the mutation frequency data at the 96-base resolution
#'  mf_96 <- readRDS(system.file("extdata", "Example_files", "mf_data_96.rds",
#'     package = "MutSeqR"))
#' # Plot the trinucleotide proportions for the control and high dose groups
#' mf_96 <- dplyr::filter(mf_96, dose_group %in% c("Control", "High"))
#'
#' # Scale y-axis the same for all groups
#'   plots <- plot_trinucleotide(
#'     mf_96 = mf_96,
#'     response = "proportion",
#'     mf_type = "min",
#'     group_col = "dose_group",
#'     indiv_y = FALSE,
#'     output_path = NULL
#'   )
#' @export
plot_trinucleotide <- function(
    mf_96,
    response = "proportion",
    mf_type = "min",
    group_col = "dose",
    indiv_y = FALSE,
    sum_totals = TRUE,
    output_path = NULL,
    output_type = "svg") {

  # Validation
  stopifnot(
    !missing(mf_96) && is.data.frame(mf_96),
    is.logical(indiv_y),
    is.logical(sum_totals)
  )
  response <- match.arg(
      response,
      choices = c("proportion", "frequency", "sum")
  )
  mf_type <- match.arg(mf_type, choices = c("min", "max"))
  output_type <- match.arg(
      output_type,
      choices = c("eps", "ps", "tex", "pdf", "jpeg",
                  "tiff", "png", "bmp", "svg", "wmf")
  )

  # Data Prep

  # Filter out Non-SNVs
  non_snv <- MutSeqR::subtype_list$type[MutSeqR::subtype_list$type != "snv"]
  non_snv_df <- dplyr::filter(
    mf_96,
    .data$normalized_context_with_mutation %in% non_snv
    )
  mf_96 <- dplyr::filter(
      mf_96,
      !.data$normalized_context_with_mutation %in% non_snv
  )
  if (response == "proportion" && nrow(non_snv_df) > 0) {
    warning(
      "Non-SNV mutation subtypes detected and removed from data.", 
      " Plotted proportions will not include these subtypes, thus will not",
      " sum to 1. Suggest: re-calculate mf_96 with variant_types = 'snv'"
    )
  }


  # Define columns
  response_col <- switch(response,
                         "proportion" = paste0("proportion_", mf_type),
                         "frequency" = paste0("mf_", mf_type),
                         "sum" = paste0("sum_", mf_type))

  sum_col <- paste0("sum_", mf_type)

  # Create Group Column
  if (length(group_col) > 1) {
      mf_96$group <- do.call(paste, c(mf_96[group_col], sep = "_"))
  } else {
      mf_96$group <- mf_96[[group_col]]
  }
  
  # Check for data granularity mismatch
  # Check if there are multiple rows per group and subtype, which would
  # indicate that the data is not at the correct granularity for plotting
  data_check <- mf_96 %>%
      dplyr::group_by(.data$group, .data$normalized_context_with_mutation) %>%
      dplyr::tally() %>%
      dplyr::filter(.data$n > 1)
    if (nrow(data_check) > 0 && response != "sum") {
        # Extract a specific example to show the user
        example_row <- data_check$group[1]
        if (response == "proportion") {
            fix <- paste0("Fix: Re-run calculate_mf() with cols_to_group = c(",
                          paste(group_col, collapse = ", "), ")")
        } else if (response == "frequency") {
            fix <- paste0(
                "Fix: Input the mean mf value across ",
                paste(group_col, collapse = ","),
                " and normalized_context_with_mutation. This can be",
                " calculated by using dplyr::group_by() and",
                " dplyr::summarize() on your sample-level mf_96 data."
            )
        }
        stop(
            "\nData granularity mismatch detected. You are attempting to plot",
            " by '", paste(group_col, collapse = ", "), "', but the",
            " input data contains multiple rows per mutation subtype for this",
            " group.\n Example: Group '", example_row, "' has ",
            data_check$n[1], " entries for the same subtype.\n Likely Cause:",
            " You calculated MF by 'sample' (or another column) but are",
            " trying to plot by '", paste(group_col, collapse = ", "), "'.\n",
            fix
        )
    } else if (nrow(data_check) > 0 && response == "sum") {
        mf_96 <- mf_96 %>%
            dplyr::group_by(.data$group, .data$normalized_context_with_mutation, .data$normalized_context) %>%
            dplyr::summarise(
                sum_min = sum(.data$sum_min, na.rm = TRUE),
                sum_max = sum(.data$sum_max, na.rm = TRUE),
                .groups = "drop"
            )
    }

  # Clean Data
  data <- mf_96 %>%
    dplyr::select(
      "group",
      subtype = "normalized_context_with_mutation",
      context = "normalized_context",
      sum = dplyr::all_of(sum_col),
      response_val = dplyr::all_of(response_col)
    ) %>%
    dplyr::mutate(
      # Vectorized Regex
      mutation = stringr::str_extract(.data$subtype, "(?<=\\[)[^\\]]+(?=\\])")
    ) %>%
    # Factorize globally
    dplyr::mutate(
      mutation = factor(.data$mutation, levels = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
    ) %>%
    # Global Arrange ensures consistent factor levels across all plots
    dplyr::arrange(.data$mutation, .data$context) %>%
    dplyr::mutate(
      subtype = factor(.data$subtype, levels = unique(.data$subtype))
    )

  # Override response for 'sum' mode
  if (response == "sum") data$response_val <- data$sum

  # -Determine Global Y-Axis Max (If not indiv_y)
  if (!indiv_y) {
    global_y_max <- max(data$response_val, na.rm = TRUE)
  }

  # Define Constants
  plotcolours <- c(
    "C>A" = "#4DB6E9", "C>G" = "#000000", "C>T" = "#E74C43",
    "T>A" = "#CCCCCC", "T>C" = "#AAC96F", "T>G" = "#F7B6B5"
  )

  n_mut <- 6
  block_len <- 16
  mut_levels <- levels(data$mutation) # Should be the standard 6

  # Split data by group for mapping
  # Splitting is efficient here as we need distinct plots
  data_list <- split(data, data$group)

  # Plotting Function (Applied via Map)
  create_plot <- function(group_name, plot_data) {

    # A. Calculate Y-Max locally
    local_max <- max(plot_data$response_val, na.rm = TRUE)
    target_max <- if (indiv_y) local_max else global_y_max

    # Y-Axis Formatting Logic
    if (response == "proportion") {
      y_lab <- "Proportion of Mutations"
      y_limit <- ceiling(target_max * 10) / 10
    } else if (response == "frequency") {
      y_lab <- "Frequency of Mutations"
      # Scientific notation rounding logic
      sci <- format(target_max, scientific = TRUE)
      parts <- strsplit(sci, "e")[[1]]
      y_limit <- as.numeric(paste0(ceiling(as.numeric(parts[1])), "e", parts[2]))
    } else {
      y_lab <- "Sum of Mutations"
      y_limit <- ceiling(target_max / 5) * 5
    }

    # B. Labels
    if (sum_totals) {
      # Vectorized summary
      mut_counts <- plot_data %>%
        dplyr::group_by(.data$mutation) %>%
        dplyr::summarise(n = sum(.data$sum), .groups = "drop")
      lbls <- setNames(paste0(mut_counts$mutation, "\n(n = ", mut_counts$n, ")"), mut_counts$mutation)
    } else {
      lbls <- setNames(mut_levels, mut_levels)
    }

    # C. Rectangles (Header Strip)
    # Calculate dimensions based on dynamic Y-limit
    gap <- 0.5
    box_gap <- 0.01 * y_limit
    box_height <- 0.02 * y_limit
    rect_ymin <- y_limit + box_gap
    rect_ymax <- rect_ymin + box_height
    text_y <- rect_ymax + 2 * box_gap

    rects <- data.frame(
      xmin = seq(0.5 + gap / 2, by = block_len, length.out = n_mut),
      xmax = seq(block_len + 0.5 - gap / 2, by = block_len, length.out = n_mut),
      mutation = mut_levels,
      ymin = rect_ymin,
      ymax = rect_ymax
    )
    rects$label <- lbls[as.character(rects$mutation)]
    rects$xcenter <- (rects$xmin + rects$xmax) / 2

    # D. Labels Map
    # Extract Context labels from the factors we set up globally
    # Just need one row per subtype to get the map
    subtype_map <- plot_data[match(levels(plot_data$subtype), plot_data$subtype), ]
    x_labels_vec <- setNames(subtype_map$context, subtype_map$subtype)

    # E. GGPlot
    p <- ggplot(
          plot_data,
          aes(x = .data$subtype, y = .data$response_val, fill = .data$mutation)
        ) +
      # Vertical guide lines
      annotate(
        "segment",
        x = 0.5,
        xend = 0.5,
        y = 0,
        yend = y_limit,
        color = "gray80",
        linewidth = 0.6
    ) +
      annotate(
        "segment",
        x = 0.5,
        xend = length(levels(plot_data$subtype)) + 0.5,
        y = 0,
        yend = 0,
        color = "gray80",
        linewidth = 0.6
    ) +
      # Data Bars
      geom_col(width = 0.5, color = NA, show.legend = FALSE) +
      # Header Rectangles
      geom_rect(
        data = rects,
        aes(xmin = xmin, xmax = xmax,
            ymin = ymin, ymax = ymax,
            fill = mutation),
        inherit.aes = FALSE,
        show.legend = FALSE
    ) +
      # Header Text
      annotate(
        "text",
        x = rects$xcenter,
        y = text_y,
        label = rects$label,
        color = "black",
        size = 4.5,
        fontface = 2,
        vjust = 0
    ) +
      # Scales
      scale_fill_manual(values = plotcolours) +
      scale_x_discrete(
        labels = x_labels_vec,
        drop = FALSE,
        expand = c(0.002, 0.002)
    ) +
      coord_cartesian(ylim = c(0, y_limit), clip = "off") +
      # Labels & Theme
      labs(x = "Trinucleotide Context", y = y_lab, title = group_name) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(
            angle = 90,
            vjust = 1,
            hjust = 1,
            family = "mono",
            size = rel(0.75),
            margin = margin(t = -14)
        ),
        axis.title.x = element_text(margin = margin(t = 5, b = 0)),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line("gray80"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 30)),
        plot.title.position = "plot",
        plot.margin = margin(72, 24, 24, 48),
        panel.border = element_blank()
      )

    # Output Saving
    if (!is.null(output_path)) {
      fname <- file.path(
          output_path,
          paste0("trinucleotide_plot_", paste(group_col, collapse = "_"),
                 "_", group_name, ".", output_type)
      )
      ggsave(
          filename = fname,
          plot = p,
          device = output_type,
          width = 12,
          height = 6
      )
    }

    return(p)
  }

  # Execution 
  plot_list <- Map(create_plot, names(data_list), data_list)

  return(plot_list)
}