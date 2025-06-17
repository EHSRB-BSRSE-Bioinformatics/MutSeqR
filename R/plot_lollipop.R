#' Plot recurrent mutations in a lollipop plot using ggplot2
#'
#' This function visualizes recurrent mutations from a data frame. It first
#' calculates the frequency of each mutation at specific genomic positions and
#' then generates a lollipop plot for each group (e.g., chromosome)
#' displaying mutations that meet a minimum recurrence threshold.
#'
#' @param mutations A data frame containing mutation data. It must contain columns
#'   for genomic start position (`start`), `variation_type`, `normalized_subtype`,
#'   and a column to group by (see `group_by_col`).
#' @param min_recurrence An integer specifying the minimum number of times a
#'   mutation must be observed at the same position to be plotted. Defaults to 2.
#' @param group_by_col A string specifying the column name to group mutations by,
#'   typically representing chromosomes or contigs (e.g., "seqnames", "chr").
#'   Defaults to "seqnames".
#' @param custom_palette A named character vector for coloring the mutation
#'   subtypes. The names should match the levels in `normalized_subtype`. If NULL
#'   (default), a default palette is used.
#'
#' @return A list of ggplot objects. Each element of the list is a lollipop
#'   plot for a specific region (e.g., a chromosome) and is named accordingly.
#'
#' @importFrom dplyr %>% group_by tally filter arrange
#' @importFrom ggplot2 ggplot aes geom_segment geom_point scale_fill_manual
#' @importFrom ggplot2 scale_x_continuous labs theme_minimal theme element_blank
#' @importFrom ggplot2 element_line element_text
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' if (requireNamespace("dplyr", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'
#' example_file <- system.file("extdata", "Example_files",
#'                             "example_mutation_data_filtered.rds",
#'                             package = "MutSeqR")
#' example_data <- readRDS(example_file)
#' example_data$dose_group <- factor(example_data$dose_group,
#'                                   levels = c("Control", "Low",
#'                                              "Medium", "High"))
#'
#'   # 2. Generate the plots
#'   plot_list <- plot_lollipop(mutations = mf, min_recurrence = 2)
#'
#'   # 3. Display a plot for a specific chromosome
#'   # print(plot_list$chr1)
#'   # print(plot_list$chr2)
#' }
plot_lollipop <- function(mutations,
                            min_recurrence = 2,
                            group_by_col = "dose_group",
                            custom_palette = NULL) {

  # --- 1. Input Validation ---
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE)) {
    stop("This function requires `ggplot2` and `dplyr`. Please install them.")
  }
  if (!is.data.frame(mutations)) {
    stop("Input `mutations` must be a data.frame.")
  }
  # Check for all required columns
  required_cols <- c("start", "variation_type", "normalized_subtype", group_by_col)
  missing_cols <- setdiff(required_cols, names(mutations))
  if (length(missing_cols) > 0) {
     stop(paste("The `mutations` data.frame is missing required columns:",
                paste(missing_cols, collapse = ", ")))
  }

  # --- 2. Data Preparation ---
  # No conversion needed; we use the input data frame directly.
  plot_data <- mutations %>%
    # Filter out non-variants if that column exists and is used
    dplyr::filter(.data$variation_type != "no_variant") %>%
    # Count occurrences by region, position, and subtype
    dplyr::group_by(.data[[group_by_col]], .data$start, .data$normalized_subtype) %>%
    dplyr::tally(name = "n") %>%
    # Keep only mutations meeting the recurrence threshold
    dplyr::filter(.data$n >= min_recurrence) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(-.data$n)

  if (nrow(plot_data) == 0) {
    warning("No mutations met the minimum recurrence threshold. Returning an empty list.")
    return(list())
  }

  # --- 3. Palette and Factor Setup ---
  if (is.null(custom_palette)) {
    # Define the default palette from your colleague's code
    palette <- c(
      RColorBrewer::brewer.pal(5, "BrBG"),
      RColorBrewer::brewer.pal(6, "Spectral")
    )
    names(palette) <- c(
      "complex", "deletion", "insertion", "mnv", "sv",
      "T>G", "T>C", "T>A", "C>T", "C>G", "C>A"
    )
  } else {
    palette <- custom_palette
  }

  # Ensure the subtype column is a factor to control legend order and colors
  all_lvls <- unique(c(names(palette), unique(plot_data$normalized_subtype)))
  plot_data$normalized_subtype <- factor(
    plot_data$normalized_subtype,
    levels = all_lvls
  )

  # --- 4. Plotting Loop ---
  plot_list <- list()
  unique_regions <- unique(as.character(plot_data[[group_by_col]]))

  for (region_name in unique_regions) {
    # Subset data for the current region
    df_region <- dplyr::filter(plot_data, .data[[group_by_col]] == region_name)

    # Create the ggplot object
    p <- ggplot(df_region, aes(x = .data$start, y = .data$n)) +
      geom_segment(
        aes(x = .data$start, xend = .data$start, y = 0, yend = .data$n),
        linewidth = 0.5,
        color = "grey50"
      ) +
      geom_point(
        aes(fill = .data$normalized_subtype),
        shape = 21,
        size = 3.5,
        stroke = 0.5
      ) +
      scale_fill_manual(values = palette, name = "Mutation Subtype", drop = FALSE) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      labs(
        title = paste("Recurrent Mutations on", region_name),
        x = "Genomic Position",
        y = paste0("Recurrence (n >= ", min_recurrence, ")")
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
      )

    plot_list[[region_name]] <- p
  }

  return(plot_list)
}