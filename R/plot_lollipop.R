#' Plot recurrent mutations in a lollipop plot using ggplot2
#'
#' @param mutation_data A data frame containing mutation data.
#' @param min_recurrence An integer specifying the minimum number of times a
#'   mutation must be observed at the same position to be plotted.
#' Default is 2.
#' @param group_col A character vector specifying the column name(s) to group
#'   mutation_data by. Default is "contig".
#' @param subtype_resolution The subtype resolution by which to group and
#'   colour the mutations. Options are "none", "type", "base_6", "base_12",
#'   "base_96", and "base_192".
#' @param custom_palette A named vector of colors to be used for the mutation
#' subtypes. The names of the vector should correspond to the mutation subtypes
#' in the data. Alternatively, you can specify a color palette from the
#' RColorBrewer package. See \code{\link[RColorBrewer]{brewer.pal}} for palette
#' options. You may visualize the palettes at the ColorBrewer website:
#' \url{https://colorbrewer2.org/}. Default is `NULL`.
#' @return A list of ggplot objects.
#' @examples
#' # For this example, we will use a subset of the example mutation data.
#' # The subset contains mutations from target chr1 in samples from the high
#' # dose group (50mg).
#'  example_data <- readRDS(system.file("extdata", "Example_files",
#'      "variants_subset_d50_chr1.rds", package = "MutSeqR"))
#'  # We will plot mutations that recoccur in at least two samples, grouped
#'  # by the "label" column, which signifies the target region (chr1).
#'  # Mutations will be grouped and coloured by their base 6 subtype (default) 
#'  plot <- plot_lollipop(
#'      mutation_data = example_data,
#'      min_recurrence = 2,
#'      group_col = "label"
#'  )
#' @importFrom dplyr %>% group_by tally filter arrange ungroup mutate select across all_of rename
#' @importFrom ggplot2 ggplot aes geom_segment geom_point scale_fill_manual
#' @importFrom ggplot2 scale_x_continuous labs theme_minimal theme element_blank
#' @importFrom ggplot2 element_line element_text
#' @importFrom rlang .data
#' @export
plot_lollipop <- function(mutation_data,
                          min_recurrence = 2,
                          group_col = "contig",
                          subtype_resolution = "base_6",
                          custom_palette = NULL) {
  stopifnot(
      "mutation_data is required." = !missing(mutation_data),
      "mutation_data must be a data.frame." = is.data.frame(mutation_data),
      "min_recurrence must be a single positive integer." =
          is.numeric(min_recurrence) && length(min_recurrence) == 1 &&
              min_recurrence >= 1 && (min_recurrence %% 1 == 0)
  )
  if ("filter_mut" %in% colnames(mutation_data)) {
    dplyr::filter(mutation_data, filter_mut==FALSE)
  }
  subtype_resolution <- match.arg(
      subtype_resolution,
      choices = names(MutSeqR::subtype_dict)
  )
  if (subtype_resolution == "none") {
    subtype_col <- NULL
  } else {
    subtype_col <- MutSeqR::subtype_dict[[subtype_resolution]]
  }

  required_cols <- c("contig", "start", "variation_type")
  if (!is.null(group_col)) required_cols <- c(required_cols, group_col)
  if (!is.null(subtype_col)) required_cols <- c(required_cols, subtype_col)

  missing_cols <- setdiff(required_cols, names(mutation_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  tally_groups <- c(group_col, "contig", "start", subtype_col)

  plot_data <- mutation_data %>%
    dplyr::filter(.data$variation_type != "no_variant") %>%
    # Group by dynamic columns
    dplyr::group_by(dplyr::across(dplyr::all_of(tally_groups))) %>%
    dplyr::tally(name = "n") %>%
    dplyr::filter(.data$n >= min_recurrence) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(-.data$n)

  if (nrow(plot_data) == 0) {
    warning("No mutations met the minimum recurrence threshold. Returning an empty list.")
    return(list())
  }

  if (subtype_resolution == "none") {
    plot_data$plotting_subtype <- "Recurrent Mutation" # Generic Label
    palette <- c("Recurrent Mutation" = "black") # Simple single color
  } else {
    plot_data <- dplyr::rename(plot_data, plotting_subtype = !!subtype_col)
    palette <- get_mutation_palette(
      custom_palette = custom_palette,
      subtype_resolution = subtype_resolution,
      num_colours = length(unique(plot_data$plotting_subtype))
    )
  }
  all_levels <- unique(c(names(palette), unique(as.character(plot_data$plotting_subtype))))

  plot_data$plotting_subtype <- factor(
    plot_data$plotting_subtype,
    levels = all_levels
  )

  if (is.null(group_col)) {
    plot_data$split_var <- "All Samples"
    plot_title_prefix <- "Recurrent Mutations"
  } else {
    if (length(group_col) > 1) {
      vals <- plot_data[, group_col, drop = FALSE]
      plot_data$split_var <- do.call(paste, c(vals, sep = " | "))
    } else {
      plot_data$split_var <- as.character(plot_data[[group_col]])
    }
    plot_title_prefix <- "Recurrent Mutations in"
  }
  data_list <- split(plot_data, plot_data$split_var)
  create_single_lollipop <- function(df, region_name) {
    final_title <- if (plot_title_prefix == "Recurrent Mutations") {
      plot_title_prefix
    } else {
      paste(plot_title_prefix, region_name)
    }

    ggplot2::ggplot(df, ggplot2::aes(x = .data$start, y = .data$n)) +
      ggplot2::geom_segment(
        ggplot2::aes(x = .data$start, xend = .data$start, y = 0, yend = .data$n),
        linewidth = 0.5,
        color = "grey50"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(fill = .data$plotting_subtype),
        shape = 21,
        size = 3.5,
        stroke = 0.5
      ) +
      ggplot2::scale_fill_manual(values = palette, name = "Mutation Subtype", drop = FALSE) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      ggplot2::labs(
        title = final_title,
        x = "Genomic Position",
        y = paste0("Recurrence (n >= ", min_recurrence, ")")
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(color = "black"),
        axis.ticks = ggplot2::element_line(color = "black"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
      )
  }

  plot_list <- Map(create_single_lollipop, data_list, names(data_list))

  return(plot_list)
}