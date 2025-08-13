#' Create a heatmap plot of mutation subtype proportions.
#'
#' @description This function creates a heatmap plot of subtype proportions for
#' a given grouping variable. The groups may be facetted by a second variable.
#' Mutation sums for each facet group and normalized subtype are calculated
#' and displayed.
#' @param mf_data A data frame containing the mutation frequency data at the
#' desired base resolution. This is obtained using the 'calculate_mf' with
#' subtype_resolution set to the desired resolution. cols_to_group
#' should be the same as 'group_col'.
#' @param group_col The variable to group by.
#' @param facet_col The variable to facet by.
#' @param mf_type The type of mutation frequency to plot. Options are "min" or
#' "max". (Default: "min")
#' @param mut_proportion_scale The scale option for the mutation proportion.
#' Options are passed to viridis::scale_fill_viridis_c.
#'  One of # inferno, magma, plasma, viridis, cividis, turbo, mako, or rocket.
#' We highly reccomend the default for its ability to disciminate hard to see
#' patterns. (Default: "turbo")
#' @param max Maximum value used for plotting the proportions.
#' Proportions that are higher will have the maximum colour. (Default: 0.2)
#' @param rescale_data Logical value indicating whether to rescale the mutation
#' proportions to increase the dynamic range of colors shown on the plot.
#' (Default: TRUE)
#' @param condensed More condensed plotting format. Default = FALSE.
#' @import ggplot2
#' @importFrom dplyr group_by summarise mutate rename all_of
#' @importFrom stringr str_length str_extract str_c
#' @importFrom magrittr %>%
#' @return A ggplot object representing the heatmap plot.
#' @export
#' @examples
#' if (requireNamespace("MutSeqRData", quietly = TRUE)) {
#' # Plot the trinucleotide proportions per sample, facetted by dose group.
#' 
#' # Example data consists of 24 mouse bone marrow DNA samples imported
#' # using import_mut_data() and filtered with filter_mut as in Example 4.
#' # Sequenced on TS Mouse Mutagenesis Panel. Example data is
#' # retrieved from MutSeqRData, an ExperimentHub data package.
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' example_data <- eh[["EH9861"]]
#' 
#' # define dose_group order
#' example_data$dose_group <- factor(example_data$dose_group,
#'                                   levels = c("Control", "Low",
#'                                              "Medium", "High"))
#' mf_96 <- calculate_mf(example_data,
#'                       cols_to_group = "sample",
#'                       variant_types = "snv",
#'                       subtype_resolution = "base_96",
#'                       retain_metadata_cols = "dose_group")
#'plot <- plot_trinucleotide_heatmap(mf_96,
#'                                   group_col = "sample",
#'                                   facet_col = "dose_group")
#'
#' }

plot_trinucleotide_heatmap <- function(mf_data,
                                       group_col = "sample",
                                       facet_col = "dose",
                                       mf_type = "min",
                                       mut_proportion_scale = "turbo",
                                       max = 0.2,
                                       rescale_data = FALSE,
                                       condensed = FALSE) {
  # Remove NA values from data
  mf_data <- mf_data %>%
    dplyr::filter(!is.na(.data[[group_col]]), !is.na(.data[[facet_col]]))

  # Make group_col a factor
  mf_data[[group_col]] <- as.factor(mf_data[[group_col]])

  if (rescale_data) {
    if (!requireNamespace("scales", quietly = TRUE)) {
      stop("Package scales is required when using the rescale_data option. Please install the package using 'install.packages('scales')'")
    }
  }

  # Check for proportion column
  proportion_col <- c(paste0("proportion_", mf_type), "proportion", "prop")
  found_prop_col <- proportion_col[proportion_col %in% tolower(colnames(mf_data))] 
  if(length(found_prop_col) == 1) {
    mf_data <- dplyr::rename(mf_data, proportion = dplyr::all_of(found_prop_col))
  } else if (length(found_prop_col) > 1) {
    stop(paste("More than one possible proportion column name found in mf_data: ",
               paste(found_prop_col, collapse = ", "),
              " Please remove columns that are not the proportion column to be plotted."))
  } else if (length(found_prop_col) == 0) {
  stop(paste0("The dataframe does not contain a proportion column. Please add a column that contains the mutation proportion to be plotted or rename column to 'proportion'."))
  }

  # Check for mutation count column
  sum_column <- paste0("sum_", mf_type)
  mf_data <- dplyr::rename(mf_data, sum_column = dplyr::all_of(sum_column))

  # Check for subtype column
  subtype_column_names <- c("normalized_context_with_mutation",
                            "context_with_mutation",
                            "normalized_subtype",
                            "subtype",
                            "variation_type")
  found_subtype_cols <- subtype_column_names[subtype_column_names %in% colnames(mf_data)] 
  if (length(found_subtype_cols) == 1) {
    # Rename the subtype column to "subtype"
    mf_data <- dplyr::rename(mf_data, subtype = dplyr::all_of(found_subtype_cols))
  } else if (length(found_subtype_cols) > 1) {
    stop(paste("More than one possible subtype column name found in mf_data: ",
               paste(found_subtype_cols, collapse = ", "),
              " Please remove columns that are not the subtype column to be plotted."))
  } else if (length(found_subtype_cols) == 0) {
    stop(paste("No subtype column name found in mf_data from the following options: ", 
               paste(subtype_column_names, collapse = ", "),
                " Please add a column that contains the mutation subtype to be plotted or rename column to one of the listed options.")) 
  }

  # Context Column synonyms
  context_column_names <- c("normalized_ref",
                            "short_ref",
                            "normalized_context",
                            "context")
  found_context_cols <- context_column_names[context_column_names %in% colnames(mf_data)]
  if (length(found_context_cols) == 1) {
    mf_data$x_variable <- mf_data[[found_context_cols]]
    mf_data <- mf_data %>%
      dplyr::mutate(x_variable = ifelse(subtype %in% MutSeqR::subtype_list$type,
                                        subtype,
                                        x_variable))
    plot_context <- TRUE
  } else if (length(found_context_cols) == 0) {
    if (any(grep(".*\\[([A-Z]?>[A-Z])\\].*", mf_data$subtype))) {
      mf_data <- mf_data %>%
        dplyr::mutate(context = paste0(substr(mf_data$subtype, 1, 1),
                                       substr(mf_data$subtype, 3, 3),
                                       substr(mf_data$subtype, 7, 7)))
      mf_data$x_variable <- mf_data$context
      mf_data <- mf_data %>%
        dplyr::mutate(x_variable = ifelse(subtype %in% MutSeqR::subtype_list$type, 
                                          subtype,
                                          x_variable))
      plot_context <- TRUE
    } else {
      print("No context column found in mf_data, plotting by subtype.")
      mf_data$x_variable <- mf_data$subtype
      plot_context <- FALSE
    }
  } else if (length(found_context_cols) > 1) {
    stop(paste("More than one possible context column name found in mf_data: ",
               paste(found_context_cols, collapse = ", "),
               " Please remove columns that are not the context column to be plotted.")) 
  }

  context_size <- max(stringr::str_length(mf_data$x_variable))
  if (context_size == 1) {
    axis_size <- 10
  } else if (context_size == 2) {
    axis_size <- 8
  } else if (context_size == 3) {
    axis_size <- 4
  } else {
    axis_size <- 4
  }

  # Change plotting parameters based on whether plot should be condensed.
  if (condensed == TRUE) {
    spacing <- 0
  } else {
    spacing <- 0.5
  }

  # Facet x
  if(plot_context) {
    # Make facet labels for subtypes.
    pattern <- "[A-Z]?>[A-Z]"
    mf_data$subtype_labels <- stringr::str_extract(mf_data$subtype, pattern)
    mf_data <- mf_data %>%
      dplyr::mutate(subtype_labels = ifelse(mf_data$subtype %in% subtype_list$type,
                                            "other",
                                            subtype_labels))
    # Count number muts per subtype
    if ("sum_column" %in% colnames(mf_data)) {
      mut_counts <- mf_data %>%
        dplyr::group_by(subtype_labels) %>%
        dplyr::summarise(nrmuts = sum(sum_column), .groups = "drop_last")
        facet_labs_x <- stringr::str_c(mut_counts$subtype_labels, " (n = ", mut_counts$nrmuts, ")")
        names(facet_labs_x) <- mut_counts$subtype_labels  
    } else {
      facet_labs_x <- unique(mf_data$subtype_labels)
      names(facet_labs_x) <- unique(mut_counts$subtype_labels)
    }
    facet_x_order <- c(MutSeqR::subtype_list$base_12, "other")
    mf_data <- mf_data %>%
      dplyr::mutate(subtype_labels = factor(subtype_labels, levels = facet_x_order))
    x_label <- "context"
  } else {
    x_label <- "Subtype"
    x_order <- c(MutSeqR::subtype_list$base_192, MutSeqR::subtype_list$base_12, MutSeqR::subtype_list$type)
    mf_data <- mf_data %>%
      dplyr::mutate(x_variable = factor(x_variable, levels = x_order))
  }

  # Facet y
  if (!is.null(facet_col)) {
    if ("sum_column" %in% colnames(mf_data)) {
      # Count number muts per sample_group
      mut_counts_groups <- mf_data %>%
        dplyr::group_by(!!ensym(facet_col)) %>%
        dplyr::summarise(nrmuts = sum(sum_column), .groups = "drop_last")
      facet_labs_y <- stringr::str_c(mut_counts_groups[[facet_col]], " (n = ", mut_counts_groups$nrmuts, ")")
      names(facet_labs_y) <- mut_counts_groups[[facet_col]]
    } else {
      facet_labs_y <- unique(mf_data[[facet_col]])
      names(facet_labs_y) <- unique(mf_data[[facet_col]])
  }
    y_label <- paste(facet_col)
    mf_data <- mf_data %>%
      dplyr::rename(Facet = !!ensym(facet_col))
    mf_data$Facet <- as.factor(mf_data$Facet)
  } else {
    y_label <- "Sample"
  }
  # If user specifies a mutation proportion max, then if value is higher than max,
  # change it to max (i.e., cut off the values at max)
  if (max < 1 && !rescale_data) {
    message(paste0("Cutting off at maximum mutation proportion value of ", max))
    df <- mf_data %>%
      dplyr::mutate(ProportionPlot = ifelse(proportion > max, max, proportion))
  } else if (max == 1 && rescale_data) {
  # If user specifies scaling to max value (the default), rescale the values to 0-1
    message(paste0("Rescaling bewteen 0 and 1"))
    df <- mf_data %>%
      dplyr::mutate(ProportionPlot = scales::rescale(proportion, to = c(0, 1)))
  } else if (rescale_data && max < 1) {
    # If user specifies both scaling and cutting off at max, then do both
    message(paste0("Rescaling and cutting off at maximum mutation proportion value of ", max))
    df <- mf_data %>%
      dplyr::mutate(ProportionPlot = scales::rescale(proportion, to = c(0, 1))) %>%
      dplyr::mutate(ProportionPlot = ifelse(proportion > max, max, proportion))
  } else {
    message(paste0("No scaling or maximum mutation proportion value applied"))
    df <- mf_data
    df$ProportionPlot <- mf_data$proportion
  }
  df <- dplyr::rename(df, Group = dplyr::all_of(group_col))
  # General figure, no facetting
  fig <- ggplot2::ggplot(df, aes(x = x_variable,
                        y = Group,
                        fill = ProportionPlot)) +
                ggplot2::geom_raster() +
                ggplot2::scale_fill_viridis_c(
                  name = "Relative proportion", limits = c(0, max),
                  option = mut_proportion_scale,
                  na.value = "white") +
                ggplot2::theme_minimal() +
                ggplot2::labs(x = x_label, y = y_label) +
                ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6),
                      axis.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5, size = axis_size, family = "mono"),
                      panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                      panel.spacing.x = unit(spacing, "lines"),
                      panel.spacing.y = unit(spacing, "lines")
  )
  # facet x : plot_context == TRUE
  # facet y : !is.null(facet_col)
  # facet x and y : plot_context == TRUE & !is.null(facet_col)
  # facet none : plot_context == FALSE & is.null(facet_col)
  if (plot_context && is.null(facet_col)) {
    figfx <- fig +
      ggplot2::facet_grid(cols = vars(subtype_labels), scales = "free_x",
                 labeller = labeller(subtype_labels = facet_labs_x)) +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
    return(figfx)
  } else if (plot_context == FALSE & !is.null(facet_col)) {
    figfy <- fig +
      ggplot2::facet_grid(rows = vars(Facet), scales = "free_y",
                 labeller = labeller(Facet = facet_labs_y)) +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
    return(figfy)
  } else if (plot_context & !is.null(facet_col)) {
    figfxy <- fig +
      ggplot2::facet_grid(Facet ~ subtype_labels, scales = "free",
                 labeller = labeller(Facet = facet_labs_y, subtype_labels = facet_labs_x)) +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
    return(figfxy)
  } else {
    return(fig)
  }

}
