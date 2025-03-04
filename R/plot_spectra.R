#' Transition-transversion plot
#' @description Given mf data, construct a plot displaying the
#' mutation subtypes observed in a cohort.
#' @param mf_data A data frame containing the mutation frequency data at the
#' desired subtype resolution. This is obtained using the 'calculate_mf'
#' function with subtype_resolution set to the desired resolution.
#' Data must include a column containing the group_col,
#' a column containing the mutation subtypes, a column containing the desired
#' response variable (mf, proportion, sum) for the desired mf_type
#' (min or max), and if applicable, a column containing the variable by which
#' to order the samples/groups.
#' @param group_col The name of the column(s) in the mf data that contains the
#' sample/group names. This will generally be the same values used for the
#' cols_to_group argument in the calculate_mf function. However, you may
#' also use groups that are at a higher level of the aggregation in mf_data.
#' @param subtype_resolution The subtype resolution of the mf data.
#' Options are `base_6`, `base_12`, `base_96`, `base_192`, or `type`.
#' Default is `base_6`.
#' @param response The desired response variable to be plotted. Options are
#' mf, proportion, or sum. Default is `proportion`. Your mf_data must contain
#' columns with the name of your desired response: `mf_min`, `mf_max`,
#' `proportion_min`, `proportion_max`, `sum_min`, and `sum_max`.
#' @param mf_type The mutation counting method to use. Options are min or max.
#' Default is `min`.
#' @param group_order The method for ordering the samples within the plot.
#' Options include:
#' \itemize{
#'   \item `none`: No ordering is performed. Default.
#'   \item `smart`: Groups are automatically ordered based on the group names
#' (alphabetical, numerical)
#'   \item `arranged`: Groups are ordered based on one or more factor column(s)
#' in mf_data. Column names are passed to the function using the
#' `group_order_input`.
#'  \item `custom`: Groups are ordered based on a custom vector of group
#' names. The custom vector is passed to the function using the
#' `group_order_input`.
#' \item `clustered`: Groups are ordered based on hierarchical clustering. The
#' dissimilarity matrix can be specified using the `dist` argument. The
#' agglomeration method can be specified using the `cluster_method` argument.
#' }
#' @param group_order_input A character vector specifying details for the
#' group order method. If `group_order` is `arranged`, `group_order_input`
#' should contain the column name(s) to be used for ordering the samples. If
#' `group_order` is `custom`, `group_order_input` should contain the custom
#' vector of group names.
#' @param dist  The dissimilarity matrix for hierarchical clustering. Options
#' are `cosine`, `euclidean`, `maximum`, `manhattan`, `canberra`, `binary` or
#' `minkowski`. The default is `cosine`. See \link[stats]{dist} for details.
#' @param cluster_method The agglomeration method for hierarchical clustering.
#' Options are `ward.D`, `ward.D2`, `single`, `complete`, `average` (= UPGMA),
#' `mcquitty` (= WPGMA), `median` (= WPGMC) or `centroid` (= UPGMC). The default
#' is `Ward.D`. See \link[stats]{hclust} for details.
#' @param palette A named vector of colors to be used for the mutation subtypes.
#' The names of the vector should correspond to the mutation subtypes in the
#' data. The default is a set of colors from the RColorBrewer package.
#' @param x_lab The label for the x-axis. Default is the value of `group_col`.
#' @param y_lab The label for the y-axis. Default is the value of `response_col`.
#' @import patchwork
#' @import ggplot2
#' @importFrom dplyr select arrange across all_of
#' @export
#' @examples 
#' # Load example data
#' example_file <- system.file("extdata", "example_mutation_data_filtered.rds", package = "MutSeqR")
#' example_data <- readRDS(example_file)
#'
#' # Example 1: plot the proportion of 6-based mutation subtypes for each sample, organized by dose group:
#'
#' # Calculate the mutation frequency data at the 6-base resolution. Retain the dose_group column to use for ordering the samples.
#' mf_data <- calculate_mf(mutation_data = example_data,
#'                         cols_to_group = "sample",
#'                         subtype_resolution = "base_6",
#'                         retain_metadata_cols = "dose_group")
#' # Set the desired order for the dose_group levels.
#' mf_data$dose_group <- factor(example_data$dose_group, levels = c("Control", "Low", "Medium", "High"))
#' # Plot the mutation spectra
#' plot <- plot_spectra(mf_data = mf_data,
#'                      group_col = "sample",
#'                      subtype_resolution = "base_6",
#'                      response = "proportion",
#'                      group_order = "arranged",
#'                      group_order_input = "dose_group")
#' 
#' # Example 2: plot the proportion of 6-based mutation subtypes for each sample, ordered by hierarchical clustering:
#' plot <- plot_spectra(mf_data = mf_data,
#'                      group_col = "sample",
#'                      subtype_resolution = "base_6",
#'                      response = "proportion",
#'                      group_order = "clustered")

plot_spectra <- function(mf_data,
                         group_col = "sample",
                         subtype_resolution = "base_6",
                         response = "proportion",
                         mf_type = "min",
                         group_order = "none",
                         group_order_input = NULL,
                         dist = "cosine",
                         cluster_method = "ward.D",
                         palette = NULL,
                         x_lab = NULL,
                         y_lab = NULL) {
  
 # check package dependencies
  if (group_order == "clustered") {
      if (!requireNamespace("ggh4x", quietly = TRUE)) {
      stop("Package ggh4x is required when using the 'clustered' group_order option. Please install the package using 'install.packages('ggh4x')'")
    }
  } else if (group_order == "smart") {
    if (!requireNamespace("gtools", quietly = TRUE)) {
      stop("Package gtools is required when using the 'smart' group_order option. Please install the package using 'install.packages('gtools')'")
    }
  }
  # Desginate the response column
  if (response == "proportion") {
    response_col <- paste0("proportion_", mf_type)
  } else if (response == "mf") {
    response_col <- paste0("mf_", mf_type)
  } else if (response == "sum") {
    response_col <- paste0("sum_", mf_type)
  } else {
    stop("response must be one of 'mf', 'proportion', or 'sum'")
  }
  # concatenate multiple group columns into one
  if (length(group_col) > 1) {
    mf_data$group <- apply(mf_data[group_col], 1, paste, collapse = "_")
  } else {
    mf_data$group <- mf_data[[group_col]]
  }
  # Select and rename required columns
  plot_data <- mf_data %>%
    dplyr::select("group",
                  dplyr::all_of(MutSeqR::subtype_dict[[subtype_resolution]]),
                  dplyr::all_of(response_col)) %>%
    dplyr::rename(
      subtype = dplyr::all_of(MutSeqR::subtype_dict[[subtype_resolution]]),
      response = dplyr::all_of(response_col)
    )

  if (group_order == "none") {
    order <- as.vector(unique(plot_data$group))
    plot_data$group <- factor(plot_data$group, levels = order)
  } else if (group_order == "smart") {
    order <- as.vector(unique(plot_data$group))
    order <- gtools::mixedsort(order)
    plot_data$group <- factor(plot_data$group, levels = order)
  } else if (group_order == "arranged") {
    plot_data <- cbind(plot_data, mf_data[group_order_input])
    plot_data <- plot_data %>%
      dplyr::arrange(!!!rlang::syms(group_order_input))
    order <- as.vector(unique(plot_data$group))
    plot_data$group <- factor(plot_data$group, levels = order)
  } else if (group_order == "custom") {
    plot_data$group <- factor(plot_data$group,
                              levels = group_order_input)
  } else if (group_order == "clustered") {
    # Cluster the samples
    hc <- cluster_spectra(mf_data = plot_data,
                          group_col = "group",
                          response_col = "response",
                          subtype_col = "subtype",
                          dist = dist,
                          cluster_method = cluster_method)
    # Reorder the samples based on hierarchical clustering
    order <- hc$labels[hc$order]
    # Reorder the levels of the sample variable in your data frame
    plot_data$group <- factor(plot_data$group,
                              levels = order)
  }
  if (subtype_resolution != "type") {
  subtype_order <- c(MutSeqR::subtype_list$type, rev(MutSeqR::subtype_list[[subtype_resolution]]))
  } else {
    subtype_order <- c(MutSeqR::subtype_list$type)
  }
  plot_data$subtype <- factor(plot_data$subtype,
                              levels = subtype_order)

  if (is.null(palette)) {
    if(subtype_resolution == "base_6") {
      # Standard palette for simple spectra
    palette <- c("C>A" = "#3288BD",
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
    } else if (subtype_resolution == "base_12") {
      # Sanger colours for 12 base spectra
      palette <- c("A>C" = "limegreen",
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
                    "symbolic" = "azure2",
                    "ambiguous" = "darkgrey",
                    "uncategorized" = "white")
    } else if (subtype_resolution == "base_96") {
      base_colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "gray", "olivedrab1", "cyan", "magenta")
      palette <- colorRampPalette(base_colors)(101)
    } else if (subtype_resolution == "base_192") {
      base_colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "gray", "olivedrab1", "cyan", "magenta")
      palette <- colorRampPalette(base_colors)(297)
    } else if(subtype_resolution == "type") {
      palette <- c("mnv" = "pink",
                   "deletion" = "black",
                   "insertion" = "grey",
                   "symbolic" = "purple",
                   "snv" = "blue")
    }
  }
  # Axis labels
  if (is.null(x_lab)) {
    x_lab <- group_col
  }
  if (is.null(y_lab)) {
    y_lab <- response_col
  }
  axis_labels <- ggplot2::labs(x = x_lab, y = y_lab)

  # bar plot
  bar <- ggplot(plot_data, aes(x = .data$group,
                               y = .data$response,
                               fill = .data$subtype,
                               add = FALSE)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = palette) +
    axis_labels +
    theme_minimal() +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 90)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01)))

  if (group_order == "clustered") {
    plot <-  bar +
      ggh4x::scale_x_dendrogram(hclust = hc, position = "top", labels = NULL,
      ) +
      theme(axis.ticks.length.x = unit(10, "pt"))
    x_axis <- ggplot(plot_data, aes(x = .data$group)) +
      theme_minimal() +
      labs(x = x_lab) +
      theme(axis.text.x = element_text(angle = 90),
            axis.ticks.length = unit(0.1, "cm"))

    layout <- c(patchwork::area(t = 1, l = 1, b = 1, r = 1),
                patchwork::area(t = 1, l = 1, b = 1, r = 1))
    p <- x_axis + plot + patchwork::plot_layout(design = layout)
    return(p)
  } else {
    return(bar)
  }
}

#' Hierarchical Clustering
#' @description perform hierarchical clustering of samples
#' based on the mutation spectra.
#' @param mf_data A data frame containing the mutation data. This data must
#' include a column containing the mutation subtypes, a column containing
#' the sample/cohort names, and a column containing the response variable.
#' @param group_col The name of the column in data that contains the
#' sample/cohort names.
#' @param response_col The name of the column in data that contains the
#' response variable. Typical response variables can be the subtype mf,
#' proportion, or count.
#' @param subtype_col The name of the column in data that contains the
#' mutation subtypes.
#' @param dist the distance measure to be used.
#' This must be one of "cosine", "euclidean", "maximum",
#' "manhattan","canberra", "binary" or "minkowski". See
#' \link[stats]{dist} for details.
#' @param cluster_method The agglomeration method to be used. See
#' \link[stats]{hclust} for details.
#' @importFrom stats hclust dist as.dist
#' @details The cosine distance measure represents the inverted cosine
#' similarity between samples:
#' 
#'\eqn{\text{Cosine Dissimilarity} = 1 - \frac{\mathbf{A} \cdot \mathbf{B}}{\| \mathbf{A} \| \cdot \| \mathbf{B} \|}}
#' 
#' This equation calculates the cosine dissimilarity between two vectors A and B. 
#' @return A dendrogram object representing the hierarchical clustering of the
#' samples.
cluster_spectra <- function(mf_data = mf_data,
                    group_col = "sample",
                    response_col = "proportion_min",
                    subtype_col = "normalized_subtype",
                    dist = "cosine",
                    cluster_method = "ward.D") {

# Get unique samples and subtypes
unique_samples <- unique(mf_data[[group_col]])
unique_subtypes <- unique(mf_data[[subtype_col]])
# Pivot Wide
mat <- matrix(0, nrow = length(unique_samples),
              ncol = length(unique_subtypes),
              dimnames = list(unique_samples, unique_subtypes))
mf_data$subtype <- as.character(mf_data$subtype)
mf_data$group <- as.character(mf_data$group)
for (i in seq_len(nrow(mf_data))) {
  mat[mf_data[[group_col]][i], mf_data[[subtype_col]][i]] <- mf_data[[response_col]][i]
}

if (dist == "cosine") {
# Calculate the cosine similarity between samples
cos_sim <- matrix(0, nrow = length(unique_samples),
                  ncol = length(unique_samples))
rownames(cos_sim) <- colnames(cos_sim) <- unique_samples
  for (i in seq_along(unique_samples)) {
  for (j in seq_along(unique_samples)) {
    cos_sim[i, j] <- sum(mat[i, ] * mat[j, ]) / (sqrt(sum(mat[i, ]^2)) * sqrt(sum(mat[j, ]^2)))
  }
}
  d <- stats::as.dist(1 - cos_sim)
} else {
    d <- stats::dist(mat, method = dist)
} 
# Perform hierarchical clustering
hc <- stats::hclust(d, method = cluster_method)
return(hc)
}
