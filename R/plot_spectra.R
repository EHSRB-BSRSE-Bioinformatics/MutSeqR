#' GenVisR Transition-transversion plot
#' @description Given a data frame construct a plot displaying the
#' mutation subtypes observed in a cohort.
#' @param mf_data A data frame containing the mutation data. This data must
#' include a column containing the mutation subtypes, a column containing
#' the sample/cohort names, and a column containing the response variable.
#' Typical response variables can be the subtype frequency, proportion, or
#' count.
#' @param sample_col The name of the column in data that contains the
#' sample/cohort names.
#' @param response_col The name of the column in data that contains the
#' response variable. Typical response variables can be the subtype frequency,
#' proportion, or count.
#' @param subtype_col The name of the column in data that contains the
#' mutation subtypes.
#' @param sample_order The method for ordering the samples within the plot.
#' Options include:
#' \itemized {
#'   \item `none`: No ordering is performed. Default.
#'   \item `smart`: Samples are ordered based on the sample names.
#'   \item `arranged`: Samples are ordered based on one or more factor column(s)
#' in mf_data. Column names are passed to the function using the
#' `sample_order_input`.
#'  \item `custom`: Samples are ordered based on a custom vector of sample
#' names. The custom vector is passed to the function using the
#' `sample_order_input`.
#' \item `clustered`: Samples are ordered based on hierarchical clustering. The
#' dissimilarity matrix can be specified using the `dist` argument. The
#' agglomeration method can be specified using the `cluster_method` argument.
#' }
#' @param sample_order_input A character vector specifying details for the
#' sample order method. If `sample_order` is `arranged`, `sample_order_input`
#' should contain the column name(s) to be used for ordering the samples. If
#' `sample_order` is `custom`, `sample_order_input` should contain the custom
#' vector of sample names.
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
#' @import patchwork
#' @importFrom ggh4x scale_x_dendrogram
#' @importFrom ggplot2 
#' @importFrom dplyr select arrange across all_of
#' @importFrom gtools mixedsort
#' @export

plot_spectra <- function(mf_data = mf_data,
                         sample_col = "sample",
                         response_col = "proportion_unique",
                         subtype_col = "normalized_subtype",
                         sample_order = "none",
                         sample_order_input = NULL,
                         dist = "cosine",
                         cluster_method = "ward.D",
                         palette = NULL,
                         x_lab = NULL,
                         y_lab = NULL) {

  plot_data <- mf_data %>%
    dplyr::select({{sample_col}}, {{subtype_col}}, {{response_col}})

  if (sample_order == "none") {
    order <- as.vector(unique(plot_data[[sample_col]]))
    mf_data[[sample_col]] <- factor(mf_data[[sample_col]])
  } else if (sample_order == "smart") {
    order <- as.vector(unique(mf_data[[sample_col]]))
    order <- gtools::mixedsort(order)
    mf_data[[sample_col]] <- factor(mf_data[[sample_col]], levels = order)
  } else if (sample_order == "arranged") {
    plot_data <- mf_data %>%
      dplyr::select({{sample_col}},
                    {{subtype_col}},
                    {{response_col}},
                    {{sample_order_input}}) %>%
      dplyr::arrange(dplyr::across(dplyr::all_of({{sample_order_input}})))
    order <- as.vector(unique(plot_data[[sample_col]]))
    plot_data[[sample_col]] <- factor(plot_data[[sample_col]], levels = order)
    plot_data <- plot_data %>%
      dplyr::select({{sample_col}},
                    {{subtype_col}},
                    {{response_col}})
  } else if (sample_order == "custom") {
    plot_data[[sample_col]] <- factor(plot_data[[sample_col]],
                                     levels = sample_order_input)
  } else if (sample_order == "clustered") {
    # Cluster the samples
    hc <- cluster(mf_data = plot_data,
                  sample_col = sample_col,
                  response_col = response_col,
                  subtype_col = subtype_col,
                  dist = dist,
                  cluster_method = cluster_method)
    # Reorder the samples based on hierarchical clustering
    order <- hc$labels[hc$order]
    # Reorder the levels of the sample variable in your data frame
    plot_data[[sample_col]] <- factor(plot_data[[sample_col]],
                                      levels = order)
  }

  subtype_order <- c("symbolic", "complex", "insertion",
                     "deletion", "mnv", "T>G", "T>C",
                     "T>A", "C>T", "C>G", "C>A")

  plot_data[[subtype_col]] <- factor(plot_data[[subtype_col]],
                                     levels = subtype_order)

  if (is.null(palette)) {
    palette <- c("C>A" = "#3288BD",
                 "C>G" = "#99D594",
                 "C>T" = "#E6F598",
                 "T>A" = "#FEE08B",
                 "T>C" = "#FC8D59",
                 "T>G" = "#D53E4F",
                 "mnv" = "pink",
                 "deletion" = "black",
                 "insertion" = "grey",
                 "symbolic" = "purple")
  }

  # Axis labels
  if (is.null(x_lab)) {
    x_lab <- sample_col
  }
  if (is.null(y_lab)) {
    y_lab <- response_col
  }
  axis_labels <- ggplot2::labs(x = x_lab, y = y_lab)

  # bar plot
  bar <- ggplot(plot_data, aes(x = .data[[sample_col]],
                               y = .data[[response_col]],
                               fill = .data[[subtype_col]],
                               add = FALSE)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = palette) +
    axis_labels +
    theme_minimal() +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 90)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01)))

  if (sample_order == "clustered") {
    plot <-  bar +
      ggh4x::scale_x_dendrogram(hclust = hc, position = "top", labels = NULL,
      ) +
      theme(axis.ticks.length.x = unit(10, "pt"))
    x_axis <- ggplot(plot_data, aes(x = .data[[sample_col]])) +
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
#' @param sample_col The name of the column in data that contains the
#' sample/cohort names.
#' @param response_col The name of the column in data that contains the
#' response variable. Typical response variables can be the subtype frequency,
#' proportion, or count.
#' @param subtype_col The name of the column in data that contains the
#' mutation subtypes.
#' @param dist the distance measure to be used.
#' This must be one of "cosine", "euclidean", "maximum",
#' "manhattan","canberra", "binary" or "minkowski". See
#' \link[stats]{dist} for details.
#' @param cluster_method The agglomeration method to be used. See
#' \link[stats]{hclust} for details.
#' @details The cosine distance measure represents the inverted cosine
#' similarity between samples:
#' 
#'\eqn{\text{Cosine Dissimilarity} = 1 - \frac{\mathbf{A} \cdot \mathbf{B}}{\| \mathbf{A} \| \cdot \| \mathbf{B} \|}}
#' 
#' This equation calculates the cosine dissimilarity between two vectors A and B. 
#' @return A dendrogram object representing the hierarchical clustering of the
#' samples.
cluster <- function(mf_data = mf_data,
                    sample_col = "sample",
                    response_col = "proportion_unique",
                    subtype_col = "normalized_subtype",
                    dist = "cosine",
                    cluster_method = "ward.D") {

# Get unique samples and subtypes
unique_samples <- unique(mf_data[[sample_col]])
unique_subtypes <- unique(mf_data[[subtype_col]])
# Pivot Wide
mat <- matrix(0, nrow = length(unique_samples),
              ncol = length(unique_subtypes),
              dimnames = list(unique_samples, unique_subtypes))
for (i in seq_len(nrow(mf_data))) {
  mat[mf_data[[sample_col]][i], mf_data[[subtype_col]][i]] <- mf_data[[response_col]][i]
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
  d <- as.dist(1 - cos_sim)
} else {
    d <- dist(mat, method = dist)
} 
# Perform hierarchical clustering
hc <- hclust(d, method = cluster_method)
return(hc)
}
