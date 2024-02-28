#' Compare the overall mutation spectra between groups
#'
#' spectra_comparison compares the mutation spectra of groups using a
#' modified contingency table approach.
#' @param mf_data The data frame containing the mutation count data.
#' Mutation counts should be summarized per group for each mutation subtype of
#' interest. This data can be obtained using `calculate_mut_freq(summary=TRUE)`.
#' @param muts The column in mf_data containing the mutation count per group.
#' @param subtype_col The column in mf_data listing the mutation subtypes.
#' @param group A vector of grouping variables: this should be the groups
#' of interest across which you want to compare mutation spectra.
#' Ex. c("dose", "tissue")
#' @param contrast_table_file a filepath to a tab-delimited `.txt` file that
#' will specify the comparisons to be made between groups. The table must
#' consist of two columns. The first column will be a level within your group
#' column and the second column must be the group level that it will be
#' compared to. All values must correspond to entries in your `group`
#' column. For more than one `group` variable, separate the levels of the
#' each group with a colon. Ensure that all `groups` are represented in each
#' entry for the table. See details for examples.
#' @param cont_sep The delimiter used in the contrast_table_file.
#' Default is tab.
#' @returns the log-likelihood statistic G2 for the specified comparisons with
#' the p-value adjusted for multiple-comparisons.
#' @export
#'
#' @details
#' Examples of `contrast_table_file`:
#' If you have `group = "dose"` with dose groups 0, 25, 50, 100. The first
#' column would contain the treated groups (25, 50, 100), while the second
#' column would be 0, thus comparing each treated group to the control group.
#'
#' 25 0
#'
#' 50 0
#'
#' 100 0
#'
#' Ex. Consider two grouping variables `group = c("dose", "tissue")`;
#' with levels dose (0, 25, 50, 100) and tissue("bone_marrow", "liver").
#' To compare the mutation spectra between tissues for each dose group,
#' the contrast table would look like:
#'
#' 0:bone_marrow	0:liver
#'
#' 25:bone_marrow	25:liver
#'
#' 50:bone_marrow	50:liver
#'
#' 100:bone_marrow 100:liver
#'
#' @importFrom dplyr select mutate
#' @importFrom stats pchisq pf r2dtable

spectra_comparison <- function(mf_data,
                               muts = "dose_sum_unique",
                               subtype_col = "normalized_subtype",
                               group = c("dose", "tissue"),
                               contrast_table_file,
                               cont_sep = "\t") {
  mut_spectra <- mf_data %>%
    dplyr::select(dplyr::all_of(c(subtype_col, muts, group)))
  mut_spectra <- mut_spectra %>%
    dplyr::mutate(group_col = do.call(paste,
                                      c(dplyr::select(mut_spectra,
                                                      dplyr::all_of(group)),
                                        sep = ":"))) %>%
    dplyr::select(-all_of(group))

  pivoted_matrix <- matrix(nrow = length(unique(mut_spectra$group_col)),
                           ncol = length(unique(mut_spectra[[subtype_col]])))
  rownames(pivoted_matrix) <- unique(mut_spectra$group_col)
  colnames(pivoted_matrix) <- unique(mut_spectra[[subtype_col]])
  # Fill in the matrix with the values from the original data frame
  for (i in seq_len(nrow(mut_spectra))) {
   row_index <- which(rownames(pivoted_matrix) == mut_spectra$group_col[i])
  col_index <- which(colnames(pivoted_matrix) == mut_spectra[[subtype_col]][i])
    pivoted_matrix[row_index, col_index] <- mut_spectra[[muts]][i]
  }
  # Convert the matrix to a data frame
  pivoted_mut <- as.data.frame(pivoted_matrix)
  pivoted_mut <- replace(pivoted_mut, is.na(pivoted_mut), 0)

  # G2 Statistic - Likelihood Ratio Statistic
  calculate_g2 <- function(x, monte_carlo = FALSE, n_sim = 10000, seed = 1234) {
    n <- sum(x) # total mut sum
    r <- apply(x, 1, sum) # mut sum of each subtype
    c <- apply(x, 2, sum) # mut sum for each group
    e <- r %*% t(c) / n # matrix multiplication/total

    g2 <- 0
    for (k in seq_len(ncol(x))) {
      # flag rows in k that are < 0; do not perform operations on these rows.
      flag <- x[, k] > 0
      g2 <- g2 + t(x[flag, k]) %*% log(x[flag, k] / e[flag, k])
    }
    g2 <- 2 * g2

    r <- nrow(x) - 1
    df <- r * (ncol(x) - 1)

    if (monte_carlo == FALSE) {
      if (n / r > 20) {
        p_value <- 1 - stats::pchisq(g2, df)
        message("Using chi-squared distribution to compute p-value")
      } else {
        p_value <- 1 - stats::pf(g2 / r, r, n - df)
        message("Using F-distribution to compute p-value")
      }
      data.frame(g2 = g2, p_value = p_value)
    } else {
      # Monte Carlo
      # Generate random rxc tables
      set.seed(seed)
      r <- apply(x, 1, sum)
      c <- apply(x, 2, sum)

      ref_dist <- rep(0, n_sim)
      for (k in seq_along(ref_dist)) {
        x <- stats::r2dtable(1, r, c)[[1]]
        n <- sum(x)
        r <- apply(x, 1, sum)
        c <- apply(x, 2, sum)
        e <- r %*% t(c) / n

        g2_t <- 0
        for (j in seq_len(ncol(x))) {
          flag <- x[, j] > 0
          g2_t <- g2_t + t(x[flag, j]) %*% log(x[flag, j] / e[flag, j])
        }

        ref_dist[k] <- 2 * g2_t
      }
      flag <- ref_dist >= g2[1, 1]
      p_value <- length(ref_dist[flag]) / 10000
      data.frame(g2 = g2_t, p_value = p_value)
    }
  }


  # Load in contrast table
  if (!is.null(contrast_table_file)) {
    contrast_table <- read.delim(file.path(contrast_table_file),
      sep = cont_sep,
      header = FALSE
    )
    if (ncol(contrast_table) <= 1) {
      stop("Your contrast_table only has one column. Make sure to set
            the proper delimiter with cont_sep.")
    }
  } else {
    stop("Please provide a contrast table file")
  }

  # Create an empty list to store the results
  results_list <- list()

  # Iterate over each row of the contrast table
  for (i in seq_len(nrow(contrast_table))) {
    # Extract the row names from the contrast table
    rowname1 <- rownames(pivoted_mut)[match(contrast_table[i, "V1"], rownames(pivoted_mut))]
    rowname2 <- rownames(pivoted_mut)[match(contrast_table[i, "V2"], rownames(pivoted_mut))]

    # Extract the corresponding rows from the dataframe
    data_subset <- pivoted_mut[c(rowname1, rowname2), ]
    data_subset_t <- t(data_subset)

    # Apply the G2 function to compare the two rows
    result <- calculate_g2(data_subset_t)

    # Store the result
    results_list[[i]] <- result
  }

  # Bind all the results together into a single dataframe
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- paste(contrast_table[, 1], "vs", contrast_table[, 2])
  # adjust p-values for multiple-comparisons
  results_df$adj_p_value <- MutSeqR::my.holm.sidak(results_df$p_value)

  return(results_df)
}
