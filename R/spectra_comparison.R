#' Compare mutation spectra between groups
#'
#' Compare mutation spectra using either a G2 contingency-table test or a
#' replicate-aware, ridge-penalized Dirichlet-multinomial model.
#'
#' @param mf_data A data frame from [calculate_mf()] at the desired subtype
#' resolution, with the experimental variable columns, a subtype column, and
#' `sum_min` or `sum_max`. For `method = "rpdmm"`, the data must be calculated
#' per sample and include a `sample` column.
#' @param exp_variable Column names of the experimental variables to compare.
#' @param mf_type Mutation count type, either `"min"` (default) or `"max"`.
#' @param contrasts A two-column data frame or path to a file specifying the
#' group levels to compare. For multiple experimental variables, join the
#' levels in each entry with a colon, in the order given by `exp_variable`.
#' @param cont_sep Delimiter for a contrasts file. Defaults to a tab.
#' @param method Comparison method: `"g2"` (default) or `"rpdmm"`.
#' @param lambda Ridge penalty strength for `"rpdmm"`. Higher values increase
#' shrinkage towards a uniform distribution. Defaults to 1.
#' @param n_boot Number of bootstrap iterations for `"rpdmm"`. Defaults to 500.
#' `lambda` and `n_boot` cannot be supplied when `method = "g2"`.
#'
#' @returns A data frame with one row per contrast. The `"g2"` method returns
#' the G2 statistic and a Sidak-adjusted p-value. The `"rpdmm"` method returns
#' hybrid and penalized bootstrap p-values, likelihood-ratio statistics, and
#' concentration estimates. Its p-values are not adjusted for multiple
#' comparisons.
#' @details
#' The G2 test pools counts across samples and assumes independent
#' observations. Large mutation counts can inflate its statistic and cause
#' false positives. The `"rpdmm"` method accounts for
#' replicate heterogeneity and overdispersion.
#'
#' @examples
#' mf_example <- readRDS(
#'   system.file("extdata", "Example_files", "mf_data_6.rds",
#'               package = "MutSeqR")
#' )
#' contrasts <- data.frame(
#'   col1 = c("Low", "Medium", "High"),
#'   col2 = rep("Control", 3)
#' )
#' spectra_comparison(mf_example, "dose_group", contrasts = contrasts)
#'
#' @export
spectra_comparison <- function(mf_data,
                               exp_variable,
                               mf_type = "min",
                               contrasts,
                               cont_sep = "\t",
                               method = c("g2", "rpdmm"),
                               lambda = 1.0,
                               n_boot = 500) {
  method <- match.arg(method)

  stopifnot(
    !missing(mf_data) && is.data.frame(mf_data),
    !missing(exp_variable) && is.character(exp_variable) &&
      length(exp_variable) > 0 && all(exp_variable %in% colnames(mf_data)),
    !missing(contrasts),
    is.numeric(lambda) && lambda >= 0,
    is.numeric(n_boot) && length(n_boot) == 1 && is.finite(n_boot) &&
      n_boot > 0 && n_boot == floor(n_boot)
  )
  if (method == "rpdmm") {
    stopifnot("sample" %in% colnames(mf_data))
  }
  mf_type <- match.arg(mf_type, choices = c("min", "max"))
  sum_col <- paste0("sum_", mf_type)
  if (!sum_col %in% colnames(mf_data)) {
    stop("Missing count column in mf_data: ", sum_col)
  }

  potential_cols <- as.vector(subtype_dict)
  potential_cols <- potential_cols[!is.na(potential_cols)]
  subtype_cols <- intersect(colnames(mf_data), potential_cols)
  if (length(subtype_cols) == 0) {
    stop("No valid subtype column found in mf_data.")
  }
  subtype_col <- subtype_cols[1]
  if (method == "rpdmm" && length(subtype_cols) > 1) {
    warning("Multiple potential subtype columns found. Using: ", subtype_col)
  }

  group <- if (length(exp_variable) > 1) {
    do.call(paste, c(mf_data[exp_variable], sep = ":"))
  } else {
    mf_data[[exp_variable]]
  }
  counts <- data.frame(
    subtype = mf_data[[subtype_col]],
    group = group,
    count = mf_data[[sum_col]],
    stringsAsFactors = FALSE
  )
  if (method == "rpdmm") {
    counts$sample <- mf_data$sample
    matrix_data <- split(counts, counts$group)
    matrix_column <- "sample"
  } else {
    counts <- dplyr::arrange(counts, .data$group, .data$subtype)
    matrix_data <- list(counts)
    matrix_column <- "group"
  }
  matrices <- lapply(matrix_data, function(group_counts) {
    wide_data <- group_counts %>%
      dplyr::select("subtype", dplyr::all_of(matrix_column), "count") %>%
      tidyr::pivot_wider(
        names_from = dplyr::all_of(matrix_column),
        values_from = "count",
        values_fill = 0,
        values_fn = sum
      )
    count_matrix <- as.matrix(wide_data[, -1, drop = FALSE])
    rownames(count_matrix) <- as.character(wide_data[[1]])
    count_matrix
  })
  if (method == "g2") {
    matrices <- matrices[[1]]
    all_groups <- colnames(matrices)
  } else {
    all_groups <- names(matrices)
  }

  if (is.data.frame(contrasts)) {
    contrast_table <- contrasts
  } else {
    contrast_table <- read.delim(contrasts, sep = cont_sep, header = FALSE)
  }
  if (ncol(contrast_table) != 2) {
    stop("Contrast table must have exactly 2 columns.")
  }
  missing_groups <- setdiff(unique(unlist(contrast_table)), all_groups)
  if (length(missing_groups) > 0) {
    stop("Groups in contrast table not found in data: ",
         paste(missing_groups, collapse = ", "))
  }

  if (method == "g2") {
    spectra_comparison_G2(matrices, contrast_table)
  } else {
    spectra_comparison_rpdmm(
      matrices, contrast_table, lambda = lambda, n_boot = n_boot
    )
  }
}

#' Compare the overall mutation spectra between groups using G2
#'
#' @param count_matrix Matrix of subtype counts, with groups as columns.
#' @param contrast_table Two-column data frame of group comparisons.
#' @returns The log-likelihood statistic G2 for the specified comparisons with
#' the p-value adjusted for multiple comparisons.
#' @keywords internal
#'
#' @details
#' This function creates an R * 2 contingency table of the subtype counts,
#' where R is the number of subtypes for the 2 groups being compared. The G2
#' likelihood ratio statistic is used to evaluate whether the proportion
#' (count/group total) of each mutation subtype equals that of the other group.
#'
#' The G2 statistic refers to a chi-squared distribution to compute the p-value
#' for large sample sizes. When N / (R-1) < 20, where N is the total mutation
#' counts across both groups, the function will use an F-distribution to
#' compute the p-value in order to reduce false positive rates.
#'
#' The comparison assumes independence among the observations; as such, it is
#' highly recommended to use `mf_type = "min"`.
#'
#' Examples of `contrasts` passed to [spectra_comparison()]:
#' For `exp_variable = "dose"` with dose groups 0, 12.5, 25, 50, compare each
#' treated dose to the control:
#'
#' 12.5 0
#'
#' 25 0
#'
#' 50 0
#'
#' For `exp_variable = c("dose", "tissue")` with dose levels 0, 12.5, 25, 50
#' and tissue levels "bone_marrow" and "liver", compare tissues at each dose:
#'
#' 0:bone_marrow  0:liver
#'
#' 12.5:bone_marrow  12.5:liver
#'
#' 25:bone_marrow  25:liver
#'
#' 50:bone_marrow  50:liver
#'
#' @examples
#' mf_example <- readRDS(
#'   system.file("extdata", "Example_files", "mf_data_6.rds",
#'               package = "MutSeqR")
#' )
#' contrasts <- data.frame(
#'   col1 = c("Low", "Medium", "High"),
#'   col2 = rep("Control", 3)
#' )
#' spectra_comparison(
#'   mf_data = mf_example, exp_variable = "dose_group",
#'   mf_type = "min", contrasts = contrasts, method = "g2"
#' )
#' @importFrom dplyr select mutate
#' @importFrom stats pchisq pf r2dtable
spectra_comparison_G2 <- function(count_matrix, contrast_table) {
    # Define G2 Function
    calculate_g2_pair <- function(g1_name, g2_name) {
        # Extract columns for the two groups
        obs <- count_matrix[, c(g1_name, g2_name), drop = FALSE]
        obs <- obs[rowSums(obs) > 0, , drop = FALSE] # Remove subtypes with 0 counts in BOTH groups
        N <- sum(obs)
        R <- nrow(obs)
        C <- ncol(obs) # Always 2 here

        # Expected Counts: (RowSum * ColSum) / Total
        row_sums <- rowSums(obs)
        col_sums <- colSums(obs)
        expected <- outer(row_sums, col_sums) / N

        # G2 Statistic = 2 * sum(Obs * log(Obs/Exp))
        valid <- obs > 0 # avoid log(0) or NaN
        term <- obs[valid] * log(obs[valid] / expected[valid])
        G2_stat <- 2 * sum(term)

        # Degrees of Freedom
        df <- (R - 1) * (C - 1)

        # P-value calculation
        if (N / (R - 1) > 20) {
            p_val <- 1 - stats::pchisq(G2_stat, df)
        } else {
            # F-distribution approximation for small samples
            p_val <- 1 - stats::pf(G2_stat / (R - 1), R - 1, N - df)
        }

        return(c(G2 = G2_stat, p.value = p_val))
    }

    # Run Comparisons
    results_matrix <- mapply(calculate_g2_pair,
                            as.character(contrast_table[, 1]),
                            as.character(contrast_table[, 2]))

    # Transpose results (mapply returns Cols=Contrasts, Rows=Stats)
    results_df <- as.data.frame(t(results_matrix))

    # Final Formatting
    results_df$contrasts <- paste(contrast_table[, 1], "vs", contrast_table[, 2])

    # Adjust P-values
    results_df$adj_p.value <- MutSeqR::sidak(results_df$p.value)$SidakP

    results_df$Significance <- dplyr::case_when(
        results_df$adj_p.value < 0.001 ~ "***",
        results_df$adj_p.value < 0.01  ~ "**",
        results_df$adj_p.value < 0.05  ~ "*",
        TRUE ~ ""
    )

    # Reorder columns
    results_df <- results_df[, c("contrasts", "G2", "p.value", "adj_p.value", "Significance")]
  message("Warning: large mutation counts can lead to inflated G2 statistics
    and false positives. For comparisons involving large mutation counts
    (> 10 mutations per sample), we recommend using method = rpdmm .")
    return(results_df)
}

#' Compare mutation spectra using a ridge-penalized Dirichlet-multinomial model
#'
#' @param list_of_matrices Named list of subtype-by-sample count matrices, one
#' per experimental group.
#' @param contrast_table Two-column data frame of group comparisons.
#' @param lambda Strength of the ridge penalty.
#' @param n_boot Number of bootstrap iterations.
#' @returns A data frame with one row per contrast, including hybrid and
#' penalized bootstrap p-values, likelihood-ratio statistics, and concentration
#' estimates for each group and the shared model.
#' @keywords internal
#' @importFrom dplyr select rename left_join
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom utils read.delim
spectra_comparison_rpdmm <- function(list_of_matrices, contrast_table,
                                     lambda = 1.0, n_boot = 500) {
  # Loop through each row of the contrast table and apply
  # run_penalized_comparison to the corresponding matrices.
  results_list <- lapply(seq_len(nrow(contrast_table)), function(i) {

    # 1. Identify the groups for this comparison
    g1_name <- as.character(contrast_table[i, 1])
    g2_name <- as.character(contrast_table[i, 2])

    # 2. Extract their matrices
    g1_mat <- list_of_matrices[[g1_name]]
    g2_mat <- list_of_matrices[[g2_name]]

    # 3. Ensure both matrices have exactly the same rows (subtypes)
    all_subtypes <- sort(unique(c(rownames(g1_mat), rownames(g2_mat))))
    align_matrix <- function(mat, rnames) {
      missing <- setdiff(rnames, rownames(mat))
      if (length(missing) > 0) {
        # Create a matrix of zeros for the missing subtypes
        pad <- matrix(0, nrow = length(missing), ncol = ncol(mat),
                      dimnames = list(missing, colnames(mat)))
        mat <- rbind(mat, pad)
      }
      # Return sorted so both matrices perfectly match row-for-row
      mat[rnames, , drop = FALSE]
    }
    g1_aligned <- align_matrix(g1_mat, all_subtypes)
    g2_aligned <- align_matrix(g2_mat, all_subtypes)

    # 4. Run the statistical comparison
    stats <- run_penalized_comparison(
      g1 = g1_aligned,
      g2 = g2_aligned,
      lambda = lambda,
      n_boot = n_boot
    )

    # 5. Format results as a 1-row data frame
    data.frame(
      group1 = g1_name,
      group2 = g2_name,
      comparison = paste0(g1_name, " vs ", g2_name),
      p_value = stats$p_value,
      lrt = stats$lrt,
      p_hybrid = stats$p_hybrid,
      p_penalized = stats$p_penalized,
      lrt_penalized = stats$lrt_penalized,
      theta_g1 = stats$theta_g1,
      theta_g2 = stats$theta_g2,
      theta_shared = stats$theta_shared,
      stringsAsFactors = FALSE,
      Significance = ifelse(stats$p_value < 0.05, "*", "")
    )
  })

  # Combine the list of 1-row dataframes into one final dataframe
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL

  return(results_df)
}
