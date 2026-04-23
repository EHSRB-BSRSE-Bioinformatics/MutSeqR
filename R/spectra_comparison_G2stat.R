#' Compare the overall mutation spectra between groups
#'
#' spectra_comparison compares the mutation spectra of groups using a
#' modified contingency table approach.
#' @param mf_data A data frame containing the MF data. This
#' is the output from calculate_mf(). MF data should be at the
#' desired subtype resolution. Required columns are the exp_variable column(s),
#' the subtype column, and sum_min or sum_max.
#' @param exp_variable The column names of the experimental variable(s) to be
#' compared.
#' @param contrasts a filepath to a file OR a dataframe that specifies the
#' comparisons to be made between levels of the exp_variable(s) The table must
#' consist of two columns, each containing a level of the exp_variable. The
#' level in the first column will compared to the level in the second column
#' for each row in contrasts. When using more than one exp_variable, separate
#' the levels of each variable with a colon. Ensure that all variables listed
#' in exp_variable are represented in each entry for the table. See details for
#' examples.
#' @param cont_sep The delimiter used to import the contrasts table.
#' Default is tab.
#' @param mf_type The type of mutation frequency to use. Default is "min"
#' (recommended).
#' @returns the log-likelihood statistic G2 for the specified comparisons with
#' the p-value adjusted for multiple-comparisons.
#' @export
#'
#' @details
#' This function creates an R * 2 contigency table of the subtype counts, where
#' R is the number of subtypes for the 2 groups being compared. The G2 likelihood
#' ratio statistic is used to evaluate whether the proportion
#' (count/group total) of each mutation subtype equals that of the other group.
#'
#' The G2 statistic refers to a chi-squared distribution to compute the p-value
#' for large sample sizes. When N / (R-1) < 20, where N is the total mutation
#' counts across both groups, the function will use a F-distribution to compute
#' the p-value in order to reduce false positive rates.
#'
#' The comparison assumes independance among the observations, as such, it is
#' highly recommended to use mf_type = "min".
#'
#' Examples of `contrasts`:
#' For 'exp_variable = "dose"` with dose groups 0, 12.5, 25, 50, compare each
#' treated dose to the control:
#'
#' 12.5 0
#'
#' 25 0
#'
#' 50 0
#'
#' Ex. Consider two 'exp_variables = c("dose", "tissue")`;
#' with levels dose (0, 12.5, 25, 50) and tissue("bone_marrow", "liver").
#' To compare the mutation spectra between tissues for each dose group,
#' the contrast table would look like:
#'
#' 0:bone_marrow	0:liver
#'
#' 12.5:bone_marrow	12.5:liver
#'
#' 25:bone_marrow	25:liver
#'
#' 50:bone_marrow	50:liver
#'
#' @examples
#' # Example data consists of 24 mouse bone marrow DNA samples imported
#' # using import_mut_data() and filtered with filter_mut. Filtered
#' # mutation data is available in the MutSeqRData ExperimentHub package:
#' # eh <- ExperimentHub::ExperimentHub()
#' # Data was summarized per sample using:
#' # calculate_mf(mutation_data = eh[["EH9861"]],
#' #              cols_to_group = "dose_group",
#' #              subtype_resolution = "base_6")
#'
#' # Example: compare 6-base mutation spectra between dose groups
#' # Load the example data
#' mf_example <- readRDS(
#'      system.file("extdata", "Example_files", "mf_data_6.rds",
#'          package = "MutSeqR"
#'      )
#' )
#' # Create the contrasts table
#' contrasts <- data.frame(
#'   col1 = c("Low", "Medium", "High"),
#'   col2 = rep("Control", 3)
#' )
#' # Run the comparison
#' spectra_comparison(
#'   mf_data = mf_example,
#'   exp_variable = "dose_group",
#'   mf_type = "min",
#'   contrasts = contrasts
#' )
#' @importFrom dplyr select mutate
#' @importFrom stats pchisq pf r2dtable
spectra_comparison <- function(mf_data,
                               exp_variable,
                               mf_type = "min",
                               contrasts,
                               cont_sep = "\t") {
    # Validation & Setup
    stopifnot(
        !missing(mf_data) && is.data.frame(mf_data),
        !missing(contrasts)
    )
    mf_type <- match.arg(mf_type, choices = c("min", "max"))

    sum_col <- paste0("sum_", mf_type)

    # Identify Subtype Column dynamically
    potential_cols <- c("variation_type", "normalized_subtype", "subtype",
                        "normalized_context_with_mutation", "context_with_mutation")
    subtype_col <- intersect(colnames(mf_data), potential_cols)[1] # Take first match

    if (is.na(subtype_col)) stop("No valid subtype column found in mf_data.")

    # Data Prep
    # Concatenate group columns
    if (length(exp_variable) > 1) {
        mf_data$group_key <- do.call(paste, c(mf_data[exp_variable], sep = ":"))
    } else {
        mf_data$group_key <- mf_data[[exp_variable]]
    }

    # Summarize sums (in case multiple rows per group/subtype exist)
    # Then Pivot Wider
    wide_data <- mf_data %>%
        dplyr::group_by(.data$group_key, dplyr::across(dplyr::all_of(subtype_col))) %>%
        dplyr::summarize(count = sum(.data[[sum_col]]), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = "group_key", values_from = "count", values_fill = 0)

    # Convert to matrix
    count_matrix <- as.matrix(wide_data[, -1])
    # Assign row names (subtypes)
    rownames(count_matrix) <- as.character(wide_data[[1]]) 

    # Process Contrasts
    if (is.data.frame(contrasts)) {
        contrast_table <- contrasts
    } else {
        contrast_table <- read.delim(contrasts, sep = cont_sep, header = FALSE)
    }

    if (ncol(contrast_table) != 2) stop("Contrast table must have exactly 2 columns.")

    # Validate groups exist
    all_groups <- colnames(count_matrix)
    requested_groups <- unique(unlist(contrast_table))
    missing_groups <- setdiff(requested_groups, all_groups)

    if (length(missing_groups) > 0) {
        stop("Groups in contrast table not found in data: ", paste(missing_groups, collapse=", "))
    }

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
  message("Warning: we have recently become aware that large mutation counts
    can lead to inflated G2 statistics and false positives. We are actively
    investigating this issue and will update the function accordingly. In the
    meantime, we recommend using this function with caution, especially for
    comparisons involving large mutation counts (> 10 mutations per sample).")
    return(results_df)
}