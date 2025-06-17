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
#' # Load the example data
#' example_file <- system.file("extdata", "Example_files",
#'                             "example_mutation_data_filtered.rds",
#'                             package = "MutSeqR")
#' example_data <- readRDS(example_file)
#'
#' # Example: compare 6-base mutation spectra between dose groups
#' # Calculate the mutation frequency data at the 6-base resolution
#' mf_data <- calculate_mf(mutation_data = example_data,
#'                         cols_to_group = "dose_group",
#'                          subtype_resolution = "base_6")
#' # Create the contrasts table
#' contrasts <- data.frame(col1 = c("Low", "Medium", "High"),
#'                         col2 = rep("Control", 3))
#' # Run the comparison
#' spectra_comparison(mf_data = mf_data,
#'                    exp_variable = "dose_group",
#'                    mf_type = "min",
#'                    contrasts = contrasts)
#' @importFrom dplyr select mutate
#' @importFrom stats pchisq pf r2dtable

spectra_comparison <- function(mf_data,
                               exp_variable,
                               mf_type = "min",
                               contrasts,
                               cont_sep = "\t") {
  # Prepare Data
  sum_col <- paste0("sum_", mf_type)
  # Find the subtype column
  subtype_col <- colnames(mf_data)[which(colnames(mf_data) %in%
      c("variation_type",
        "normalized_subtype",
        "subtype",
        "normalized_context_with_mutation",
        "context_with_mutation"))]
  if (length(subtype_col) == 0) {
    stop("No subtype column found in the mf_data")
  }
  # Select the necessary columns
  mut_spectra <- mf_data %>%
    dplyr::select(dplyr::all_of(c(exp_variable, subtype_col, sum_col)))
  # Create a single group column
  mut_spectra <- mut_spectra %>%
    dplyr::mutate(group_col = do.call(paste,
      c(dplyr::select(mut_spectra, dplyr::all_of(exp_variable)), sep = ":")
    )) %>%
    dplyr::select(-dplyr::all_of(exp_variable))
  # Sum across groups, in case of higher level grouping
  mut_spectra <- mut_spectra %>%
    dplyr::group_by(.data$group_col,
                    dplyr::across(dplyr::all_of(subtype_col))) %>%
    dplyr::summarize(sum = sum(dplyr::across(dplyr::all_of(sum_col))))
  # All groups
  groups <- unique(mut_spectra$group_col)
  # Extract data for each group
  filtered_data <- list()
  for (i in seq_along(groups)) {
    group <- groups[i]
    df_i <- mut_spectra %>%
      dplyr::filter(.data$group_col == group)
    filtered_data[[i]] <- df_i
  }

  # G2 Statistic - Likelihood Ratio Statistic
  ## Piegorsch and Bailer 1994 doi: 10.1093/genetics/136.1.403.
  G2 <- function(x, monte.carlo = FALSE, n.sim = 10000, seed = 1234) {
    N <- sum(x)
    r <- apply(x, 1, sum)
    c <- apply(x, 2, sum)
    e <- r %*% t(c) / N
    G2 <- 0
    for (k in 1:ncol(x)){
      flag <- x[, k] > 0
      G2 <- G2 + t(x[flag, k]) %*% log(x[flag, k] / e[flag, k])
    }
    G2 <- 2 * G2
    R <- nrow(x) - 1
    df <- R * (ncol(x) - 1)
    if (monte.carlo == FALSE) {
      if (N / R > 20) {
        p.value <- 1 - pchisq(G2, df)
        message("Using chi-squared distribution to compute p-value")
      } else {
        p.value <- 1 - pf(G2 / R, R, N - df)
        message("Using F-distribution to compute p-value")
      }
    } else {
      #Monte Carlo
      #Generate random rxc tables
      set.seed(seed)
      r <- apply(x, 1, sum)
      c <- apply(x, 2, sum)
      rtbl <- r2dtable(1, r, c)
      ref.dist <- rep(0, n.sim)
      for (k in 1:length(ref.dist)){
        x <- r2dtable(1, r, c)[[1]]
        N <- sum(x)
        r <- apply(x, 1, sum)
        c <- apply(x, 2, sum)
        e <- r %*% t(c) / N
        G2.t <- 0
        for (j in 1:ncol(x)){
          flag <- x[, j] > 0
          G2.t <- G2.t + t(x[flag, j]) %*% log(x[flag, j] / e[flag, j])
        }
        ref.dist[k] <- 2 * G2.t
      }
      flag <- ref.dist >= G2[1, 1]
      p.value <- length(ref.dist[flag]) / 10000
    }
    data.frame(G2 = G2, p.value = p.value)
  }

  # Load in contrast table
  if (is.null(contrasts)) {
    stop("Please provide a contrasts table")
  }
  if (is.data.frame(contrasts)) {
    contrast_table <- contrasts
  } else {
    contrast_table <- read.delim(file.path(contrasts), sep = cont_sep, header = F)
    if (ncol(contrast_table) <= 1) {
      stop("Your contrast_table only has one column. Make sure to set the proper delimiter with cont_sep.")
    }
  }

  contrast_data <- list()
  # Loop over the rows of the contrasts table and select
  # the corresponding data frames from the 'filtered_data' list.
  # Combine the mut counts of both dataframes into one.
  for (i in seq_len(nrow(contrast_table))) {
    indices <- match(contrast_table[i, ], groups)
    dfs <- filtered_data[indices]
    df_combined <- do.call(cbind, lapply(dfs, function(df) df[, 3])) # here column 3 is the sum data
    contrast_data[[i]] <- df_combined
  }

  # Run the G2 function for all contrasts
  results <- data.frame(contrasts = character(),
                        G2 = numeric(),
                        p.value = numeric(),
                        stringsAsFactors = FALSE)
  for (i in seq_len(length(contrast_data))) {
    result <- G2(contrast_data[[i]])
    contrast_str <- paste(contrast_table[i, 1], "vs", contrast_table[i, 2])
    results <- rbind(results, data.frame(contrasts = contrast_str,
                                         G2 = result$G2,
                                         p.value = result$p.value,
                                         stringsAsFactors = FALSE))
  }

  # Apply the Holm-Sidak correction for multiple comparisons
  results$adj_p.value <- MutSeqR::sidak(results$p.value)$SidakP
  results$Significance <- ""
  results$Significance[results$adj_p.value < 0.05] <- "***"

  return(results)
}
