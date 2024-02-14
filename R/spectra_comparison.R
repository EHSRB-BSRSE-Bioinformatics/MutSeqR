#' Compare the overall mutation spectra between groups
#' 
#' spectra_comparison compares the mutation spectra of groups using a modified contingency table approach
#' @param mf_data The data frame containing the mutation count data.
#' Mutation counts should be summarized per group for each mutation subtype of 
#' interest. This data can be obtained using `calculate_mut_freq(summary=TRUE)`. 
#' @param muts The column in mf_data containing the mutation count per group.
#' @param subtype_col The column in mf_data listing the mutation subtypes.
#' @param group A vector of grouping variables: this should be the groups
#' of interest across which you want to compare mutation spectra. 
#' Ex. c("dose", "tissue")
#' @param contrast_table_file a filepath to a tab-delimited `.txt` file that will 
#' specify the comparisons to be made between groups. The table must consist of
#' two columns. The first column will be a level within your group column and the 
#' second column must be the group level that it will be compared to. All values 
#' must correspond to entries in your `group` column. For more than one `group` 
#' variable, separate the levels of the each group with a colon. Ensure that all 
#' `groups` are represented in each entry for the table. See details
#' for examples.
#' @param cont_sep The delimiter used in the contrast_table_file. Default is tab.
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
#' Ex. Consider two grouping variables `group = c("dose", "tissue")`; with levels
#' dose(0, 25, 50, 100) and tissue("bone_marrow", "liver"). To compare
#' the mutation spectra between tissues for each dose group, the 
#' contrast table would look like:
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

# Piergorsch & Bailer, 1993
  # R x T
  # Ri number of subtypes (rows)
  # Tj number of treatment groups (columns)
  # N mut total
  # p proportion; reps the response probability in the(i,j)th cell
   # p1j + p2j + ... + prj = 1.0 (ie the proportion of subtypes (R) in each group (Tj) should equal 1)
    # H0 = pi1 = pi2 = ... = piT (ie the proportion of each subtype (Ri) is the same across groups T)
# To test the H0, use the likelihood ratio statistic:
  # G2 = 2 Ep(R; i=1)Ep(T; j=1)Yij * log{Yij/Eij}
  # df = (R - 1)(T - 1)
# When the test-stat > a chi2 table value w df, reject H0
# When ALL treatment groups have a subtype count as 0, you can remove it to reduce the df and increase the sensitivity of the test stat.
  # Yi1 + Yi2 + ... + YiT = 0
# Monte Carlo's approximation to Fisher's exact test:
  # When TRUE; randomly generate a large number of R x T tables with the 
  # same row and column totals as observed in R x T table. Each random table is compared
  # to the original observed table to determine if the random table exhibits greater
  # departure from H0. The proportion of randomly generated tables that exhibits this
  # departure is an estimate of the "exact" P-value. 
  # recomended when R x T table becomes too large (when N > 5df)
# G2 exhibits high false positive rates in small samples when referred to with a chi2 distribution. 

# Gabriel 1966
# Test which subtype is driving the difference between groups
# simultaneous inference of R x T sub-tables (ie which combination of subtypes exhibits H0?)
# Sub-table of form A x B
  # 2 <= A <= I (non-0 rows)
  # 2 <= B <= T columns
# Compute the G2 for any/all A x B sub-tables.
# df = (i - 1)(T - 1) Note it uses I and T, not A and B; restricts false positive rate
# Although G2 isn't the best test for the global spectra, it is the only one that can be used for Gabriel's method
# For large values of N/(R-1) G2 is stable.. (tot mut)/(# subtypes - 1)... so lots of mutants = good test!
  # When N/(R-1)<20, do some adjustments, switch to an F-distribution... but we are waaay above that. 
# So we run G2ab on lots of the sub-tables, until we find significance - do we need to run multiple comparisons?


spectra_comparison <- function(mf_data,
                               muts = "dose_sum_unique",
                               subtype_col = "normalized_subtype",
                               group = c("dose", "tissue"),
                               contrast_table_file,
                               cont_sep = "\t"){
  
  mut_spectra <- mf_data %>%
    dplyr::select(dplyr::all_of(c(subtype_col, muts, group)))
  mut_spectra <- mut_spectra %>%
    dplyr::mutate(group_col = do.call(paste, c(dplyr::select(mut_spectra, dplyr::all_of(group)), sep = ":"))) %>%
    dplyr::select(-all_of(group))
 # pivot wider
  pivoted_matrix <- matrix(nrow = length(unique(mut_spectra$group_col)), ncol = length(unique(mut_spectra[[subtype_col]])))
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
  
  #G2 Statistic - Likelihood Ratio Statistic
  g2 <- function(x, monte.carlo = FALSE, n.sim = 10000, seed = 1234) {
    n <- sum(x) # total mut sum
    r <- apply(x, 1, sum) # mut sum of each subtype
    c <- apply(x, 2, sum) # mut sum for each group
    e <- r %*% t(c) / n # matrix multiplication/total
    
    g2 <- 0
    for (k in seq_len(ncol(x))) {
      flag <- x[, k] > 0 # flag rows in k that are < 0; do not perform operations on these rows.
      g2 <- g2 + t(x[flag, k]) %*% log(x[flag, k] / e[flag, k])
    }
    g2 <- 2 * g2
    
    r <- nrow(x) - 1
    df <- r * (ncol(x) - 1)
    
    if (monte.carlo == FALSE) {
      if (n / r > 20) {
        p.value <- 1 - stats::pchisq(g2, df)
        message("Using chi-squared distribution to compute p-value")
      } else {
        p.value <- 1 - stats::pf(g2 / r, r, n - df)
        message("Using F-distribution to compute p-value")
      }
    data.frame(g2 = g2, p.value = p.value)
    } else {
      # Monte Carlo
      # Generate random rxc tables
      set.seed(seed)
      r <- apply(x, 1, sum)
      c <- apply(x, 2, sum)

      ref.dist <- rep(0, n.sim)
      for (k in seq_along(ref.dist)) {
        x <- stats::r2dtable(1, r, c)[[1]]
        n <- sum(x)
        r <- apply(x, 1, sum)
        c <- apply(x, 2, sum)
        e <- r %*% t(c) / n
        
        g2.t <- 0
        for (j in seq_len(ncol(x))) {
          flag <- x[, j] > 0
          g2.t <- g2.t + t(x[flag, j]) %*% log(x[flag, j] / e[flag, j])
        }
        
        ref.dist[k] <- 2 * g2.t
      }
      flag <- ref.dist >= g2[1, 1]
      p.value <- length(ref.dist[flag]) / 10000
      data.frame(g2 = g2.t, p.value = p.value)
    }
  }

 
  # Load in contrast table 
  if (!is.null(contrast_table_file)) {
    contrast_table <- read.delim(file.path(contrast_table_file), sep = cont_sep,
                                 header = F)
    if (ncol(contrast_table) <= 1) {
      stop("Your contrast_table only has one column. Make sure to set the proper delimiter with cont_sep.")
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
    result <- g2(data_subset_t)
    
    # Store the result
    results_list[[i]] <- result
  }
  
  # Bind all the results together into a single dataframe
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- paste(contrast_table[, 1], "vs", contrast_table[, 2])
 # adjust p-values for multiple-comparisons 
  results_df$adj_p.value <- MutSeqR::my.holm.sidak(results_df$p.value)
  
  return(results_df)
 }
##############################
  # G2ab
   # x = mut count
  # R = rows of counts w each row representing a different mutant site (subtype)
  # T = number of treatment groups under study (dose)
  # N = column totals? Total sample size?
  # G2ab <- function(x, R, T, N, monte.carlo = FALSE, n.sim = 10000, seed = 1234){
  # 
  #   flag <- x > 0
  #   S1 <- t(x[flag]) %*% log(x[flag])
  #   flag <- apply(x, 1, sum) > 0 # I added flags here - not sure if needed...
  #   S2 <- sum(apply(x, 1, sum)[flag] * log(apply(x, 1, sum)[flag])) # sum of each row
  #   flag <- apply(x, 2, sum) > 0
  #   S3 <- sum(apply(x, 2, sum)[flag] * log(apply(x, 2, sum)[flag])) # sum of each column
  #   flag <- x > 0
  #   S4 <- sum(x) * log(sum(x)) # total 
  # 
  #   G2ab <- 2*(S1 - S2 - S3 + S4)
  # 
  #   df <- (R-1)*(T-1)
  # 
  #   if(monte.carlo == FALSE){
  # 
  #     if(N/R > 20){
  #       p.value <- 1-stats::pchisq(G2ab, df)
  #     } else{
  #       p.value <- 1-stats::pf(G2ab/df, df, N-df)
  #     }
  #   } else {
  #     #Monte Carlo
  #     #Generate random rxc tables
  #     set.seed(seed)
  #     r <- apply(x, 1, sum)
  #     c <- apply(x, 2, sum)
  #     rtbl <- stats::r2dtable(1, r, c)
  # 
  #     ref.dist <- rep(0, n.sim)
  #     for(k in 1:length(ref.dist)){
  #       x <- stats::r2dtable(1, r, c)[[1]]
  #       flag <- x > 0
  #       S1 <- t(x[flag]) %*% log(x[flag])
  #       S2 <- sum(apply(x, 1, sum) * log(apply(x, 1, sum)))
  #       S3 <- sum(apply(x, 2, sum) * log(apply(x, 2, sum)))
  #       S4 <- sum(x) * log(sum(x))
  # 
  #       ref.dist[k] <- 2*(S1 - S2 - S3 + S4)
  # 
  #     }
  #     flag <- ref.dist >= G2ab[1,1]
  #     p.value <- length(ref.dist[flag])/10000
  #   }
  # 
  #   data.frame(G2ab = G2ab, p.value = p.value)
  # }