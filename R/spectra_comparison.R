#' Compare the overall mutation spectra between groups
#'
#' spectra_comparison compares the mutation spectra of groups using a
#' modified contingency table approach.
#' @param mutation_data A data frame containing the mutation data. This
#' is the output from import_mut_data or import_vcf_data.
#' @param cols_to_group A character vector of the column names in the mutation
#' data that you want to group by. Ex. c("dose", "tissue"). This
#' function will sum the mutations across groups before running the comparison.
#' @param subtype_resolution The resolution of the mutation spectra to be
#' compared. Options include "base_6", "base_12", "base_96", and "base_192" and "type".
#' See calculate_mf for more details.
#' @param variant_types A character vector of the mutation types to include
#' in the comparison. Default is all types of mutations.
#' @param contrasts a filepath to a tab-delimited `.txt` file OR a dataframe that
#' will specify the comparisons to be made between groups. The table must
#' consist of two columns. The first column will be a level within your group
#' column and the second column must be the group level that it will be
#' compared to. All values must correspond to entries in your `cols_to_group`
#' column. For more than one `group` variable, separate the levels of
#' each group with a colon. Ensure that all `groups` listed in cols_to_group
#' are represented in each entry for the table. See details for examples.
#' @param cont_sep The delimiter used to import the contrasts table.
#' Default is tab.
#' @param mf_type The type of mutation frequency to use. Default is "min".
#' @returns the log-likelihood statistic G2 for the specified comparisons with
#' the p-value adjusted for multiple-comparisons.
#' @export
#'
#' @details
#' Examples of `contrasts`:
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

spectra_comparison <- function(mutation_data,
                               cols_to_group,
                               subtype_resolution = "base_6",
                               variant_types = c("snv",
                                                 "deletion",
                                                 "insertion",
                                                 "complex",
                                                 "mnv",
                                                 "sv",
                                                 "ambiguous",
                                                 "uncategorized"),
                               mf_type = "min",
                               contrasts,
                               cont_sep = "\t") {
  
  mf_data <- MutSeqR::calculate_mf(mutation_data = mutation_data,
                                         cols_to_group = cols_to_group,
                                         subtype_resolution = subtype_resolution,
                                         variant_types = variant_types,
                                         summary = TRUE)
  # Prepare Data
  sum_col <- paste0("sum_", mf_type)
  ## Find the subtype column
  subtype_col <- MutSeqR::subtype_dict[[subtype_resolution]]
  ## Select the necessary columns
  mut_spectra <- mf_data %>%
    dplyr::select(dplyr::all_of(c(cols_to_group, subtype_col, sum_col)))
  ## Create a single group column
  mut_spectra <- mut_spectra %>%
    dplyr::mutate(group_col = do.call(paste,
                                      c(dplyr::select(mut_spectra,
                                                      dplyr::all_of(cols_to_group)),
                                        sep = ":"))) %>%
    dplyr::select(-all_of(cols_to_group))
  
  # All groups
  groups <- unique(mut_spectra$group_col)
  # Extract data for each group
  filtered_data <- list()
  for (i in seq_along(groups)) {
    group = groups[i]
    df_i <- mut_spectra %>%
      filter(group_col == group)
    filtered_data[[i]] <- df_i
  }
 
  # G2 Statistic - Likelihood Ratio Statistic
  ## Piegorsch and Bailer 1994 doi: 10.1093/genetics/136.1.403.
  G2 <- function(x, monte.carlo = FALSE, n.sim = 10000, seed = 1234){
    N <- sum(x)
    r <- apply(x, 1, sum)
    c <- apply(x, 2, sum)
    e <- r %*% t(c)/N
    G2 <- 0
    for(k in 1:ncol(x)){
      flag <- x[,k] > 0
      G2 <- G2 + t(x[flag,k]) %*% log(x[flag,k]/e[flag,k])	
    }
    G2 <- 2*G2
    R <- nrow(x)-1
    df <- R * (ncol(x)-1)
    if(monte.carlo == FALSE){
      if(N/R > 20){
        p.value <- 1-pchisq(G2, df)
        message("Using chi-squared distribution to compute p-value")
      } else{
        p.value <- 1-pf(G2/R, R, N-df)
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
      for(k in 1:length(ref.dist)){
        x <- r2dtable(1, r, c)[[1]]
        N <- sum(x)
        r <- apply(x, 1, sum)
        c <- apply(x, 2, sum)
        e <- r %*% t(c)/N
        G2.t <- 0
        for(j in 1:ncol(x)){
          flag <- x[,j] > 0
          G2.t <- G2.t + t(x[flag,j]) %*% log(x[flag,j]/e[flag,j])	
        }
        ref.dist[k] <- 2*G2.t
      } 
      flag <- ref.dist >= G2[1,1]
      p.value <- length(ref.dist[flag])/10000
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
    df_combined <- do.call(cbind, lapply(dfs, function(df) df[,2]))
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
    results <- rbind(results, data.frame(contrasts = contrast_str, G2 = result$G2, p.value = result$p.value, stringsAsFactors = FALSE))
  }

  # Apply the Holm-Sidak correction for multiple comparisons
  results$adjP <- MutSeqR::sidak(results$p.value)$SidakP
  results$sign <- ""
 results$sign[results$adjP < 0.05] <- "***"
  
  return(results)
}
