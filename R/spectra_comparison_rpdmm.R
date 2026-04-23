#' @title spectra_comparison_rpdmm
#' @description Ridge-penalized Dirichlet-Multinomial test for comparing
#' mutation spectra between two groups. Robust for large mutation counts.
#' @param mf_data A data frame containing the MF data. This
#' is the output from calculate_mf(). MF data should be calculate at the
#' *sample level* and  at the desired subtype resolution. Required columns are
#' sample, the exp_variable column(s), the subtype column, and sum_min or
#' sum_max.
#' @param exp_variable The column names of the experimental variable(s) to be
#' compared.
#' @param contrasts A filepath (character) OR a data frame specifying the
#' comparisons to be made. Must consist of exactly two columns. The level in
#' the first column will be compared to the level in the second column for
#' each row. If using multiple exp_variables, separate levels with a colon
#' (e.g., "Drug:High").
#' @param cont_sep Character. The delimiter used to import the contrasts
#' table if a filepath is provided. Default is tab.
#' @param mf_type Character. The type of mutation frequency count to use.
#' Choices  are "min" or "max". Default is "min" (recommended).
#' @param lambda Numeric. The strength of the ridge penalty. Default is 1.0.
#' Higher values increase shrinkage towards a uniform distribution,
#' stabilizing estimation in sparse datasets.
#' @param n_boot Integer. The number of bootstrap iterations to perform for
#' p-value calculation. Default is 500.
#' @returns A data frame containing one row per specified contrast. Columns
#' include the group names, comparison string, bootstrap p-value, observed
#' Likelihood Ratio Test (LRT) statistic, and robust overdispersion metrics
#' (theta) for group 1, group 2, and the shared null model.
#' @export
#' @importFrom dplyr select rename left_join
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames

spectra_comparison_rpdmm <- function(
  mf_data,
  exp_variable,
  contrasts,
  cont_sep = "\t",
  mf_type = "min",
  lambda = 1.0,
  n_boot = 500
) {

  # Validation & Setup
  stopifnot(
    !missing(mf_data) && is.data.frame(mf_data),
    "sample" %in% colnames(mf_data),
    !missing(exp_variable) && all(exp_variable %in% colnames(mf_data)),
    !missing(contrasts),
    is.numeric(lambda) && lambda >= 0,
    is.numeric(n_boot) && n_boot > 0
  )
  # Identify the appropriate count column based on mf_type
  mf_type <- match.arg(mf_type, choices = c("min", "max"))
  sum_col <- paste0("sum_", mf_type)
  mf_data$count <- mf_data[[sum_col]]

  # Identify Subtype Column
  potential_cols <- as.vector(subtype_dict)[!is.na(as.vector(subtype_dict))]
  subtype_cols <- intersect(colnames(mf_data), potential_cols)
  subtype_col <- subtype_cols[1]
  if (is.na(subtype_cols)) {
    stop("No valid subtype column found in mf_data.")
  } else if (length(subtype_cols) > 1) {
    warning("Multiple potential subtype columns found. Using: ", subtype_col)
  } else {
    mf_data <- dplyr::rename(mf_data, subtype = !!subtype_col)
  }

  # Concatenate Exp Variable columns into one.
  if (length(exp_variable) > 1) {
    mf_data$exp_var <- do.call(paste, c(mf_data[exp_variable], sep = ":"))
  } else {
    mf_data$exp_var <- mf_data[[exp_variable]]
  }
  # Exp variables = exp_var
  # Subtype column = subtype
  # Count column = count
  split_data <- split(mf_data, mf_data$exp_var)
  list_of_matrices <- lapply(split_data, function(sub_df) {
    wide_df <- sub_df %>%
      dplyr::select("subtype", "sample", "count") %>%
      tidyr::pivot_wider(
        names_from = "sample",
        values_from = "count",
        values_fill = 0, # Fills 0 if a sample lacks a specific subtype
        values_fn = sum  # Sums counts if a sample has multiples of a subtype
      ) %>%
      tibble::column_to_rownames(var = "subtype") # add row names
    as.matrix(wide_df)
  })

  # Load Contrasts
  if (is.data.frame(contrasts)) {
    contrast_table <- contrasts
  } else {
    contrast_table <- read.delim(contrasts, sep = cont_sep, header = FALSE)
  }
  # Validate the Contrasts Table
  if (ncol(contrast_table) != 2) {
    stop("Contrast table must have exactly 2 columns.")
  }
  all_groups <- names(list_of_matrices)
  requested_groups <- unique(unlist(contrast_table))
  missing_groups <- setdiff(requested_groups, all_groups)
  if (length(missing_groups) > 0) {
    stop("Groups in contrast table not found in data: ",
      paste(missing_groups, collapse = ", ")
    )
  }
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


#' @title dm_loglik_ridge_stable
#' @description Computes the robust, ridge-penalized log-likelihood for a
#' Dirichlet-Multinomial distribution. It safely handles massive counts
#' utilizing  lgamma components and matrix operations. The ridge penalty
#' shrinks probability vectors toward a uniform distribution to stabilize
#' optimization for sparse or heavily overdispersed data.
#' @param x A numeric matrix of mutational counts (subtypes as rows, samples
#' as columns).
#' @param p A numeric vector representing the expected probability of each
#' subtype.
#' @param theta Numeric scalar. The overdispersion parameter (concentration
#' parameter). Larger values of theta indicate lower overdispersion
#' (approaching multinomial).
#' @param lambda Numeric scalar. The strength of the ridge penalty.
#' Default is 1.0. Set to 0 for unpenalized (raw) likelihood calculation.
#' @returns A numeric scalar representing the computed log-likelihood value,
#' adjusted by the penalty term.
dm_loglik_ridge_stable <- function(x, p, theta, lambda = 1.0) {
  # --- 1. Force Numeric Matrix ---
  # This fixes the 'list' error by ensuring x is a pure numeric matrix
  x <- as.matrix(x)
  if (!is.numeric(x)) {
    # If the matrix contains strings, try to convert them
    class(x) <- "numeric"
  }
  # --- 2. Safety Constraints ---
  theta <- min(max(theta, 1e-4), 400000)
  p <- pmax(p, 1e-12)
  p <- p / sum(p)

  alpha <- p * theta
  N_j <- colSums(x)
  K <- length(p)

  # --- 3. Log-Likelihood Calculation ---
  # Term 1: Denominator (Sample Depth Component)
  term_denominator <- sum(lgamma(theta) - lgamma(N_j + theta))

  # Term 2: Numerator (Counts Component)
  alpha_mat <- matrix(alpha, nrow = K, ncol = ncol(x))

  # We use as.matrix here to double-ensure we don't return a list structure
  numerator_mat <- as.matrix(lgamma(x + alpha_mat) - lgamma(alpha_mat))

  # Check for finiteness (safety fallback for optimization steps)
  if (!all(is.finite(numerator_mat))) {
    return(-1e10)
  }

  ll_dm <- term_denominator + sum(numerator_mat)

  # --- 4. Scaled Ridge Penalty ---
  N_total <- sum(x)
  penalty <- (lambda * N_total) * sum((p - (1 / K))^2)

  return(ll_dm - penalty)
}

#' @title rdm_fast
#' @description Rapidly simulates count matrices drawn from a
#' Dirichlet-Multinomial distribution. This function is highly optimized for
#' parametric bootstrapping, utilizing a Gamma-Multinomial approximation for
#' speed.
#' @param depths A numeric vector representing the total sequencing depth
#' (total mutation counts) for each sample to be simulated.
#' @param p A numeric vector representing the expected probability of each
#' subtype.
#' @param theta Numeric scalar. The overdispersion parameter (concentration
#' parameter).
#' @returns A numeric matrix of simulated mutation counts, with subtypes as
#' rows and simulated samples as columns.
#' @importFrom stats rgamma rmultinom
rdm_fast <- function(depths, p, theta) {
  alpha <- p * theta
  K <- length(p)
  gammas <- matrix(
    stats::rgamma(K * length(depths), shape = alpha, rate = 1),
    nrow = K
  )
  probs <- sweep(gammas, 2, colSums(gammas), "/")
  return(
    sapply(
      seq_along(depths),
      function(j) stats::rmultinom(1, size = depths[j], prob = probs[, j])
    )
  )
}
#' @title run_penalized_comparison
#' @description Core statistical engine for the RP-DMM test. Evaluates whether
#' two count matrices represent significantly different mutational spectra. It
#' utilizes a "hybrid logic": robust parameters (thetas) are estimated utilizing
#' a ridge penalty (QA phase), while the test statistic (LRT) is evaluated using
#' the  unpenalized likelihood. Significance is evaluated via parametric
#' bootstrap.
#' @param g1 A numeric count matrix for Group 1 (subtypes as rows, samples as
#' columns).
#' @param g2 A numeric count matrix for Group 2 (subtypes as rows, samples as
#' columns).
#' @param lambda Numeric. The ridge penalty strength used during parameter
#' estimation.
#' @param n_boot Integer. The number of parametric bootstrap iterations to
#' perform.
#' @returns A list containing 5 items:
#' \itemize{
#'   \item \code{p_value}: Numeric. The bootstrapped p-value.
#'   \item \code{lrt}: Numeric. The observed unpenalized Likelihood Ratio Test
#' statistic.
#'   \item \code{theta_g1}: Numeric. The optimized robust theta
#' (overdispersion) for Group 1.
#'   \item \code{theta_g2}: Numeric. The optimized robust theta for Group 2.
#'   \item \code{theta_shared}: Numeric. The optimized robust theta under the
#' Null hypothesis.
#' }
#' @importFrom stats optimize
run_penalized_comparison <- function(g1, g2, lambda = 1.0, n_boot = 500) {
  # 1. Data Setup
  all_data <- cbind(g1, g2)
  depths1 <- colSums(g1)
  depths2 <- colSums(g2)
  p_null <- rowSums(all_data) / sum(all_data)
  p_g1 <- rowSums(g1) / sum(g1)
  p_g2 <- rowSums(g2) / sum(g2)
  safe_bounds <- c(1e-2, 1000000)

  # --- 2. ROBUST ESTIMATION (The QA Phase) ---
  # We estimate the individual group thetas first.
  # These serve as our QA metrics.
  opt_n_pen <- stats::optimize(
    f = function(th) dm_loglik_ridge_stable(all_data, p_null, th, lambda),
    interval = safe_bounds,
    maximum = TRUE
  )
  opt_1_pen <- stats::optimize(
    f = function(th) dm_loglik_ridge_stable(g1, p_g1, th, lambda),
    interval = safe_bounds,
    maximum = TRUE
  )
  opt_2_pen <- stats::optimize(
    f = function(th) dm_loglik_ridge_stable(g2, p_g2, th, lambda),
    interval = safe_bounds,
    maximum = TRUE
  )

  # Capture the individual group metrics
  theta_g1 <- opt_1_pen$maximum
  theta_g2 <- opt_2_pen$maximum
  theta_null <- opt_n_pen$maximum

  # --- 3. OBSERVED STATISTIC (The Testing Phase) ---
  # Evaluate using UNPENALIZED likelihood (Hybrid Logic)
  ll_null_raw <- dm_loglik_ridge_stable(
    all_data, p_null,
    theta_null,
    lambda = 0
  )
  ll_alt_raw  <- dm_loglik_ridge_stable(g1, p_g1, theta_g1, lambda = 0) +
    dm_loglik_ridge_stable(g2, p_g2, theta_g2, lambda = 0)

  obs_lrt <- max(0, 2 * (ll_alt_raw - ll_null_raw))

  # --- 4. BOOTSTRAP LOOP ---
  boot_stats <- numeric(n_boot)
  for (i in 1:n_boot) {
    # Simulate Null Data using the shared robust theta
    s1 <- rdm_fast(depths1, p_null, theta_null)
    s2 <- rdm_fast(depths2, p_null, theta_null)
    s_all <- cbind(s1, s2)
    ps1 <- rowSums(s1) / sum(s1)
    ps2 <- rowSums(s2) / sum(s2)
    psN <- rowSums(s_all) / sum(s_all)

    # Estimate bootstrapped parameters with penalty
    sn_p <- stats::optimize(
      f = function(th) dm_loglik_ridge_stable(s_all, psN, th, lambda),
      interval = safe_bounds,
      maximum = TRUE
    )
    s1_p  <- stats::optimize(
      f = function(th) dm_loglik_ridge_stable(s1, ps1, th, lambda),
      interval = safe_bounds,
      maximum = TRUE
    )
    s2_p  <- stats::optimize(
      f = function(th) dm_loglik_ridge_stable(s2, ps2, th, lambda),
      interval = safe_bounds,
      maximum = TRUE
    )

    # Calculate bootstrap LRT (Hybrid logic)
    b_ll_n   <- dm_loglik_ridge_stable(s_all, psN, sn_p$maximum, lambda = 0)
    b_ll_alt <- dm_loglik_ridge_stable(s1, ps1, s1_p$maximum, lambda = 0) +
      dm_loglik_ridge_stable(s2, ps2, s2_p$maximum, lambda = 0)

    boot_stats[i] <- max(0, 2 * (b_ll_alt - b_ll_n))
  }

  # --- 5. RETURN COMPREHENSIVE RESULTS ---
  return(list(
    p_value = mean(boot_stats >= (obs_lrt - 1e-8)),
    lrt = obs_lrt,
    theta_g1 = theta_g1,     # QA Metric for Group 1
    theta_g2 = theta_g2,     # QA Metric for Group 2
    theta_shared = theta_null # Overall consistency under Null
  ))
}

#' @title select_optimal_lambda
#' @description Performs leave-one-out cross-validation (LOOCV) to empirically
#' select the optimal ridge penalty (lambda) for a given count matrix. Models
#' are fit on $N-1$ samples using a specified lambda, and unpenalized
#' likelihood is evaluated on the held-out sample.
#' @param count_matrix A numeric matrix of mutation counts with subtypes
#' as rows and samples as columns.
#' @param lambdas A numeric vector of lambda penalty values to evaluate.
#' Default is c(0.1, 0.5, 1, 5, 10, 50, 100).
#' @returns A list containing two items:
#' \itemize{
#'   \item \code{best_lambda}: Numeric scalar. The lambda value from the
#' provided vector that yielded the highest overall cross-validation
#' log-likelihood.
#'   \item \code{scores}: A named numeric vector detailing the cumulative
#' log-likelihood score achieved by each tested lambda.
#' }
#' @export
#' @importFrom stats optimize
select_optimal_lambda <- function(
  count_matrix,
  lambdas = c(0.1, 0.5, 1, 5, 10, 50, 100)
) {
  n_samples <- ncol(count_matrix)
  cv_scores <- numeric(length(lambdas))
  names(cv_scores) <- as.character(lambdas)
  K <- nrow(count_matrix)

  for (i in seq_along(lambdas)) {
    lam <- lambdas[i]
    loglik_sum <- 0
    for (j in 1:n_samples) {
      test_sample <- count_matrix[, j, drop = FALSE]
      train_data <- count_matrix[, -j, drop = FALSE]
      # Estimate proportions with shrinkage
      raw_counts <- rowSums(train_data)
      p_shrunk <- (raw_counts + lam) / (sum(raw_counts) + (K * lam))
      # Optimize theta on training set
      opt <- stats::optimize(
        f = function(th) dm_loglik_ridge_stable(train_data, p_shrunk, th, lambda = lam),
        interval = c(1, 100000),
        maximum = TRUE
      )
      # Evaluate unpenalized likelihood on held-out sample
      loglik_sum <- loglik_sum +
        dm_loglik_ridge_stable(test_sample, p_shrunk, opt$maximum, lambda = 0)
    }
    cv_scores[i] <- loglik_sum
  }
  return(list(best_lambda = lambdas[which.max(cv_scores)], scores = cv_scores))
}