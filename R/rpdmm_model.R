#' @title dm_loglik_ridge_stable
#' @keywords internal
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

  ll_dm <- term_denominator + sum(numerator_mat)

  # --- 4. Scaled Ridge Penalty ---
  N_total <- sum(x)
  penalty <- (lambda * N_total) * sum((p - (1 / K))^2)

  return(ll_dm - penalty)
}

#' @title softmax_stable
#' @description Converts logits to probabilities using a numerically stable
#' softmax calculation.
#' @param phi Numeric vector of category logits.
#' @returns A numeric probability vector that sums to one.
#' @keywords internal
softmax_stable <- function(phi) {
  shifted <- phi - max(phi)
  exponentials <- exp(shifted)
  exponentials / sum(exponentials)
}

#' @title extract_params
#' @description Converts optimizer parameters to category probabilities and
#' the concentration parameter.
#' @param pars Numeric optimizer vector containing log-theta followed by the
#' free category logits.
#' @param category_names Optional names to assign to the probabilities.
#' @returns A list containing the named probability vector \code{p} and
#' concentration parameter \code{theta}.
#' @keywords internal
extract_params <- function(pars, category_names = NULL) {
  phi <- c(pars[-1], 0)
  p <- softmax_stable(phi)
  names(p) <- category_names
  list(p = p, theta = exp(pars[1]))
}

#' @title joint_dm_target
#' @description Returns the negative ridge-penalized Dirichlet-Multinomial
#' log-likelihood for joint parameter optimization.
#' @param pars Numeric optimizer vector containing log-theta and free logits.
#' @param x Numeric count matrix with categories in rows and samples in columns.
#' @param lambda Non-negative ridge penalty strength.
#' @returns A numeric negative penalized log-likelihood.
#' @keywords internal
joint_dm_target <- function(pars, x, lambda) {
  params <- extract_params(pars)
  -dm_loglik_ridge_stable(x, params$p, params$theta, lambda)
}

#' @title joint_dm_gradient
#' @description Computes the gradient of the negative ridge-penalized
#' Dirichlet-Multinomial log-likelihood.
#' @param pars Numeric optimizer vector containing log-theta and free logits.
#' @param x Numeric count matrix with categories in rows and samples in columns.
#' @param lambda Non-negative ridge penalty strength.
#' @returns A numeric gradient matching the order of \code{pars}.
#' @keywords internal
joint_dm_gradient <- function(pars, x, lambda) {
  K <- nrow(x)
  theta <- exp(pars[1])
  phi <- c(pars[-1], 0)
  p <- softmax_stable(phi)
  alpha <- p * theta
  N_j <- colSums(x)

  g_alpha <- rowSums(digamma(x + alpha) - digamma(alpha))

  d_ll_dtheta_direct <- sum(digamma(theta) - digamma(N_j + theta))
  d_ll_dtheta_alpha <- sum(p * g_alpha)
  d_ll_dtheta <- d_ll_dtheta_direct + d_ll_dtheta_alpha
  d_ll_d_logtheta <- theta * d_ll_dtheta # d_ll_d_logtheta

  weighted_g <- sum(p * g_alpha)
  d_ll_d_phi <- theta * p * (g_alpha - weighted_g)
  d_ll_d_phi <- d_ll_d_phi[-K]

  total_counts <- sum(x)
  u <- 1 / K
  weighted_penalty <- sum(p * (p - u))
  d_penalty_d_phi <- 2 * lambda * total_counts * p * ((p - u) - weighted_penalty)
  d_penalty_d_phi <- d_penalty_d_phi[-K]

  c(-d_ll_d_logtheta, -d_ll_d_phi + d_penalty_d_phi)
}

#' @title fit_dmm_fast
#' @description Fits a ridge-penalized Dirichlet-Multinomial model using
#' bounded L-BFGS-B optimization.
#' @param x Numeric count matrix with categories in rows and samples in columns.
#' @param init_pars Numeric starting values for log-theta and free logits.
#' @param lambda Non-negative ridge penalty strength.
#' @returns The optimization result, or \code{NULL} if optimization fails.
#' @keywords internal
#' @importFrom stats optim
fit_dmm_fast <- function(x, init_pars, lambda) {
  lower <- c(log(1e-4), rep(-20, nrow(x) - 1))
  upper <- c(log(400000), rep(20, nrow(x) - 1))
  init_pars <- pmax(lower, pmin(upper, init_pars))
  result <- tryCatch(
    stats::optim(
      par = init_pars,
      fn = joint_dm_target,
      gr = joint_dm_gradient,
      x = x,
      lambda = lambda,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = 2000, factr = 1e7, pgtol = 1e-7)
    ),
    error = function(error) {
      message("Optimization error: ", conditionMessage(error))
      NULL
    }
  )
  if (is.null(result)) return(NULL)
  if (result$convergence != 0 || !is.finite(result$value)) {
    message("Optimization did not converge. Code: ", result$convergence,
            "; message: ", result$message)
    return(NULL)
  }
  result
}

#' @title rdm_fast
#' @keywords internal
#' @description Rapidly simulates count matrices drawn from a
#' Dirichlet-Multinomial distribution. This function is highly optimized for
#' parametric bootstrapping, utilizing log-space Gamma draws followed by
#' multinomial sampling to remain stable for tiny concentration parameters.
#' @param depths A numeric vector representing the total sequencing depth
#' (total mutation counts) for each sample to be simulated.
#' @param p A numeric vector representing the expected probability of each
#' subtype.
#' @param theta Numeric scalar. The overdispersion parameter (concentration
#' parameter).
#' @returns A numeric matrix of simulated mutation counts, with subtypes as
#' rows and simulated samples as columns.
#' @importFrom stats rgamma rmultinom runif
rdm_fast <- function(depths, p, theta) {
  categories <- length(p)
  samples <- length(depths)
  alpha <- pmax(p * theta, 1e-300)
  log_gamma <- log(matrix(
    stats::rgamma(categories * samples, shape = alpha + 1), categories, samples
  )) + log(matrix(stats::runif(categories * samples), categories, samples)) / alpha
  log_gamma <- sweep(log_gamma, 2, apply(log_gamma, 2, max), "-")
  p <- exp(log_gamma)
  p <- sweep(p, 2, colSums(p), "/")
  matrix(
    vapply(seq_len(samples), function(sample_index) {
      as.numeric(stats::rmultinom(1, depths[sample_index],
                                 p[, sample_index]))
    }, numeric(categories)),
    nrow = categories, ncol = samples
  )
}
#' @title run_penalized_comparison
#' @keywords internal
#' @description Core statistical engine for the RP-DMM test. Evaluates whether
#' two count matrices represent significantly different mutational spectra. It
#' jointly estimates category probabilities and concentration parameters using
#' a ridge penalty. Both unpenalized (hybrid) and penalized LRT statistics
#' are calibrated by parametric bootstrap. Failed bootstrap fits are excluded;
#' at least half the requested replicates must succeed.
#' @param g1 A numeric count matrix for Group 1 (subtypes as rows, samples as
#' columns).
#' @param g2 A numeric count matrix for Group 2 (subtypes as rows, samples as
#' columns).
#' @param lambda Numeric. The ridge penalty strength used during parameter
#' estimation.
#' @param n_boot Integer. The number of parametric bootstrap iterations to
#' perform.
#' @returns A list containing:
#' \itemize{
#'   \item \code{p_value}: Numeric. The hybrid bootstrap p-value.
#'   \item \code{lrt}: Numeric. The observed unpenalized Likelihood Ratio Test
#' statistic.
#'   \item \code{p_penalized}: Numeric. The penalized bootstrap p-value.
#'   \item \code{lrt_penalized}: Numeric. The observed penalized LRT statistic.
#'   \item \code{theta_g1}: Numeric. The optimized robust theta
#' (overdispersion) for Group 1.
#'   \item \code{theta_g2}: Numeric. The optimized robust theta for Group 2.
#'   \item \code{theta_shared}: Numeric. The optimized robust theta under the
#' Null hypothesis.
#' }
run_penalized_comparison <- function(g1, g2, lambda = 1.0, n_boot = 500) {
  stopifnot(is.numeric(n_boot), length(n_boot) == 1, is.finite(n_boot),
            n_boot > 0, n_boot == floor(n_boot))
  all_data <- cbind(g1, g2)
  categories <- nrow(all_data)
  category_names <- rownames(all_data)
  depths1 <- colSums(g1)
  depths2 <- colSums(g2)
  failed_result <- list(
    p_value = NA_real_, p_penalized = NA_real_,
    lrt = NA_real_, lrt_penalized = NA_real_,
    theta_shared = NA_real_, theta_g1 = NA_real_, theta_g2 = NA_real_
  )
  if (sum(g1) == 0 || sum(g2) == 0) return(failed_result)
  observed <- pmax(rowSums(all_data) / sum(all_data), 1e-6)
  observed <- observed / sum(observed)
  init_pars <- c(log(100), log(observed[-categories] / observed[categories]))
  fit_null <- fit_dmm_fast(all_data, init_pars, lambda)
  fit_g1 <- fit_dmm_fast(g1, init_pars, lambda)
  fit_g2 <- fit_dmm_fast(g2, init_pars, lambda)
  if (is.null(fit_null) || is.null(fit_g1) || is.null(fit_g2)) {
    return(failed_result)
  }
  null_parameters <- extract_params(fit_null$par, category_names)
  parameters_g1 <- extract_params(fit_g1$par, category_names)
  parameters_g2 <- extract_params(fit_g2$par, category_names)
  raw_null <- dm_loglik_ridge_stable(
    all_data, null_parameters$p, null_parameters$theta, 0
  )
  raw_alt <- dm_loglik_ridge_stable(g1, parameters_g1$p, parameters_g1$theta, 0) +
    dm_loglik_ridge_stable(g2, parameters_g2$p, parameters_g2$theta, 0)
  observed_hybrid <- max(0, 2 * (raw_alt - raw_null))
  observed_penalized <- max(0, 2 * (fit_null$value - fit_g1$value - fit_g2$value))
  boot_hybrid <- boot_penalized <- numeric(n_boot)
  valid_boots <- 0L
  for (iteration in seq_len(n_boot)) {
    simulated_g1 <- rdm_fast(depths1, null_parameters$p, null_parameters$theta)
    simulated_g2 <- rdm_fast(depths2, null_parameters$p, null_parameters$theta)
    simulated_all <- cbind(simulated_g1, simulated_g2)
    rownames(simulated_g1) <- rownames(simulated_g2) <-
      rownames(simulated_all) <- category_names
    boot_null <- fit_dmm_fast(simulated_all, fit_null$par, lambda)
    boot_g1 <- fit_dmm_fast(simulated_g1, fit_g1$par, lambda)
    boot_g2 <- fit_dmm_fast(simulated_g2, fit_g2$par, lambda)
    if (is.null(boot_null) || is.null(boot_g1) || is.null(boot_g2)) next
    valid_boots <- valid_boots + 1L
    boot_parameters_null <- extract_params(boot_null$par, category_names)
    boot_parameters_g1 <- extract_params(boot_g1$par, category_names)
    boot_parameters_g2 <- extract_params(boot_g2$par, category_names)
    boot_raw_null <- dm_loglik_ridge_stable(
      simulated_all, boot_parameters_null$p, boot_parameters_null$theta, 0
    )
    boot_raw_alt <- dm_loglik_ridge_stable(
      simulated_g1,
      boot_parameters_g1$p,
      boot_parameters_g1$theta,
      0
    ) + dm_loglik_ridge_stable(
          simulated_g2,
          boot_parameters_g2$p,
          boot_parameters_g2$theta,
          0
    )
    boot_hybrid[valid_boots] <- max(0, 2 * (boot_raw_alt - boot_raw_null))
    # 2(l alt - l null) = 2*(-g1 + -g2 - -null) (please keep this comment for clarity during review)
    boot_penalized[valid_boots] <- max(
      0, 2 * (boot_null$value - boot_g1$value - boot_g2$value)
    )
  }
  if (valid_boots < n_boot * 0.5) return(failed_result)
  result <- list(
    p_value = (1 + sum(boot_hybrid[seq_len(valid_boots)] >= observed_hybrid)) /
      (valid_boots + 1),
    p_penalized = (1 + sum(boot_penalized[seq_len(valid_boots)] >=
                           observed_penalized)) / (valid_boots + 1),
    lrt = observed_hybrid, lrt_penalized = observed_penalized,
    theta_shared = null_parameters$theta,
    theta_g1 = parameters_g1$theta, theta_g2 = parameters_g2$theta
  )
  result
}

#' @title select_optimal_lambda
#' @description Performs leave-one-out cross-validation (LOOCV) to empirically
#' select the optimal ridge penalty (lambda) for a given count matrix. Models
#' are fit on $N-1$ samples using a specified lambda, and unpenalized
#' likelihood is evaluated on the held-out sample.
#' @param mf_data A data frame of sample-level mutation-frequency data,
#' typically returned by \code{calculate_mf()}. It must contain a sample
#' column, a recognized subtype column, and the selected sum-count column.
#' @param mf_type Character. The sum-count column to use: \code{"min"} or
#' \code{"max"}. Defaults to \code{"min"}.
#' @param lambdas A numeric vector of lambda penalty values to evaluate.
#' Default is c(0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 5).
#' @returns A list containing three items:
#' \itemize{
#'   \item \code{best_lambda}: Numeric scalar. The lambda value from the
#' provided vector that yielded the highest overall cross-validation
#' log-likelihood.
#'   \item \code{scores}: A named numeric vector detailing the cumulative
#' log-likelihood score achieved by each tested lambda.
#'   \item \code{lambdas}: The evaluated candidate penalties. Failed fits
#' produce \code{NA} scores; if all candidates fail, \code{best_lambda}
#' is \code{NA_real_}.
#' }
#' @export
select_optimal_lambda <- function(
  mf_data,
  mf_type = "min",
  lambdas = c(0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 5)
) {
  stopifnot(is.data.frame(mf_data), "sample" %in% colnames(mf_data))
  mf_type <- match.arg(mf_type, choices = c("min", "max"))
  sum_col <- paste0("sum_", mf_type)
  if (!sum_col %in% colnames(mf_data)) {
    stop("mf_data must contain the selected count column: ", sum_col)
  }
  potential_cols <- as.vector(subtype_dict)[!is.na(as.vector(subtype_dict))]
  subtype_cols <- intersect(colnames(mf_data), potential_cols)
  if (length(subtype_cols) == 0) {
    stop("No valid subtype column found in mf_data.")
  }
  if (length(subtype_cols) > 1) {
    warning("Multiple potential subtype columns found. Using: ", subtype_cols[1])
  }
  count_matrix <- data.frame(
    subtype = mf_data[[subtype_cols[1]]],
    sample = mf_data$sample,
    count = mf_data[[sum_col]]
  ) %>%
    dplyr::select("subtype", "sample", "count") %>%
    tidyr::pivot_wider(
      names_from = "sample",
      values_from = "count",
      values_fill = 0,
      values_fn = sum
    ) %>%
    tibble::column_to_rownames(var = "subtype") %>%
    as.matrix()
  stopifnot(is.numeric(count_matrix),
            nrow(count_matrix) >= 1, ncol(count_matrix) >= 2,
            all(is.finite(count_matrix)), all(count_matrix >= 0),
            all(colSums(count_matrix) > 0),
            is.numeric(lambdas), length(lambdas) > 0,
            all(is.finite(lambdas)), all(lambdas >= 0))
  n_samples <- ncol(count_matrix)
  cv_scores <- numeric(length(lambdas))
  names(cv_scores) <- as.character(lambdas)
  K <- nrow(count_matrix)

  for (i in seq_along(lambdas)) {
    lam <- lambdas[i]
    loglik_sum <- 0
    for (j in seq_len(n_samples)) {
      test_sample <- count_matrix[, j, drop = FALSE]
      train_data <- count_matrix[, -j, drop = FALSE]
      # Initial guesses for training data
      observed <- pmax(rowSums(train_data) / sum(train_data), 1e-6)
      observed <- observed / sum(observed)
      init_pars <- c(log(100), log(observed[-K] / observed[K]))
      # Joint optimization on TRAINING data under current lambda
      fit_train <- fit_dmm_fast(train_data, init_pars, lam)
      if (is.null(fit_train)) {
        loglik_sum <- NA_real_
        break
      }
      # Extract trained parameters
      params <- extract_params(fit_train$par, rownames(count_matrix))
      # Evaluate UNPENALIZED log-likelihood on TEST sample
      # We want to know how well the model predicts held-out data
      loglik_sum <- loglik_sum + dm_loglik_ridge_stable(
        test_sample, params$p, params$theta, lambda = 0
      ) # unpenalized evaluation on test data
    }
    cv_scores[i] <- loglik_sum
  }
  best <- which.max(replace(cv_scores, !is.finite(cv_scores), -Inf))
  best_lambda <- if (any(is.finite(cv_scores))) lambdas[best] else NA_real_
  list(best_lambda = best_lambda, scores = cv_scores, lambdas = lambdas)
}

