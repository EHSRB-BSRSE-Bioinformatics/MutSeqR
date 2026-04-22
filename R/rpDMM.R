
run_rank_permutation_test <- function(g1, g2, n_perm = 1000) {
  combined_data <- cbind(g1, g2)
  n1 <- ncol(g1); n_total <- ncol(combined_data)
  ranks_matrix <- apply(combined_data, 2, function(x) rank(x, ties.method = "average"))
  
  get_rank_dist <- function(r_mat, idx1) {
    mean_r1 <- rowMeans(r_mat[, idx1, drop = FALSE])
    mean_r2 <- rowMeans(r_mat[, -idx1, drop = FALSE])
    return(sum((mean_r1 - mean_r2)^2))
  }
  
  obs_dist <- get_rank_dist(ranks_matrix, 1:n1)
  perm_stats <- replicate(n_perm, get_rank_dist(ranks_matrix, sample(1:n_total, n1)))
  return(list(p_value = mean(perm_stats >= (obs_dist - 1e-8)), obs_dist = obs_dist))
}

# ------------------------------------------------------------------------------
# 1. Stabilized Ridge-Penalized Dirichlet-Multinomial Log-Likelihood
# ------------------------------------------------------------------------------
dm_loglik_ridge_stable <- function(x, p, theta, lambda = 1.0) {
  # --- 1. Force Numeric Matrix ---
  # This fixes the 'list' error by ensuring x is a pure numeric matrix
  x <- as.matrix(x)
  if(!is.numeric(x)) {
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
  penalty <- (lambda * N_total) * sum((p - (1/K))^2)
  
  return(ll_dm - penalty)
}

# ==============================================================================
# Hybrid RP-DMM Test Function with QA/QC Metrics
# ==============================================================================
run_penalized_comparison <- function(g1, g2, lambda = 1.0, n_boot = 500) {
  
  # 1. Data Setup
  all_data <- cbind(g1, g2)
  depths1 <- colSums(g1); depths2 <- colSums(g2)
  p_null <- rowSums(all_data) / sum(all_data)
  p_g1 <- rowSums(g1) / sum(g1); p_g2 <- rowSums(g2) / sum(g2)
  safe_bounds <- c(1e-2, 1000000)
  
  # --- 2. ROBUST ESTIMATION (The QA Phase) ---
  # We estimate the individual group thetas first. 
  # These serve as our QA metrics.
  opt_n_pen <- optimize(f = function(th) dm_loglik_ridge_stable(all_data, p_null, th, lambda), interval = safe_bounds, maximum = TRUE)
  opt_1_pen <- optimize(f = function(th) dm_loglik_ridge_stable(g1, p_g1, th, lambda), interval = safe_bounds, maximum = TRUE)
  opt_2_pen <- optimize(f = function(th) dm_loglik_ridge_stable(g2, p_g2, th, lambda), interval = safe_bounds, maximum = TRUE)
  
  # Capture the individual group metrics
  theta_g1 <- opt_1_pen$maximum
  theta_g2 <- opt_2_pen$maximum
  theta_null <- opt_n_pen$maximum
  
  # --- 3. OBSERVED STATISTIC (The Testing Phase) ---
  # Evaluate using UNPENALIZED likelihood (Hybrid Logic)
  ll_null_raw <- dm_loglik_ridge_stable(all_data, p_null, theta_null, lambda = 0)
  ll_alt_raw  <- dm_loglik_ridge_stable(g1, p_g1, theta_g1, lambda = 0) + 
    dm_loglik_ridge_stable(g2, p_g2, theta_g2, lambda = 0)
  
  obs_lrt <- max(0, 2 * (ll_alt_raw - ll_null_raw))
  
  # --- 4. BOOTSTRAP LOOP ---
  boot_stats <- numeric(n_boot)
  for(i in 1:n_boot) {
    # Simulate Null Data using the shared robust theta
    s1 <- rdm_fast(depths1, p_null, theta_null)
    s2 <- rdm_fast(depths2, p_null, theta_null)
    s_all <- cbind(s1, s2)
    ps1 <- rowSums(s1)/sum(s1); ps2 <- rowSums(s2)/sum(s2); psN <- rowSums(s_all)/sum(s_all)
    
    # Estimate bootstrapped parameters with penalty
    sn_p  <- optimize(f=function(th) dm_loglik_ridge_stable(s_all, psN, th, lambda), interval=safe_bounds, maximum=TRUE)
    s1_p  <- optimize(f=function(th) dm_loglik_ridge_stable(s1, ps1, th, lambda), interval=safe_bounds, maximum=TRUE)
    s2_p  <- optimize(f=function(th) dm_loglik_ridge_stable(s2, ps2, th, lambda), interval=safe_bounds, maximum=TRUE)
    
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

# --- Helper 2: Simulation for Bootstrap ---
rdm_fast <- function(depths, p, theta) {
  alpha <- p * theta
  K <- length(p)
  gammas <- matrix(rgamma(K * length(depths), shape = alpha, rate = 1), nrow = K)
  probs <- sweep(gammas, 2, colSums(gammas), "/")
  return(sapply(seq_along(depths), function(j) rmultinom(1, size = depths[j], prob = probs[, j])))
}

select_optimal_lambda <- function(count_matrix, lambdas = c(0.1, 0.5, 1, 5, 10, 50, 100)) {
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
      opt <- optimize(f = function(th) dm_loglik_ridge_stable(train_data, p_shrunk, th, lambda = lam), 
                      interval = c(1, 100000), maximum = TRUE)
      # Evaluate unpenalized likelihood on held-out sample
      loglik_sum <- loglik_sum + dm_loglik_ridge_stable(test_sample, p_shrunk, opt$maximum, lambda = 0)
    }
    cv_scores[i] <- loglik_sum
  }
  return(list(best_lambda = lambdas[which.max(cv_scores)], scores = cv_scores))
}

