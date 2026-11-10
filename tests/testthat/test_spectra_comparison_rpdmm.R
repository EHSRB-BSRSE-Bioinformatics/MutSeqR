test_that("parameter extraction infers categories from the parameter vector", {
  pars <- c(log(12), log(2), log(3))
  parameters <- extract_params(pars)
  expect_equal(parameters$p, c(2, 3, 1) / 6)
  expect_null(names(parameters$p))
  expect_equal(parameters$theta, 12)
  named_parameters <- extract_params(pars, c("A", "B", "C"))
  expect_named(named_parameters$p, c("A", "B", "C"))
  expect_equal(unname(named_parameters$p), parameters$p)
  expect_equal(extract_params(log(12))$p, 1)
})

test_that("stable softmax handles large logits", {
  probabilities <- softmax_stable(c(1000, 1001))
  expect_true(all(is.finite(probabilities)))
  expect_equal(probabilities, exp(c(0, 1)) / sum(exp(c(0, 1))))
})

test_that("DM log-likelihood applies the expected ridge penalty", {
  counts <- matrix(c(2, 1, 1, 3), nrow = 2)
  probabilities <- c(0.25, 0.75)
  theta <- 4
  lambda <- 0.5
  alpha <- probabilities * theta
  depths <- colSums(counts)
  expected <- sum(lgamma(theta) - lgamma(depths + theta)) +
    sum(lgamma(counts + alpha) - lgamma(alpha)) -
    lambda * sum(counts) * sum((probabilities - 1 / 2)^2)
  expect_equal(dm_loglik_ridge_stable(counts, probabilities, theta, lambda),
               expected)
  expect_equal(joint_dm_target(c(log(theta), log(1 / 3)), counts, lambda),
               -expected)
})

test_that("joint DM gradient agrees with finite differences", {
  counts <- matrix(c(20, 5, 2, 15, 8, 3, 25, 4, 1), nrow = 3)
  for (pars in list(c(log(12), 0.6, -0.3), c(log(12), 20, -20))) {
    for (lambda in c(0, 0.5)) {
      step <- 1e-5
      numerical <- vapply(seq_along(pars), function(index) {
        offset <- rep(0, length(pars))
        offset[index] <- step
        (joint_dm_target(pars + offset, counts, lambda) -
           joint_dm_target(pars - offset, counts, lambda)) / (2 * step)
      }, numeric(1))
      expect_equal(joint_dm_gradient(pars, counts, lambda), numerical,
                   tolerance = 1e-5)
    }
  }
})

test_that("joint fitting shrinks probabilities and preserves category names", {
  counts <- matrix(c(90, 8, 2, 80, 15, 5, 95, 4, 1, 85, 10, 5), 3)
  category_names <- c("A", "B", "C")
  initial <- c(log(100), log(c(0.9, 0.08) / 0.02))
  raw <- fit_dmm_fast(counts, initial, 0)
  penalized <- fit_dmm_fast(counts, initial, 2)
  expect_false(is.null(raw))
  expect_false(is.null(penalized))
  raw_parameters <- extract_params(raw$par, category_names)
  penalized_parameters <- extract_params(penalized$par, category_names)
  expect_named(penalized_parameters$p, category_names)
  expect_equal(sum(penalized_parameters$p), 1)
  expect_lt(sum((penalized_parameters$p - 1 / 3)^2),
            sum((raw_parameters$p - 1 / 3)^2))
  expect_true(penalized_parameters$theta >= 1e-4)
  expect_true(penalized_parameters$theta <= 400000)
})

test_that("log-space DM simulation preserves dimensions and depths", {
  set.seed(241)
  depths <- c(10, 0, 100)
  draws <- rdm_fast(depths, c(1 - 2e-12, 1e-12, 1e-12), 1e-4)
  expect_equal(dim(draws), c(3L, 3L))
  expect_equal(colSums(draws), depths)
  expect_true(all(is.finite(draws)))
  expect_true(all(draws >= 0))
  expect_equal(dim(rdm_fast(10, c(0.5, 0.5), 1)), c(2L, 1L))
  expect_equal(dim(rdm_fast(c(10, 20), 1, 1)), c(1L, 2L))
})

test_that("comparison returns corrected hybrid and penalized bootstrap results", {
  counts <- matrix(c(20, 5, 2, 15, 8, 3, 25, 4, 1), 3)
  set.seed(242)
  result <- run_penalized_comparison(counts, counts[, 3:1], n_boot = 4)
  expect_equal(result$p_value, result$p_hybrid)
  expect_true(all(is.finite(unlist(result))))
  expect_true(result$p_hybrid >= 1 / 5 && result$p_hybrid <= 1)
  expect_true(result$p_penalized >= 1 / 5 && result$p_penalized <= 1)
  expect_gte(result$lrt, 0)
  expect_gte(result$lrt_penalized, 0)
  expect_error(run_penalized_comparison(counts, counts, n_boot = 0.5))
})

test_that("comparison failures return complete missing results", {
  counts <- matrix(c(20, 5, 15, 8), 2)
  local_mocked_bindings(fit_dmm_fast = function(...) NULL)
  result <- run_penalized_comparison(counts, counts, n_boot = 2)
  expect_named(result, c("p_hybrid", "p_penalized", "lrt", "lrt_penalized",
                         "theta_shared", "theta_g1", "theta_g2", "p_value"))
  expect_true(all(is.na(unlist(result))))
})

test_that("failed bootstrap fits are excluded and require half to succeed", {
  counts <- matrix(c(20, 5, 15, 8), 2)
  calls <- 0L
  fail_after <- 6L
  local_mocked_bindings(fit_dmm_fast = function(x, init_pars, lambda) {
    calls <<- calls + 1L
    if (calls > fail_after) return(NULL)
    list(par = init_pars, value = joint_dm_target(init_pars, x, lambda))
  })
  set.seed(243)
  result <- run_penalized_comparison(counts, counts, n_boot = 2)
  expect_true(result$p_hybrid %in% c(0.5, 1))
  expect_true(result$p_penalized %in% c(0.5, 1))
  calls <- 0L
  fail_after <- 3L
  expect_true(all(is.na(unlist(run_penalized_comparison(counts, counts, n_boot = 2)))))
})

test_that("joint LOOCV scores held-out samples without a penalty", {
  counts <- matrix(c(20, 5, 2, 15, 8, 3, 25, 4, 1,
                     18, 6, 4, 22, 7, 2), 3)
  rownames(counts) <- c("C>A", "C>G", "C>T")
  colnames(counts) <- paste0("sample", seq_len(ncol(counts)))
  lambda <- 0.1
  expected <- sum(vapply(seq_len(ncol(counts)), function(index) {
    training <- counts[, -index, drop = FALSE]
    observed <- pmax(rowSums(training) / sum(training), 1e-6)
    observed <- observed / sum(observed)
    fit <- fit_dmm_fast(training, c(log(100), log(observed[-3] / observed[3])),
                        lambda)
    parameters <- extract_params(fit$par)
    dm_loglik_ridge_stable(counts[, index, drop = FALSE], parameters$p,
                          parameters$theta, lambda = 0)
  }, numeric(1)))
  mf_data <- expand.grid(
    normalized_subtype = rownames(counts),
    sample = colnames(counts),
    stringsAsFactors = FALSE
  )
  mf_data$sum_min <- as.vector(counts)
  mf_data$sum_max <- 2 * mf_data$sum_min
  result <- select_optimal_lambda(mf_data, lambdas = lambda)
  expect_equal(result$best_lambda, lambda)
  expect_equal(unname(result$scores), expected)
  expect_equal(result$lambdas, lambda)
  result_max <- select_optimal_lambda(mf_data, mf_type = "max", lambdas = lambda)
  expect_false(isTRUE(all.equal(result_max$scores, result$scores)))
  expect_error(select_optimal_lambda(mf_data[, "sum_min", drop = FALSE]))
})

test_that("joint LOOCV reports all failed candidates without a best lambda", {
  local_mocked_bindings(fit_dmm_fast = function(...) NULL)
  mf_data <- data.frame(
    normalized_subtype = c("C>A", "C>G", "C>A", "C>G"),
    sample = c("s1", "s1", "s2", "s2"),
    sum_min = c(20, 5, 15, 8)
  )
  result <- select_optimal_lambda(mf_data, lambdas = c(0, 1))
  expect_identical(result$best_lambda, NA_real_)
  expect_true(all(is.na(result$scores)))
})

test_that("public comparison aligns subtype rows across groups", {
  observed <- NULL
  local_mocked_bindings(run_penalized_comparison = function(g1, g2, ...) {
    observed <<- list(g1 = g1, g2 = g2)
    list(p_value = 0.5, p_hybrid = 0.5, p_penalized = 0.5,
         lrt = 0, lrt_penalized = 0, theta_g1 = 1, theta_g2 = 1,
         theta_shared = 1)
  })
  mf_data <- data.frame(
    normalized_subtype = c("A", "B", "B", "C"),
    sample = paste0("s", 1:4),
    dose = c("g1", "g1", "g2", "g2"),
    sum_min = c(10, 5, 8, 7)
  )
  result <- spectra_comparison_rpdmm(
    mf_data, "dose", data.frame("g1", "g2")
  )
  expect_equal(rownames(observed$g1), c("A", "B", "C"))
  expect_equal(rownames(observed$g2), c("A", "B", "C"))
  expect_equal(observed$g1["C", ], 0)
  expect_equal(observed$g2["A", ], 0)
  expect_equal(nrow(result), 1L)
})

test_that("public spectra comparison retains old columns and adds both tests", {
  mf_data <- expand.grid(normalized_subtype = c("C>A", "C>G", "C>T"),
                         sample = paste0("sample", 1:6), stringsAsFactors = FALSE)
  mf_data$dose <- rep(c("control", "treated"), each = 9)
  mf_data$sum_min <- rep(c(20, 5, 2, 15, 8, 3, 25, 4, 1), 2)
  set.seed(244)
  result <- spectra_comparison_rpdmm(
    mf_data, "dose", data.frame("control", "treated"), n_boot = 2
  )
  expect_equal(nrow(result), 1L)
  expect_equal(result$p_value, result$p_hybrid)
  expect_true(all(c("p_penalized", "lrt_penalized", "Significance",
                    "theta_shared", "theta_g1", "theta_g2") %in% names(result)))
})