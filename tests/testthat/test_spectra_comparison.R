# Define test cases for spectra_comparison function
test_that("spectra_comparison works as expected", {

  # Create a sample mf_data dataframe
  mf_data <- data.frame(
    dose = c(0, 25, 50, 100, 0, 25, 50, 100),
    tissue = c(rep("BM", 4), rep("LV", 4)),
    normalized_subtype = rep(c("subtype1", "subtype2"), 4),
    sum_min = c(10, 20, 30, 40, 50, 60, 70, 80)
  )

  # Create a contrast table file
  contrast_df <- data.frame(
    V1 = c("0:BM", "25:BM", "50:BM", "100:BM"),
    V2 = c("0:LV", "25:LV", "50:LV", "100:LV")
  )

  # Call the spectra_comparison function on the sample data
  result <- spectra_comparison(mf_data = mf_data,
                               mf_type = "min",
                               exp_variable = c("dose", "tissue"),
                               contrasts = contrast_df)

  # Check if the result is a data frame
  expect_equal(class(result), "data.frame",
               info = "Check if the result is a data frame")

  # Check if the result has the expected columns
  expect_true(all(c("contrasts", "G2", "p.value", "adj_p.value", "Significance") %in% names(result)),
              info = "Check if the result has the expected columns")

  # Check if the result has the correct number of rows
  expect_equal(nrow(result), 4,
               info = "Check if the result has the correct number of rows")

  # Check if the p-value is within the expected range
  expect_true(all(result$adjP >= 0 & result$adjP <= 1),
              info = "Check if the p-value is within the expected range")
})

test_that("G2 dispatch preserves existing calls and results", {
  mf_data <- data.frame(
    dose = rep(c("control", "treated"), each = 2),
    normalized_subtype = rep(c("C>A", "C>G"), 2),
    sum_min = c(10, 20, 15, 25)
  )
  contrasts <- data.frame("treated", "control")
  count_matrix <- matrix(
    c(10, 20, 15, 25), nrow = 2,
    dimnames = list(c("C>A", "C>G"), c("control", "treated"))
  )
  direct <- MutSeqR:::spectra_comparison_G2(count_matrix, contrasts)
  expect_equal(spectra_comparison(mf_data, "dose", "min", contrasts), direct)
  expect_equal(
    spectra_comparison(mf_data, "dose", "min", contrasts, method = "g2"),
    direct
  )
  expect_error(
    spectra_comparison(mf_data, "dose", contrasts = contrasts, method = "other"),
    "'arg' should be"
  )
  expect_error(
    spectra_comparison(mf_data, "dose", contrasts = contrasts, n_boot = 2),
    "only be supplied with method"
  )
  expect_error(
    spectra_comparison(mf_data, "dose", contrasts = contrasts,
                       method = "g2", lambda = 1),
    "only be supplied with method"
  )
})

test_that("rpdmm dispatch forwards method-specific and common arguments", {
  arguments <- names(formals(spectra_comparison))
  expect_true(all(c("lambda", "n_boot") %in% arguments))
  expect_false("..." %in% arguments)

  mf_data <- expand.grid(
    normalized_subtype = c("C>A", "C>G", "C>T"),
    sample = paste0("sample", 1:6), stringsAsFactors = FALSE
  )
  mf_data$dose <- rep(c("control", "treated"), each = 9)
  mf_data$sum_min <- rep(c(20, 5, 2, 15, 8, 3, 25, 4, 1), 2)
  mf_data$sum_max <- mf_data$sum_min + 1
  contrast_file <- tempfile(fileext = ".csv")
  on.exit(unlink(contrast_file), add = TRUE)
  writeLines("control,treated", contrast_file)

  count_matrix <- matrix(
    mf_data$sum_max, nrow = 3,
    dimnames = list(c("C>A", "C>G", "C>T"), paste0("sample", 1:6))
  )
  set.seed(244)
  direct <- MutSeqR:::spectra_comparison_rpdmm(
    list(control = count_matrix[, 1:3], treated = count_matrix[, 4:6]),
    data.frame("control", "treated"),
    lambda = 0.1, n_boot = 2
  )
  set.seed(244)
  dispatched <- spectra_comparison(
    mf_data, "dose", "max", contrast_file, ",", method = "rpdmm",
    lambda = 0.1, n_boot = 2
  )
  expect_equal(dispatched, direct)
  expect_true(all(c("p_value", "p_hybrid", "p_penalized", "lrt",
                    "lrt_penalized") %in% names(dispatched)))
})

test_that("shared preparation preserves pooled and replicate counts", {
  observed <- NULL
  local_mocked_bindings(
    run_penalized_comparison = function(g1, g2, lambda, n_boot) {
      observed <<- list(g1 = g1, g2 = g2)
      list(p_value = 0.5, p_hybrid = 0.5, p_penalized = 0.5,
           lrt = 0, lrt_penalized = 0, theta_g1 = 1, theta_g2 = 1,
           theta_shared = 1)
    }
  )
  mf_data <- data.frame(
    dose = c(rep("control", 5), rep("treated", 4)),
    tissue = "BM",
    sample = c("s1", "s1", "s1", "s2", "s2", "s3", "s3", "s4", "s4"),
    normalized_subtype = c("A", "A", "B", "B", "C", "A", "C", "B", "C"),
    sum_min = c(3, 7, 5, 4, 8, 6, 9, 10, 4)
  )
  mf_data$sum_max <- 2 * mf_data$sum_min
  contrast_file <- tempfile(fileext = ".csv")
  on.exit(unlink(contrast_file), add = TRUE)
  writeLines("control:BM,treated:BM", contrast_file)

  g2 <- spectra_comparison(
    mf_data, c("dose", "tissue"), mf_type = "max",
    contrasts = contrast_file, cont_sep = ",", method = "g2"
  )
  pooled <- 2 * matrix(c(10, 9, 8, 6, 10, 13), nrow = 3)
  expected <- outer(rowSums(pooled), colSums(pooled)) / sum(pooled)
  expected_g2 <- 2 * sum(pooled * log(pooled / expected))
  expect_equal(g2$G2, expected_g2)
  expect_equal(g2$p.value, stats::pchisq(expected_g2, 2, lower.tail = FALSE))
  expect_equal(g2$contrasts, "control:BM vs treated:BM")

  spectra_comparison(
    mf_data, c("dose", "tissue"), mf_type = "max",
    contrasts = contrast_file, cont_sep = ",", method = "rpdmm"
  )
  expect_equal(
    observed$g1,
    2 * matrix(c(10, 5, 0, 0, 4, 8), nrow = 3,
               dimnames = list(c("A", "B", "C"), c("s1", "s2")))
  )
  expect_equal(
    observed$g2,
    2 * matrix(c(6, 0, 9, 0, 10, 4), nrow = 3,
               dimnames = list(c("A", "B", "C"), c("s3", "s4")))
  )
})

test_that("both comparison methods validate shared inputs", {
  local_mocked_bindings(
    spectra_comparison_G2 = function(...) stop("G2 method called"),
    spectra_comparison_rpdmm = function(...) stop("RPDMM method called")
  )
  mf_data <- data.frame(
    dose = c("control", "treated"),
    sample = c("s1", "s2"),
    normalized_subtype = c("A", "B"),
    sum_min = c(10, 20)
  )
  contrasts <- data.frame("control", "treated")
  for (method in c("g2", "rpdmm")) {
    expect_error(spectra_comparison(mf_data, "missing", contrasts = contrasts,
                                    method = method),
                 "exp_variable")
    expect_error(spectra_comparison(mf_data, "dose", mf_type = "max",
                                    contrasts = contrasts, method = method),
                 "Missing count column in mf_data: sum_max")
    expect_error(
      spectra_comparison(mf_data[, c("dose", "sample", "sum_min")], "dose",
                         contrasts = contrasts, method = method),
      "No valid subtype column"
    )
    expect_error(
      spectra_comparison(mf_data, "dose", contrasts = data.frame("control"),
                         method = method),
      "exactly 2 columns"
    )
    expect_error(
      spectra_comparison(mf_data, "dose",
                         contrasts = data.frame("control", "missing"),
                         method = method),
      "Groups in contrast table not found in data: missing"
    )
  }
  expect_error(
    spectra_comparison(
      mf_data[, c("dose", "normalized_subtype", "sum_min")], "dose",
      contrasts = contrasts, method = "rpdmm"
    ),
    "sample"
  )
  expect_error(
    spectra_comparison(mf_data, "dose", contrasts = contrasts,
                       method = "rpdmm", lambda = -1),
    "lambda"
  )
  expect_error(
    spectra_comparison(mf_data, "dose", contrasts = contrasts,
                       method = "rpdmm", n_boot = 0),
    "n_boot"
  )
})
