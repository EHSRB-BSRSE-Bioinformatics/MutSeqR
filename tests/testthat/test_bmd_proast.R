test_that("bmd_proast optionally retains model-averaged bootstrap curves", {
  skip_if_not_installed("withr")

  mf_file <- system.file(
    "extdata/Example_files/mf_data_global.rds",
    package = "MutSeqR"
  )
  expect_true(nzchar(mf_file))
  mf_data <- readRDS(mf_file)

  suppressWarnings(
    suppressMessages(
      invisible(capture.output(
        result <- bmd_proast(
          mf_data = mf_data,
          dose_col = "dose",
          response_col = "mf_min",
          bmr = 0.5,
          model_averaging = TRUE,
          num_bootstraps = 2,
          raw_results = TRUE,
          return_bootstrap_curves = TRUE,
          seed = 125
        )
      ))
    )
  )

  ma_name <- grep(
    "model_averaging$",
    names(result$raw_results),
    value = TRUE
  )
  expect_length(ma_name, 1)

  curves <- result$raw_results[[ma_name]]$MA$bootstrap_curves
  expect_s3_class(curves, "data.frame")
  expect_named(
    curves,
    c("bootstrap", "model", "covariate", "dose", "response")
  )
  expect_setequal(unique(curves$bootstrap), 1:2)
  expect_true(all(curves$model == "model_average"))
  expect_true(all(is.finite(curves$dose)))
  expect_true(all(is.finite(curves$response)))
  expect_true(all(curves$dose > 0))
  expect_equal(
    as.integer(table(curves$bootstrap)),
    rep(length(unique(curves$dose)), 2)
  )
  expect_false(formals(bmd_proast)$return_bootstrap_curves)
})
