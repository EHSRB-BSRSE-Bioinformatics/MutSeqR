# Define test cases for model_mf function
test_that("model_mf works as expected ", {
  # Create a sample mf_data dataframe
  mf_data <- data.frame(
    dose = c(0, 25, 50, 100, 0, 25, 50, 100),
    sum_min = c(10, 20, 30, 40, 50, 60, 70, 80),
    group_depth = c(100, 200, 300, 400, 500, 600, 700, 800)
  )

  # Call the model_mf function on the sample data
  model_results <- suppressWarnings(model_mf(mf_data = mf_data,
                                             fixed_effects = "dose",
                                             reference_level = 25))

  # Check if the specified columns are converted to factors
  expect_equal(class(model_results$model_data$dose), "factor",
               info = "Check if the dose column is converted to factor")
  expect_equal(levels(model_results$model_data$dose),
               c("25", "0", "50", "100"),
               info = "Check if the reference level for 'dose' is correctly set")
  expect_true(all(c("model_data", "model_formula", "summary",
                    "residuals_histogram",
                    "residuals_qq_plot",
                    "point_estimates_matrix",
                    "point_estimates")
                  %in% names(model_results)),
              info = "All objects are created on model_results")
})
