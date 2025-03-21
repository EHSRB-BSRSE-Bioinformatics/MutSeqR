library(testthat)

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
                               cols_to_group = c("dose", "tissue"),
                               contrasts = contrast_df)

  # Check if the result is a data frame
  expect_equal(class(result), "data.frame",
               info = "Check if the result is a data frame")

  # Check if the result has the expected columns
  expect_true(all(c("G2", "p.value", "adjP", "sign") %in% names(result)),
              info = "Check if the result has the expected columns")

  # Check if the result has the correct number of rows
  expect_equal(nrow(result), 4,
               info = "Check if the result has the correct number of rows")

  # Check if the p-value is within the expected range
  expect_true(all(result$adjP >= 0 & result$adjP <= 1),
              info = "Check if the p-value is within the expected range")
})
