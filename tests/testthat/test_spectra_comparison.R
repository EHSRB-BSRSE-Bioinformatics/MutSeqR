library(testthat)

# Define test cases for spectra_comparison function
test_that("spectra_comparison works as expected", {
  # Create a sample mf_data dataframe
  mf_data <- data.frame(
    dose = c(0, 25, 50, 100, 0, 25, 50, 100),
    tissue = c("bone_marrow", "bone_marrow", "bone_marrow", "bone_marrow", "liver", "liver", "liver", "liver"),
    normalized_subtype = c("subtype1", "subtype2", "subtype1", "subtype2", "subtype1", "subtype2", "subtype1", "subtype2"),
    dose_sum_min = c(10, 20, 30, 40, 50, 60, 70, 80)
  )
  
  # Create a contrast table file
  contrast_table_file <- tempfile(fileext = ".txt")
  write.table(
    data.frame(
      V1 = c("0:bone_marrow", "25:bone_marrow", "50:bone_marrow", "100:bone_marrow"),
      V2 = c("0:liver", "25:liver", "50:liver", "100:liver")
    ),
    file = contrast_table_file,
    sep = "\t", row.names = FALSE, col.names = FALSE  # Remove column names
  )

  
  # Call the spectra_comparison function on the sample data
  result <- spectra_comparison(mf_data = mf_data,
                               muts = "dose_sum_min",
                               subtype_col = "normalized_subtype",
                               group = c("dose", "tissue"),
                               contrast_table_file = contrast_table_file)
  
  # Check if the result is a data frame
  expect_equal(class(result), "data.frame", info = "Check if the result is a data frame")
  
  # Check if the result has the expected columns
  expect_true(all(c("g2", "p_value", "adj_p_value") %in% names(result)), info = "Check if the result has the expected columns")
  
  # Check if the result has the correct number of rows
  expect_equal(nrow(result), 4, info = "Check if the result has the correct number of rows")
  
  # Check if the p-value is within the expected range
  for(i in seq_len(nrow(result))) {
    # Check if the p.value is within the expected range
    expect_true(result$p_value[i] >= 0 && result$p_value[i] <= 1,
                info = paste("Check if the p-value in row", i, "is within the expected range"))
  }
  
  for(i in seq_len(nrow(result))) {
    # Check if the p.value is within the expected range
    expect_true(result$adj_p_value[i] >= 0 && result$adj_p_value[i] <= 1,
                info = paste("Check if the p-value in row", i, "is within the expected range"))
  }})
