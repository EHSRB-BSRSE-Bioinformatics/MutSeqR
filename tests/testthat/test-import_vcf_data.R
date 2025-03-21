library(testthat)

# Define a test case for import_mut_data function
test_that("import_vcf_datafunction correctly imports vcf files", {
  # Create temporary test file with example mutation data

  test_file <- system.file("extdata",
                           "example_import_vcf_data_cleaned.vcf.bgz",
                           package = "MutSeqR")


  # Call the import_mut_data function on the test data
  mut_data <- suppressWarnings(import_vcf_data(vcf_file = test_file,
                              regions = "TSpanel_mouse",
                              species = "mouse",
                              genome = "mm10",
                              output_granges = FALSE))

  expect_true(is(mut_data, "data.frame"),
              info = "Check if the resulting object is a data frame")
  expect_true(all(c("short_ref", "normalized_ref", "context",
                    "normalized_context", "variation_type", "subtype",
                    "normalized_subtype", "context_with_mutation",
                    "normalized_context_with_mutation", "nchar_ref",
                    "nchar_alt", "varlen", "ref_depth", "vaf", "gc_content",
                    "row_has_duplicate") %in% colnames(mut_data)),
              info = "Check if the resulting object has the correct columns")

})
