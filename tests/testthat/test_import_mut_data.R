library(testthat)
library(fs)

# Define test cases for import_mut_data function
test_that("import_mut_data function correctly imports mutation data", {
  # Create temporary test file with example mutation data
  tmpfile <- tempfile(fileext = ".mut")
  write.table(
    data.frame(
      sample = c("mouse1", "mouse2", "mouse1", "mouse2"),
      contig = c("chr1", "chr1", "chr2", "chr2"),
      start = c(69304225, 69304240, 50833424, 50833439),
      end = c(69304226, 69304241, 50833425, 50833440),
      total_depth = c(50, 100, 75, 150),
      alt_depth = c(1, 20, 30, 50),
      ref = c("C", "G", "T", "AA"),
      alt = c("T", "A", ".", "A")
    ),
    file = tmpfile,
    sep = "\t", row.names = FALSE
  )

  #create a temporary custom regions file
  tmpfile2 <- tempfile(fileext = ".mut")
  write.table(
    data.frame(
      contig = c("chr1", "chr2"),
      start = c(69304217, 50833175),
      end = c(69306617, 50835575),
      description = c("region_330", "region_4547"),
      location_relative_to_genes = c("intergenic", "intergenic")
    ),
    file = tmpfile2,
    sep = "\t", row.names = FALSE
  )

  # Call the import_mut_data function on the test data
  mut_data <- import_mut_data(mut_file = tmpfile,
                              regions = tmpfile2,
                              species = "mouse",
                              genome = "mm10",
                              output_granges = FALSE)

  expect_true(is(mut_data, "data.frame"),
              info = "Check if the resulting object is a data frame")
  expect_equal(NROW(mut_data), 4,
               info = "Check if the resulting object has the correct number of rows")
  expect_true(all(c("short_ref", "normalized_ref", "context",
                    "normalized_context", "variation_type", "subtype",
                    "normalized_subtype", "context_with_mutation",
                    "normalized_context_with_mutation", "nchar_ref",
                    "nchar_alt", "varlen", "ref_depth", "vaf", "gc_content",
                    "row_has_duplicate") %in% colnames(mut_data)),
              info = "Check if the resulting object has the correct columns")

  # Clean up temporary file
  unlink(tmpfile)
  unlink(tmpfile2)
})
