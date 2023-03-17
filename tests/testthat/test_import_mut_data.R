library(testthat)

# Define a test case for import_mut_data function
test_that("import_mut_data function correctly imports mutation data", {
  # Create temporary test file with example mutation data
  tmpfile <- tempfile(fileext = ".mut")
  write.table(
    data.frame(
      contig = c("chr1", "chr1", "chr2", "chr2"),
      start = c(100, 200, 300, 400),
      end = c(101, 201, 301, 401),
      context = c("GCA", "GGC", "ATC", "ACC"),
      subtype = c("G>T", "G>A", "A>G", "A>C"),
      total_depth = c(50, 100, 75, 150),
      alt_depth = c(10, 20, 30, 50)
    ),
    file = tmpfile,
    sep = "\t", row.names = FALSE
  )
  
  # Call the import_mut_data function on the test data
  mut_data <- import_mut_data(mut_file = tmpfile)
  
  # Check if the resulting object has the correct number of rows and columns
  expect_equal(dim(mut_data), c(4, 13))
  # Check if the resulting object has the correct column names
  expect_equal(colnames(mut_data), c(
    "contig", "start", "end", "context", "subtype", "total_depth",
    "alt_depth", "ref_depth", "context_with_mutation",
    "normalized_context", "normalized_subtype",
    "normalized_context_with_mutation", "gc_content"
  ))
  # Check if the resulting object has the correct data type for each column
  expect_equal(
    sapply(mut_data, class),
    c("Rle", "integer", "integer", "character", "character", "integer", "integer", "integer",
      "character", "character", "character", "character", "numeric")
  )
  
  # Clean up temporary file
  unlink(tmpfile)
})
