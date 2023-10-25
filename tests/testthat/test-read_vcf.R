library(testthat)

# Define a test case for import_mut_data function
test_that("read_vcf function correctly vcf files", {
  # Create temporary test file with example mutation data
 
test_file <- system.file("extdata", "vcf_sample.vcf", package = "DupSeqR")


  # Call the import_mut_data function on the test data
  mut_data <- read_vcf(vcf_file = test_file, regions_file = "mouse", assembly = "GRCm38")
  
  
  expect_true(is(mut_data, "GRanges"), info = "Check if the resulting object is a granges object")
  
  expect_equal(NROW(mut_data), 4, info = "Check if the resulting object has the correct number of rows")
  
  # Ranges: start+1 because we changed from 0based to 1based. 
  expect_equal(GenomicRanges::start(mut_data), c(69304225, 69304240, 50833424, 50833439) +1, info = "Check if the start ranges are correct")
  expect_equal(GenomicRanges::end(mut_data), c(69304226, 69304241, 50833425, 50833440), info = "Check if the end ranges are correct")
  
  expect_equal(names(mcols(mut_data)), c(
    "sample", "context", "subtype", "variation_type", "total_depth",
    "alt_depth", "ref", "alt", "ref_depth", "context_with_mutation",
    "normalized_context", "normalized_subtype", "short_ref", "normalized_ref", 
    "normalized_context_with_mutation", "gc_content", "VAF", 
    "description", "location_relative_to_genes"
  ), info = "Check if the resulting object has the correct meta data column names")
  
  expect_equal(sapply(mcols(mut_data), class),
               c(sample = "character", context = "character", subtype = "character", variation_type = "character", 
                 total_depth = "integer", alt_depth = "integer", ref = "character", alt = "character", ref_depth = "integer",
                 context_with_mutation = "character", normalized_context = "character", normalized_subtype = "character", 
                 short_ref = "character", normalized_ref = "character", normalized_context_with_mutation = "character", 
                 gc_content = "numeric", VAF = "numeric", description = "character", location_relative_to_genes = "character" ),
               info = " Check if the resulting object has the correct data type for each metadata column" )
  
  # Clean up temporary file
  unlink(tmpfile)
  unlink(tmpfile2)
})