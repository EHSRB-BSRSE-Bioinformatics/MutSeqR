library(testthat)

# Define a test case for import_mut_data function
test_that("read_vcf function correctly vcf files", {
  # Create temporary test file with example mutation data
 
test_file <- system.file("extdata", "vcf_sample.vcf", package = "DupSeqR")


  # Call the import_mut_data function on the test data
  mut_data <- read_vcf(vcf_file = test_file, regions_file = "mouse", assembly = "GRCm38")
  
  expected_message <-"Expected ' alt ' but found ' alt.value ', matching columns in input data\nExpected ' variation_type ' but found ' TYPE ', matching columns in input data"
  captured_output <- capture_output(read_vcf(vcf_file = test_file, regions_file = "mouse", assembly = "GRCm38"))
  expect_equal(expected_message, captured_output, info = "Check that the function reports messages correctly")
  
  expect_true(is(mut_data, "GRanges"), info = "Check if the resulting object is a granges object")
  
  expect_equal(NROW(mut_data), 13, info = "Check if the resulting object has the correct number of rows")
  
  expect_equal(GenomicRanges::start(mut_data), c(69304218, 69304219, 69304220,
                                                 69304221, 69306180, 50833809, 
                                                 50833809, 50835177, 50835330, 
                                                 50835330, 50835421, 96824726,
                                                 21443246), info = "Check that the start ranges are correct")
  expect_equal(GenomicRanges::end(mut_data), c(69304218, 69304219, 69304220,
                                               69304221, 69306180, 50833810, 
                                               50833809, 50835177, 50835331, 
                                               50835330, 50835421, 96826154,
                                               21444398), info = "Check that the end ranges are correct")
  
  expect_that(all(op$processed_required_mut_cols %in% colnames(GenomicRanges::mcols(mut_data))), is_true(), 
              info() = "Check that the resuting object contains the required columns to calculate mutation frequency")
  
  expect_that(all(c("sample", "nchar_ref", "nchar_alt", "VARLEN", "subtype", "context",
                    "total_depth", "no_calls", "sequence", "ext_start", "ext_end",
                    "context_with_mutation", "normalized_context", "normalized_subtype",
                    "short_ref", "normalized_ref", "normalized_context_with_mutation", 
                    "gc_content", "VAF", ) %in% colnames(GenomicRanges::mcols(mut_data))), is_true(), 
              info() = "Check that the function created all the new columns")
  
  expect_equal(unique(GenomicRanges::mcols(mut_data)$sample), "PRC_ST_115.1", 
               info = "Check that the sample column is created correctly")
  expect_equal(GenomicRanges::mcols(mut_data)$variation_type, c("no_variant", "no_variant", 
                                                                "no_variant", "no_variant", "snv",
                                                                "deletion", "no_variant", 
                                                                "insertion", "mnv", 
                                                                "snv", "snv", "symbolic", 
                                                                "symbolic"), 
               info = "Check that the variation_type column is changes to default names")
  expect_equal(GenomicRanges::mcols(mut_data)$nchar_ref, c(1, 1, 1, 1, 1, 2, 1,
                                                           1, 2, 1, 1, 1, 1), 
               info = "Check that the nchar_ref column is created correctly")
  expect_equal(GenomicRanges::mcols(mut_data)$nchar_alt, c(0, 0, 0, 0, 1, 1, 0, 
                                                           2, 2, 1, 1, NA, NA), 
               info = "Check that the nchar_alt column is created correctly")
  expect_equal(GenomicRanges::mcols(mut_data)$VARLEN, c(NA, NA, NA, NA, 1, -1,
                                                        NA, 1, 2, 1, 1, NA, NA), 
               info = "Check that the VARLEN column is created correctly")
  expect_equal(GenomicRanges::mcols(mut_data)$subtype, c("no_variant", "no_variant", 
                                                         "no_variant", "no_variant", "G>A",
                                                         "deletion", "no_variant", 
                                                         "insertion", "mnv", 
                                                         "T>C", "T>C", "symbolic", 
                                                         "symbolic"), 
               info = "Check that the subtype column is created correctly")
  expect_equal(GenomicRanges::mcols(mut_data)$context, c("ACA", "CAA", "AAT", "ATC",
                                                         "CGT", "AGA", "ACT", "TTA",
                                                         "TTA", "CTA", ".", "TTG"), 
               info = "Check that the subtype column is created correctly")

})