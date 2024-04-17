library(testthat)

# Define a test case for import_mut_data function
test_that("import_vcf_datafunction correctly imports vcf files", {
  # Create temporary test file with example mutation data
 
test_file <- system.file("extdata", "vcf_sample.vcf", package = "MutSeqR")


  # Call the import_mut_data function on the test data
  mut_data <- import_vcf_data(vcf_file = test_file, vaf_cutoff = 0.1, regions = "mouse", output_granges = TRUE)
  
  expect_true(is(mut_data, "GRanges"), info = "Check if the resulting object is a granges object")
  
  expect_equal(NROW(mut_data), 12, info = "Check if the resulting object has the correct number of rows; Remove 1 INV outside of region ranges")
  
  expect_equal(GenomicRanges::start(mut_data), c(69304218, 69304219, 69304220,
                                                 69304221, 69306180, 50833809, 
                                                 50833809, 50835177, 50835330, 
                                                 50835330, 50835421, 21443246), 
               info = "Check that the start ranges are correct")
  expect_equal(GenomicRanges::end(mut_data), c(69304218, 69304219, 69304220,
                                               69304221, 69306180, 50833810, 
                                               50833809, 50835177, 50835331, 
                                               50835330, 50835421, 21444398), 
               info = "Check that the end ranges are correct")
  
  expect_true(all(op$processed_required_mut_cols %in% colnames(GenomicRanges::mcols(mut_data))),
              info = "Check that the resulting object contains the required columns to calculate mutation frequency")
  
  expect_true(all(c("sample", "nchar_ref", "nchar_alt", "varlen", "subtype", 
                    "context", "total_depth", "no_calls", "context_with_mutation", 
                    "normalized_context", "normalized_subtype", "short_ref", 
                    "normalized_ref", "normalized_context_with_mutation", 
                    "gc_content", "vaf", "is_germline") %in% colnames(GenomicRanges::mcols(mut_data))),
              info = "Check that the function created all the new columns")
  
  expect_equal(unique(GenomicRanges::mcols(mut_data)$sample), "PRC_ST_115.1", 
               info = "Check that the sample column is created correctly")
  expect_equal(GenomicRanges::mcols(mut_data)$variation_type, c("no_variant", "no_variant", 
                                                                "no_variant", "no_variant", "snv",
                                                                "deletion", "no_variant", 
                                                                "insertion", "mnv", 
                                                                "snv", "snv", "symbolic"), 
               info = "Check that the variation_type column is changes to default names")
  expect_equal(GenomicRanges::mcols(mut_data)$nchar_ref, c(1, 1, 1, 1, 1, 2, 1,
                                                           1, 2, 1, 1, 1), 
               info = "Check that the nchar_ref column is created correctly")
  expect_equal(GenomicRanges::mcols(mut_data)$nchar_alt, c(0, 0, 0, 0, 1, 1, 0, 
                                                           2, 2, 1, 1, NA), 
               info = "Check that the nchar_alt column is created correctly")
  expect_equal(GenomicRanges::mcols(mut_data)$varlen, c(NA, NA, NA, NA, 1, -1,
                                                        NA, 1, 2, 1, 1, NA), 
               info = "Check that the VARLEN column is created correctly")
  expect_equal(GenomicRanges::mcols(mut_data)$subtype, c("no_variant", "no_variant", 
                                                         "no_variant", "no_variant", "G>A",
                                                         "deletion", "no_variant", 
                                                         "insertion", "mnv", 
                                                         "T>C", "T>C", "symbolic"), 
               info = "Check that the subtype column is created correctly")
  expect_equal(GenomicRanges::mcols(mut_data)$context, c("ACA", "CAA", "AAT", "ATC",
                                                         "CGT", "AGA", "AGA", "ACT", "TTA",
                                                         "TTA", "CTA", "TTG"), 
               info = "Check that the context column is created correctly")
  expect_equal(GenomicRanges::mcols(mut_data)$total_depth, c(4123, 3491, 3524, 
                                                             4268, 21306, 22157, 
                                                             22157, 11601, 11584, 11584, 
                                                             11564, 21569), 
               info = "Check that the total_depth column is created correctly using the take_del method") 
  expect_equal(GenomicRanges::mcols(mut_data)$no_calls, c(384, 1066, 1113, 405, 
                                                          720, 106, 48, 0, 222,
                                                          222, 188, 170), 
               info = "Check that the no_calls column is created correctly using the take_del method")
  expect_equal(GenomicRanges::mcols(mut_data)$vaf, c(0, 0, 0, 0, 1/21306, 1/22157, 
                                                     0, 10814/11601, 11180/11584, 
                                                     403/11584, 11561/11564, 1/21569), 
               info = "Check that the vaf column is created correctly using the take_del method")
   mut_data_take_mean <- import_vcf_data(vcf_file = test_file, vaf_cutoff = 0.1, regions = "mouse", depth_calc = "take_mean", output_granges = TRUE)
  expect_equal(GenomicRanges::mcols(mut_data_take_mean)$total_depth, c(4123, 3491, 3524, 
                                                                      4268, 21306, 22110, 
                                                                      22110, 11601, 11584, 11584, 
                                                                      11564, 21569), 
               info = "Check that the total_depth column is created correctly using the take_mean method") 
  expect_equal(GenomicRanges::mcols(mut_data_take_mean)$no_calls, c(384, 1066, 1113, 405, 
                                                                    720, 153, 95, 0, 222,
                                                                    222, 188, 170), 
               info = "Check that the no_calls column is created correctly using the take_mean method") 
  expect_equal(GenomicRanges::mcols(mut_data_take_mean)$vaf, c(0, 0, 0, 0, 1/21306, 1/22110, 
                                                     0, 10814/11601, 11180/11584, 
                                                     403/11584, 11561/11564,
                                                     1/21569), 
               info = "Check that the vaf column is created correctly using the take_mean method") 
  expect_equal(GenomicRanges::mcols(mut_data)$context_with_mutation, 
               c("no_variant", "no_variant", "no_variant", "no_variant", "C[G>A]T",
                 "deletion", "no_variant", "insertion", "mnv", "T[T>C]A", "C[T>C]A","symbolic"), 
               info = "Check that the context_with_mutation column is created correctly")
  expect_equal(GenomicRanges::mcols(mut_data)$normalized_context, 
               c("ACA", "TTG", "ATT", "ATC", "ACG", "TCT", "TCT", "ACT", "TTA",
                  "TTA", "CTA", "TTG"), 
               info = "Check that the normalized_context column is created correctly")
  expect_equal(GenomicRanges::mcols(mut_data)$normalized_subtype, c("no_variant", "no_variant", 
                                                         "no_variant", "no_variant", "C>T",
                                                         "deletion", "no_variant", 
                                                         "insertion", "mnv", 
                                                         "T>C", "T>C", "symbolic"), 
               info = "Check that the normalized_subtype column is created correctly") 
  expect_equal(GenomicRanges::mcols(mut_data)$short_ref, c("C", "A", "A", "T",
                                                           "G", "G", "G", "C", 
                                                           "T", "T", "T", "T"), 
               info = "Check that the short_ref column is created correctly") 
  expect_equal(GenomicRanges::mcols(mut_data)$normalized_ref, c("C", "T", "T", "T",
                                                           "C", "C", "C", "C", 
                                                           "T", "T", "T", "T"), 
               info = "Check that the normalized_ref column is created correctly") 
  expect_equal(GenomicRanges::mcols(mut_data)$normalized_context_with_mutation, 
               c("no_variant", "no_variant", "no_variant", "no_variant", "A[C>T]G",
                 "deletion", "no_variant", "insertion", "mnv", "T[T>C]A", "C[T>C]A",
                  "symbolic"), 
               info = "Check that the normalized_context_with_mutation column is created correctly") 
  expect_equal(GenomicRanges::mcols(mut_data)$gc_content, c(1/3, 1/3, 0/3, 1/3,
                                                            2/3, 1/3, 1/3, 1/3,
                                                            0/3, 0/3, 1/3, 1/3), 
                 info = "Check that the gc_content column is created correctly")
      })
