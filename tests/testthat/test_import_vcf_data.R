library(testthat)

# Define a test case for import_mut_data function
test_that("import_vcf_datafunction correctly imports vcf files", {
  # Create temporary test file with example mutation data
  file <- file.path("./testdata/simple_vcf_data.vcf")

  # Call the import_mut_data function on the test data
  mut_data <- import_vcf_data(vcf_file = file,
    regions = NULL,
    species = "mouse",
    genome = "mm10",
    output_granges = FALSE
  )
  colnames <- c(
    MutSeqR::op$base_required_mut_cols,
    MutSeqR::op$processed_required_mut_cols, # subtype/context cols
    "total_depth", "ref_depth", "vaf", # depth cols
    "nchar_ref", "nchar_alt", "varlen",
    "gc_content", "row_has_duplicate",
    "strand", "width", # added by GRanges
    "alt.group", "alt.group_name", "AD_1", "AD_2" # vcf cols
  )
  expect_named(mut_data, colnames, ignore.order = TRUE) # check columns
  expect_equal(nrow(mut_data), 10)
  expect_equal(
    mut_data$variation_type,
    c("no_variant", "snv", "no_variant", "insertion", "snv",
      "no_variant", "mnv", "snv", "deletion", "no_variant")
  )
  expect_equal(mut_data$vaf, mut_data$alt_depth / mut_data$total_depth)
})
