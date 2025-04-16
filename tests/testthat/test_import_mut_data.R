test_that("import simple_mut_import.txt", {
  file <- file.path("./testdata/simple_mut_import.txt")
  #file <- system.file("testdata", "simple_mut_import.txt", package = "MutSeqR")
  regions <- data.frame(contig = c("chr1", "chr2"),
                        start = c(101, 201),
                        end = c(110, 210),
                        rg_metadata = c("R1", "R2"))
  sampledata <- data.frame(sample = c("sample1", "sample2", "sample3"),
                           sd_metadata = c("A", "B", "C"))
  mut_data <- import_mut_data(mut_file = file,
                              sample_data = sampledata,
                              regions = "custom",
                              custom_regions = regions)

  expect_s3_class(mut_data, "data.frame") # check class
  colnames <- c(MutSeqR::op$base_required_mut_cols,
                MutSeqR::op$processed_required_mut_cols, # subtype/context cols
                "total_depth", "ref_depth", "vaf", # depth cols
                "nchar_ref", "nchar_alt", "varlen",
                "gc_content", "row_has_duplicate",
                "rg_metadata", "in_regions", "sd_metadata", # metadata cols
                "strand", "width") # added by GRanges
                # is_known not added without ID col
  expect_named(mut_data, colnames, ignore.order = TRUE) # check columns
  expect_true(nrow(mut_data) == 36) # check row #
  expect_equal(sum(mut_data$in_regions), 35) # 1 row outside regions
  expect_equal(sum(mut_data$row_has_duplicate), 10) # 10 overlaping positions
  # check classify_variation
  expect_equal(mut_data$variation_type,
               c("snv", "no_variant", "snv", "snv", "no_variant",
                 "no_variant", "deletion", "no_variant", "snv", "insertion",
                 "sv", "no_variant", "insertion", "mnv", "mnv", "mnv", "snv",
                 "no_variant", "snv", "snv", "no_variant", "snv", "complex",
                 "no_variant", "snv", "sv", "ambiguous", "no_variant",
                 "insertion", "snv", "no_variant", "no_variant", "deletion",
                 "snv", "no_variant", "mnv"))
  # expect_warning("1 rows were outside of the specified regions. To remove these rows, use the filter_mut() function")
  # expect_warning(mut_data, "10 rows were found whose position was the same as that of at least one other row for the same sample\\.")
  # expect_warning(mut_data, "The total_depth may be double-counted in some instances due to overlapping positions\\. Use the filter_mut\\(\\) function to correct the total_depth for these instances\\.")
})