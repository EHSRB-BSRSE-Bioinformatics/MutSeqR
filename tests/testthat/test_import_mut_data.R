test_that("import simple_mut_import.txt produces all expected warnings", {
  file <- file.path("./testdata/simple_mut_import.txt")
  regions <- data.frame(
    contig = c("chr1", "chr2"),
    start = c(101, 201),
    end = c(110, 210),
    rg_metadata = c("R1", "R2")
  )
  sampledata <- data.frame(
    sample = c("sample1", "sample2", "sample3"),
    sd_metadata = c("A", "B", "C")
  )
  warn_msgs <- testthat::capture_warnings(
    mut_data <- import_mut_data(
      mut_file = file,
      sample_data = sampledata,
      regions = regions
    )
  )
  expect_true(any(grepl("outside of the specified regions", warn_msgs)))
  expect_true(any(grepl(
    "position was the same as that of at least one other row",
    warn_msgs
  )))
  expect_true(any(grepl("total_depth may be double-counted", warn_msgs)))
  expect_s3_class(mut_data, "data.frame") # check class

  colnames <- c(
    MutSeqR::op$base_required_mut_cols,
    MutSeqR::op$processed_required_mut_cols, # subtype/context cols
    "total_depth",
    "ref_depth",
    "vaf", # depth cols
    "nchar_ref",
    "nchar_alt",
    "varlen",
    "gc_content",
    "row_has_duplicate",
    "rg_metadata",
    "in_regions",
    "sd_metadata", # metadata cols
    "strand",
    "width"
  ) # added by GRanges
  expect_named(mut_data, colnames, ignore.order = TRUE) # check columns
  expect_true(nrow(mut_data) == 36) # check row #
  expect_equal(sum(mut_data$in_regions), 35) # 1 row outside regions
  expect_equal(sum(mut_data$row_has_duplicate), 10) # 10 overlaping positions

  # check classify_variation
  expect_equal(
    mut_data$variation_type,
    c(
      "snv",
      "no_variant",
      "snv",
      "snv",
      "no_variant",
      "no_variant",
      "deletion",
      "no_variant",
      "snv",
      "insertion",
      "sv",
      "no_variant",
      "insertion",
      "mnv",
      "mnv",
      "mnv",
      "snv",
      "no_variant",
      "snv",
      "snv",
      "no_variant",
      "snv",
      "complex",
      "no_variant",
      "snv",
      "sv",
      "ambiguous",
      "no_variant",
      "insertion",
      "snv",
      "no_variant",
      "no_variant",
      "deletion",
      "snv",
      "no_variant",
      "mnv"
    )
  )
  expect_equal(mut_data$vaf, mut_data$alt_depth / mut_data$total_depth)
})

test_that("import_mut_data leaves BS_genome as NULL when context is already present", {
  dat <- data.frame(
    contig = c("chr1", "chr1"),
    start = c(101, 102),
    end = c(102, 103),
    sample = c("sample1", "sample1"),
    ref = c("A", "C"),
    alt = c("T", "G"),
    context = c("TAC", "GCG"),
    alt_depth = c(1, 1),
    total_depth = c(100, 100)
  )

  expect_no_warning(
    expect_no_error(
      out <- import_mut_data(mut_file = dat, BS_genome = NULL)
    )
  )

  expect_s3_class(out, "data.frame")
  expect_true("context" %in% names(out))
})
