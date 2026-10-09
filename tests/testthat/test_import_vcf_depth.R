test_that("single-file depth precedence is preserved", {
  directory <- withr::local_tempdir()
  fields <- list(VD = c("1", "2"), AD = c("99,1", "198,2"),
                 DP = c("500", "600"), no_calls = c("20", "30"))
  file <- write_depth_vcf(directory, "explicit.vcf", format_fields = fields,
    info_fields = list(total_depth = c("150", "250")))
  expect_depth_values(import_vcf_data(file), c(150, 250))

  file <- write_depth_vcf(directory, "nocalls.vcf", format_fields = fields)
  expect_depth_values(import_vcf_data(file), c(480, 570))

  file <- write_depth_vcf(directory, "ad.vcf", format_fields = fields[c("VD", "AD", "DP")])
  expect_depth_values(import_vcf_data(file), c(100, 200))

  file <- write_depth_vcf(directory, "dp.vcf", format_fields = fields[c("VD", "DP")])
  expect_warning(dat <- import_vcf_data(file), "set to DP")
  expect_depth_values(dat, c(500, 600))
})

test_that("sample and region metadata can still supply depth", {
  directory <- withr::local_tempdir()
  file <- write_depth_vcf(directory, "metadata.vcf", format_fields = list(VD = c("1", "2")))
  metadata <- data.frame(sample = "sample1", total_depth = 150)
  expect_depth_values(import_vcf_data(file, sample_data = metadata), c(150, 150))

  regions <- data.frame(contig = "chr1", start = 1, end = 100, total_depth = 175)
  expect_depth_values(import_vcf_data(file, regions = regions), c(175, 175))
})

test_that("annotated tables retain the parent import schema and numeric types", {
  directory <- withr::local_tempdir()
  file <- write_depth_vcf(directory, "annotated.vcf")
  metadata <- data.frame(sample = "sample1", dose = 10)
  regions <- data.frame(contig = "chr1", start = 1, end = 100, label = "target")
  dat <- import_vcf_data(file, sample_data = metadata, regions = regions)
  # Verified against the pre-refactor importer, not inferred from new helpers.
  expect_named(dat, c(
    "contig", "start", "end", "width", "strand", "ref", "alt",
    "alt_depth", "AD_1", "AD_2", "context", "sample", "dose", "label",
    "in_regions", "variation_type", "nchar_ref", "nchar_alt", "varlen",
    "short_ref", "normalized_ref", "subtype", "normalized_subtype",
    "normalized_context", "context_with_mutation", "normalized_context_with_mutation",
    "gc_content", "filter_mut", "total_depth", "row_has_duplicate", "vaf", "ref_depth"
  ))
  expect_type(dat$alt_depth, "integer")
  expect_type(dat$total_depth, "double")
  expect_depth_values(dat, c(100, 200))
})

test_that("AD arithmetic and absent alternate depth defaults stay unchanged", {
  directory <- withr::local_tempdir()
  file <- write_depth_vcf(directory, "partial.vcf",
    format_fields = list(VD = c("1", "0"), AD = c("99,.", ".,.")))
  dat <- import_vcf_data(file)
  expect_depth_values(dat, c(99, 0), c(1, 0))
  expect_true(is.nan(dat$vaf[2]))

  file <- write_depth_vcf(directory, "default.vcf", format_fields = list(AD = c("99,1", "198,2")))
  expect_depth_values(import_vcf_data(file), c(100, 200), c(1, 1))
})

test_that("homogeneous directory depths match individual imports", {
  directory <- withr::local_tempdir()
  first <- write_depth_vcf(directory, "a.vcf", sample = "sampleA")
  second <- write_depth_vcf(directory, "b.vcf", sample = "sampleB",
    format_fields = list(VD = c("3", "4"), AD = c("297,3", "396,4")))
  expected <- dplyr::bind_rows(import_vcf_data(first), import_vcf_data(second))
  dat <- import_vcf_data(directory)
  expect_equal(dat, expected)
  mf <- calculate_mf(dat)
  expect_equal(mf$group_depth[match(c("sampleA", "sampleB"), mf$sample)], c(300, 700))
})

test_that("global site depth correction stays unchanged across files", {
  directory <- withr::local_tempdir()
  write_depth_vcf(directory, "a.vcf", starts = c(10, 20))
  write_depth_vcf(directory, "b.vcf", starts = c(10, 30),
    refs = c("CC", "T"), alts = c("C", "C"),
    format_fields = list(VD = c("3", "4"), AD = c("297,3", "396,4")))
  dat <- suppressWarnings(import_vcf_data(directory))
  expect_equal(dat$total_depth, c(100, 200, 300, 400))
  expect_equal(dat$row_has_duplicate, c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(calculate_mf(dat, correct_depth = FALSE)$group_depth, 1000)
  expect_equal(calculate_mf(dat, correct_depth = TRUE)$group_depth, 700)
  expect_equal(calculate_mf(dat, correct_depth_by_indel_priority = TRUE)$group_depth, 900)
})

test_that("entirely depth-free files still import without derived depths", {
  directory <- withr::local_tempdir()
  write_depth_vcf(directory, "a.vcf", sample = "sampleA", format_fields = list(VD = c("1", "2")))
  write_depth_vcf(directory, "b.vcf", sample = "sampleB", format_fields = list(VD = c("3", "4")))
  expect_warning(dat <- import_vcf_data(directory), "appropriate depth column")
  expect_false(any(c("total_depth", "ref_depth", "vaf") %in% names(dat)))
  expect_setequal(dat$sample, c("sampleA", "sampleB"))
})
