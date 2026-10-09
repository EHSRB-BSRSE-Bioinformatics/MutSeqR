test_that("heterogeneous files resolve depth before combining fields", {
  directory <- withr::local_tempdir()
  write_depth_vcf(directory, "a.vcf", sample = "sampleA",
    info_fields = list(tag = c("one", "two")))
  write_depth_vcf(directory, "b.vcf", sample = "sampleB",
    format_fields = list(VD = c("3", "4"), DP = c("300", "400")))
  expect_warning(dat <- import_vcf_data(directory), "set to DP")
  expect_equal(dat$sample, c("sampleA", "sampleA", "sampleB", "sampleB"))
  expect_depth_values(dat, c(100, 200, 300, 400), c(1, 2, 3, 4))
  expect_equal(dat$tag, c("one", "two", NA, NA))
  expect_false(any(dat$row_has_duplicate))
  mf <- calculate_mf(dat)
  expect_equal(mf$group_depth[match(c("sampleA", "sampleB"), mf$sample)], c(300, 700))
})

test_that("depth precedence is local to each file's fields", {
  directory <- withr::local_tempdir()
  write_depth_vcf(directory, "a.vcf", sample = "sampleA",
    info_fields = list(total_depth = c("150", "250")))
  write_depth_vcf(directory, "b.vcf", sample = "sampleB")
  write_depth_vcf(directory, "c.vcf", sample = "sampleC",
    format_fields = list(VD = c("3", "4"), DP = c("300", "400"),
                         no_calls = c("10", "20")))
  dat <- import_vcf_data(directory)
  expect_depth_values(dat, c(150, 250, 100, 200, 290, 380), c(1, 2, 1, 2, 3, 4))
})

test_that("metadata supplies missing file-level sources before depth resolution", {
  directory <- withr::local_tempdir()
  write_depth_vcf(directory, "a.vcf", sample = "sampleA",
    info_fields = list(total_depth = c("150", "250")))
  write_depth_vcf(directory, "b.vcf", sample = "sampleB",
    format_fields = list(VD = c("3", "4")))
  metadata <- data.frame(sample = c("sampleA", "sampleB"), total_depth = c(900, 700))
  dat <- import_vcf_data(directory, sample_data = metadata)
  expect_depth_values(dat, c(150, 250, 700, 700), c(1, 2, 3, 4))
})

test_that("AD expansion padding cannot change another file's depth", {
  directory <- withr::local_tempdir()
  write_depth_vcf(directory, "a.vcf", sample = "sampleA")
  write_depth_vcf(directory, "b.vcf", sample = "sampleB",
    alts = c("A,G", "C,G"),
    format_fields = list(VD = c("3", "4"), AD = c("97,1,2", "196,2,2")))
  dat <- import_vcf_data(directory)
  expect_depth_values(dat, c(100, 200, 100, 200), c(1, 2, 3, 4))
  expect_equal(dat$alt, c("A", "C", "A,G", "C,G"))
  expect_equal(dat$AD_3, c(NA, NA, 2, 2))
  expect_equal(nrow(dat), 4)
})

test_that("depth-bearing and depth-free files cannot be silently mixed", {
  directory <- withr::local_tempdir()
  write_depth_vcf(directory, "a.vcf", sample = "sampleA")
  write_depth_vcf(directory, "b.vcf", sample = "sampleB",
    format_fields = list(VD = c("3", "4")))
  expect_error(suppressWarnings(import_vcf_data(directory)), "depth")
})

test_that("missing chosen depth sources fail without switching to AD", {
  directory <- withr::local_tempdir()
  file <- write_depth_vcf(directory, "missing.vcf", sample = "missingSample",
    info_fields = list(total_depth = c(".", "250")))
  expect_error(import_vcf_data(file), "depth")
  expect_error(import_vcf_data(file), "missingSample")
  expect_error(import_vcf_data(file), "missing.vcf")

  file <- write_depth_vcf(directory, "missing_dp.vcf",
    format_fields = list(VD = c("1", "2"), AD = c("99,1", "198,2"),
                         DP = c(".", "400"), no_calls = c("10", "20")))
  expect_error(import_vcf_data(file), "depth")

  file <- write_depth_vcf(directory, "missing_vd.vcf",
    format_fields = list(VD = c(".", "2"), AD = c("99,1", "198,2")))
  expect_error(import_vcf_data(file), "alt_depth")
})

test_that("GRanges and region expansion retain resolved depth", {
  directory <- withr::local_tempdir()
  write_depth_vcf(directory, "a.vcf", sample = "sampleA")
  write_depth_vcf(directory, "b.vcf", sample = "sampleB",
    format_fields = list(VD = c("3", "4"), DP = c("300", "400")))
  regions <- data.frame(contig = "chr1", start = c(1, 5), end = c(100, 100),
                        label = c("one", "two"))
  gr <- suppressWarnings(import_vcf_data(directory, regions = regions, output_granges = TRUE))
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 8)
  expect_equal(sort(gr$total_depth), c(100, 100, 200, 200, 300, 300, 400, 400))
  expect_true(all(gr$in_regions))
  expect_true(all(gr$row_has_duplicate))
})

test_that("mixed INFO and header sample identification stays local", {
  directory <- withr::local_tempdir()
  write_depth_vcf(directory, "a.vcf", sample = "genericHeader",
    info_fields = list(sample = c("sampleA", "sampleA")))
  write_depth_vcf(directory, "b.vcf", sample = "sampleB")
  dat <- import_vcf_data(directory)
  expect_equal(dat$sample, c("sampleA", "sampleA", "sampleB", "sampleB"))
  expect_depth_values(dat, c(100, 200, 100, 200), c(1, 2, 1, 2))
})

test_that("header-only files are skipped but all-empty imports fail clearly", {
  directory <- withr::local_tempdir()
  empty <- write_depth_vcf(directory, "a_empty.vcf", sample = "emptySample")
  lines <- readLines(empty)
  writeLines(lines[startsWith(lines, "#")], empty)
  expect_error(import_vcf_data(directory), "[Ee]mpty|[Nn]o variant")
  expect_error(import_vcf_data(empty), "[Ee]mpty|[Nn]o variant")
  write_depth_vcf(directory, "b.vcf", sample = "sampleB")
  dat <- import_vcf_data(directory)
  expect_equal(nrow(dat), 2)
  expect_equal(dat$sample, c("sampleB", "sampleB"))
})

test_that("empty sites-only VCFs need no sample header", {
  directory <- withr::local_tempdir()
  empty <- file.path(directory, "a_sites_only.vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=1000>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  ), empty)

  expect_error(import_vcf_data(empty), "VCF file contains no variant records")
  expect_error(import_vcf_data(directory), "No variant records")

  named_empty <- write_depth_vcf(directory, "b_empty.vcf")
  lines <- readLines(named_empty)
  writeLines(lines[startsWith(lines, "#")], named_empty)
  expect_error(import_vcf_data(directory), "No variant records")

  populated <- write_depth_vcf(directory, "c_populated.vcf", sample = "sampleC")
  expect_identical(import_vcf_data(directory), import_vcf_data(populated))
})

test_that("nonempty sites-only VCFs still require a named sample", {
  directory <- withr::local_tempdir()
  file <- file.path(directory, "sites_only.vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=1000>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t10\t.\tC\tA\t.\t.\t."
  ), file)

  expect_error(import_vcf_data(file),
    "Expected one named sample in VCF file: sites_only.vcf", fixed = TRUE)
  expect_error(import_vcf_data(directory),
    "Expected one named sample in VCF file: sites_only.vcf", fixed = TRUE)
})

test_that("incompatible INFO field types identify the column and files", {
  directory <- withr::local_tempdir()
  first <- write_depth_vcf(directory, "a.vcf", sample = "sampleA",
    info_fields = list(tag = c("one", "two")))
  second <- write_depth_vcf(directory, "b.vcf", sample = "sampleB",
    info_fields = list(tag = c("1", "2")))
  lines <- readLines(second)
  lines <- sub('ID=tag,Number=1,Type=String', 'ID=tag,Number=1,Type=Integer', lines)
  writeLines(lines, second)
  expect_error(import_vcf_data(directory), "tag")
  expect_error(import_vcf_data(directory), "a.vcf")
  expect_error(import_vcf_data(directory), "b.vcf")
})

test_that("invalid metadata depth text cannot silently become zero AD", {
  directory <- withr::local_tempdir()
  file <- write_depth_vcf(directory, "invalid.vcf", format_fields = list(VD = c("1", "2")))
  metadata <- data.frame(sample = "sample1", AD_1 = "unknown", AD_2 = 1)
  expect_error(import_vcf_data(file, sample_data = metadata), "AD_1")
  expect_error(import_vcf_data(file, sample_data = metadata), "invalid.vcf")
})

test_that("filtering and subtype depth remain usable after heterogeneous import", {
  directory <- withr::local_tempdir()
  write_depth_vcf(directory, "a.vcf", sample = "sampleA")
  write_depth_vcf(directory, "b.vcf", sample = "sampleB",
    format_fields = list(VD = c("3", "4"), DP = c("300", "400")))
  dat <- suppressWarnings(import_vcf_data(directory))
  filtered <- filter_mut(dat, vaf_cutoff = 0.5)
  expect_equal(filtered$total_depth, dat$total_depth)
  mf <- calculate_mf(filtered, subtype_resolution = "base_6")
  expect_true(all(mf$group_depth[mf$sample == "sampleA"] == 300))
  expect_true(all(mf$group_depth[mf$sample == "sampleB"] == 700))
  expect_false(anyNA(mf$group_depth))
})
