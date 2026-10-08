# Context is supplied in INFO so these tests need no reference genome package.
# Filenames deliberately differ from sample IDs; the files have unequal sizes
# and share a position to expose both mapping errors and false duplicates.
test_that("directory imports preserve header samples without metadata", {
  dat <- import_vcf_data(test_path("testdata", "vcf_header_samples"))

  expect_equal(dat$sample, c(rep("S-MGT112-01", 2), rep("S-MGT112-02", 3)))
  expect_equal(dat$start, c(10, 20, 10, 30, 40))
  expect_false(any(dat$row_has_duplicate))
  expect_equal(dat$variation_type, c("snv", "snv", "snv", "no_variant", "sv"))
  expect_equal(dat$end, c(10, 20, 10, 30, 45))

  mf <- calculate_mf(dat, cols_to_group = "sample")
  expect_setequal(mf$sample, c("S-MGT112-01", "S-MGT112-02"))
  expect_equal(mf$sum_min[match(c("S-MGT112-01", "S-MGT112-02"), mf$sample)], c(2, 2))
})

test_that("directory header samples join the correct metadata rows", {
  # Reverse metadata order to ensure the join uses names rather than row order.
  metadata <- data.frame(
    sample = c("S-MGT112-02", "S-MGT112-01"),
    dose = c(10, 0),
    group = c("treated", "control")
  )
  dat <- import_vcf_data(
    test_path("testdata", "vcf_header_samples"), sample_data = metadata
  )

  expect_equal(dat$sample, c(rep("S-MGT112-01", 2), rep("S-MGT112-02", 3)))
  expect_equal(dat$dose, c(0, 0, 10, 10, 10))
  expect_equal(dat$group, c("control", "control", "treated", "treated", "treated"))
})

test_that("single-file header fallback is preserved", {
  dat <- import_vcf_data(test_path("testdata", "vcf_header_samples", "01_variants.vcf"))

  expect_equal(dat$sample, rep("S-MGT112-01", 2))
  expect_equal(dat$start, c(10, 20))
})

test_that("directory header samples are retained in GRanges output", {
  gr <- import_vcf_data(
    test_path("testdata", "vcf_header_samples"), output_granges = TRUE
  )

  expect_s4_class(gr, "GRanges")
  expect_equal(gr$sample, c(rep("S-MGT112-01", 2), rep("S-MGT112-02", 3)))
  expect_false(any(gr$row_has_duplicate))
})

test_that("compressed directory VCFs preserve sample identities", {
  directory <- withr::local_tempdir()
  files <- list.files(test_path("testdata", "vcf_header_samples"), full.names = TRUE)
  for (file in files) {
    connection <- gzfile(file.path(directory, paste0(basename(file), ".gz")), "wt")
    writeLines(readLines(file), connection)
    close(connection)
  }
  dat <- import_vcf_data(directory)

  expect_equal(dat$sample, c(rep("S-MGT112-01", 2), rep("S-MGT112-02", 3)))
  expect_false(any(dat$row_has_duplicate))
})

test_that("INFO sample identifiers take precedence over directory headers", {
  directory <- withr::local_tempdir()
  files <- list.files(test_path("testdata", "vcf_header_samples"), full.names = TRUE)
  for (i in seq_along(files)) {
    lines <- readLines(files[i])
    header_row <- which(startsWith(lines, "#CHROM"))
    lines <- append(lines,
      '##INFO=<ID=sample,Number=1,Type=String,Description="Sample identifier">',
      after = header_row - 1L
    )
    record_rows <- which(!startsWith(lines, "#"))
    lines[record_rows] <- sub("context=", paste0("sample=INFO-", i, ";context="), lines[record_rows])
    writeLines(lines, file.path(directory, basename(files[i])))
  }
  dat <- import_vcf_data(directory)

  expect_equal(as.character(dat$sample), c(rep("INFO-1", 2), rep("INFO-2", 3)))
})

test_that("directory imports reject multi-sample VCF headers", {
  directory <- withr::local_tempdir()
  lines <- readLines(test_path("testdata", "vcf_header_samples", "01_variants.vcf"))
  header_row <- which(startsWith(lines, "#CHROM"))
  lines[header_row] <- paste(lines[header_row], "second_sample", sep = "\t")
  record_rows <- which(!startsWith(lines, "#"))
  lines[record_rows] <- paste(lines[record_rows], "1:99,1", sep = "\t")
  writeLines(lines, file.path(directory, "multi_sample.vcf"))

  expect_error(import_vcf_data(directory), "Expected one named sample in VCF file: multi_sample.vcf")
})
