library(testthat)

test_that("VCF ALT values are converted without expanding records", {
  alt <- IRanges::CharacterList(
    "T",
    c("A", "C"),
    "<DEL>",
    character()
  )

  expect_identical(
    MutSeqR:::vcf_alt_to_character(alt),
    c("T", "A,C", "<DEL>", NA_character_)
  )
})

test_that("import_vcf_data preserves alleles and depths from the VCF records", {
  file <- file.path("./testdata/simple_vcf_data.vcf")

  mut_data <- import_vcf_data(vcf_file = file,
    BS_genome = "BSgenome.Mmusculus.UCSC.mm10",
    add_chr = TRUE
  )
  colnames <- c(
    MutSeqR::op$base_required_mut_cols,
    MutSeqR::op$processed_required_mut_cols, # subtype/context cols
    "total_depth", "ref_depth", "vaf", # depth cols
    "nchar_ref", "nchar_alt", "varlen",
    "gc_content", "row_has_duplicate",
    "strand", "width", # added by GRanges
    "AD_1", "AD_2" # vcf cols
  )
  expect_named(mut_data, colnames, ignore.order = TRUE) # check columns
  expect_equal(nrow(mut_data), 10)
  expect_type(mut_data$ref, "character")
  expect_type(mut_data$alt, "character")
  expect_equal(as.character(mut_data$contig), rep(c("chr1", "chr2"), each = 5))
  expect_equal(
    mut_data$start,
    c(5819109, 5819110, 5819111, 5819112, 5819113,
      148627878, 148627879, 148627880, 148627881, 148627882)
  )
  expect_equal(
    mut_data$ref,
    c("G", "A", "T", "G", "C", "A", "AGA", "G", "AA", "A")
  )
  expect_equal(
    mut_data$alt,
    c("G", "T", "T", "GA", "A", "A", "TCA", "A", "A", "A")
  )
  expect_equal(mut_data$AD_1, c(100, 149, 145, 297, 2, 300, 558, 119, 349, 540))
  expect_equal(mut_data$AD_2, c(0, 1, 0, 3, 158, 0, 2, 1, 1, 0))
  # alt_depth comes from VD when present (when VD and AD differ - row 1 alt_depth=5)
  expect_equal(mut_data$alt_depth, c(5, 1, 0, 3, 1, 0, 2, 1, 1, 0))
  expect_equal(
    mut_data$total_depth,
    c(100, 150, 145, 300, 160, 300, 560, 120, 350, 540)
  )
  expect_equal(
    mut_data$ref_depth,
    c(95, 149, 145, 297, 159, 300, 558, 119, 349, 540)
  )
  expect_equal(
    mut_data$variation_type,
    c("no_variant", "snv", "no_variant", "insertion", "snv",
      "no_variant", "mnv", "snv", "deletion", "no_variant")
  )
  expect_equal(mut_data$end[mut_data$start == 5819112], 5819113)
  expect_equal(
    mut_data$vaf,
    c(5/100, 1 / 150, 0, 3 / 300, 1 / 160, 0, 2 / 560, 1 / 120, 1 / 350, 0)
  )
})

test_that("import_vcf_data respects INFO END and derives SV end from SVLEN", {
  file <- file.path("./testdata/structural_vcf_data.vcf")

  mut_data <- import_vcf_data(
    vcf_file = file,
    BS_genome = "BSgenome.Hsapiens.UCSC.hg38",
  )

  expect_equal(mut_data$variation_type, "sv")
  expect_equal(mut_data$start, 23665136)
  expect_equal(mut_data$end, 23666093)
  expect_equal(mut_data$alt, "<DEL>")
  expect_equal(mut_data$variation_type, "sv")
})

test_that("import_vcf_data does not coerce a missing BS_genome into an installed check", {
  file <- file.path("./testdata/simple_vcf_data.vcf")

  expect_error(
    import_vcf_data(vcf_file = file),
    "no BS_genome was provided"
  )
})

test_that("alt_depth and total_depth come from AD when VD is missing", {
  expect_warning(
    dat <- import_vcf_data(test_path("testdata", "vcf_ad_depth.vcf")),
    "1 records have an incomplete AD; their total_depth is set to DP."
  )

  expect_equal(dat$variation_type, c("no_variant", "snv", "mnv", "deletion"))
  expect_equal(dat$alt_depth, c(0, 2, 3, 1))
  # The MNV's AD lacks its REF depth, so its total_depth is DP.
  expect_equal(dat$total_depth, c(100, 100, 97, 96))
  expect_equal(dat$vaf, dat$alt_depth / dat$total_depth)
})

test_that("alt_depth comes from the ALT depths of AD", {
  dat <- data.frame(
    variation_type = c("no_variant", "snv", "snv"),
    AD_1 = c(100, 98, 50),
    AD_2 = c(NA, 2, NA)
  )

  expect_equal(
    MutSeqR:::alt_depth_from_ad(dat, c("AD_1", "AD_2")),
    c(0, 2, NA)
  )
  expect_null(MutSeqR:::alt_depth_from_ad(dat, "AD_1"))
})

test_that("an AD missing its REF depth does not undercount total_depth", {
  dat <- data.frame(
    variation_type = c("snv", "mnv"),
    depth = c(100, 60),
    AD_1 = c(95, NA),
    AD_2 = c(5, 5)
  )

  expect_warning(
    total_depth <- MutSeqR:::estimate_total_depth(dat),
    "1 records have an incomplete AD"
  )
  expect_equal(total_depth, c(100, 60))
})

test_that("FORMAT SITE_DP is read as total_depth, and its gaps are filled from AD or DP", {
  expect_warning(
    dat <- import_vcf_data(test_path("testdata", "vcf_site_dp.vcf")),
    "1 records have an incomplete AD; their total_depth is set to DP."
  )

  expect_equal(dat$start, c(10, 11, 12, 13, 14, 900, 950))
  # SITE_DP where present; otherwise the AD sum, or DP where AD lacks REF.
  expect_equal(dat$total_depth, c(100, 100, 99, 3, 96, 50, 60))
  expect_equal(dat$vaf, dat$alt_depth / dat$total_depth)
})

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
  # Use BGZF rather than ordinary gzip for portable htslib header rewinding.
  # Exercise both supported compressed suffixes.
  suffixes <- c(".gz", ".bgz")
  for (i in seq_along(files)) {
    Rsamtools::bgzip(files[i],
      dest = file.path(directory, paste0(basename(files[i]), suffixes[i]))
    )
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

