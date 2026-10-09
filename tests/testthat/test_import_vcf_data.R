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

# Define a test case for import_mut_data function
test_that("import_vcf_datafunction correctly imports vcf files", {
  # Create temporary test file with example mutation data
  file <- file.path("./testdata/simple_vcf_data.vcf")

  # Call the import_mut_data function on the test data
  mut_data <- import_vcf_data(vcf_file = file,
    regions = NULL,
    BS_genome = "BSgenome.Mmusculus.UCSC.mm10",
    output_granges = FALSE,
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
  expect_equal(mut_data$ref[mut_data$start == 5819110], "A")
  expect_equal(mut_data$alt[mut_data$start == 5819110], "T")
  expect_equal(
    mut_data$variation_type,
    c("no_variant", "snv", "no_variant", "insertion", "snv",
      "no_variant", "mnv", "snv", "deletion", "no_variant")
  )
  expect_equal(mut_data$end[mut_data$start == 5819112], 5819113)
  expect_equal(mut_data$vaf, mut_data$alt_depth / mut_data$total_depth)
})

test_that("import_vcf_data respects INFO END and derives SV end from SVLEN", {
  file <- file.path("./testdata/structural_vcf_data.vcf")

  mut_data <- import_vcf_data(
    vcf_file = file,
    regions = NULL,
    BS_genome = "BSgenome.Hsapiens.UCSC.hg38",
    output_granges = FALSE
  )

  expect_equal(mut_data$variation_type, "sv")
  expect_equal(mut_data$start, 23665136)
  expect_equal(mut_data$end, 23666093)
  expect_type(mut_data$ref, "character")
  expect_type(mut_data$alt, "character")
  expect_equal(mut_data$alt, "<DEL>")
})

test_that("import_vcf_data does not coerce a missing BS_genome into an installed check", {
  file <- file.path("./testdata/simple_vcf_data.vcf")

  expect_error(
    import_vcf_data(vcf_file = file, regions = NULL, output_granges = FALSE),
    "no BS_genome was provided"
  )
})

test_that("import_vcf_data warns when no appropriate depth column is supplied", {
  input_file <- file.path("./testdata/simple_vcf_data.vcf")
  no_depth_file <- tempfile(fileext = ".vcf")
  on.exit(unlink(no_depth_file), add = TRUE)
  lines <- readLines(input_file)
  lines <- lines[!grepl("^##FORMAT=<ID=AD,", lines)]
  record_lines <- !grepl("^#", lines)
  records <- strsplit(lines[record_lines], "\t", fixed = TRUE)
  records <- lapply(records, function(record) {
    record[[9]] <- "VD"
    record[[10]] <- sub(":.*$", "", record[[10]])
    paste(record, collapse = "\t")
  })
  lines[record_lines] <- unlist(records)
  writeLines(lines, no_depth_file)

  warning_messages <- testthat::capture_warnings(
    import_vcf_data(
      vcf_file = no_depth_file,
      BS_genome = "BSgenome.Mmusculus.UCSC.mm10",
      output_granges = FALSE
    )
  )
  expect_true(any(grepl(
    "Could not find an appropriate depth column. Some package functionality may be limited.",
    warning_messages,
    fixed = TRUE
  )), info = "Summary_report.Rmd depends on this warning to use precalculated depth.")
})

test_that("import_vcf_data warns when duplicate positions can double-count depth", {
  input_file <- file.path("./testdata/simple_vcf_data.vcf")
  duplicate_file <- tempfile(fileext = ".vcf")
  on.exit(unlink(duplicate_file), add = TRUE)
  lines <- readLines(input_file)
  record_lines <- lines[!grepl("^#", lines)]
  writeLines(c(lines, record_lines[[1]]), duplicate_file)

  warning_messages <- testthat::capture_warnings(
    import_vcf_data(
      vcf_file = duplicate_file,
      BS_genome = "BSgenome.Mmusculus.UCSC.mm10",
      output_granges = FALSE
    )
  )
  expect_true(any(grepl(
    "The total_depth may be double-counted in some instances due to overlapping positions. Set the correct_depth parameter in calculate_mf() to correct the total_depth for these instances.",
    warning_messages,
    fixed = TRUE
  )), info = "Summary_report.Rmd depends on this warning to enable depth correction.")
})
