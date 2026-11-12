test_that("import_vcf_data keeps FILTER in a filter column", {
  dat <- import_vcf_data(test_path("testdata", "vcf_filter.vcf"))

  expect_equal(dat$filter, c(".", "PASS", "NEAR_READ_END"))
})

test_that("filter_mut can flag records by their FILTER", {
  dat <- import_vcf_data(test_path("testdata", "vcf_filter.vcf"))
  filtered <- filter_mut(
    dat,
    custom_filter_col = "filter",
    custom_filter_val = "NEAR_READ_END"
  )

  expect_equal(filtered$filter_mut, c(FALSE, FALSE, TRUE))
})

test_that("filter_mut flags a FILTER label among several", {
  dat <- import_vcf_data(test_path("testdata", "vcf_filter.vcf"))
  dat$filter[3] <- "PASS;NEAR_READ_END"
  filtered <- filter_mut(
    dat,
    custom_filter_col = "filter",
    custom_filter_val = "(^|;)NEAR_READ_END(;|$)",
    custom_filter_match = "regex"
  )

  expect_equal(filtered$filter_mut, c(FALSE, FALSE, TRUE))
})

test_that("filter_mut warns of a FILTER label no record carries", {
  dat <- import_vcf_data(test_path("testdata", "vcf_filter.vcf"))
  expect_warning(
    filtered <- filter_mut(
      dat,
      custom_filter_col = "filter",
      custom_filter_val = c("NEAR_READ_END", "STRAND_BIAS")
    ),
    "No rows in <filter> match: STRAND_BIAS"
  )

  expect_equal(filtered$filter_mut, c(FALSE, FALSE, TRUE))
})

test_that("filter_mut asks for the values of a custom filter column", {
  dat <- import_vcf_data(test_path("testdata", "vcf_filter.vcf"))

  expect_error(
    filter_mut(dat, custom_filter_col = "filter"),
    "did not specify the filter value"
  )
})

test_that("filter_mut treats a dot as a literal custom filter value", {
  dat <- import_vcf_data(test_path("testdata", "vcf_filter.vcf"))
  dat$filter <- c("PASS", ".", "LOW_HICNT;NEAR_READ_END")

  filtered <- filter_mut(
    dat,
    custom_filter_col = "filter",
    custom_filter_val = "."
  )

  expect_equal(filtered$filter_mut, c(FALSE, TRUE, FALSE))
})

test_that("filter_mut matches multiple values in compound FILTER entries", {
  dat <- import_vcf_data(test_path("testdata", "vcf_filter.vcf"))
  dat$filter[2] <- "PASS;NEAR_READ_END"
  dat$filter[3] <- "NEAR_READ_END;PASS"

  filtered <- filter_mut(
    dat,
    custom_filter_col = "filter",
    custom_filter_val = c(
      "(^|;)NEAR_READ_END(;|$)",
      "(^|;)PASS(;|$)"
    ),
    custom_filter_match = "regex"
  )

  expect_equal(filtered$filter_mut, c(FALSE, TRUE, TRUE))
})

test_that("filter_mut rejects an empty custom filter value vector", {
  dat <- import_vcf_data(test_path("testdata", "vcf_filter.vcf"))

  expect_error(
    filter_mut(
      dat,
      custom_filter_col = "filter",
      custom_filter_val = character()
    ),
    "must contain at least one value"
  )
})

test_that("filter_mut exact matching distinguishes real allele and depth values", {
  dat <- utils::read.delim(
    test_path("testdata", "simple_mut_import.txt"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  alt_filtered <- filter_mut(
    dat,
    custom_filter_col = "alt",
    custom_filter_val = "A"
  )
  depth_filtered <- filter_mut(
    dat,
    custom_filter_col = "alt_depth",
    custom_filter_val = "1"
  )

  expect_equal(alt_filtered$filter_mut, dat$alt == "A")
  expect_equal(
    depth_filtered$filter_mut,
    as.character(dat$alt_depth) == "1"
  )
  expect_false(any(alt_filtered$filter_mut[dat$alt == "ACTG"]))
  expect_false(any(depth_filtered$filter_mut[dat$alt_depth == 10]))
})

test_that("filter_mut preserves substring matching by default", {
  dat <- utils::read.delim(
    test_path("testdata", "simple_mut_import.txt"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  filtered <- filter_mut(
    dat,
    custom_filter_col = "alt_depth",
    custom_filter_val = "1"
  )

  expect_equal(
    filtered$filter_mut,
    grepl("1", as.character(dat$alt_depth), fixed = TRUE)
  )
  expect_true(any(filtered$filter_mut[dat$alt_depth == 10]))
})
