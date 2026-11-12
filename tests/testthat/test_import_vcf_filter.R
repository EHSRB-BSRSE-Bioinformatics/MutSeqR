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
  dat$filter[3] <- "LOW_HICNT;NEAR_READ_END"
  filtered <- filter_mut(
    dat,
    custom_filter_col = "filter",
    custom_filter_val = "NEAR_READ_END"
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
    "No rows in <filter> contain: STRAND_BIAS"
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
  dat$filter <- c("PASS", "LOW_HICNT;NEAR_READ_END", "STRAND_BIAS;PASS")

  filtered <- filter_mut(
    dat,
    custom_filter_col = "filter",
    custom_filter_val = c("NEAR_READ_END", "STRAND_BIAS")
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
