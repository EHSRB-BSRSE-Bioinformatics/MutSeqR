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
