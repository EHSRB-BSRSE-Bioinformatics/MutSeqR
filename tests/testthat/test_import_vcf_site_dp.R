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
