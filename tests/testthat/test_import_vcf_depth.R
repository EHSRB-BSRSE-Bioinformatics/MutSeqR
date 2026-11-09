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
