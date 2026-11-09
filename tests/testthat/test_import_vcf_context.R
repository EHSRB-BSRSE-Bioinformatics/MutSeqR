test_that("INFO CONTEXT_3mer is read as the context without a BSgenome", {
  dat <- import_vcf_data(test_path("testdata", "vcf_context_3mer.vcf"))

  expect_equal(dat$context, c("ACA", "CAT"))
  expect_equal(dat$normalized_context, c("ACA", "ATG"))
})
