test_that("get_seq retrieves sequences and creates GRanges object", {
  tmpfile <- tempfile(fileext = ".txt")
  write.table(
    data.frame(
      contig = c("chr11", "chr13"),
      start = c(108510788, 75803913),
      end = c(108510791, 75803916),
      gene = c("GeneA", "GeneB"),
      transcription_status = c("genic", "intergenic")
    ),
  file = tmpfile,
  sep = "\t", row.names = FALSE
  )

  gr <- get_seq(regions = tmpfile,
                species = "mouse",
                genome = "mm10",
                is_0_based_rg = FALSE)

  # Check if the result is a GRanges object
  expect_true(is(gr, "GRanges"))

  # Check if the extra metadata columns are retained
  expect_equal(gr$gene, c("GeneA", "GeneB"))
  expect_equal(gr$transcription_status, c("genic", "intergenic"))

  # Check if the sequence data is added to the GRanges object
  expect_equal(as.character(gr$sequence), c("GGTT", "ACAC"))

   unlink(tmpfile)
})
