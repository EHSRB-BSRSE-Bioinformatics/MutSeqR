test_that("get_seq retrieves sequences and creates GRanges object", {
  
  regions_df <- data.frame(
    contig = c("chr11", "chr13"),
    start = c(108510788, 75803913),
    end = c(108510791, 75803916),
    gene = c("GeneA", "GeneB"),
    transcription_status = c("genic", "intergenic")
  )
  
  species <- "human"
  gr <- get_seq(species, regions_df)
  
  # Check if the result is a GRanges object
  expect_true(is(gr, "GRanges"))
  
  # Check if the extra metadata columns are retained
  expect_identical(gr$gene, regions_df$gene)
  expect_identical(gr$transcription_status, regions_df$transcription_status)
  
  # Check if the sequence data is added to the GRanges object
  expect_identical(as.character(gr$sequence), c("GTTT", "AGAA"))
})
