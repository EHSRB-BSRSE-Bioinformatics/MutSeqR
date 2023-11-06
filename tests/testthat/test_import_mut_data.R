library(testthat)

# Define a test case for import_mut_data function
test_that("import_mut_data function correctly imports mutation data", {
  # Create temporary test file with example mutation data
  tmpfile <- tempfile(fileext = ".mut")
  write.table(
    data.frame(
      sample = c("mouse1", "mouse2", "mouse1", "mouse2"),
      chromosome = c("chr1", "chr1", "chr2", "chr2"),
      start = c(69304225, 69304240, 50833424, 50833439),
      end = c(69304226, 69304241, 50833425, 50833440),
      context = c("GCA", "GGC", "ATC", "AAC"),
      subtype = c("C>T", "G>A", "T>G", "."),
      variation_type = c("snv", "snv", "snv", "indel"),
      total_depth = c(50, 100, 75, 150),
      alt_depth = c(10, 20, 30, 50),
      reference = c("C", "G", "T", "AA"),
      alt = c("T", "A", "G", "A" )
    ),
    file = tmpfile,
    sep = "\t", row.names = FALSE
  )
 
#create a temporary custom regions file
  tmpfile2 <- tempfile(fileext = ".mut")
  write.table(
    data.frame(
      contig = c("chr1", "chr2"),
      start = c(69304217, 50833175),
      end = c(69306617, 50835575),
      description = c("region_330", "region_4547"), 
      location_relative_to_genes = c("intergenic", "intergenic")
    ), 
    file = tmpfile2,
    sep = "\t", row.names = FALSE
    )
  
  # Call the import_mut_data function on the test data
  mut_data <- import_mut_data(mut_file = tmpfile, 
                              regions = "custom", 
                              custom_regions_file = tmpfile2,
                              vaf_cutoff = 0.1)
  
    
  expect_true(is(mut_data, "GRanges"), info = "Check if the resulting object is a granges object")

  expect_equal(NROW(mut_data), 4, info = "Check if the resulting object has the correct number of rows")
   
   # Ranges: start+1 because we changed from 0based to 1based. 
  expect_equal(GenomicRanges::start(mut_data), c(69304225, 69304240, 50833424, 50833439) +1, info = "Check if the start ranges are correct")
  expect_equal(GenomicRanges::end(mut_data), c(69304226, 69304241, 50833425, 50833440), info = "Check if the end ranges are correct")
   
  expect_equal(names(mcols(mut_data)), c(
    "sample", "context", "subtype", "variation_type", "total_depth",
    "alt_depth", "ref", "alt", "nchar_ref", "nchar_alt", "VARLEN", "ref_depth", 
    "context_with_mutation", "normalized_context", "normalized_subtype", 
    "short_ref", "normalized_ref", "normalized_context_with_mutation", 
    "gc_content", "VAF", "is_germline", "description", "location_relative_to_genes"
  ), info = "Check if the resulting object has the correct meta data column names")
  
  expect_equal(sapply(mcols(mut_data), class),
    c(sample = "character", context = "character", subtype = "character", 
      variation_type = "character", total_depth = "integer", alt_depth = "integer", 
      ref = "character", alt = "character", nchar_ref = "integer", nchar_alt = "integer",
      VARLEN = "integer", ref_depth = "integer", context_with_mutation = "character", 
      normalized_context = "character", normalized_subtype = "character", 
      short_ref = "character", normalized_ref = "character", 
      normalized_context_with_mutation = "character", gc_content = "numeric", 
      VAF = "numeric", is_germline = "logical", description = "character", 
      location_relative_to_genes = "character" ),
    info = " Check if the resulting object has the correct data type for each metadata column" )
  
  # Clean up temporary file
  unlink(tmpfile)
  unlink(tmpfile2)
})
