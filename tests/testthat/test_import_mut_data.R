library(testthat)

# Define test cases for import_mut_data function
test_that("import_mut_data function correctly imports mutation data from default complete file", {
  # Create temporary test file with example mutation data
  tmpfile <- tempfile(fileext = ".mut")
  write.table(
    data.frame(
      sample = c("mouse1", "mouse2", "mouse1", "mouse2"),
      contig = c("chr1", "chr1", "chr2", "chr2"),
      start = c(69304225, 69304240, 50833424, 50833439),
      end = c(69304226, 69304241, 50833425, 50833440),
      context = c("GCA", "GGC", "ATC", "AAC"),
      subtype = c("C>T", "G>A", ".", "."),
      variation_type = c("snv", "snv", "no_variant", "indel"),
      total_depth = c(50, 100, 75, 150),
      alt_depth = c(10, 20, 30, 50),
      ref = c("C", "G", "T", "AA"),
      alt = c("T", "A", ".", "A" )
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
  mut_data <- import_mut_data(mut_file = tmpfile, regions_file = "custom", custom_regions_file = tmpfile2)
  
  expect_true(is(mut_data, "GRanges"), info = "Check if the resulting object is a granges object")
  expect_equal(NROW(mut_data), 4, info = "Check if the resulting object has the correct number of rows")
   
   # Ranges: start+1 because we changed from 0based to 1based. 
  expect_equal(GenomicRanges::start(mut_data), c(69304225, 69304240, 50833424, 50833439) +1, info = "Check if the start ranges are correct")
  expect_equal(GenomicRanges::end(mut_data), c(69304226, 69304241, 50833425, 50833440), info = "Check if the end ranges are correct")
   
  expect_equal(names(mcols(mut_data)), c(
    "sample", "context", "subtype", "variation_type", "total_depth",
    "alt_depth", "ref", "alt", "ref_depth", "context_with_mutation",
    "normalized_context", "normalized_subtype", "short_ref", "normalized_ref", 
    "normalized_context_with_mutation", "gc_content", "VAF", 
    "description", "location_relative_to_genes"
  ), info = "Check if the resulting object has the correct meta data column names")
  
  expect_equal(sapply(mcols(mut_data), class),
    c(sample = "character", context = "character", subtype = "character", variation_type = "character", 
                        total_depth = "integer", alt_depth = "integer", ref = "character", alt = "character", ref_depth = "integer",
                        context_with_mutation = "character", normalized_context = "character", normalized_subtype = "character", 
                        short_ref = "character", normalized_ref = "character", normalized_context_with_mutation = "character", 
                        gc_content = "numeric", VAF = "numeric", description = "character", location_relative_to_genes = "character" ),
    info = " Check if the resulting object has the correct data type for each metadata column" )
  
  # Check that the output is correct
  expect_equal(mut_data$ref_depth, c(40, 80, 45, 100), info = "Check if the ref_depth values are correct")
  expect_equal(mut_data$context_with_mutation, c("G[C>T]A", "G[G>A]C", "no_variant", "indel"), info = "Check if the context_with_mutation values are correct")
  expect_equal(mut_data$normalized_context, c("GCA", "GCC", "ATC", "GTT"), info = "Check if the normalized_context values are correct")
  expect_equal(mut_data$normalized_subtype, c("C>T", "C>T", "no_variant", "indel"), info = "Check if the normalized_subtype values are correct")
  expect_equal(mut_data$short_ref, c("C", "G", "T", "A"), info = "Check if the short_ref values are correct")
  expect_equal(mut_data$normalized_ref, c("C", "C", "T", "T"), info = "Check if the normalized_ref values are correct")
  expect_equal(mut_data$normalized_context_with_mutation, c("G[C>T]A", "G[C>T]C", "no_variant", "indel"), info = "Check if the normalized_context_with_mutation values are correct")
  expect_equal(mut_data$gc_content, c(2/3, 1, 1/3, 1/3), info = "Check if the gc_content values are correct")
  expect_equal(mut_data$VAF, c(10/50, 20/100, 30/75, 50/150), info = "Check if the VAF values are correct")
  expect_equal(mut_data$description, c("region_330", "region_330", "region_4547", "region_4547"), info = "Check if the description values are correct")
  expect_equal(mut_data$location_relative_to_genes, c("intergenic", "intergenic", "intergenic", "intergenic"), info = "Check if the location_relative_to_genes values are correct")
  
  # Clean up temporary file
  unlink(tmpfile)
  unlink(tmpfile2)
})

test_that("import_mut_data function fails to import mutation data from an empty file", {
  # Create temporary empty test file
  tmpfile <- tempfile(fileext = ".mut")
  
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
  expect_error(import_mut_data(mut_file = tmpfile, regions_file = "custom", custom_regions_file = tmpfile2),
               "Error: You are trying to import an empty file/folder.", 
               info = "Check if we get an error message when imported file is empty")
  
  # Clean up temporary file
  unlink(tmpfile)
  unlink(tmpfile2)
})

test_that("import_mut_data function fails to import mutation data from an incomplete file", {
  # Create temporary test file with incomplete mutation data - missing 'end' column
  tmpfileA <- tempfile(fileext = ".mut")
  write.table(
    data.frame(
      #contig = c("chr1", "chr1", "chr2", "chr2"),
      start = c(69304225, 69304240, 50833424, 50833439),
      end = c(69304226, 69304241, 50833425, 50833440),
      context = c("GCA", "GGC", "ATC", "AAC"),
      subtype = c("C>T", "G>A", ".", "."),
      variation_type = c("snv", "snv", "no_variant", "indel"),
      total_depth = c(50, 100, 75, 150),
      alt_depth = c(10, 20, 30, 50),
      ref = c("C", "G", "T", "AA")
    ),
    file = tmpfileA,
    sep = "\t", row.names = FALSE
  )
  
  tmpfileB <- tempfile(fileext = ".mut")
  write.table(
    data.frame(
      #contig = c("chr1", "chr1", "chr2", "chr2"),
      start = c(69304225, 69304240, 50833424, 50833439),
      #end = c(69304226, 69304241, 50833425, 50833440),
      context = c("GCA", "GGC", "ATC", "AAC"),
      subtype = c("C>T", "G>A", ".", "."),
      variation_type = c("snv", "snv", "no_variant", "indel"),
      total_depth = c(50, 100, 75, 150),
      alt_depth = c(10, 20, 30, 50),
      ref = c("C", "G", "T", "AA")
    ),
    file = tmpfileB,
    sep = "\t", row.names = FALSE
  )
  
  tmpfileC <- tempfile(fileext = ".mut")
  write.table(
    data.frame(
      #contig = c("chr1", "chr1", "chr2", "chr2"),
      start = c(69304225, 69304240, 50833424, 50833439),
      #end = c(69304226, 69304241, 50833425, 50833440),
      context = c("GCA", "GGC", "ATC", "AAC"),
      subtype = c("C>T", "G>A", ".", "."),
      variation_type = c("snv", "snv", "no_variant", "indel"),
      #total_depth = c(50, 100, 75, 150),
      alt_depth = c(10, 20, 30, 50),
      ref = c("C", "G", "T", "AA")
    ),
    file = tmpfileC,
    sep = "\t", row.names = FALSE
  )
  
  tmpfileD <- tempfile(fileext = ".mut")
  write.table(
    data.frame(
      #contig = c("chr1", "chr1", "chr2", "chr2"),
      start = c(69304225, 69304240, 50833424, 50833439),
      #end = c(69304226, 69304241, 50833425, 50833440),
      context = c("GCA", "GGC", "ATC", "AAC"),
      subtype = c("C>T", "G>A", ".", "."),
      variation_type = c("snv", "snv", "no_variant", "indel"),
      depth = c(50, 100, 75, 150),
      alt_depth = c(10, 20, 30, 50),
      ref = c("C", "G", "T", "AA")
    ),
    file = tmpfileD,
    sep = "\t", row.names = FALSE
  )
  
  tmpfileE <- tempfile(fileext = ".mut")
  write.table(
    data.frame(
      contig = c("chr1", "chr1", "chr2", "chr2"),
      start = c(69304225, 69304240, 50833424, 50833439),
      end = c(69304226, 69304241, 50833425, 50833440),
      context = c("GCA", "GGC", "ATC", "AAC"),
      subtype = c("C>T", "G>A", ".", "."),
      variation_type = c("snv", "snv", "no_variant", "indel"),
      total_depth = c(50, 100, 75, 150),
      alt_depth = c(10, 20, 30, 50),
      ref = c("C", "G", "T", "AA")
    ),
    file = tmpfileE,
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
  expect_error(import_mut_data(mut_file = tmpfileA, regions_file = "custom", custom_regions_file = tmpfile2), 
               "Required column(s) missing: contig", fixed=TRUE,
               info = "Check if we get an error message when imported file has 1 missing column")
  
  expect_error(import_mut_data(mut_file = tmpfileB, regions_file = "custom", custom_regions_file = tmpfile2), 
               "Required column(s) missing: contig, end", fixed=TRUE,
               info = "Check if we get an error message when imported file has multiple missing columns")
  
  expect_error(import_mut_data(mut_file = tmpfileC, regions_file = "custom", custom_regions_file = tmpfile2), 
               "Required column(s) missing: contig, end, (depth and no_calls) OR total_depth", fixed=TRUE,
               info = "Check if we get an error message when imported file has multiple missing columns, including one of the depth columns")
  
  expect_error(import_mut_data(mut_file = tmpfileD, regions_file = "custom", custom_regions_file = tmpfile2), 
               "Required column(s) missing: contig, end, no_calls", fixed=TRUE,
               info = "Check if we get an error message when imported file has multiple missing columns, including a no_calls column in the presence of a depth column")
  expect_silent(import_mut_data(mut_file = tmpfileE, regions_file = "custom", custom_regions_file = tmpfile2)) #info = "Check if that we get no error if all the required columns are present"
  
  # Clean up temporary file
  unlink(tmpfileA)
  unlink(tmpfileB)
  unlink(tmpfileC)
  unlink(tmpfileD)
  unlink(tmpfileE)
  unlink(tmpfile2)
})
