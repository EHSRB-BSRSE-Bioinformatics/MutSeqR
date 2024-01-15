library(testthat)
library(fs)

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
  mut_data <- import_mut_data(mut_file = tmpfile, 
                              regions = "custom", 
                              custom_regions_file = tmpfile2,
                              vaf_cutoff = 0.1,
                              output_granges = TRUE)
  
  expect_true(is(mut_data, "GRanges"), info = "Check if the resulting object is a granges object")
  expect_equal(NROW(mut_data), 4, info = "Check if the resulting object has the correct number of rows")
   
   # Ranges: start+1 because we changed from 0based to 1based. 
  expect_equal(GenomicRanges::start(mut_data), c(69304225, 69304240, 50833424, 50833439) +1, info = "Check if the start ranges are correct")
  expect_equal(GenomicRanges::end(mut_data), c(69304226, 69304241, 50833425, 50833440), info = "Check if the end ranges are correct")
   
  expect_equal(names(mcols(mut_data)), c(
    "sample", "context", "subtype", "variation_type",
    "alt_depth", "ref", "alt", "nchar_ref", "nchar_alt", "varlen", "total_depth", 
    "VAF", "is_germline",  "ref_depth", "context_with_mutation", 
    "normalized_context", "normalized_subtype", "short_ref", "normalized_ref", 
    "normalized_context_with_mutation", "gc_content", "description", 
    "location_relative_to_genes"
  ), info = "Check if the resulting object has the correct meta data column names")
  
  expect_equal(sapply(mcols(mut_data), class),
    c(sample = "character", context = "character", subtype = "character", 
      variation_type = "character", alt_depth = "integer", 
      ref = "character", alt = "character", nchar_ref = "integer", 
      nchar_alt = "integer", varlen = "integer", total_depth = "integer", 
      VAF = "numeric", is_germline = "logical", ref_depth = "integer", 
      context_with_mutation = "character", normalized_context = "character", 
      normalized_subtype = "character", short_ref = "character", 
      normalized_ref = "character", normalized_context_with_mutation = "character", 
      gc_content = "numeric", description = "character", location_relative_to_genes = "character" ),
    info = " Check if the resulting object has the correct data type for each metadata column" )
  
  # Check that the output is correct
  expect_equal(mut_data$ref_depth, c(40, 80, 45, 100), info = "Check if the ref_depth values are correct")
  expect_equal(mut_data$context_with_mutation, c("G[C>T]A", "G[G>A]C", "no_variant", "deletion"), info = "Check if the context_with_mutation values are correct")
  expect_equal(mut_data$normalized_context, c("GCA", "GCC", "ATC", "GTT"), info = "Check if the normalized_context values are correct")
  expect_equal(mut_data$normalized_subtype, c("C>T", "C>T", "no_variant", "deletion"), info = "Check if the normalized_subtype values are correct")
  expect_equal(mut_data$short_ref, c("C", "G", "T", "A"), info = "Check if the short_ref values are correct")
  expect_equal(mut_data$normalized_ref, c("C", "C", "T", "T"), info = "Check if the normalized_ref values are correct")
  expect_equal(mut_data$normalized_context_with_mutation, c("G[C>T]A", "G[C>T]C", "no_variant", "deletion"), info = "Check if the normalized_context_with_mutation values are correct")
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
  empty_file <- "empty.txt"
  file.create(empty_file)
  
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
  expect_error(import_mut_data(mut_file = empty_file, regions = "custom", custom_regions_file = tmpfile2),
               "Error: You are trying to import an empty file", fixed=TRUE,
               info = "Check if we get an error message when imported file is empty")

  # Clean up temporary file
  unlink(empty_file)
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
      #sample = c("DNA1", "DNA2", "DNA3", "DNA4"),
      alt = c("T", "A", ".", "A"),
      #context = c("GCA", "GGC", "ATC", "AAC"),
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
      contig = c("chr1", "chr1", "chr2", "chr2"),
      start = c(NA, 69304240, 50833424, 50833439),
      end = c(69304226, 69304241, 50833425, 50833440),
      sample = c("DNA1", "DNA2", "DNA3", "DNA4"),
      alt = c("T", "A", ".", "A"),
      context = c("GCA", "GGC", "ATC", "AAC"),
      subtype = c("C>T", "G>A", ".", "."),
      variation_type = c("snv", NA, "no_variant", "indel"),
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
      contig = c("chr1", "chr1", "chr2", "chr2"),
      start = c(69304225, 69304240, 50833424, 50833439),
      end = c(69304226, 69304241, 50833425, 50833440),
      sample = c("DNA1", "DNA2", "DNA3", "DNA4"),
      alt = c("T", "A", ".", "A"),
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
      contig = c("chr1", "chr1", "chr2", "chr2"),
      start = c(69304225, 69304240, 50833424, 50833439),
      end = c(69304226, 69304241, 50833425, 50833440),
      sample = c("DNA1", "DNA2", "DNA3", "DNA4"),
      alt = c("T", "A", ".", "A"),
      context = c("GCA", "GGC", "ATC", "AAC"),
      subtype = c("C>T", "G>A", ".", "."),
      variation_type = c("snv", "snv", "no_variant", "indel"),
      no_calls = c(1, 2, 3, 4),
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
      sample = c("DNA1", "DNA2", "DNA3", "DNA4"),
      alt = c("T", "A", ".", "A"),
      context = c("GCA", "GGC", "ATC", "AAC"),
      subtype = c("C>T", "G>A", ".", "."),
      variation_type = c("snv", "snv", "no_variant", "indel"),
      depth = c(50, 100, 75, 150),
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
  
# Test that the function throws the correct errors
  # Use TRyCatch to catch the errors thrown by the function
  # Define the expected error message
  # Use expect equal to compare the expected error with the actual error. 
  
  error <- tryCatch({
    import_mut_data(mut_file = tmpfileA, vaf_cutoff = 0.1, regions = "custom", custom_regions_file = tmpfile2)},
    error = function(e) e$message)
  expected_error <- "Some required columns are missing or their synonyms are not found:  contig, sample, context"
  expect_equal(error, expected_error, 
               info = "Check if we get an error message when imported file has missing columns")
  expect_error(import_mut_data(mut_file = tmpfileB,vaf_cutoff = 0.1, regions = "custom", custom_regions_file = tmpfile2), 
               "Function stopped: NA values were found within the following required column(s):  start, variation_type . Please confirm that your data is complete before proceeding.", fixed=TRUE,
               info = "Check if we get an error message when imported file NA values in required columns")
  expect_error(import_mut_data(mut_file = tmpfileC, vaf_cutoff = 0.1, regions = "custom", custom_regions_file = tmpfile2), 
               "Required columns are missing or could not be determined: depth column ('depth' and 'no_calls' OR 'total_depth')", fixed=TRUE,
               info = "Check if we get an error message when imported file has multiple missing columns, including one of the depth columns")
  expect_error(import_mut_data(mut_file = tmpfileD,  vaf_cutoff = 0.1, regions = "custom", custom_regions_file = tmpfile2), 
               "Required columns are missing or could not be determined: depth column ('depth' OR 'total_depth')", fixed=TRUE,
               info = "Check if we get an error message when imported file has missing depth column but does have a no_calls column")
  warning <- tryCatch({import_mut_data(mut_file = tmpfileE,  vaf_cutoff = 0.1, regions = "custom", custom_regions_file = tmpfile2)},
                      warning = function(w) w$message)
  expected_warning <- " 'total_depth' column was not found and cannot be calculated without the 'no_calls' column. \n            'Depth_col' will be set as 'depth'. Please review the definitions of each column ~here~ (README)\n            before proceeding"
  expect_equal(warning, expected_warning, 
               info = "Check if we get an warning message when we cannot calculate total_depth")

  # Clean up temporary file
  unlink(tmpfileA)
  unlink(tmpfileB)
  unlink(tmpfileC)
  unlink(tmpfileD)
  unlink(tmpfileE)
  unlink(tmpfileF)
  unlink(tmpfile2)
})

test_that("import_mut_data function fails to import mutation data if an invalid file path is specified", {
  # Create temporary non-empty test file with mutation data
  tmpfile <- tempfile(fileext = ".mut")
  write.table(
    data.frame(
      contig = c("chr1", "chr1", "chr2", "chr2")
    ),
    file = tmpfile,
    sep = "\t", row.names = FALSE
  )
  
  invalid_file <- paste0(file.path(tmpfile), "(1)")
  
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
  expect_error(import_mut_data(mut_file = invalid_file, vaf_cutoff = 0.1,  regions = "custom", custom_regions_file = tmpfile2),
               "Error: The file path you've specified is invalid",
               info = "Check that we get an error message when an incorrect file path is specified")
  
  expect_error(import_mut_data(mut_file = "",  vaf_cutoff = 0.1, regions = "custom", custom_regions_file = tmpfile2),
               "Error: The file path you've specified is invalid",
               info = "Check that we get an error message when a blank path is specified")
  
  expect_error(import_mut_data(mut_file = is.null(), vaf_cutoff = 0.1,  regions = "custom", custom_regions_file = tmpfile2),
               info = "Check if we get an error message when input is NULL")
  
  expect_error(import_mut_data(mut_file = is.na(), vaf_cutoff = 0.1,  regions = "custom", custom_regions_file = tmpfile2),
               info = "Check if we get an error message when input is NA")
  

  # Clean up temporary file
  unlink(tmpfile)
  unlink(tmpfile2)
})

test_that("import_mut_data function correctly imports mutation data when file contains synonymous headings", {
  # Create temporary test file with mutation data
  tmpfile <- tempfile(fileext = ".mut")
  write.table(
    data.frame(
      "SaMpLE_ID " = c("mouse1", "mouse2", "mouse1", "mouse2"),
      "cHromosOMe" = c("chr1", "chr1", "chr2", "chr2"),
      " POs" = c(69304225, 69304240, 50833424, 50833439),
      "  ENd" = c(69304226, 69304241, 50833425, 50833440),
      "fLaNkinG.seQUeNCE  " = c("GCA", "GGC", "ATC", "AAC"),
      "mUtaTiOn.sUBtYpE" = c("C>T", "G>A", ".", "."),
      "muTatioN_tYpe" = c("snv", "snv", "no_variant", "indel"),
      "   DeptH   " = c(50, 100, 75, 150),
      "No_Depth" = c(5, 5, 5, 5),
      "alt_rEAD_DEPtH" = c(10, 20, 30, 50),
      "REFEREnCe" = c("C", "G", "T", "AA"),
      "alTERnATE" = c("T", "A", ".", "A" )
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
  mut_data <- import_mut_data(mut_file = tmpfile, vaf_cutoff = 0.1,  regions = "custom", custom_regions_file = tmpfile2)
  
  expect_equal(names(mut_data), c(
    "contig", "start", "end", "width", "strand", "sample", "context", "subtype", 
    "variation_type", "depth", "no_calls", "alt_depth", "ref", "alt", "nchar_ref", 
    "nchar_alt", "varlen", "total_depth", "VAF", "is_germline", "ref_depth", 
    "context_with_mutation", "normalized_context", "normalized_subtype", 
    "short_ref", "normalized_ref", "normalized_context_with_mutation", 
    "gc_content","description", "location_relative_to_genes"
  ), info = "Check if the resulting object has the correct meta data column names, showing that lower-case/trimming and the synonym dictionary work")
  
  # Check that the output is correct
  expect_equal(mut_data$total_depth, c(45, 95, 70, 145), info = "Check if the total_depth values are correct")
  expect_equal(mut_data$VAF, c(10/45, 20/95, 30/70, 50/145), info = "Check if the VAF values are correct")

  # Clean up temporary file
  unlink(tmpfile)
  unlink(tmpfile2)
})

test_that("import_mut_data function correctly imports mutation data from a folder with complete default files", {
  # Create an empty directory
  test_folder <- file.path(tempdir(), "Temp_test_folder")
  dir.create(test_folder, recursive = TRUE, showWarnings = FALSE, mode = "0777")
  
  # Populate folder with complete files
  file1 <- file.path(test_folder, "file1.txt")
  write.table(
    data.frame(
      sample = c("mouse1", "mouse2"),
      contig = c("chr1", "chr1"),
      start = c(69304225, 69304240),
      end = c(69304226, 69304241),
      context = c("GCA", "GGC"),
      subtype = c("C>T", "G>A"),
      variation_type = c("snv", "snv"),
      total_depth = c(50, 100),
      alt_depth = c(10, 20),
      ref = c("C", "G"),
      alt = c("T", "A")
    ),
    file = file1,
    sep = "\t", row.names = FALSE
  )
  file2 <- file.path(test_folder, "file2.txt")
  write.table(
    data.frame(
      sample = c("mouse1", "mouse2"),
      contig = c("chr2", "chr2"),
      start = c(50833424, 50833439),
      end = c(50833425, 50833440),
      context = c("ATC", "AAC"),
      subtype = c(".", "."),
      variation_type = c("no_variant", "indel"),
      total_depth = c(75, 150),
      alt_depth = c(30, 50),
      ref = c("T", "AA"),
      alt = c(".", "A" )
    ),
    file = file2,
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
  mut_data <- import_mut_data(mut_file = test_folder, vaf_cutoff = 0.1,  regions = "custom", custom_regions_file = tmpfile2, output_granges = TRUE)

  expect_true(is(mut_data, "GRanges"), info = "Check if the resulting object is a granges object")
  expect_equal(NROW(mut_data), 4, info = "Check if the resulting object has the correct number of rows")  
  
  # Clean up temporary file
  unlink(test_folder, recursive = TRUE)
  unlink(tmpfile2)
})

test_that("import_mut_data function correctly throws and error when it imports mutation data from an empty folder", {
  # Create an empty directory
  test_folder <- file.path(tempdir(), "Temp_test_folder1")
  dir.create(test_folder, recursive = TRUE, showWarnings = FALSE, mode = "0777")

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
  expect_error(import_mut_data(mut_file = test_folder, vaf_cutoff = 0.1, regions = "custom", custom_regions_file = tmpfile2),
               "Error: The folder you've specified is empty", fixed=TRUE,
               info = "Check if we get an error message when imported folder is empty")

  # Clean up temporary file
  unlink(test_folder, recursive = TRUE)
  unlink(tmpfile2)
})

test_that("import_mut_data function fails to import mutation data if an invalid folder path is specified", {
  # Create temporary non-empty folder with mutation data
  test_folder <- file.path(tempdir(), "Temp_test_folder1")
  dir.create(test_folder, recursive = TRUE, showWarnings = FALSE, mode = "0777")
  
  invalid_folder <- paste0(file.path(test_folder), "(1)")
  
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
  expect_error(import_mut_data(mut_file = invalid_folder, vaf_cutoff = 0.1, regions = "custom", custom_regions_file = tmpfile2),
               "Error: The file path you've specified is invalid", fixed=TRUE,
               info = "Check that we get an error message when an incorrect folder path is specified")
  
  # Clean up temporary file
  unlink(test_folder, recursive = TRUE)
  unlink(tmpfile2)
})

test_that("import_mut_data function correctly imports mutation data from a folder with some empty files", {
  # Create a test directory
  test_folder <- file.path(tempdir(), "Temp_test_folder")
  dir.create(test_folder, recursive = TRUE, showWarnings = FALSE, mode = "0777")
  
  # Populate folder with 2 empty files
  empty_file1 <- "empty(1).txt"
  file.create(file.path(test_folder, empty_file1))
  
  empty_file2 <- "empty(2).txt"
  file.create(file.path(test_folder, empty_file2))
  
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
  expect_error(import_mut_data(mut_file = test_folder, vaf_cutoff = 0.1,  regions = "custom", custom_regions_file = tmpfile2), 
               "Error: All the files in the specified directory are empty", fixed=TRUE,
               info = "Check that we are unable to process mut data if all of the files in the specified folder are empty")
  
  # Create a non-empty file with mutation data
  empty_file3 <- file.path(test_folder, "non-empty(3).txt")
  write.table(
    data.frame(
      sample = c("mouse1", "mouse2"),
      contig = c("chr1", "chr1"),
      start = c(69304225, 69304240),
      end = c(69304226, 69304241),
      context = c("GCA", "GGC"),
      subtype = c("C>T", "G>A"),
      variation_type = c("snv", "snv"),
      total_depth = c(50, 100),
      alt_depth = c(10, 20),
      ref = c("C", "G"),
      alt = c("T", "A")
    ), 
    file = empty_file3,
    sep = "\t", row.names = FALSE
  )
  
  expect_warning(import_mut_data(mut_file = test_folder, vaf_cutoff = 0.1,  regions = "custom", custom_regions_file = tmpfile2),
                 "Warning: The following files in the specified directory are empty: empty(1).txt, empty(2).txt", fixed=TRUE,
                 info = "Check if we have an empty file in the directory, it will trigger a warning")
  
  # Clean up temporary file
  unlink(test_folder, recursive = TRUE)
  unlink(tmpfile2)
})
