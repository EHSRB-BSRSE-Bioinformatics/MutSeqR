#' Import a vcf file
#' 
#' Imports a .vcf file into the R environment and converts it into a dataframe
#' TO DO: vcf files need to first be modified because they contain different column names...
#' The last column name is the sample name. This must be removed in order to bind them all together. 
#' @param vcf_file The .vcf file containing mutation data to be imported.
#' @param sample_data_file An optional file containing additional sample metadata (dose, timepoint, etc.)
#' @param sd_sep The delimiter for importing sample metadata tables
#' @param regions_file "human", "mouse", or "custom". The argument refers to the TS Mutagenesis panel of the specified species, or to a custom panel. If custom, provide file path in custom_regions_file. TO DO: add rat.
#' @param custom_regions_file "filepath". If regions_file is set to custom, provide the file path for the tab-delimited file containing regions metadata. Required columns are "contig", "start", and "end"
#' @param rg_sep The delimiter for importing the custom_regions_file
#' @returns A table where each row is a mutation, and columns indicate the location, type, and other data.
#' @importFrom  VariantAnnotation alt info geno readVcf ref  
#' @importFrom dplyr mutate select rename
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stringr str_sub str_count
#' @importFrom SummarizedExperiment colData
#' @export
#' 
# C:/Users/ADODGE/OneDrive - HC-SC PHAC-ASPC/Documents/DupSeq R Package Building/Test Data/vcf files/Small test
# /PRC_ST_113.1.consensus.variant-calls.genome.vcf

read_vcf <- function(
                      vcf_file,
                      sample_data_file = NULL,
                      sd_sep = "\t",
                      regions_file = c("human", "mouse", "custom"),
                      custom_regions_file = NULL,
                      rg_sep = "\t") {
  
  vcf_file <- file.path(vcf_file)
  if (file.info(vcf_file)$isdir == T) {
    vcf_files <- list.files(path = vcf_file, pattern = "\\.vcf$", full.names = T)
    # Initialize an empty VCF object to store the combined data
    vcf <- NULL
    # Read and combine VCF files
    for (file in vcf_files) {
      vcf_list <- VariantAnnotation::readVcf(file)
      
      # Ensure consistent sample names and INFO columns here if needed
      rownames(SummarizedExperiment::colData(vcf_list)) <- "sample"
      # Combine the VCF data 
      if (is.null(vcf)) {
        vcf <- vcf_list
      } else {
        vcf <- VariantAnnotation::rbind(vcf, vcf_list)
      }
    }
  } else {
    vcf <- VariantAnnotation::readVcf(vcf_file)
  }
  
  # Extract mutation data into a dataframe
  dat <-data.frame(
    sample = VariantAnnotation::info(vcf)$SAMPLE,
    contig = seqnames(vcf),
    start = start(vcf)-1,
    end = VariantAnnotation::info(vcf)$END,
    ref = VariantAnnotation::ref(vcf),
    alt = VariantAnnotation::alt(vcf),
    LSEQ = VariantAnnotation::info(vcf)$LSEQ,
    RSEQ = VariantAnnotation::info(vcf)$RSEQ,
    qual = VariantAnnotation::info(vcf)$QUAL, #is this the right one?
    depth = VariantAnnotation::geno(vcf)$DP[, c(1)],
    alt_depth = VariantAnnotation::geno(vcf)$VD[, c(1)],
    variation_type = VariantAnnotation::info(vcf)$TYPE,
    SVTYPE = VariantAnnotation::info(vcf)$SVTYPE,
    SVLEN = VariantAnnotation::info(vcf)$SVLEN
  )
  
  # Read in sample data if it's provided
  if (!is.null(sample_data_file)) {
    sampledata <- read.delim(file.path(sample_data_file),
                            sep = sd_sep,
                            header = T
                            
    )
    dat <- left_join(dat, sampledata, suffix = c("", ".sampledata"))
  }
  
  # Clean data
  #Clean up variation_type column
  dat <- dat %>%
    dplyr::mutate(
      variation_type = tolower(dat$variation_type),
      variation_type = 
        ifelse(.data$variation_type == "ref", "no_variant", 
          ifelse(.data$variation_type %in% c("inv", "dup", "del", "ins", "fus"), "sv",
            ifelse(.data$variation_type == "insertion" | .data$variation_type == "deletion", "indel",
               .data$variation_type))
          ),
      nchar_ref = nchar(ref),
      nchar_alt = nchar(alt.value),
      variation_type = 
        ifelse( .data$variation_type == "complex" & .data$nchar_ref != .data$nchar_alt, "indel",
           ifelse(.data$variation_type == "complex" & .data$nchar_ref == .data$nchar_alt, "mnv",
             .data$variation_type)
      ),
      INDELTYPE = 
        ifelse(.data$variation_type == "indel" & .data$nchar_ref < .data$nchar_alt, "insertion",
          ifelse(.data$variation_type == "indel" & .data$nchar_ref > .data$nchar_alt, "deletion",
                 ".")),
      INDELLEN = 
        ifelse(.data$variation_type == "indel", abs(.data$nchar_ref - .data$nchar_alt), 
               ".")
    )
    
  #create ref_depth, short_ref, subtype, context 
  dat <- dat %>%
    dplyr::mutate(
      ref_depth = .data$depth - .data$alt_depth,
      subtype = 
        ifelse(.data$variation_type == "snv",
               paste0(.data$ref, ">", .data$alt.value),
               "."),
      short_ref = substr(.data$ref, 1, 1),
      context = 
        ifelse(.data$variation_type %in% c("snv", "sv") | (.data$INDELTYPE == "insertion" & nchar_ref == 1),
               paste0(stringr::str_sub(.data$LSEQ, -1, -1), ref, stringr::str_sub(.data$RSEQ, 1, 1)), 
               ifelse(.data$variation_type =="mnv" | .data$INDELTYPE == "deletion" | .data$INDELTYPE == "insertion" & nchar_ref > 1 , 
                      paste0(stringr::str_sub(.data$LSEQ, -1, -1), stringr::str_sub(.data$ref, 1, 2)), 
                      "."))) %>%
    dplyr::select(-.data$LSEQ, -.data$RSEQ)
  
  # Define substitution dictionary to normalize to pyrimidine context
  sub_dict <- c(
    "G>T" = "C>A", "G>A" = "C>T", "G>C" = "C>G",
    "A>G" = "T>C", "A>C" = "T>G", "A>T" = "T>A"
  )
  # Clean up data:
  # Get reverse complement of sequence context where mutation is listed on purine context
  # Change all purine substitutions to pyrimidine substitutions
  # Make new column with COSMIC-style 96 base context
  # Get GC Content
  # TODO - describe better
  
  # Context with mutation
  dat <- dat %>%
    dplyr::mutate(
      context_with_mutation =
        ifelse(.data$subtype != ".",
               paste0(
                 stringr::str_sub(.data$context, 1, 1),
                 "[", .data$subtype, "]",
                 stringr::str_sub(.data$context, 3, 3)
               ),
               .data$variation_type
        ),
      normalized_context = ifelse(
        stringr::str_sub(.data$context, 2, 2) %in% c("G", "A", "g", "a"),
        mapply(function(x) reverseComplement.default(x, case = "upper"), .data$context),
        .data$context
      ),
      normalized_subtype = 
        ifelse(.data$subtype %in% names(sub_dict),
               sub_dict[.data$subtype],
               .data$subtype
        ),
      normalized_ref = dplyr::case_when(
        substr(.data$ref, 1, 1) == "A" ~ "T",
        substr(.data$ref, 1, 1) == "G" ~ "C",
        substr(.data$ref, 1, 1) == "C" ~ "C",
        substr(.data$ref, 1, 1) == "T" ~ "T"
      ),
      normalized_context_with_mutation =
        ifelse(.data$normalized_subtype != ".",
               paste0(
                 stringr::str_sub(.data$normalized_context, 1, 1),
                 "[", .data$normalized_subtype, "]",
                 stringr::str_sub(.data$normalized_context, 3, 3)
               ),
               .data$variation_type
        ),
      gc_content = (stringr::str_count(string = .data$context, pattern = "G") +
                      stringr::str_count(string = .data$context, pattern = "C"))
      / stringr::str_count(.data$context)
    ) %>%
    mutate(
      normalized_subtype = ifelse(
        .data$normalized_subtype == ".",
        .data$variation_type,
        .data$normalized_subtype
      ),
      subtype = ifelse(
        .data$subtype == ".",
        .data$variation_type,
        .data$subtype
      )
    )
  
  # We don't have a total_depth column yet. 
  # When we do - we will have to add depth_col 
  dat <- dat %>% mutate(VAF = .data$alt_depth / .data$depth)
  
  mut_ranges <- makeGRangesFromDataFrame(
    df = as.data.frame(dat),
    keep.extra.columns = T,
    seqnames.field = "contig",
    start.field = "start",
    end.field = "end",
    starts.in.df.are.0based = TRUE
  )
  
  # Annotate the mut file with additional information about genomic regions in the file
  if (regions_file == "human") {
    genic_regions <- read.table(system.file("extdata", "genic_regions_hg38.txt", package = "DupSeqR"), header = TRUE)
  } else if (regions_file == "mouse") {
    genic_regions <- read.table(system.file("extdata", "genic_regions_mm10.txt", package = "DupSeqR"), header = TRUE)
  } else if (regions_file == "custom") {
    if (!is.null(custom_regions_file)) {
      genic_regions <- read.table(custom_regions_file, header = TRUE, sep = rg_sep)
    } else {
      warning("You must provide a file path to custom_regions_file when regions_file is set to 'custom'.")
    }
  } else {
    warning("Invalid regions_file parameter. Choose from 'human', 'mouse', or 'custom'.")
  }
  
  
  region_ranges <- makeGRangesFromDataFrame(
    df = genic_regions,
    keep.extra.columns = T,
    seqnames.field = "contig",
    start.field = "start",
    end.field = "end",
    starts.in.df.are.0based = TRUE
  )
  
  ranges_joined <- plyranges::join_overlap_left(mut_ranges, region_ranges)
  return(ranges_joined)
}







