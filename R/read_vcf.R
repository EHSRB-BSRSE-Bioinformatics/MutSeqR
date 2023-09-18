#' Import a vcf file
#' 
#' Imports a .vcf file into the R environment and converts it into a dataframe
#' TO DO: vcf files need to first be modified because they contain different column names...
#' The last column name is the sample name. This must be removed in order to bind them all together. 
#' @param vcf_file The .vcf file containing mutation data to be imported.
#' TO DO: get it to import everything in a folder and bind it together.
#' #Add sample data
#' #Add regions file
#' @returns A table where each row is a mutation, and columns indicate the location, type, and other data.
#
<<<<<<< HEAD
#' @importFrom  VariantAnnotations read_vcf info() geno()
#' @importFrom dplyr mutate select rename
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stringr tolower
#' 
# "C:/Users/ADODGE/OneDrive - HC-SC PHAC-ASPC/Documents/DupSeq R Package Building/Test Data/vcf files/PRC_ST_113.1.consensus.variant-calls.genome.vcf"
=======
#' @importFrom VariantAnnotation read_vcf info() geno()
#' @importFrom dplyr mutate select rename
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stringr str_sub str_count
#' 
# "C:/Users/ADODGE/OneDrive - HC-SC PHAC-ASPC/Documents/DupSeq R Package Building/Test Data/vcf files/Small test/PRC_ST_113.1.consensus.variant-calls.genome.vcf"
>>>>>>> 5ba6989ce12106cc390e94eafb559cf554e5f7d8

read_vcf <- function(vcf_file) {
  vcf_file <- file.path(vcf_file)
  if (file.info(vcf_file)$isdir == T) {
    vcf_files <- list.files(path = vcf_file, pattern = "\\.vcf$", full.names = T)
    # Initialize an empty VCF object to store the combined data
    vcf <- NULL
    # Read and combine VCF files
    for (file in vcf_files) {
      vcf_list <- readVcf(file)
      # Ensure consistent sample names and INFO columns here if needed
      if (is.null(vcf)) {
        vcf <- vcf_list
      } else {
        vcf <- rbind(vcf, vcf_list)
      }
    }
  } else {
    vcf <- VariantAnnotation::readVcf(vcf_file)
  }

# Extract mutation data into a dataframe
dat <-data.frame(
  sample = info(vcf)$SAMPLE,
  contig = seqnames(vcf),
  start = start(vcf)-1,
  end = info(vcf)$END,
  ref = ref(vcf),
  alt = alt(vcf),
  LSEQ = info(vcf)$LSEQ,
  RSEQ = info(vcf)$RSEQ,
<<<<<<< HEAD
  qual = info(vcf)$QUAL,
  depth = geno(vcf)$DP[, c(1)],
  alt_depth = geno(vcf)$VD[, c(1)],
  variation_type = info(vcf)$TYPE,
  SVTYPE = info(vcf)$SVTYPE
 )

# Clean data
=======
  qual = info(vcf)$QUAL, #is this the right one?
  depth = geno(vcf)$DP[, c(1)],
  alt_depth = geno(vcf)$VD[, c(1)],
  variation_type = info(vcf)$TYPE,
  SVTYPE = info(vcf)$SVTYPE,
  SVLEN = info(vcf)$SVLEN
 )

# Clean data
#Clean up variation_type column
>>>>>>> 5ba6989ce12106cc390e94eafb559cf554e5f7d8
dat <- dat %>%
  dplyr::mutate(
    variation_type2 = 
      ifelse (.data$variation_type == "Complex", "mnv",
              ifelse(.data$variation_type == "REF", "no_variant", 
                     ifelse(.data$variation_type %in% c("inv", "dup", "del", "ins", "fus"), "sv", 
              .data$variation_type)))
  ) %>%
  dplyr::select(-.data$variation_type) %>%
  dplyr::rename(variation_type = .data$variation_type2)

<<<<<<< HEAD
dat$variation_type <- stringr::tolower(dat$variation_type)
 
=======
dat$variation_type <- tolower(dat$variation_type)

#create ref_depth, subtype, context 
>>>>>>> 5ba6989ce12106cc390e94eafb559cf554e5f7d8
dat <- dat %>%
  dplyr::mutate(
    ref_depth = .data$depth - .data$alt_depth,
    subtype = 
<<<<<<< HEAD
      ifelse(.data$variation_type == "SNV",
             paste0(.data$ref, ">", .data$alt.value),
             "."),
    context = 
      ifelse(.data$variation_type != "no_variant", 
             paste0(substr(.data$LSEQ, 20, 20), ref, substr(.data$RSEQ, 1, 1)),
             ".")) %>%
  select(-.data$LSEQ, -.data$Rseq)
# change the variation type to lower case to match mut file

  
=======
      ifelse(.data$variation_type == "snv",
             paste0(.data$ref, ">", .data$alt.value),
             "."),
    short_ref = substr(.data$ref, 1, 1),
    context = 
      ifelse(.data$variation_type != "no_variant", 
             paste0(substr(.data$LSEQ, 20, 20), ref, substr(.data$RSEQ, 1, 1)), # can do stringr::str_sub(.data$LSEQ, -1, -1), ref, 
             ".")) %>%
  select(-.data$LSEQ, -.data$RSEQ)
# change the variation type to lower case to match mut file

# Define substitution dictionary to normalize to pyrimidine context
sub_dict <- c(
  "G>T" = "C>A", "G>A" = "C>T", "G>C" = "C>G",
  "A>G" = "T>C", "A>C" = "T>G", "A>T" = "T>A"
)
# Clean up data:
# Get reverse complement of sequence context where mutation is listed on purine context
# Change all purine substitutions to pyrimidine substitutions
# Make new column with COSMIC-style 96 base context
# TODO - describe better

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
  mapply(function(x) reverseComplement(x, case = "upper"), .data$context),
  .data$context
),
normalized_subtype = ifelse(
  .data$subtype %in% names(sub_dict),
  sub_dict[.data$subtype],
  .data$subtype
),
normalized_ref = dplyr::case_when(
  substr(.data$ref, 1, 1) == "A" ~ "T",
  substr(.data$ref, 1, 1) == "G" ~ "C",
  substr(.data$ref, 1, 1) == "C" ~ "C",
  substr(.data$ref, 1, 1) == "T" ~ "T"
)
) %>%
  mutate(
    normalized_context_with_mutation =
      ifelse(.data$subtype != ".",
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
    genic_regions <- read.table(custom_regions_file, header = TRUE)
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

>>>>>>> 5ba6989ce12106cc390e94eafb559cf554e5f7d8
  
  

  
  

