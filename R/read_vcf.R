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
#' @importFrom  VariantAnnotations read_vcf info() geno()
#' @importFrom dplyr mutate select rename
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stringr tolower
#' 
# "C:/Users/ADODGE/OneDrive - HC-SC PHAC-ASPC/Documents/DupSeq R Package Building/Test Data/vcf files/PRC_ST_113.1.consensus.variant-calls.genome.vcf"

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
  qual = info(vcf)$QUAL,
  depth = geno(vcf)$DP[, c(1)],
  alt_depth = geno(vcf)$VD[, c(1)],
  variation_type = info(vcf)$TYPE,
  SVTYPE = info(vcf)$SVTYPE
 )

# Clean data
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

dat$variation_type <- stringr::tolower(dat$variation_type)
 
dat <- dat %>%
  dplyr::mutate(
    ref_depth = .data$depth - .data$alt_depth,
    subtype = 
      ifelse(.data$variation_type == "SNV",
             paste0(.data$ref, ">", .data$alt.value),
             "."),
    context = 
      ifelse(.data$variation_type != "no_variant", 
             paste0(substr(.data$LSEQ, 20, 20), ref, substr(.data$RSEQ, 1, 1)),
             ".")) %>%
  select(-.data$LSEQ, -.data$Rseq)
# change the variation type to lower case to match mut file

  
  
  

  
  

