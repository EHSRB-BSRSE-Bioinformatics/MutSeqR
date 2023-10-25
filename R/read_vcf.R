#' Import a vcf file
#' 
#' Imports a .vcf file into the R environment and converts it into a dataframe
#' @param vcf_file The path to the .vcf file  to be imported. If you
#' specify a folder, the function will attempt to read all files in the folder and
#' combine them into a single data frame. Multisample vcf files are not supported.
#'  vcf files must contain one sample each.  
#' @param sample_data_file An optional file containing additional sample metadata 
#' (dose, timepoint, etc.)
#' @param sd_sep The delimiter for importing sample metadata tables. Default is "\t".
#' @param regions_file "human", "mouse", or "custom". The argument refers to the
#'  TS Mutagenesis panel of the specified species, or to a custom panel. 
#'  If custom, provide file path in custom_regions_file. TO DO: add rat.
#' @param custom_regions_file "filepath". If regions_file is set to custom,
#'  provide the file path for the file containing regions metadata. 
#'  Required columns are "contig", "start", and "end"
#' @param rg_sep The delimiter for importing the custom_regions_file. Default is "\t".
#' @param assembly The genome assembly. Accepted values: "GRCh37", "GRCh38", "GRCm38", "GRCm39" TO DO: MAKE GENERAL
#' @param depth_calc In the instance when there are two or more calls at the 
#' same location within a sample, and the depths differ, this parameter chooses 
#' the method of calculation for the total_depth. take_mean calculates the 
#' total_depth by taking the mean reference depth and then adding all the alt depths. 
#' take_del calculates the total_depth by choosing only the reference depth of 
#' the deletion in the group, or if no deletion is present, the complex variant,
#'  then adding all alt depths, if there is no deletion or complex variant, 
#'  then it takes the mean of the reference depths. Default is "take_del".
#' @returns A GRanges object where each row is a mutation, and columns indicate the location, type, and other data.
#' @importFrom  VariantAnnotation alt info geno readVcf ref rbind 
#' @importFrom dplyr filter group_by left_join mutate rename select summarize ungroup
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stringr str_sub str_count
#' @importFrom SummarizedExperiment colData
#' @export
#' 
# C:/Users/ADODGE/OneDrive - HC-SC PHAC-ASPC/Documents/DupSeq R Package Building/Test Data/vcf files/Small test
# /PRC_ST_208.1.consensus.variant-calls.genome.vcf
# C:/Users/ADODGE/OneDrive - HC-SC PHAC-ASPC/Documents/DupSeq R Package Building/Test Data/PRC_ST_sample_data.txt
read_vcf <- function(
    vcf_file,
    sample_data_file = NULL,
    sd_sep = "\t",
    regions_file = c("human", "mouse", "custom"),
    custom_regions_file = NULL,
    rg_sep = "\t",
    assembly = c("GRCh37", "GRCh38", "GRCm38", "GRCm39"),
    depth_calc = "take_del"
) {
  vcf_file <- file.path(vcf_file)
  
  # Check if a sample identifier is already present in the INFO field
  check_and_rename_sample <- function(vcf) {
  
    # Define possible variations of sample identifier names
    possible_sample_names <- c("sample","sample_name", "sample_id")
    
    # Initialize the sample_info column
    sample_info <- NULL
    
    # Search for variations of sample identifier names
    for (sample_name_var in possible_sample_names) {
      sample_name_var <- tolower(sample_name_var)  # Make the comparison case-insensitive
      if (sample_name_var %in% tolower(names(info(vcf)))) {
        # If found, rename it to "sample" and break the loop
        names(info(vcf))[tolower(names(info(vcf))) == sample_name_var] <- "sample"
        break
      }
    }
    
    # If "sample" still doesn't exist, create it using colData
    if (!"sample" %in% colnames(info(vcf))) {
      sample_info <- rownames(SummarizedExperiment::colData(vcf))
      info(vcf)$sample <- sample_info
    }
    
    return(vcf)
  }

  # Read and bind vcfs from folder  
  if (file.info(vcf_file)$isdir == TRUE) {
    vcf_files <- list.files(path = vcf_file, pattern = "\\.vcf$", full.names = TRUE)
    
    # Initialize an empty VCF object to store the combined data
    vcf <- NULL
    
    # Read and combine VCF files
    for (file in vcf_files) {
      vcf_list <- readVcf(file)
      
      # Rename or create the "sample" column in the INFO field
      vcf_list <- check_and_rename_sample(vcf_list)
      # Ensure consistent column names
      rownames(SummarizedExperiment::colData(vcf_list)) <- "sample_info" 
      # Combine the VCF data
      if (is.null(vcf)) {
        vcf <- vcf_list
      } else {
        vcf <- rbind(vcf, vcf_list)
      }
    }
  } else {

  # Read a single vcf file
        vcf <- readVcf(vcf_file)
    # Rename or create the "sample" column in the INFO field
    vcf <- check_and_rename_sample(vcf)
  }
  

   # Extract mutation data into a dataframe
  dat <- data.frame(
    sample = VariantAnnotation::info(vcf)$sample,
    contig = SummarizedExperiment::seqnames(vcf),
    start = SummarizedExperiment::start(vcf),
    ref = VariantAnnotation::ref(vcf),
    alt = VariantAnnotation::alt(vcf),
    depth = VariantAnnotation::geno(vcf)$DP[, c(1)],
    alt_depth = VariantAnnotation::geno(vcf)$VD[, c(1)]
    )  
  # Retain all INFO fields
  info <- as.data.frame(info(vcf))
  dat <- cbind(dat, info)
  row.names(dat) <- NULL 
  # Rename columns
  names(dat)[names(dat) == "TYPE"] <- "variation_type"
  names(dat)[names(dat) == "END"] <- "end"
 
  
  
  # Read in sample data if it's provided
  if (!is.null(sample_data_file)) {
    sampledata <- read.delim(file.path(sample_data_file),
                            sep = sd_sep,
                            header = T
                            
    )
    dat <- left_join(dat, sampledata, suffix = c("", ".sampledata"))
  }
  
# Clean up variation_type column to match .mut
  # REF --> no_variant
  # sv subtypes and IUPAC symbols --> symbolic
  # n(ref) = n(alt) > 1 --> mnv
  # n(ref) > n(alt) == 1 --> deletion
  # n(ref) == 1 < n(alt) --> insertion
  # n(ref) != n(alt) & neither == 1 --> "complex"
# Create an VARLEN column
  dat <- dat %>%
    dplyr::mutate(
       nchar_ref = nchar(ref),
      nchar_alt = nchar(alt.value),
      variation_type = tolower(dat$variation_type),
      #TO DO: fix sv - it doesn't change it to symbolic. 
      variation_type = 
        ifelse(.data$variation_type == "ref", "no_variant", 
          ifelse(.data$variation_type %in% c("inv", "dup", "del", "ins", "fus", "cnv",
                                             "cnv:tr", "dup:tandem", "del:me", "ins:me",
                                      #These ambiguity codes may need to be defined from the alt column instead
                                             "r", "k", "s", "y", "m", "w", "b", "h", "n", "d", "v"),"symbolic",
                 .data$variation_type))) %>%
    dplyr::mutate(                
           ifelse(.data$variation_type != "symbolic" & .data$nchar_ref == .data$nchar_alt & .data$nchar_ref > 1 , "mnv",
            ifelse(.data$variation_type != "symbolic" & .data$nchar_ref > .data$nchar_alt & .data$nchar_alt == 1, "deletion",
             ifelse(.data$variation_type != "symbolic" & .data$nchar_ref < .data$nchar_alt & .data$nchar_ref == 1, "insertion", 
              ifelse(.data$variation_type != "symbolic" & .data$nchar_ref != .data$nchar_alt & .data$nchar_alt > 1 &  .data$nchar_ref > 1, "complex",
                    .data$variation_type)))),
      VARLEN = 
        ifelse(.data$variation_type %in% c("insertion", "deletion", "complex"), .data$nchar_alt - .data$nchar_ref,
         ifelse(.data$variation_type %in% c("snv", "mnv"), .data$nchar_ref,
               "."))
    )
    
  #create ref_depth, short_ref, subtype
  dat <- dat %>%
    dplyr::mutate(
      ref_depth = .data$depth - .data$alt_depth,
      subtype = 
        ifelse(.data$variation_type == "snv",
               paste0(.data$ref, ">", .data$alt.value),
               "."),
      short_ref = substr(.data$ref, 1, 1))
  

# Create total_depth and no_calls columns based on set parameter depth_calc.
  # Requires AD field in FORMAT of vcf. If this field is missing, we use depth instead of total_depth
AD <- VariantAnnotation::geno(vcf)$AD  

if (length(AD) == 0) {
  # Handle the case where AD is missing
  dat <- dat %>%
    dplyr::mutate(
      no_calls = 0,  # Since AD is missing, no calls can't be calculated
      VAF = .data$alt_depth / .data$depth  # Calculate VAF using depth
    )
  cat("Warning: no_calls cannot be calculated because there is no AD field.\n")
  cat("VAF calculated with depth (DP; includes N-calls) because AD field is missing.\n")
  
} else {  
AD <- as.data.frame(do.call(rbind, AD))
colnames(AD) <- paste0("depth", 1:ncol(AD))

dat$ref_depth <- AD$depth1
dat$var_depth <- AD$depth2
dat <- dat %>%
  dplyr::mutate(
        var_depth = ifelse(is.na(.data$var_depth), 0, .data$var_depth))

# Filter the rows where "sample," "contig," and "start" are duplicated
duplicated_rows <- dat %>%
  dplyr::group_by(.data$sample, .data$contig, .data$start) %>%
  dplyr::filter(n() > 1)

# Calculate total_depth for the duplicated rows
if (depth_calc == "take_del") {
  total_depth_duplicated <- duplicated_rows %>%
    dplyr::group_by(.data$sample, .data$contig, .data$start) %>%
    dplyr::summarize(
      total_depth = case_when(
        any(.data$variation_type == "deletion") & length(unique(.data$ref_depth[!is.na(.data$ref_depth)])) > 1 ~
          sum(ifelse(.data$variation_type == "deletion", .data$ref_depth, 0), na.rm = TRUE) + sum(.data$var_depth, na.rm = TRUE),
        any(.data$variation_type == "complex" & !any(.data$variation_type == "deletion")) & length(unique(.data$ref_depth[!is.na(.data$ref_depth)])) > 1 ~
          sum(ifelse(.data$variation_type == "complex", .data$ref_depth, 0), na.rm = TRUE + sum(.data$var_depth, na.rm = TRUE)),
        TRUE ~ round(mean(.data$ref_depth, na.rm = TRUE) + sum(.data$var_depth, na.rm = TRUE))
      )
    ) %>%
    dplyr::ungroup()
  
} else if (depth_calc == "take_mean") {
  total_depth_duplicated <- duplicated_rows %>%
    dplyr::group_by(.data$sample, .data$contig, .data$start) %>%
    dplyr::summarize(
      total_depth = round(mean(.data$ref_depth, na.rm = TRUE) + sum(.data$var_depth, na.rm = TRUE))
    ) %>%
    dplyr::ungroup()

  } else {
  stop("Invalid depth_calc input. Please choose 'take_mean' or 'take_del'.")
}

# Merge the total_depth values for duplicated rows back into the original data frame
dat <- dat %>%
  dplyr::left_join(total_depth_duplicated, by = c("sample", "contig", "start"))

# Calculate total_depth for non-duplicated rows (ref_depth + var_depth)
dat <- dat %>%
  dplyr::mutate(
    total_depth = ifelse(!duplicated(dat[, c("sample", "contig", "start")]), .data$ref_depth + .data$var_depth, .data$total_depth),
    no_calls = .data$depth - .data$total_depth,
    VAF = .data$alt_depth / .data$total_depth
  ) %>%
  dplyr::select(-var_depth, ignore.case = FALSE)
}

##############################################################
# Create Context Column using target Sequences
# Turn dat into a GRanges object.  
mut_ranges <- makeGRangesFromDataFrame(
    df = as.data.frame(dat),
    keep.extra.columns = T,
    seqnames.field = "contig",
    start.field = "start",
    end.field = "end"
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
#####################################################################################################
# Retrieve reference sequences

  #get the species based on the genome  
  get_species_param <- function(assembly) {
    if (assembly %in% c("GRCm38", "GRCm39")) {
      return("mouse")
    } else if (assembly %in% c("GRCh37", "GRCh38")) {
      return("human")
    } else {
      stop("Invalid assembly parameter")
    }
  }
  
#retrieve the sequences (GRanges object)
region_ranges <- DupSeqR::get_seq(regions_df = genic_regions, species = get_species_param(assembly), genome_version = assembly)
  
# Join with mutation data
 ranges_joined <- plyranges::join_overlap_left(mut_ranges, region_ranges, suffix = c("_mut", "_regions"))
  return(ranges_joined)

# Get no_variant context
dat <- ranges_joined %>%
  plyranges::mutate(
                    start_string = start - ext_start +1,
                    context = substr(sequence, start_string - 1, start_string + 1)) %>%
  plyranges::select(-start_string)

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
    plyranges::mutate(
      context_with_mutation =
        ifelse(subtype != ".",
               paste0(
                 stringr::str_sub(context, 1, 1),
                 "[", subtype, "]",
                 stringr::str_sub(context, 3, 3)
               ),
               variation_type
        ),
       normalized_context = ifelse(
        stringr::str_sub(context, 2, 2) %in% c("G", "A", "g", "a"),
        mapply(function(x) reverseComplement.default(x, case = "upper"), context),
        context
      ),
      normalized_subtype = 
        ifelse(subtype %in% names(sub_dict),
               sub_dict[subtype],
               subtype
        ),
      normalized_ref = dplyr::case_when(
        substr(ref, 1, 1) == "A" ~ "T",
        substr(ref, 1, 1) == "G" ~ "C",
        substr(ref, 1, 1) == "C" ~ "C",
        substr(ref, 1, 1) == "T" ~ "T"
      ),
      normalized_context_with_mutation =
        ifelse(normalized_subtype != ".",
               paste0(
                 stringr::str_sub(normalized_context, 1, 1),
                 "[", normalized_subtype, "]",
                 stringr::str_sub(normalized_context, 3, 3)
               ),
               variation_type
        ),
      gc_content = (stringr::str_count(string = context, pattern = "G") +
                      stringr::str_count(string = context, pattern = "C"))
      / stringr::str_count(context)
    ) %>%
    plyranges::mutate(
      normalized_subtype = ifelse(
        normalized_subtype == ".",
        variation_type,
        normalized_subtype
      ),
      subtype = ifelse(
        subtype == ".",
        variation_type,
        subtype
      )
    )


}







