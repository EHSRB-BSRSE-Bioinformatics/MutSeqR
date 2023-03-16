#' Import a .mut file
#'
#' Imports a .mut file into the local R environment.
#' @param mut_file The .mut file containing mutation data to be imported. Columns required are... (fill in one day)
#' @param rsids TRUE or FALSE; whether or not the .mut file contains rsID information (existing SNPs)
#' @param sample_data_file An optional file containing additional sample metadata (dose, timepoint, etc.)
#' @param sd_sep The delimiter for importing sample metadata tables
#' @param mut_sep The delimiter for importing the .mut file
#' @returns A table where each row is a mutation, and columns indicate the location, type, and other data.
#' @import tidyverse
#' @import plyranges
#' @export
import_mut_data <- function(mut_file = "../../data/Jonatan_Mutations_in_blood_and_sperm_samples_221021_MM.txt",
                            rsids = F,
                            sample_data_file = NULL,
                            sd_sep = "\t",
                            mut_sep = "\t",
                            regions_file = "../../inst/extdata/genic_regions_hg38.txt") {
  if (!require(tidyverse)) {
    stop("tidyverse not installed")
  }
  if (!require(GenomicRanges)) {
    stop("GenomicRanges not installed")
  }
  if (!require(plyranges)) {
    stop("plyranges not installed")
  }
  # Read in mut file
  # Note: col names
  # mut_depth = final_somatic_alt_depth
  # total_depth_ = informative_total_depth
  dat <- read.table(mut_file, header = T, sep = "\t", fileEncoding = "UTF-8-BOM")
  dat <- read.table(mut_file, header = T, sep = mut_sep, fileEncoding = "UTF-8-BOM")
  if (ncol(dat)<=1) { stop("Your imported data only has one column.
                           You may want to set mut_sep to properly reflect
                           the delimiter used for the data you are importing.")}
  if (rsids == T) {
    if(!"id" %in% colnames(dat)) {
      stop("Error: you have set rsids to TRUE,
      but there is no id column in the mut file!")
    }
    # If we have rs IDs, add a column indicating whether the mutation is a known SNP
    dat <- dat %>% mutate(is_known = ifelse(!id == ".", "Y", "N"))
  }

  # Read in sample data if it's provided
  if (!is.null(sample_data_file)) {
    sampledata <- read.delim(file.path(sample_data_file), sep = sd_sep,
                             header = T)
    dat <- left_join(dat, sampledata, suffix = c("", ".sampledata"))
  }

  ################
  # Clean up data:
  # Get reverse complement of sequence context where mutation is listed on purine context
  # Change all purine substitutions to pyrimidine substitutions
  # Make new column with COSMIC-style 96 base contex

  # Define substitution dictionary to normalize to pyrimidine context
  sub_dict <- c("G>T" = "C>A", "G>A" = "C>T", "G>C" = "C>G",
                "A>G" = "T>C", "A>C" = "T>G", "A>T" = "T>A")
  
  dat <- dat %>%
   mutate(
      ref_depth = total_depth - alt_depth,
      context_with_mutation = paste0(
        str_sub(context, 1, 1),
        "[", subtype, "]",
        str_sub(context, 3, 3)),
      normalized_context = ifelse(
        test = subtype %in% names(sub_dict),
        yes = mapply(function(x) spgs::reverseComplement(x, case = "upper"), context),
        no = context),
      normalized_subtype = ifelse(
        test = subtype %in% names(sub_dict),
        yes = sub_dict[subtype],
        no = subtype)) %>%
    mutate(normalized_context_with_mutation = paste0(
      str_sub(normalized_context, 1, 1),
      "[", normalized_subtype, "]",
      str_sub(normalized_context, 3, 3)
    )) %>%
    dplyr::mutate(gc_content = (str_count(string = context, pattern = "G") +
                                  str_count(string = context, pattern = "C"))
                  / str_count(context))

  mut_ranges <- makeGRangesFromDataFrame(
    df = dat,
    keep.extra.columns = T,
    seqnames.field = "contig",
    start.field = "start",
    end.field = "end",
    starts.in.df.are.0based = TRUE
  )
  
  # Annotate the mut file with additional information about genomic regions in the file
  genic_regions <- read.delim(regions_file)

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


