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
#' @import GenomicRanges
#' @export
import_mut_data <- function(mut_file = "../../data/Jonatan_Mutations_in_blood_and_sperm_samples_221021_MM.txt",
                            rsids = F,
                            sample_data_file = NULL,
                            sd_sep = "\t",
                            mut_sep = "\t",
                            regions_file = "../../inst/extdata/genic_regions_hg38.txt") {

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
  # Make new column with COSMIC-style 96 base context

  # Define substitution dictionary to normalize to pyrimidine context
  sub_dict <- c("G>T" = "C>A", "G>A" = "C>T", "G>C" = "C>G",
                "A>G" = "T>C", "A>C" = "T>G", "A>T" = "T>A")
  
  # The column that represents depth might vary
  depth_col <- ifelse("total_depth" %in% colnames(dat),
                      "total_depth",
                      ifelse("depth" %in% colnames(dat),
                             "depth",  stop("Error: I'm not sure which column
                             specifies depth.")))
  
  dat <- dat %>%
    mutate(
      ref_depth = !!sym(depth_col) - alt_depth,
      context_with_mutation =
        ifelse(subtype != ".",
               paste0(str_sub(context, 1, 1),
                      "[", subtype, "]",
                      str_sub(context, 3, 3)),
               "."),
      normalized_context = ifelse(
        subtype %in% names(sub_dict),
        mapply(function(x) spgs::reverseComplement(x, case = "upper"), context),
        context),
      normalized_subtype = ifelse(
        subtype %in% names(sub_dict),
        sub_dict[subtype],
        subtype),
      short_ref = substr(ref, 1, 1),
      normalized_ref = case_when(
        substr(ref, 1, 1) == "A" ~ "T",
        substr(ref, 1, 1) == "G" ~ "C",
        substr(ref, 1, 1) == "C" ~ "C",
        substr(ref, 1, 1) == "T" ~ "T"
      )
     ) %>%
    mutate(normalized_context_with_mutation =
             ifelse(subtype != ".",
                    paste0(str_sub(normalized_context, 1, 1),
                           "[", normalized_subtype, "]",
                           str_sub(normalized_context, 3, 3)),
                    "."),
           gc_content = (str_count(string = context, pattern = "G") +
                           str_count(string = context, pattern = "C"))
           / str_count(context)) %>%
    { if ("depth" %in% names(.))
      mutate(., total_depth = depth - no_calls)
      else
      .
    }
  
  dat <- dat %>% mutate(VAF = alt_depth / total_depth)
  
  mut_ranges <- makeGRangesFromDataFrame(
    df = as.data.frame(dat),
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
