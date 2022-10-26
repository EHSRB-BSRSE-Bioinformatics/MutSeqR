#!/usr/bin/R
library(tidyverse)
library(fuzzyjoin)
library(GenomicRanges)
library(plyranges)

#' Import a .mut file
#'
#' Imports a .mut file into the local R environment.
#' @param mut_file The .mut file containing mutation data to be imported. Columns required are... (fill in one day)
#' @param rsids TRUE or FALSE; whether or not the .mut file contains rsID information (existing SNPs)
#' @param sample_data_file An optional file containing additional sample metadata (dose, timepoint, etc.)
#' @param sd_sep The delimiter for importing sample metadata tables
#' @param mut_sep The delimiter for importing the .mut file
#' @param grouping_variable Analagous to the DESeq2 nomenclature for "interesting groups", this should be a column in your data for which you would like to build comparisons upon, e.g., dose, or tissue, or sex
#' @returns A table where each row is a mutation, and columns indicate the location, type, and other data.
#' @export
import_mut_data <- function(mut_file = "../../data/Jonatan_Mutations_in_blood_and_sperm_samples_221021_MM.txt",
                           rsids = F,
                           sample_data_file = NULL,
                           sd_sep = "\t",
                           mut_sep = "\t",
                           regions_file = "../../inst/genic_regions_hg38.txt",
                           grouping_variable = "dose") {

  # Read in mut file
  # Note: col names
  # mut_depth = final_somatic_alt_depth
  # total_depth_ = informative_total_depth
  dat <- read.table(mut_file, header = T, sep = "\t", fileEncoding = "UTF-8-BOM")
  if (rsids == T) {
    # If we have rs IDs, add a column indicating whether the mutation is a known SNP
    dat <- dat %>% mutate(is_known = ifelse(!id==".", "Y","N"))
  }
  
  # Read in sample data if it's provided
  if (!is.null(sample_data_file)) {
    sampledata <- read.delim(file.path(sample_data_file),
                             sep = sd_sep,
                             header = T)
    dat <- dplyr::left_join(dat, sampledata, suffix = c("",".sampledata"))
  }
  
  ################
  # Clean up data:
  # Select only SNVs
  # Remove sites where mut_depth (final_somatic_alt_depth) is zero
  # Get reverse complement of sequence context where mutation is listed on purine context
  # Change all purine substitutions to pyrimidine substitutions
  # Make new column with COSMIC-style 96 base context
  # Calculate depth for each sequence context and dose group
  # Calculate frequency for each mouse within each 96 trinucleotide mutation
  
  dat <- dat %>%
    dplyr::mutate(normalized_context = ifelse(
      test = subtype %in% c("G>T","G>A","G>C","A>T","A>C","A>G"),
      yes = mapply(function(x) spgs::reverseComplement(x, case="upper"), context),
      no = context)) %>%
    dplyr::mutate(normalized_subtype = subtype) %>%
    dplyr::mutate(normalized_subtype = str_replace(normalized_subtype, "G>T", "C>A")) %>%
    dplyr::mutate(normalized_subtype = str_replace(normalized_subtype, "G>T", "C>A")) %>%
    dplyr::mutate(normalized_subtype = str_replace(normalized_subtype, "G>A", "C>T")) %>%
    dplyr::mutate(normalized_subtype = str_replace(normalized_subtype, "G>C", "C>G")) %>%
    dplyr::mutate(normalized_subtype = str_replace(normalized_subtype, "A>T", "T>A")) %>%
    dplyr::mutate(normalized_subtype = str_replace(normalized_subtype, "A>C", "T>G")) %>%
    dplyr::mutate(normalized_subtype = str_replace(normalized_subtype, "A>G", "T>C")) %>%
    dplyr::mutate(context_with_mutation = paste0(str_sub(normalized_context, 1, 1),
                                          "[",normalized_subtype,"]",
                                          str_sub(normalized_context, 3, 3)) ) %>%
    dplyr::group_by(normalized_context, !!sym(grouping_variable)) %>%
    dplyr::mutate(group_depth = sum(total_depth)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(!!sym(grouping_variable)) %>%
    dplyr::mutate(group_mut_count = sum(mut_depth)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(context_with_mutation, !!sym(grouping_variable)) %>%
    dplyr::mutate(group_mut_count_by_type = sum(mut_depth)) %>%
    dplyr::mutate(group_frequency = group_mut_count_by_type/group_mut_count/group_depth) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(normalized_context, sample) %>%
    dplyr::mutate(sample_depth = sum(total_depth)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sample_frequency = (sum(mut_depth)/sample_depth) ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(variation_type == "snv") %>%
    dplyr::filter(!mut_depth == 0) %>%
    dplyr::mutate(gc_content = (str_count(string = context, pattern = "G") +
                                  str_count(string = context, pattern = "C"))
                  /str_count(context))
  
  mut_ranges <- makeGRangesFromDataFrame(df = dat,
                           keep.extra.columns = T,
                           seqnames.field = "contig",
                           start.field = "start",
                           end.field = "end",
                           starts.in.df.are.0based = TRUE)

  # Annotate the mut file with additional information about genomic regions in the file
  genic_regions <- read.delim(regions_file)
  
  region_ranges <- makeGRangesFromDataFrame(df = genic_regions,
                                         keep.extra.columns = T,
                                         seqnames.field = "contig",
                                         start.field = "start",
                                         end.field = "end",
                                         starts.in.df.are.0based = TRUE)
  
  ranges_joined <- plyranges::join_overlap_left(mut_ranges, region_ranges)
  return(ranges_joined)
}
