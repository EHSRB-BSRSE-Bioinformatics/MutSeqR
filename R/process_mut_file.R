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
#' @returns A table where each row is a mutation, and columns indicate the location, type, and other data.
#' @export
import_ds_data <- function(mut_file = "../../data/Jonatan_Mutations_in_blood_and_sperm_samples_221021_MM.txt",
                           rsids = F,
                           sample_data_file = NULL,
                           sd_sep = "\t",
                           mut_sep = "\t",
                           regions_file = "../../inst/genic_regions_hg38.txt") {

  # Read in mut file
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
  # ranges_joined <- findOverlaps(mut_ranges, region_ranges)
  
  # dat <- fuzzyjoin::genome_join(dat, genic_regions, by=c("contig", "start", "end"), mode="inner") %>%
  #   dplyr::rename("start" = "start.x",
  #                 "end" = "end.x",
  #                 "contig" = "contig.x")
  
  return(ranges_joined)
}
