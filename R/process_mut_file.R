#!/usr/bin/R
library(tidyverse)
library(fuzzyjoin)

import_ds_data <- function(mut_file = "./data/Jonatan_Mutations_in_blood_and_sperm_samples_221021_MM.txt",
                           rsids = F,
                           sample_data_file = NULL,
                           sd_sep = "\t",
                           mut_sep = "\t",
                           regions_file = "./inst/genic_regions_hg38.txt") {
  
  
  
  # Read in mut file
  dat <- read.table(mut_file, header = T, sep = "\t", fileEncoding = "UTF-8-BOM")
  if (rsids == T) {
    # If we have rs IDs, add a column indicating whether the mutation is a known SNP
    dat <- dat %>% mutate(is_known = ifelse(!id==".", "Y","N"))
  }
  
  # Annotate the mut file with additional information about genomic regions in the file
  genic_regions <- read.delim(regions_file)
  #dat <- fuzzyjoin::genome_join(dat, genic_regions, by=c("contig", "start", "end"), mode="inner")
  
  
  # Read in sample data if it's provided
  if (!is.null(sample_data_file)) {
    sampledata <- read.delim(file.path(sample_data_file),
                             sep = sd_sep,
                             header = T)
    dat <- dplyr::left_join(dat, sampledata)
  }
  
  return(dat)
}
