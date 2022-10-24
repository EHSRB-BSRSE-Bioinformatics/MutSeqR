library(tidyverse)
library(GenomicRanges)


get_regions <- function(species) {
  # Step 1 - Import the regions of interest
  # Load species database
  if (species == "human") {
    db <- "BSgenome.Hsapiens.UCSC.hg38"
    # library(org.Hs.eg.db) # May be useful...
    regions_file = "./inst/genic_regions_hg38.txt"
  } else if (species == "mouse") {
    db <- "BSgenome.Mmusculus.UCSC.mm10"
    #library(org.Mm.eg.db)
    regions_file = "./inst/genic_regions_mm10.txt"
  }
  
  regions <- read.delim(regions_file)
  regions_ranges <- makeGRangesFromDataFrame(df = regions,
                                            keep.extra.columns = T,
                                            seqnames.field = "contig",
                                            start.field = "start",
                                            end.field = "end",
                                            starts.in.df.are.0based = TRUE)
  
  # Step 2 - Get reference sequences for regions of interest
  seqs <- getSeq(get(db), names = regions_ranges)
  mcols(regions_ranges)$sequence <- seqs
  return(regions_ranges)
}


# Step 3 - enumerate all the CpG sites within those regions identified

Biostrings::matchPattern(pattern = "CG", subject = regions_ranges[i]$sequence[[1]])
