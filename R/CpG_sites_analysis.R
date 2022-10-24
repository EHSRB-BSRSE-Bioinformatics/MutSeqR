library(tidyverse)
library(GenomicRanges)

# Step 1 - Import the regions of interest
get_regions <- function(species = "human") {
  # Load species database
  if (species == "human") {
    db <- "BSgenome.Hsapiens.UCSC.hg38"
    # library(org.Hs.eg.db) # May be useful...
    regions_file = "./inst/genic_regions_hg38.txt"
  } else if (species == "mouse") {
    library(org.Mm.eg.db)
    regions_file = "./inst/genic_regions_mm10.txt"
  }
  
  regions <- read.delim(regions_file)
  regions_ranges <- makeGRangesFromDataFrame(regions,
                                            keep.extra.columns = T,
                                            seqnames.field = "contig",
                                            start.field = "start",
                                            end.field = "end")
  
  seqs <- getSeq(get(db), names = regions_ranges )
  return(regions)
}

# Step 2 - Get FASTA reference sequences for regions of interest

# Step 3 - enumerate all the CpG sites within those regions identified