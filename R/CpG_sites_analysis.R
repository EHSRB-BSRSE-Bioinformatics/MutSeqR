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
annotate_CpG_sites <- function(regions, mut_data) {
  all_CpGs <- list()
  for (i in seq_along(regions_ranges)) {
    CpG_sites <- Biostrings::matchPattern(pattern = "CG",
                                          subject = regions_ranges[i]$sequence[[1]])
    CpG_sites <- GRanges(seqnames = seqnames(regions_ranges[i]),
                         ranges = IRanges(start = start(ranges(CpG_sites)) + start(regions_ranges[i])-1,
                                          end = end(ranges(CpG_sites)) + start(regions_ranges[i])-1))
    all_CpGs[[i]] <- CpG_sites
  }
  CpGs_combined <- do.call("c", all_CpGs)
  # Step 4 - join mutation data with CpG sites
  CpGs_in_data <- plyranges::find_overlaps(mut_data, CpGs_combined)
  return(CpGs_in_data)
}
