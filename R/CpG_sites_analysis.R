library(tidyverse)
library(GenomicRanges)
library(Biostrings)

get_region_seqs <- function(species) {
  # Step 1 - Import the regions of interest
  # Load species database
  if (species == "human") {
    db <- "BSgenome.Hsapiens.UCSC.hg38"
    if (!require(BSgenome.Hsapiens.UCSC.hg38)) {
      stop("BSgenome.Hsapiens.UCSC.hg38 not installed")
    }
    # library(org.Hs.eg.db) # May be useful...
    regions_file = "../../inst/genic_regions_hg38.txt"
  } else if (species == "mouse") {
    db <- "BSgenome.Mmusculus.UCSC.mm10"
    if (!require(BSgenome.Mmusculus.UCSC.mm10)) {
      stop("BSgenome.Mmusculus.UCSC.mm10 not installed")
    }
    #library(org.Mm.eg.db) # May be useful...
    regions_file = "../../inst/genic_regions_mm10.txt"
  }
  
  regions <- read.delim(regions_file)
  regions_ranges <- makeGRangesFromDataFrame(df = regions,
                                            keep.extra.columns = T,
                                            seqnames.field = "contig",
                                            start.field = "start",
                                            end.field = "end",
                                            starts.in.df.are.0based = TRUE)
  
  # Step 2 - Get reference sequences for regions of interest
  seqs <- Biostrings::getSeq(get(db), names = regions_ranges)
  mcols(regions_ranges)$sequence <- seqs
  return(regions_ranges)
}

get_CpG_mutations <- function(regions, mut_data) {
  # Step 3 - find all the CpG sites within those regions identified
  all_CpGs <- list()
  for (i in seq_along(regions)) {
    CpG_sites <- Biostrings::matchPattern(pattern = "CG",
                                          subject = regions[i]$sequence[[1]])
    CpG_sites <- GRanges(seqnames = seqnames(regions[i]),
                         ranges = IRanges(start = start(ranges(CpG_sites)) + start(regions[i])-1,
                                          end = end(ranges(CpG_sites)) + start(regions[i])-1))
    all_CpGs[[i]] <- CpG_sites
  }
  CpGs_combined <- do.call("c", all_CpGs)
  # Step 4 - join mutation data with CpG sites
  CpGs_in_data <- plyranges::find_overlaps(mut_data, CpGs_combined)
  return(CpGs_in_data)
}
