#' Get sequence of Duplex Sequencing target regions
#'
#' This is mostly a helper function. It imports package data to find target
#' regions and some associated information, and further extends the table by
#' getting raw nucleotide sequences for each region of the genome. Note that the
#'  way this is written, currently, the default genomes are hg38 and mm10 for
#'  human and mouse, respectively.
#' @param species One of "mouse" or "human", to determine which regions to return.
#' @returns A GRanges object where each range is a target region.
#' @export
get_region_seqs <- function(species) {
  if (!require(tidyverse)) {
    stop("tidyverse not installed")
  }
  if (!require(GenomicRanges)) {
    stop("GenomicRanges not installed")
  }
  if (!require(Biostrings)) {
    stop("Biostrings not installed")
  }
  # Step 1 - Import the regions of interest
  # Load species database
  if (species == "human") {
    db <- "BSgenome.Hsapiens.UCSC.hg38"
    if (!require(BSgenome.Hsapiens.UCSC.hg38)) {
      stop("BSgenome.Hsapiens.UCSC.hg38 not installed")
    }
    # library(org.Hs.eg.db) # May be useful...
    regions_file <- system.file("extdata",
                                "genic_regions_hg38.txt",
                                package="DupSeqR")
  } else if (species == "mouse") {
    db <- "BSgenome.Mmusculus.UCSC.mm10"
    if (!require(BSgenome.Mmusculus.UCSC.mm10)) {
      stop("BSgenome.Mmusculus.UCSC.mm10 not installed")
    }
    # library(org.Mm.eg.db) # May be useful...
    regions_file <- system.file("extdata",
                                "genic_regions_mm10.txt",
                                package="DupSeqR")
  }
  
  regions <- read.delim(regions_file)
  regions_ranges <- makeGRangesFromDataFrame(
    df = regions,
    keep.extra.columns = T,
    seqnames.field = "contig",
    start.field = "start",
    end.field = "end",
    starts.in.df.are.0based = TRUE
  )
  
  # Step 2 - Get reference sequences for regions of interest
  seqs <- Biostrings::getSeq(get(db), names = regions_ranges)
  mcols(regions_ranges)$sequence <- seqs
  return(regions_ranges)
}
