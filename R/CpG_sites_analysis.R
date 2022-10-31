#' Get sequence of Duplex Sequencing target regions
#'
#' Imports package data to find target regions and some associated information, and further extends the table by getting raw nucleotide sequences for each region of the genome. Note that the way this is written, currently, the default genomes are hg38 and mm10 for human and mouse, respectively.
#' @param species One of "mouse" or "human", to determine which regions to return.
#' @returns A GRanges object where each range is a target region
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
    regions_file <- "../../inst/genic_regions_hg38.txt"
  } else if (species == "mouse") {
    db <- "BSgenome.Mmusculus.UCSC.mm10"
    if (!require(BSgenome.Mmusculus.UCSC.mm10)) {
      stop("BSgenome.Mmusculus.UCSC.mm10 not installed")
    }
    # library(org.Mm.eg.db) # May be useful...
    regions_file <- "../../inst/genic_regions_mm10.txt"
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

#' Get mutations at CpG sites
#'
#' Subset the mutations provided and return only mutations that are found at CpG sites.
#' @param regions A GRanges object containing the genomic regions of interest in which to look for CpG sites. Must have the metadata column "sequence" populated with the raw nucleotide sequence to search for CpGs.
#' @param mut_data A GRanges object containing the mutation data to be interrogated.
#' @returns A GRanges object where each range is a mutation at a CpG site (a subset of mutations from the larger object provided to the function).
#' @export
get_CpG_mutations <- function(regions, mut_data) {
  if (!require(tidyverse)) {
    stop("tidyverse not installed")
  }
  if (!require(GenomicRanges)) {
    stop("GenomicRanges not installed")
  }
  if (!require(Biostrings)) {
    stop("Biostrings not installed")
  }
  # Step 3 - find all the CpG sites within those regions identified
  all_CpGs <- list()
  for (i in seq_along(regions)) {
    CpG_sites <- Biostrings::matchPattern(
      pattern = "CG",
      subject = regions[i]$sequence[[1]]
    )
    CpG_sites <- GRanges(
      seqnames = seqnames(regions[i]),
      ranges = IRanges(
        start = start(ranges(CpG_sites)) + start(regions[i]) - 1,
        end = end(ranges(CpG_sites)) + start(regions[i]) - 1
      )
    )
    all_CpGs[[i]] <- CpG_sites
  }
  CpGs_combined <- do.call("c", all_CpGs)
  # Step 4 - join mutation data with CpG sites
  CpGs_in_data <- plyranges::find_overlaps(mut_data, CpGs_combined)
  return(CpGs_in_data)
}

#' Get the coordinates of CpG sites
#'
#' Imports package data to find target regions and some associated information, and further extends the table by getting raw nucleotide sequences for each region of the genome. Note that the way this is written, currently, the default genomes are hg38 and mm10 for human and mouse, respectively.
#' @param regions A GRanges object containing the genomic regions of interest in which to look for CpG sites. Must have the metadata column "sequence" populated with the raw nucleotide sequence to search for CpGs.
#' @returns A GRanges object where each range is a mutation at a CpG site (a subset of mutations from the larger object provided to the function).
#' @export
get_CpG_regions <- function(regions) {
  if (!require(tidyverse)) {
    stop("tidyverse not installed")
  }
  if (!require(GenomicRanges)) {
    stop("GenomicRanges not installed")
  }
  if (!require(Biostrings)) {
    stop("Biostrings not installed")
  }
  # Similar to the above function but instead returns all the sites where CpGs are found in the reference (instead of the mutation data)
  all_CpGs <- list()
  for (i in seq_along(regions)) {
    CpG_sites <- Biostrings::matchPattern(
      pattern = "CG",
      subject = regions[i]$sequence[[1]]
    )
    CpG_sites <- GRanges(
      seqnames = seqnames(regions[i]),
      ranges = IRanges(
        start = start(ranges(CpG_sites)) + start(regions[i]) - 1,
        end = end(ranges(CpG_sites)) + start(regions[i]) - 1
      )
    )
    all_CpGs[[i]] <- CpG_sites
  }
  CpGs_combined <- do.call("c", all_CpGs)
  return(CpGs_combined)
}

#' Annotate CpG sites
#'
#' A simple method to test whether your trinucleotide context contains a CpG site. Vectorized version of Biostrings::vcountPattern is used.
#' @param mut_data A GRanges object containing the genomic regions of interest in which to look for CpG sites. Must have the metadata column "sequence" populated with the raw nucleotide sequence to search for CpGs.
#' @returns A data frame with the same number of rows as there were ranges in the input, but with an additional metadata column indicating CpG sites in the target sequence of the mutation.
#' @export
annotate_CpG_sites <- function(mut_data) {
  if (!require(tidyverse)) {
    stop("tidyverse not installed")
  }
  if (!require(Biostrings)) {
    stop("Biostrings not installed")
  }
  annotated_data <- as.data.frame(mut_data) %>%
    dplyr::mutate(CpG_site = Biostrings::vcountPattern(pattern = "CG", context)) %>%
    dplyr::mutate(CpG_site = ifelse(CpG_site == 0, F, T))
  return(annotated_data)
}
