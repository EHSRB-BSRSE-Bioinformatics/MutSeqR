#' Get mutations at CpG sites
#'
#' Subset the mutations provided and return only mutations that are found at CpG sites.
#' @param regions A GRanges object containing the genomic regions of interest in which to look for CpG sites. Must have the metadata column "sequence" populated with the raw nucleotide sequence to search for CpGs.
#' @param mut_data A GRanges object containing the mutation data to be interrogated.
#' @returns A GRanges object where each range is a mutation at a CpG site (a subset of mutations from the larger object provided to the function).
#' @importFrom Biostrings matchPattern
#' @importFrom GenomicRanges GRanges
#' @importFrom plyranges find_overlaps
#' @export
get_CpG_mutations <- function(regions, mut_data) {
  # Step 3 - find all the CpG sites within those regions identified
  all_CpGs <- list()
  for (i in seq_along(regions)) {
    CpG_sites <- Biostrings::matchPattern(
      pattern = "CG",
      subject = regions[i]$sequence[[1]]
    )
    CpG_sites <- GenomicRanges::GRanges(
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
#' @importFrom Biostrings matchPattern
#' @importFrom GenomicRanges GRanges
#' @export
get_CpG_regions <- function(regions) {
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
#' @importFrom Biostrings vcountPattern
#' @export
annotate_CpG_sites <- function(mut_data) {
  annotated_data <- as.data.frame(mut_data) %>%
    dplyr::mutate(CpG_site = Biostrings::vcountPattern(pattern = "CG", context)) %>%
    dplyr::mutate(CpG_site = ifelse(CpG_site == 0, F, T))
  return(annotated_data)
}
