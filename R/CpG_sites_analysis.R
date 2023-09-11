#' Get mutations at CpG sites
#'
#' Subset the mutations provided and return only mutations that are found at CpG sites.
#' @param regions A GRanges object containing the genomic regions of interest in which to look for CpG sites. Must have the metadata column "sequence" populated with the raw nucleotide sequence to search for CpGs.
#' @param mut_data A GRanges object containing the mutation data to be interrogated.
#' @param variant_types TODO Copy in info
#' @param include_no_variants TRUE or FALSE to indicate whether the table should
#' include CpG sites with no variants. Useful if you want to know how many of 
#' the potential sites were mutated.
#' @param motif Default "CG", which returns CpG sites. You could in theory use
#' an arbitrary string to look at different motifs. Use with caution.
#' @returns A GRanges object where each range is a mutation at a CpG site (a subset of mutations from the larger object provided to the function).
#' @importFrom Biostrings matchPattern
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges ranges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom dplyr filter mutate
#' @importFrom plyranges find_overlaps
#' @export
get_CpG_mutations <- function(regions, mut_data,
                              variant_types = c("snv","indel","mnv","sv"),
                              include_no_variants = T,
                              motif = "CG") {
  # Step 3 - find all the CpG sites within those regions identified
  all_CpGs <- list()
  for (i in seq_along(regions)) {
    CpG_sites <- Biostrings::matchPattern(
      pattern = motif,
      subject = regions[i]$sequence[[1]]
    )
    CpG_sites <- GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(regions[i]),
      ranges = IRanges::IRanges(
        start = start(IRanges::ranges(CpG_sites)) + start(regions[i]) - 1,
        end = end(IRanges::ranges(CpG_sites)) + start(regions[i]) - 1
      )
    )
    all_CpGs[[i]] <- CpG_sites
  }
  CpGs_combined <- do.call("c", all_CpGs)
  # Step 4 - join mutation data with CpG sites
  CpGs_in_data <- plyranges::find_overlaps(mut_data, CpGs_combined) |>
    dplyr::filter(variation_type %in% variant_types)
  if (include_no_variants == T) {
    return(CpGs_in_data)
  } else {
    return(CpGs_in_data |> dplyr::filter(!variation_type == "no_variant"))
  }
}

#' Get the coordinates of CpG sites
#'
#' Imports package data to find target regions and some associated information, and further extends the table by getting raw nucleotide sequences for each region of the genome. Note that the way this is written, currently, the default genomes are hg38 and mm10 for human and mouse, respectively.
#' @param regions A GRanges object containing the genomic regions of interest in which to look for CpG sites. Must have the metadata column "sequence" populated with the raw nucleotide sequence to search for CpGs.
#' @param motif Default "CG", which returns CpG sites. You could in theory use
#' an arbitrary string to look at different motifs. Use with caution.
#' @returns A GRanges object where each range is a mutation at a CpG site (a subset of mutations from the larger object provided to the function).
#' @importFrom Biostrings matchPattern
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges IRanges ranges
#' @export
get_CpG_regions <- function(regions, motif = "CG") {
  # Similar to the above function but instead returns all the sites where CpGs are found in the reference (instead of the mutation data)
  all_CpGs <- list()
  for (i in seq_along(regions)) {
    CpG_sites <- Biostrings::matchPattern(
      pattern = motif,
      subject = regions[i]$sequence[[1]]
    )
    CpG_sites <- GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(regions[i]),
      ranges = IRanges::IRanges(
        start = start(IRanges::ranges(CpG_sites)) + start(regions[i]) - 1,
        end = end(IRanges::ranges(CpG_sites)) + start(regions[i]) - 1
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
#' @param motif Default "CG", which returns CpG sites. You could in theory use
#' an arbitrary string to look at different motifs. Use with caution. In this 
#' case the pattern being searched must be a column in the mutation data.
#' @param column_query Default "context" but can be any column  in the mutation
#' data that you wish to look for a motif in.
#' @param ... Additional arguments to vcountPattern()
#' @returns A data frame with the same number of rows as there were ranges in the input, but with an additional metadata column indicating CpG sites in the target sequence of the mutation.
#' @importFrom Biostrings vcountPattern
#' @export
annotate_CpG_sites <- function(mut_data,
                               motif = "CG",
                               column_query = "context",
                               ...) {
  annotated_data <- as.data.frame(mut_data) |>
    dplyr::mutate(CpG_site = Biostrings::vcountPattern(
      pattern = motif,
      subject = .data[[column_query]],
      ...)) |>
    dplyr::mutate(CpG_site = ifelse(CpG_site == 0, F, T))
  return(annotated_data)
}

#' Summarize CpG sites
#' 
#' Creates a summary table of CpG sites based on groupings of interest. This is
#' basically a convenience function that wraps `calculate_mut_freqs()` over CpG
#' data (or any data). See the documentation for that function for parameters. 
#' It is up to the user to supply proper data to the function.
#' @param cpg_muts A data frame containing CpG mutations
#' @param ... Additional arguments to calculate_mut_freqs()
#' @export
make_CpG_summary_table <- function(cpg_muts = cpg_mutations, ...) {
    # Verify that cpg_muts has a 
  # Do we even need a special function for this? Might get by with 
  # calculate_mut_freqs().
  
}
