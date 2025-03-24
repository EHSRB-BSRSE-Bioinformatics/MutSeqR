#' Get mutations at CpG sites
#'
#' Subset the mutation data provided and return only mutations that are found 
#' at CpG sites.
#' @param regions A GRanges object containing the genomic regions of interest 
#' in which to look for CpG sites. Must have the metadata column "sequence" 
#' populated with the raw nucleotide sequence to search for CpGs. This object can be 
#' obtained using the get_seq.R function.
#' @param mut_data A dataframe or GRanges object containing the mutation data to 
#' be interrogated.
#' @param variant_types Include these variant types. A vector of one or more
#'  "snv", "complex", "deletion", "insertion", "mnv", "symbolic", "no_variant".
#'  Default includes all variants. 
#' @param include_no_variants TRUE or FALSE to indicate whether the table should
#' include CpG sites with no variants. Useful if you want to know how many of 
#' the potential sites were mutated.
#' @param motif Default "CG", which returns CpG sites. You could in theory use
#' an arbitrary string to look at different motifs. Use with caution.
#' @returns A GRanges object where each range is a mutation at a CpG site 
#' (a subset of mutations from the larger object provided to the function).
#' @importFrom Biostrings matchPattern
#' @importFrom GenomicRanges GRanges start end
#' @importFrom IRanges IRanges ranges
#' @importFrom dplyr filter mutate
#' @importFrom plyranges find_overlaps
#' @importFrom rlang .data
#' @export
get_CpG_mutations <- function(regions, mut_data,
                              variant_types = c("snv","insertion", "deletion", "mnv","symbolic"),
                              include_no_variants = TRUE,
                              motif = "CG") {

   if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    stop("Package GenomeInfoDb is required. Please install from Bioconductor.")
  }
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
        start = GenomicRanges::start(IRanges::ranges(CpG_sites)) + GenomicRanges::start(regions[i]) - 1,
        end = GenomicRanges::end(IRanges::ranges(CpG_sites)) + GenomicRanges::start(regions[i]) - 1
      )
    )
    all_CpGs[[i]] <- CpG_sites
  }
  CpGs_combined <- do.call("c", all_CpGs)
  # Step 4 - join mutation data with CpG sites

  if (inherits(mut_data, "data.frame")) { 
    mut_data <- GenomicRanges::makeGRangesFromDataFrame(
      df = mut_data,
      keep.extra.columns = T,
      seqnames.field = "contig",
      start.field = "start",
      end.field = "end"
    )}
  
  CpGs_in_data <- plyranges::find_overlaps(mut_data, CpGs_combined) %>%
    as.data.frame %>%
    dplyr::filter(.data$variation_type %in% variant_types)
  if (include_no_variants == T) {
    return(CpGs_in_data)
  } else {
    return(CpGs_in_data %>% dplyr::filter(!.data$variation_type == "no_variant"))
  }
}

#' Get the coordinates of the CpG sites within your genomic regions
#'
#' Filters the ranges of your genomic regions to find a positions with a 
#' specific motif. The default is CpG sites, but can be customizable. 
#' @param regions A GRanges object containing the genomic regions of interest in 
#' which to look for CpG sites. Must have the metadata column "sequence" populated 
#' with the raw nucleotide sequence to search for CpGs. This object can be 
#' obtained using the get_seq.R function.
#' @param motif Default "CG", which returns CpG sites. You could in theory use
#' an arbitrary string to look at different motifs. Use with caution.
#' @returns A GRanges object where each range is a CpG site (a subset of ranges 
#' from the larger object provided to the function).
#' @importFrom Biostrings matchPattern
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges ranges
#' @export
get_CpG_regions <- function(regions, motif = "CG") {
  
  if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    stop("Package GenomeInfoDb is required. Please install from Bioconductor.")
  }
  all_CpGs <- list()
  for (i in seq_along(regions)) {
    CpG_sites <- Biostrings::matchPattern(
      pattern = motif,
      subject = regions[i]$sequence[[1]]
    )
    CpG_sites <- GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(regions[i]),
      ranges = IRanges::IRanges(
        start = GenomicRanges::start(IRanges::ranges(CpG_sites)) + GenomicRanges::start(regions[i]) - 1,
        end = GenomicRanges::end(IRanges::ranges(CpG_sites)) + GenomicRanges::start(regions[i]) - 1
      )
    )
    all_CpGs[[i]] <- CpG_sites
  }
  CpGs_combined <- do.call("c", all_CpGs)
  return(CpGs_combined)
}

#' Annotate CpG sites
#'
#' A simple method to test whether your trinucleotide context contains a CpG site. 
#' Vectorized version of Biostrings::vcountPattern is used.
#' @param mut_data A dataframe or GRanges object containing the genomic regions 
#' of interest in which to look for CpG sites. 
#' @param motif Default "CG", which returns CpG sites. You could in theory use
#' an arbitrary string to look at different motifs. Use with caution. In this 
#' case the pattern being searched must be a column in the mutation data.
#' @param column_query Default "context" but can be any column  in the mutation
#' data that you wish to look for a motif in.
#' @param ... Additional arguments to vcountPattern()
#' @returns A data frame with the same number of rows as there were ranges in the 
#' input, but with an additional metadata column indicating CpG sites in the 
#' target sequence of the mutation.
#' @importFrom Biostrings vcountPattern
#' @importFrom rlang .data
#' @export
annotate_CpG_sites <- function(mut_data,
                               motif = "CG",
                               column_query = "context",
                               ...) {
  annotated_data <- as.data.frame(mut_data) %>%
    dplyr::mutate(CpG_site = Biostrings::vcountPattern(
      pattern = motif,
      subject = .data[[column_query]],
      ...)) %>%
    dplyr::mutate(CpG_site = ifelse(.data$CpG_site == 0, F, T))
  return(annotated_data)
}

#' Summarize CpG sites
#' 
#' Creates a summary table of CpG sites based on groupings of interest. This is
#' basically a convenience function that wraps `calculate_mfs()` over CpG
#' data (or any data). See the documentation for that function for parameters. 
#' It is up to the user to supply proper data to the function.
#' @param cpg_muts A data frame containing CpG mutations 
#' TO DO: cpg_muts = df "cpg_mutations" is created in the .Rmd file, but is not created by any other function. 
#' @param ... Additional arguments to calculate_mfs()
#' @export
make_CpG_summary_table <- function(cpg_muts, ...) {
    # Verify that cpg_muts has a 
  # Do we even need a special function for this? Might get by with 
  # calculate_mfs().
  
}
