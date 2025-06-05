#' Get mutations at CpG sites.
#'
#' *Needs to be reworked for variants >1bp*. Subset the mutation data and return only mutations that are found
#' at positions with a specific motif. The default is CpG sites, but can be
#' customizable.
#' @param mutation_data A dataframe or GRanges object containing the mutation
#' data to be interrogated. If supplying a data frame, the genomic coordinates
#' must be 1-based (true for mutation data imported using import_mut_data or
#' import_vcf_data).
#' @param regions A GRanges object containing the genomic regions of interest
#' in which to look for CpG sites. Must have the metadata column "sequence"
#' populated with the raw nucleotide sequence to search for CpGs. This object
#' can be obtained using the get_seq.R function.
#' @param variant_types Use this parameter to choose which variation_types
#' to include in the output. Provide a character vector of the variation _types
#' that you want to include. Options are "ambiguous", "complex", "deletion",
#' "insertion", "mnv", "no_variant", "snv", "sv", "uncategorized".
#' Alternatively, provide a character vector of the variation_types that you
#' want to exclude preceded by "-". All variation_types except those excluded
#' will be returned. Ex. inclusion: variant_types = "snv", will return only
#' rows with variation_type == "snv". Ex. exclusion:
#' variant_types = "-no_variant" will return all rows, except those with
#' variation_type == "no_variant" (default).
#' @param motif Default "CG", which returns CpG sites. You could in theory use
#' an arbitrary string to look at different motifs. Use with caution.
#' @param filter_mut A logical value indicating whether the function should
#' exclude rows flagged in the filter_mut column from the output. Default
#' is TRUE.
#' @returns A GRanges object where each range is a mutation at a CpG site
#' (a subset of mutations from the larger object provided to the function).
#' @importFrom Biostrings matchPattern
#' @importFrom GenomicRanges GRanges start end
#' @importFrom IRanges IRanges ranges
#' @importFrom dplyr filter mutate
#' @importFrom plyranges find_overlaps
#' @importFrom rlang .data
#' @importFrom GenomeInfoDb seqnames
#' @export
get_cpg_mutations <- function(mutation_data,
                              regions,
                              variant_types = c("-no_variant"),
                              motif = "CG",
                              filter_mut = TRUE) {
  # Set the list of variation_types to return
  all_variant_types <- MutSeqR::subtype_list$type
  filter_variants <- function(selected_types, all_variant_types) {
    if (is.character(selected_types) && length(selected_types) == 1) {
      selected_types <- unlist(strsplit(selected_types, ","))
    }
    exclusions <- grepl("^-", selected_types)
    if (any(exclusions)) {
      excluded_variants <- sub("^-", "", selected_types[exclusions])
      selected_variants <- setdiff(all_variant_types, excluded_variants)
    } else {
      selected_variants <- intersect(all_variant_types, selected_types)
    }
    return(selected_variants)
  }
  variant_types <- filter_variants(variant_types, all_variant_types)

  # Find all the CpG sites in specified regions
  CpGs_rg_sites <- get_cpg_regions(regions = regions, motif = motif)

  # Join mutation data with CpG sites
  if (inherits(mutation_data, "data.frame")) {
    mutation_data <- GenomicRanges::makeGRangesFromDataFrame(
      df = mutation_data,
      keep.extra.columns = TRUE,
      seqnames.field = "contig",
      start.field = "start",
      end.field = "end"
    )
  }
  # variants > 1bp that overlap with multiple cpg sites are added multiple times.
  CpGs_in_data <- plyranges::find_overlaps(mutation_data, CpGs_rg_sites) %>%
    as.data.frame %>%
    unique() %>%
    dplyr::filter(.data$variation_type %in% variant_types)

## Fix Overlaps for > 1bp: only consider start position
## need to also consider the end position (for the other strand)
  if (filter_mut) {
    CpGs_in_data <- dplyr::filter(CpGs_in_data, filter_mut == FALSE)
  }
  return(CpGs_in_data)
}

#' Get the coordinates of the CpG sites within your genomic regions
#'
#' Filters the ranges of your genomic regions to find all positions with a
#' specific motif. The default is CpG sites, but can be customizable.
#' @param regions A GRanges object containing the genomic regions of interest
#' in  which to look for CpG sites. Must have the metadata column "sequence"
#' populated  with the raw nucleotide sequence to search for CpGs. This object
#' can be obtained using the get_seq() function.
#' @param motif Default "CG", which returns CpG sites. You could in theory use
#' an arbitrary string to look at different motifs. Use with caution.
#' @returns A GRanges object where each range is a CpG site (a subset of ranges
#' from the larger object provided to the function).
#' @importFrom Biostrings matchPattern
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges ranges
#' @importFrom GenomeInfoDb seqnames
#' @export
get_cpg_regions <- function(regions, motif = "CG") {

  all_CpGs_rg <- list()
  for (i in seq_along(regions)) {
    CpG_sites_rg <- Biostrings::matchPattern(
      pattern = motif,
      subject = regions[i]$sequence[[1]]
    )
    CpG_sites_rg <- GenomicRanges::GRanges(
      seqnames = GenomeInfoDb::seqnames(regions[i]),
      ranges = IRanges::IRanges(
        start = GenomicRanges::start(IRanges::ranges(CpG_sites_rg)) + GenomicRanges::start(regions[i]) - 1,
        end = GenomicRanges::end(IRanges::ranges(CpG_sites_rg)) + GenomicRanges::start(regions[i]) - 1
      )
    )
    all_CpGs_rg[[i]] <- CpG_sites_rg
  }
  CpGs_combined <- do.call("c", all_CpGs_rg)
  return(CpGs_combined)
}

#' Annotate CpG sites
#'
#' A simple method to test whether your trinucleotide context contains a CpG
#' site.  Vectorized version of Biostrings::vcountPattern is used.
#' @param mutation_data A dataframe or GRanges object containing the genomic
#' regions of interest in which to look for CpG sites.
#' @param motif Default "CG", which returns CpG sites. You could in theory use
#' an arbitrary string to look at different motifs. Use with caution. In this
#' case the pattern being searched must be a column in the mutation data.
#' @param column_query Default "context" but can be any column  in the mutation
#' data that you wish to look for a motif in.
#' @param ... Additional arguments to vcountPattern()
#' @returns A data frame with the same number of rows as there were ranges in
#' the input, but with an additional metadata column indicating CpG sites in
#' the target sequence of the mutation.
#' @importFrom Biostrings vcountPattern
#' @importFrom rlang .data
#' @export
annotate_cpg_sites <- function(mutation_data,
                               motif = "CG",
                               column_query = "context",
                               ...) {
  annotated_data <- as.data.frame(mutation_data) %>%
    dplyr::mutate(
      CpG_site = Biostrings::vcountPattern(
        pattern = motif,
        subject = .data[[column_query]],
        ...
      )
    ) %>%
    dplyr::mutate(CpG_site = ifelse(.data$CpG_site == 0, FALSE, TRUE))
  return(annotated_data)
}
