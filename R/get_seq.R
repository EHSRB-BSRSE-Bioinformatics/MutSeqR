#' Get sequence of Duplex Sequencing target regions
#'
#' To replace get_region_seqs.R for its reliance on importing the entire genomes
#' This will create a granges object from the target metadata and import raw nucleotide sequences from ensemble
#' Current defaults are GRCh38 and GRCm39 for human and mouse. Will add to specify genome
#' @param species "human", "mouse", or "rat"
#' @param genome_version "Genome version", ex. "GRCm38". Default = NULL (no version specified; human = GRCh38, mouse = GRCm39, rat = mRatBN7)
#' @param regions_df data frame with target locations. Contains columns: contig, start, and end
#' @param is_0_based TRUE or FALSE. Are the target region coordinates 0 based (TRUE) or 1 based (FALSE)
#' @return a GRanges object with sequences and metadata of targeted regions
#' @examples
#' regions_df <- data.frame(
#'   contig = c("chr11", "chr13"),
#'   start = c(108510788, 75803913),
#'   end = c(108513187, 75806312),
#'   gene = c("GeneA", "GeneB"),
#'   transcription_status = c("genic", "intergenic")
#' )
#' t <- get_seq(species = "human", genome_version = "GRCh37", regions_df = regions_df)
#' t$sequence
#' @export
get_seq <- function(species, genome_version = NULL, regions_df, is_0_based = TRUE) {
  process_region <- function(contig, start, end) {
    if (is_0_based) {
      start <- start + 1
    }

    ext <- paste0("https://rest.ensembl.org/sequence/region/", species, "/", contig, ":", start, "..", end, ifelse(!is.null(genome_version), paste0("?coord_system_version=", genome_version), ""))
    r <- httr::GET(paste(ext, sep = ""), httr::content_type("text/plain"))
    return(httr::content(r))
  }

  seq_list <- lapply(1:nrow(regions_df), function(i) {
    process_region(regions_df$contig[i], regions_df$start[i], regions_df$end[i])
  })

  seqs <- unlist(seq_list)
  seqs_df <- data.frame(sequence = seqs)

  gr <- GenomicRanges::makeGRangesFromDataFrame(
    df = regions_df,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqnames.field = "contig",
    start.field = "start",
    end.field = "end",
    starts.in.df.are.0based = is_0_based
  )

  gr$sequence <- seqs_df$sequence
  return(gr)
}


# https://rest.ensembl.org/sequence/region/human/chr11:108510787-108513187
# dat <- read.delim("~/DupSeq R Package Building/duplex-sequencing/inst/extdata/genic_regions_mm10.txt")
