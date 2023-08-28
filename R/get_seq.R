#' Get sequence of Duplex Sequencing target regions
#'
#' To replace get_region_seqs.R for its relaince on importing the entire genomes
#' This will create a granges object from the target metadata and import raw nucleotide sequences from ensemble 
#' Current defaults are GRCh38 and GRCm39 for human and mouse. Will add to specify genome
#' The code is functional, but inconsistent. It will work, but occasionally it will not be able to retrieve the sequences from ensembl
#'   Perhaps an issue w network connectivity or the ensemble api server load

#' @param species "human", "mouse", or "rat"
#' @param genome_version "Genome version", ex. "GRCm38". Default = NULL (no version specified; human = GRCh38, mouse = GRCm39, rat = mRatBN7)
#' @param regions_df data frame with target locations. Contains columns: contig, start, and end
#' @param is_0_based TRUE or FALSE. Are the target region coordinates 0 based (TRUE) or 1 based (FALSE)
#' @return a GRanges object with sequences and metadata of targeted regions
#' @examples
#' regions_df <- data.frame(contig = c("chr11", "chr13"),
#'                      start = c(108510788, 75803913),
#'                      end = c(108513187, 75806312),
#'                      gene = c("GeneA", "GeneB"),
#'                      transcription_status = c("genic", "intergenic"))
#' t <- get_seq(speceis = "human", genome_version = "GRCh37", regions_df = regions_df)
#' t$sequence
#' @export
get_seq <- function(species, genome_version = NULL, regions_df, is_0_based = TRUE) {
  seq_list <- list()
  for (i in 1:nrow(regions_df)) {
    contig <- regions_df$contig[i]
    start <- regions_df$start[i]
    end <- regions_df$end[i]
    
    # Adjust start based on is_0_based parameter
    if (is_0_based) {
      start <- start + 1
    }
    
    ext <- paste0("https://rest.ensembl.org/sequence/region/", species, "/", contig, ":", start, "..", end, ifelse(!is.null(genome_version), paste0("?coord_system_version=", genome_version), ""))
    r <- httr::GET(paste(ext, sep = ""), httr::content_type("text/plain"))
    seq_list[[i]] <- httr::content(r)
  }
  
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


#https://rest.ensembl.org/sequence/region/human/chr11:108510787-108513187
