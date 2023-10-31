# TO DO: Refine the get seq function now that we changed the regions_df import
# Specify the species and the genome only if you are using a custom regions file. 
# Otherwise, use TS defaults
# Fix custom and ensure genome version works as expected with NULL default. 
#' Get sequence of Duplex Sequencing target regions
#'
#' To replace get_region_seqs.R for its reliance on importing the entire genomes
#' This will create a granges object from the target metadata and import raw nucleotide sequences from ensemble
#' Current defaults are GRCh38 and GRCm39 for human and mouse. Will add to specify genome
#' @param regions_file "human", "mouse", or "custom". The argument refers to the TS Mutagenesis panel of the specified species, or to a custom panel. If custom, provide file path in custom_regions_file. TO DO: add rat.
#' @param custom_regions_file "filepath". If regions_file is set to custom, provide the file path for the tab-delimited file containing regions metadata. Required columns are "contig", "start", and "end".
#' @param rg_sep The delimiter for importing the custom_regions_file. The default is tab-delimited.
#' @param species If a custom regions file is provided, indicate species: "human", "mouse", or "rat"
#' @param genome_version If a custom regions file is provided, indicate the genome assembly version, ex. "GRCm38". Default is human = GRCh38, mouse = GRCm39, rat = mRatBN7
##' @param is_0_based TRUE or FALSE. Are the target region coordinates 0 based (TRUE) or 1 based (FALSE)
#' @param padding An interger value by which the function will extend the range of the target sequence on both sides. Modified region ranges will be reported in ext_start and ext_end. 
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
get_seq <- function( 
                    regions_file = c("human", "mouse", "custom"),
                    custom_regions_file = NULL,
                    species = NULL, 
                    genome_version = NULL, 
                    rg_sep = "\t",
                    is_0_based = TRUE,
                    padding = 1) {
  
  process_region <- function(contig, start, end) {
    if (is_0_based) {
      start <- start + 1
    }
 # Add padding to start and subtract padding from end
  start <- max(1, start - padding)
  end <- end + padding
  
  # State species and genome version defaults for Mutagenesis Panel
  if (regions_file == "human") {
   species_param <- "human"
   genome_version_param <- "GRCh38"
  } else if (regions_file == "mouse") {
     species_param <- "mouse"
     genome_version_param <- "GRCm38"
   } else if (regions_file == "custom") {
    if (!is.null(species)) {
      species_param <- species
      genome_version_param <- genome_version
  } else {
      warning("You must provide a species when using a custom regions file")
  } }
   
    ext <- paste0("https://rest.ensembl.org/sequence/region/", species_param, "/", contig, ":", start, "..", end, ifelse(!is.null(genome_version_param), paste0("?coord_system_version=", genome_version_param), ""))
    r <- httr::GET(paste(ext, sep = ""), httr::content_type("text/plain"))
    return(httr::content(r))
  }
   regions_df <- load_regions_file(regions_file, custom_regions_file, rg_sep)

  seq_list <- lapply(1:nrow(regions_df), function(i) {
    process_region(regions_df$contig[i], regions_df$start[i], regions_df$end[i])
  })

  seqs <- unlist(seq_list)
  ext_df <- data.frame(
                      sequence = seqs,
                      ext_start =  
                        if (is_0_based) {
                        regions_df$start + 1 - padding
                      } else {
                        regions_df$start - padding},
                      ext_end = regions_df$end + padding)


  gr <- GenomicRanges::makeGRangesFromDataFrame(
    df = regions_df,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqnames.field = "contig",
    start.field = "start",
    end.field = "end",
    starts.in.df.are.0based = is_0_based
  )

  gr$sequence <- ext_df$sequence
  gr$ext_start <-ext_df$ext_start
  gr$ext_end <-ext_df$ext_end
  return(gr)
}


# https://rest.ensembl.org/sequence/region/human/chr11:108510787-108513187
# dat <- read.delim("~/DupSeq R Package Building/duplex-sequencing/inst/extdata/genic_regions_mm10.txt")
