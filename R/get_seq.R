#' Get sequence of genomic target regions
#'
#' Create a GRanges object from the genomic target ranges and import raw
#' nucleotide sequences.
#' @param regions The file containing the genomic regions. Regions TwinStrand
#' Mutagenesis Panels are stored in the package files and can be accessed using
#' the values `TSpanel_human`, `TSpanel_mouse`, and `TSpanel_rat`. If you have
#' a custom range of genomic regions, set the value to `custom` and provide the
#' regions file using the `custom_regions` argument.
#' @param custom_regions If `regions = "custom"`, provide  the genomic ranges.
#' Can be a file path or a data frame. Required columns are `contig`, `start`,
#' and `end`.
#' @param rg_sep The delimiter for importing the custom_regions. The default is
#' tab-delimited.
#' @param is_0_based A logical variable. Indicates whether the position
#' coordinates in the `custom_regions` are 0 based (TRUE) or 1 based (FALSE).
#' If TRUE, positions will be converted to 1-based. Need not be supplied for
#' TSpanels.
#' @param species The species for which to retrieve the sequences.
#' Species may be given as the scientific name or the common name.
#' Ex. "Human", "Homo sapien". Used to choose the appropriate
#' BS genome. Need not be supplied for TSpanels.
#' @param genome The genome assembly version for which to retrieve the
#' sequencies. Used to choose the appropriate genome (BS genome or UCSC).
#' Ex. hg38, hg19, mm10, mm39, rn6, rn7. Need not be supplied for TSpanels.
#' @param masked A logical value indicating whether to use the masked version
#' of the BS genome when retrieving sequences. Default is FALSE.
#' @param padding An integer value by which the function will extend the range
#' of the target sequence on both sides. Start and end coordinates will be
#' adjusted accordingly. Default is 0.
#' @param ucsc A logical value. If TRUE, the function will retrieve the
#' sequences from the UCSC genome browser using an API. If FALSE, the function
#' will retrieve sequencies using the appropriate BSgenome package, which will
#' be installed as needed.
#' @return a GRanges object with sequences of targeted regions. Region ranges
#' coordinates will become 1-based.
#' @details Consult
#' \code{available.genomes(splitNameParts=FALSE, type=getOption("pkgType"))}
#' for a full list of the available BS genomes and their associated
#' species/genome/masked values. The BSgenome package will be installed if
#' not already available. If using the UCSC API, the function will retrieve
#' the sequences from the UCSC genome browser using the DAS API. See the
#' UCSC website for available genomes: \url{https://genome.ucsc.edu}.
#' @examples
#' # Example 1: Retrieve the sequences for TwinStrand Mouse Mutagenesis Panel
#' regions_seq <- get_seq(regions = "TSpanel_mouse")
#' 
#' # Example 2: Retrieve the sequences for custom regions
#' # We will load the TSpanel_human regions file as an example
#' human <- load_regions_file("TSpanel_human")
#' regions_seq <- get_seq(regions = "custom",
#'                        custom_regions = human,
#'                        is_0_based = TRUE,
#'                        species = "human",
#'                        genome = "hg38",
#'                        masked = FALSE,
#'                        padding = 0)
#' @importFrom httr content content_type GET
#' @importFrom xml2 read_xml xml_text xml_find_first
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export
get_seq <- function(regions,
                    custom_regions = NULL,
                    rg_sep = "\t",
                    is_0_based = TRUE,
                    species = NULL,
                    genome = NULL,
                    masked = FALSE,
                    padding = 0,
                    ucsc = FALSE) {

  regions_df <- MutSeqR::load_regions_file(regions = regions,
                                           custom_regions = custom_regions,
                                           rg_sep = rg_sep)
  if (regions == "TSpanel_human") {
    species <- "human"
    genome <- "hg38"
  }
  if (regions == "TSpanel_mouse") {
    species <- "mouse"
    genome <- "mm10"
  }
  if (regions == "TSpanel_rat") {
    species <- "rat"
    genome <- "rn6"
  }
  if (regions %in% c("TSpanel_mouse", "TSpanel_human", "TSpanel_rat")) {
    is_0_based <- TRUE
    masked <- FALSE
  }
  # Change to 1-based
  regions_df$start <- regions_df$start + as.numeric(is_0_based)
  # Add padding to the region
  regions_df$start <- regions_df$start - padding
  regions_df$end <- regions_df$end + padding

  if (ucsc) {
    # Define the API base URL
    get_sequence_for_region <- function(contig, start, end) {
      base_url <- paste0("https://genome.ucsc.edu/cgi-bin/das/", genome, "/dna")
      params <- list(segment = paste(contig, ":", start, ",", end, sep = ""))
      response <- httr::GET(url = base_url, query = params)
      parsed_xml <- xml2::read_xml(content(response, "text"))
      sequence <- xml2::xml_text(xml2::xml_find_first(parsed_xml, "//DASDNA/SEQUENCE/DNA"))
      cleaned_sequence <- gsub("\n", "", sequence)  # Remove newline characters
      return(toupper(cleaned_sequence))
    }

    # Apply the function to each row of the dataframe
    regions_df$sequence <- mapply(get_sequence_for_region,
                                  regions_df$contig,
                                  regions_df$start,
                                  regions_df$end)
    regions_gr <- GenomicRanges::makeGRangesFromDataFrame(df = as.data.frame(regions_df),
                                                          keep.extra.columns = TRUE,
                                                          seqnames.field = "contig",
                                                          start.field = "start",
                                                          end.field = "end")
  } else {
    ref_genome <- install_ref_genome(organism = species,
                                     genome = genome,
                                     masked = masked)
    regions_gr <- GenomicRanges::makeGRangesFromDataFrame(
      df = regions_df,
      keep.extra.columns = TRUE,
      seqnames.field = "contig",
      start.field = "start",
      end.field = "end"
    )
    seqs <- Biostrings::getSeq(ref_genome, names = regions_gr)
    S4Vectors::mcols(regions_gr)$sequence <- seqs
  }
  return(regions_gr)
}