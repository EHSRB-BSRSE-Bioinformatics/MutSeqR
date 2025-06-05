#' Get sequence of genomic target regions
#'
#' Create a GRanges object from the genomic target ranges and import raw
#' nucleotide sequences.
#' @param regions The regions metadata file to import. Can be either a file
#' path, a data frame, or a GRanges object. File paths will be read using
#' the rg_sep. Users can also choose from the built-in TwinStrand's Mutagenesis
#' Panels by inputting "TSpanel_human",  "TSpanel_mouse", or "TSpanel_rat".
#' Required columns for the regions file are "contig", "start", and "end".
#' In a GRanges object, the required columns are "seqnames", "start", and
#' "end".
#' @param rg_sep The delimiter for importing the regions file. The default is
#' tab-delimited ("\t").
#' @param is_0_based_rg A logical variable. Indicates whether the position
#' coordinates in `regions` are 0 based (TRUE) or 1 based (FALSE).
#' If TRUE, positions will be converted to 1-based (start + 1). Need not be
#' supplied for TSpanels. Default is TRUE.
#' @param species The species for which to retrieve the sequences.
#' Species may be given as the scientific name or the common name.
#' Ex. "Human", "Homo sapien". Used to choose the appropriate
#' BS genome. Need not be supplied for TSpanels.
#' @param genome The genome assembly version for which to retrieve the
#' sequences. Used to choose the appropriate genome (BS genome or UCSC).
#' Ex. hg38, hg19, mm10, mm39, rn6, rn7. Need not be supplied for TSpanels.
#' @param masked A logical value indicating whether to use the masked version
#' of the BS genome when retrieving sequences. Default is FALSE.
#' @param padding An integer value by which the function will extend the range
#' of the target sequence on both sides. Start and end coordinates will be
#' adjusted accordingly. Default is 0.
#' @param ucsc A logical value. If TRUE, the function will retrieve the
#' sequences from the UCSC genome browser using an API. If FALSE, the function
#' will retrieve sequences using the appropriate BSgenome package, which will
#' be installed as needed. Default is FALSE.
#' @return a GRanges object with sequences of targeted regions.
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
#' # and supply it to the function as a GRanges object.
#' human <- load_regions_file("TSpanel_human")
#' regions_seq <- get_seq(regions = human,
#'                        is_0_based_rg = FALSE,
#'                        species = "human",
#'                        genome = "hg38",
#'                        masked = FALSE,
#'                        padding = 0)
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom Biostrings getSeq
#' @importFrom BiocGenerics start end
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb seqnames
#'
#' @export
get_seq <- function(regions,
                    rg_sep = "\t",
                    is_0_based_rg = TRUE,
                    species = NULL,
                    genome = NULL,
                    masked = FALSE,
                    padding = 0,
                    ucsc = FALSE) {
  if (ucsc && !requireNamespace("xml2", quietly = TRUE)) {
    stop("The 'xml2' package is required for UCSC API access.")
  }
  if (ucsc && !requireNamespace("httr", quietly = TRUE)) {
    stop("The 'httr' package is required for UCSC API access.")
  }
  regions_gr <- MutSeqR::load_regions_file(regions = regions,
                                           rg_sep = rg_sep,
                                           is_0_based_rg = is_0_based_rg)
  if (is.character(regions)) {
    if (regions == "TSpanel_human") {
      species <- "human"
      genome <- "hg38"
      masked <- FALSE
    }
    if (regions == "TSpanel_mouse") {
      species <- "mouse"
      genome <- "mm10"
      masked <- FALSE
    }
    if (regions == "TSpanel_rat") {
      species <- "rat"
      genome <- "rn6"
      masked <- FALSE
    }
  }

  # Add padding to the region
  BiocGenerics::start(regions_gr) <- pmax(BiocGenerics::start(regions_gr) - padding, 1)
  BiocGenerics::end(regions_gr) <- BiocGenerics::end(regions_gr) + padding

  if (ucsc) {
    # Define the API base URL

    get_sequence_for_region <- function(contig, start, end) {
      base_url <- paste0("https://genome.ucsc.edu/cgi-bin/das/", genome, "/dna")
      params <- list(segment = paste(contig, ":", start, ",", end, sep = ""))
      response <- httr::GET(url = base_url, query = params)
      parsed_xml <- xml2::read_xml(httr::content(response, "text"))
      sequence <- xml2::xml_text(xml2::xml_find_first(parsed_xml, "//DASDNA/SEQUENCE/DNA"))
      cleaned_sequence <- gsub("\n", "", sequence)  # Remove newline characters
      return(toupper(cleaned_sequence))
    }

    # Apply the function to each row of the GR
    seqs <- mapply(get_sequence_for_region,
                   as.vector(GenomeInfoDb::seqnames(regions_gr)),
                   BiocGenerics::start(regions_gr),
                   BiocGenerics::end(regions_gr))
    S4Vectors::mcols(regions_gr)$sequence <- seqs
  } else {
    ref_genome <- install_ref_genome(organism = species,
                                     genome = genome,
                                     masked = masked)
    seqs <- Biostrings::getSeq(ref_genome, names = regions_gr)
    S4Vectors::mcols(regions_gr)$sequence <- seqs
  }
  return(regions_gr)
}