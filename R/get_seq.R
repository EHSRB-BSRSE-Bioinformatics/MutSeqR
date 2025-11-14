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
#' tab-delimited ("\\t").
#' @param is_0_based_rg A logical variable. Indicates whether the position
#' coordinates in `regions` are 0 based (TRUE) or 1 based (FALSE).
#' If TRUE, positions will be converted to 1-based (start + 1). Need not be
#' supplied for TSpanels. Default is TRUE.
#' @param BS_genome The name of the appropriate BSgenome package to use
#' for sequence retrieval. Ex. "BSgenome.Hsapiens.UCSC.hg38",
#' "BSgenome.Mmusculus.UCSC.mm10", "BSgenome.Rnorvegicus.UCSC.rn6".
#' Use the function find_BS_genome() to help identify the appropriate BSgenome
#' package if needed. Need not be supplied for TSpanels.
#' BS_genome must be installed if using this method.
#' @param ucsc A logical value. If TRUE, the function will retrieve the
#' sequences from the UCSC genome browser using an API. If FALSE, the function
#' will retrieve sequences using the appropriate BSgenome package, which will
#' be installed as needed. Default is FALSE.
#' @param species The species for which to retrieve the sequences.
#' Only required if using the UCSC method.
#' Species may be given as the scientific name or the common name.
#' Ex. "Human", "Homo sapien". Used to choose the appropriate
#' BS genome. Need not be supplied for TSpanels.
#' @param genome The genome assembly version for which to retrieve the
#' sequences. Only required if using the UCSC method.
#' Ex. hg38, hg19, mm10, mm39, rn6, rn7. Need not be supplied for TSpanels.
#' @param padding An integer value by which the function will extend the range
#' of the target sequence on both sides. Start and end coordinates will be
#' adjusted accordingly. Default is 0.
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
#' regions_seq <- get_seq(
#'   regions = human,
#'   is_0_based_rg = FALSE,
#'   BS_genome = "BSgenome.Hsapiens.UCSC.hg38",
#'   padding = 0
#' )
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom Biostrings getSeq
#' @importFrom BiocGenerics start end
#' @importFrom S4Vectors mcols
#' @importFrom Seqinfo seqnames
#' @importFrom BSgenome getBSgenome installed.genomes
#' @export

get_seq <- function(regions,
                    rg_sep = "\t",
                    is_0_based_rg = TRUE,
                    padding = 0,
                    BS_genome = NULL,
                    ucsc = FALSE,
                    species = NULL,
                    genome = NULL) {
    stopifnot(
        "regions must be a file path, data frame, or GRanges object." =
            is.character(regions) ||
            is.data.frame(regions) ||
            methods::is(regions, "GRanges"),
        "rg_sep must be a character" = is.character(rg_sep),
        "is_0_based_rg must be a logical" = is.logical(is_0_based_rg),
        "padding must be a non-negative integer" =
            is.numeric(padding) && padding >= 0 && padding == floor(padding),
        "BS_genome must be a character or NULL" =
            is.null(BS_genome) || is.character(BS_genome),
        "ucsc must be a logical" = is.logical(ucsc),
        "species must be a character or NULL" =
            is.null(species) || is.character(species),
        "genome must be a character or NULL" =
            is.null(genome) || is.character(genome)
    )

    if (ucsc && !requireNamespace("xml2", quietly = TRUE)) {
        stop("The 'xml2' package is required for UCSC API access.")
    }
    if (ucsc && !requireNamespace("httr", quietly = TRUE)) {
        stop("The 'httr' package is required for UCSC API access.")
    }
    regions_gr <- MutSeqR::load_regions_file(
        regions = regions,
        rg_sep = rg_sep,
        is_0_based_rg = is_0_based_rg
    )
    if (is.character(regions)) {
        if (regions == "TSpanel_human") {
            BS_genome <- "BSgenome.Hsapiens.UCSC.hg38"
            species <- "human"
            genome <- "hg38"
        }
        if (regions == "TSpanel_mouse") {
            BS_genome <- "BSgenome.Mmusculus.UCSC.mm10"
            species <- "mouse"
            genome <- "mm10"
        }
        if (regions == "TSpanel_rat") {
            BS_genome <- "BSgenome.Rnorvegicus.UCSC.rn6"
            species <- "rat"
            genome <- "rn6"
        }
    }

    # Add padding to the region
    BiocGenerics::start(regions_gr) <- pmax(
        BiocGenerics::start(regions_gr) - padding, 1
    )
    BiocGenerics::end(regions_gr) <- BiocGenerics::end(regions_gr) + padding

    if (ucsc) {
        # Define the API base URL
        get_sequence_for_region <- function(contig, start, end) {
            base_url <- paste0(
                "https://genome.ucsc.edu/cgi-bin/das/",
                genome, "/dna"
            )
            params <- list(
                segment = paste(contig, ":", start, ",", end, sep = "")
            )
            response <- httr::GET(url = base_url, query = params)
            parsed_xml <- xml2::read_xml(httr::content(response, "text"))
            sequence <- xml2::xml_text(
                xml2::xml_find_first(parsed_xml, "//DASDNA/SEQUENCE/DNA")
            )
        cleaned_sequence <- gsub("\n", "", sequence)
        return(toupper(cleaned_sequence))
    }

    # Apply the function to each row of the GR
    seqs <- mapply(
        get_sequence_for_region,
        as.vector(Seqinfo::seqnames(regions_gr)),
        BiocGenerics::start(regions_gr),
        BiocGenerics::end(regions_gr)
    )
    S4Vectors::mcols(regions_gr)$sequence <- seqs
    } else {
        if (is.null(BS_genome)) {
            stop(
                "If not using the UCSC method, please indicate the",
                " appropriate BS genome (must be installed). If you",
                " are not sure which BS genome to use, please provide",
                " the species and reference genome to find_BS_genome()."
            )
        }
        installed_BS_genomes <- BSgenome::installed.genomes()
        if (!(BS_genome %in% installed_BS_genomes)) {
            stop(
                "The specified BS genome is not installed. Please install the",
                " appropriate BS genome using BiocManager::install('pkgname')",
                " where pkgname is the name of the BSgenome package. If you",
                " are not sure which BS genome to use, please provide the",
                " species and reference genome to find_BS_genome()."
            )
        }
        message("Loading reference genome: ", BS_genome, ".")
        ref_genome <- BSgenome::getBSgenome(BS_genome)
        seqs <- Biostrings::getSeq(ref_genome, names = regions_gr)
        S4Vectors::mcols(regions_gr)$sequence <- seqs
    }
    return(regions_gr)
}
