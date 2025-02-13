#' Get sequence of genomic target regions
#'
#' Create a GRanges object from the target metadata and import raw nucleotide 
#' sequences from the UCSC database. https://genome.ucsc.edu
#' @param regions "TSpanel_human", "TSpanel_mouse", "TSpanel_rat, or "custom". The argument refers to the 
#' TS Mutagenesis panel of the specified species, or to a custom panel. 
#' If custom, provide custom region intervals in custom_regions.
#' @param custom_regions If regions is set to custom,
#' provide regions metadata. Can be a file path or a data frame.
#' Required columns are "contig", "start", and "end".
#' @param rg_sep The delimiter for importing the custom_regions. 
#' The default is tab-delimited.
#' @param genome If a custom regions file is provided, indicate the genome 
#' assembly for the reference genome. Please refer to the UCSC genomes. 
#' Ex.Human GRCh38 = hg38 | Human GRCh37 = hg19 | Mouse GRCm38 = mm10 | 
#' Mouse GRCm39 = mm39 | Rat RGSC 6.0 = rn6 | Rat mRatBN7.2 = rn7
#' @param is_0_based TRUE or FALSE. Indicates whether the target region 
#' coordinates are 0 based (TRUE) or 1 based (FALSE). If TRUE, ranges will be converted
#' to 1-based.
#' @param padding An integer value by which the function will extend the range 
#' of the target sequence on both sides. Modified region ranges will be reported 
#' in seq_start and seq_end. Default is 0.
#' @return a GRanges object with sequences and metadata of targeted regions. 
#' Region ranges coordinates will become 1-based.
#' @importFrom httr content content_type GET
#' @importFrom xml2 read_xml xml_text xml_find_first
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export
get_seq <- function(regions,
                    custom_regions = NULL,
                    rg_sep = "\t",
                    genome = NULL,
                    is_0_based = TRUE,
                    padding = 0) {

  regions_df <- MutSeqR::load_regions_file(regions = regions,
                                           custom_regions = custom_regions,
                                           rg_sep = rg_sep)

  if (is_0_based) {
    regions_df$start <- regions_df$start + 1
  }

  regions_df$seq_start <- regions_df$start - padding
  regions_df$seq_end <- regions_df$end + padding

  # Define the API base URL
  # Function to retrieve sequence for a given region
  get_sequence_for_region <- function(contig, start, end) {
    # Specify the genome to be searched
    if (regions == "TSpanel_human") {
      genome <- "hg38"
    } else if (regions == "TSpanel_mouse") {
      genome <- "mm10"
    } else if (regions == "TSpanel_rat") {
      genome <- "rn6"
    } else if (regions == "custom") {
      genome <- genome
      }
    
    base_url <- paste0("https://genome.ucsc.edu/cgi-bin/das/", genome, "/dna")
    params <- list(segment = paste(contig, ":", start, ",", end, sep = ""))
    response <- httr::GET(url = base_url, query = params)
      parsed_xml <- xml2::read_xml(content(response, "text"))
      sequence <- xml2::xml_text(xml2::xml_find_first(parsed_xml, "//DASDNA/SEQUENCE/DNA"))
      cleaned_sequence <- gsub("\n", "", sequence)  # Remove newline characters
      return(toupper(cleaned_sequence))
  }

  # Apply the function to each row of the dataframe
  regions_df$sequence <- mapply(get_sequence_for_region, regions_df$contig, regions_df$seq_start, regions_df$seq_end)

  regions_gr <- GenomicRanges::makeGRangesFromDataFrame(df = as.data.frame(regions_df),
                                                        keep.extra.columns = TRUE,
                                                        seqnames.field = "contig",
                                                        start.field = "start",
                                                        end.field = "end")
  return(regions_gr)
}