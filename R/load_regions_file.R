#'Imports the regions file
#'
#' @description A helper function to import the regions metadata file and
#' return a GRanges object.
#' @param regions The regions metadata file to import. Can be either a file
#' path, a data frame, or a GRanges object. File paths will be read using
#' the rg_sep. Users can also choose from the built-in TwinStrand's Mutagenesis
#' Panels by inputting "TSpanel_human",  "TSpanel_mouse", or "TSpanel_rat".
#' Required columns for the regions file are "contig", "start", and "end".
#' In a GRanges object, the required columns are "seqnames", "start", and "end".
#' @param rg_sep The delimiter for importing the custom_regions. The default is
#' tab-delimited "\t".
#' @param is_0_based_rg A logical variable. Indicates whether the position
#' coordinates in `regions` are 0 based (TRUE) or 1 based (FALSE).
#' If TRUE, positions will be converted to 1-based (start + 1).
#' Need not be supplied for TSpanels. Default is TRUE.
#' @returns a GRanges object of the imported regions metadata file.
#' @export

load_regions_file <- function(regions,
                              rg_sep = "\t",
                              is_0_based_rg = TRUE) {
  if (is(regions, "GRanges")) {
    return(regions)
  } else if (is.data.frame(regions)) {
    regions_df <- regions
  } else if (is.character(regions)) {
    if (regions == "TSpanel_human") {
      regions_df <- read.table(system.file("extdata",
                                           "inputs",
                                           "metadata",
                                           "human_mutagenesis_panel_hg38.txt",
                                           package = "MutSeqR"),
                               header = TRUE)
      is_0_based_rg <- TRUE
    } else if (regions == "TSpanel_mouse") {
      regions_df <- read.table(system.file("extdata",
                                           "inputs",
                                           "metadata",
                                           "mouse_mutagenesis_panel_mm10.txt",
                                           package = "MutSeqR"),
                               header = TRUE)
      is_0_based_rg <- TRUE
    } else if (regions == "TSpanel_rat") {
      regions_df <- read.table(system.file("extdata",
                                           "inputs",
                                           "metadata",
                                           "rat_mutagenesis_panel_rn6.txt",
                                           package = "MutSeqR"),
                               header = TRUE)
      is_0_based_rg <- TRUE
    } else {
      regions_file <- file.path(regions)
      # Check if the file exists
      if (!file.exists(regions_file)) {
        stop("Error: could not load your regions file because the file path is invalid.")
      }
      regions_df <- read.table(regions_file, header = TRUE, sep = rg_sep)
      if (nrow(regions_df) == 0) {
        stop("Error: your imported regions file is empty.")
      }
      if (ncol(regions_df) == 1) {
        stop("Error: your imported regions file has only one column. Please check the delimiter in rg_sep.")
      }
    }
  } else {
    stop("Invalid regions parameter.")
  }
  if (!all(c("contig", "start", "end") %in% colnames(regions_df))) {
    stop("Error: your regions file is missing the required columns 'contig', 'start', and 'end'.")
  }
  # Turn regions_df into a GRanges object
  regions_gr <- GenomicRanges::makeGRangesFromDataFrame(
      df = regions_df,
      keep.extra.columns = TRUE,
      seqnames.field = "contig",
      start.field = "start",
      end.field = "end",
      starts.in.df.are.0based = is_0_based_rg
  )
  return(regions_gr)
}