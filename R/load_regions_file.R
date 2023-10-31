#'Imports the regions file
#'
#' @param regions_file "human", "mouse", or "custom". The argument refers to the 
#' TS Mutagenesis panel of the specified species, or to a custom panel. If custom,
#'  provide file path in custom_regions_file. TO DO: add rat.
#' @param custom_regions_file "filepath". If regions_file is set to custom, 
#' provide the file path for the tab-delimited file containing regions metadata.
#'  Required columns are "contig", "start", and "end".
#'  @param rg_sep The delimiter for importing the custom_regions_file
#' @export

load_regions_file <- function(regions_file, 
                              custom_regions_file = NULL,
                              rg_sep = "\t") {
  # Annotate the mut file with additional information about genomic regions in the file
  if (regions_file == "human") {
    genic_regions <- read.table(system.file("extdata", "genic_regions_hg38.txt", package = "DupSeqR"), header = TRUE)
  } else if (regions_file == "mouse") {
    genic_regions <- read.table(system.file("extdata", "genic_regions_mm10.txt", package = "DupSeqR"), header = TRUE)
  } else if (regions_file == "custom") {
    if (!is.null(custom_regions_file)) {
      genic_regions <- read.table(custom_regions_file, header = TRUE, sep = rg_sep)
    } else {
      warning("You must provide a file path to custom_regions_file when regions_file is set to 'custom'.")
    }
  } else {
    warning("Invalid regions_file parameter. Choose from 'human', 'mouse', or 'custom'.")
  }
  return(genic_regions)
}