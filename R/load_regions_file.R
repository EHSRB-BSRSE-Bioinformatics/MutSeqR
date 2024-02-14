#'Imports the regions file
#'
#' A helper function to import the regions metadata file. It is used in import_mut_data and get_seq. 
#' @param regions "human", "mouse", 'rat', or "custom". The argument refers to the 
#' TS Mutagenesis panel of the specified species, or to a custom panel. If custom,
#' provide file path in custom_regions_file.
#' @param custom_regions_file "filepath". If regions is set to custom, 
#' provide the file path for the tab-delimited file containing regions metadata. 
#' Required columns are "contig", "start", and "end".
#' @param rg_sep The delimiter for importing the custom_regions_file. The default is tab-delimited
#' @export

load_regions_file <- function(regions, 
                              custom_regions_file = NULL,
                              rg_sep = "\t") {
  # Annotate the mut file with additional information about genomic regions in the file
  if (regions == "human") {
    regions_df <- read.table(system.file("extdata", "human_mutagenesis_panel_hg38.txt", package = "MutSeqR"), header = TRUE)
  } else if (regions == "mouse") {
    regions_df <- read.table(system.file("extdata", "mouse_mutagenesis_panel_mm10.txt", package = "MutSeqR"), header = TRUE)
  } else if (regions == "rat") {
    regions_df <- read.table(system.file("extdata", "rat_mutagenesis_panel_rn6.txt", package = "MutSeqR"), header = TRUE)
  } else if (regions == "custom") {
    if (!is.null(custom_regions_file)) {
      regions_df <- read.table(custom_regions_file, header = TRUE, sep = rg_sep)
    } else {
      warning("You must provide a file path to custom_regions_file when regions is set to 'custom'.")
    }
  } else {
    warning("Invalid regions parameter. Choose from 'human', 'mouse', 'rat', or 'custom'.")
  }
  return(regions_df)
}