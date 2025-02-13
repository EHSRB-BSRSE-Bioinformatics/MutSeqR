#'Imports the regions file
#'
#' A helper function to import the regions metadata file. It is used in import_mut_data and get_seq. 
#' @param regions "TSpanel_human", "TSpanel_mouse", 'TSpanel_rat', or "custom". The argument refers to the 
#' TS Mutagenesis panel of the specified species, or to a custom panel. If custom,
#' provide file path in custom_regions.
#' @param custom_regions If regions is set to 'custom',
#' provide the regions metadata. Can be either a file path or a data frame.
#' Required columns are "contig", "start", and "end".
#' @param rg_sep The delimiter for importing the custom_regions. The default is tab-delimited
#' @export

load_regions_file <- function(regions,
                              custom_regions = NULL,
                              rg_sep = "\t") {
  # Annotate the mut file with additional information about genomic regions in the file
  if (regions == "TSpanel_human") {
    regions_df <- read.table(system.file("extdata", "human_mutagenesis_panel_hg38.txt", package = "MutSeqR"), header = TRUE)
  } else if (regions == "TSpanel_mouse") {
    regions_df <- read.table(system.file("extdata", "mouse_mutagenesis_panel_mm10.txt", package = "MutSeqR"), header = TRUE)
  } else if (regions == "TSpanel_rat") {
    regions_df <- read.table(system.file("extdata", "rat_mutagenesis_panel_rn6.txt", package = "MutSeqR"), header = TRUE)
  } else if (regions == "custom") {
    
    if (!is.null(custom_regions)) {
      if (is.data.frame(custom_regions)) {
        regions_df <- custom_regions
      } else {
        file <- file.path(custom_regions)
        if (!file.exists(file)) {
          stop("Error: could not load your custom_regions because the file
          path is invalid")
        }
        regions_df <- read.table(file, header = TRUE, sep = rg_sep)
        if (nrow(regions_df) == 0) {
          stop("Error: your imported regions file is empty")
        }
      }
    } else {
      stop("You must provide either a data frame or a filepath to
      custom_regions when regions is set to 'custom'.")
    }
  } else {
    stop("Invalid regions parameter. Choose from 'TSpanel_human',
    'TSpanel_mouse','TSpanel_rat', or 'custom'.")
  }
  return(regions_df)
}