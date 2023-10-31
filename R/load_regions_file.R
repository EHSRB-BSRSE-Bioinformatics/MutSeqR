#' @export
load_regions_file <- function(regions_file, custom_regions_file = NULL) {
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