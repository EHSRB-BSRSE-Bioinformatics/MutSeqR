.onLoad <- function(libname, pkgname) {
  reticulate::configure_environment(pkgname, force = TRUE)
  
  print_ascii_art()

  }

#' This function prints ASCII art when the package is loaded
#' @export
print_ascii_art <- function() {
  art_path <- system.file("extdata", "ASCII_art_MutSeqR.txt", package = "MutSeqR")
  art <- readLines(art_path)
  cat(art, sep = "\n")
}
