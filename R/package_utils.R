# .onLoad <- function(libname, pkgname) {
#   reticulate::configure_environment(pkgname, force = TRUE)
# }
#' Set up Python environment for MutSeqR
#'
#' This function initializes the Python environment used by MutSeqR.
#' It is not run automatically to avoid issues during installation and checks.
#'
#' @param force Logical. Whether to force reconfiguration even if an environment already exists.
#' @export
setup_mutseqr_python <- function(force = FALSE) {
  reticulate::configure_environment("MutSeqR", force = force)
}

.onAttach <- function(libname, pkgname) {
  print_ascii_art()
}
#' This function prints ASCII art when the package is loaded
#' @export
print_ascii_art <- function() {
  art_path <- system.file("extdata", "ASCII_art_MutSeqR.txt", package = "MutSeqR")
  art <- readLines(art_path)
  packageStartupMessage(paste(art, collapse = "\n"))
}
