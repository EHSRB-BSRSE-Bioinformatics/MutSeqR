.onLoad <- function(libname, pkgname) {
  reticulate::configure_environment(pkgname, force = TRUE)
}

utils::globalVariables(c("mutation_data"))