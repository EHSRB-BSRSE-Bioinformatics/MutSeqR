#' @import BiocManager

.onAttach <- function(libname, pkgname) {
  auto_install_options <- c("Yes, install Bioconductor packages", "No, do not install")
  
  user_choice <- menu("Would you like to automatically install required Bioconductor packages?",
                      choices = auto_install_options, title = "Package Installation")
  
  if (user_choice == 1) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    
    # List of Bioconductor packages your package depends on
    bioconductor_packages <- c("Biostrings", "VariantAnnotation", "trackViewer", "GenVisR", "plyranges", "GenomicRanges", "BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Mmusculus.UCSC.mm10,")
    
    # Install missing Bioconductor packages
    missing_packages <- bioconductor_packages[!sapply(bioconductor_packages, requireNamespace, quietly = TRUE)]
    
    if (length(missing_packages) > 0) {
      message("Installing missing Bioconductor packages...")
      BiocManager::install(missing_packages, dependencies = TRUE)
      message("Bioconductor packages installed.")
    } else {
      message("All required Bioconductor packages are already installed.")
    }
  }
}
