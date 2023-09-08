.onLoad <- function(libname, pkgname) {
  reticulate::configure_environment(pkgname, force = TRUE)
}

# Set an environment variable to control package installation
Sys.setenv(DUPSEQR_INSTALL_PACKAGES = "0")

.onLoad <- function(libname, pkgname) {
  install_required_packages <- Sys.getenv("DUPSEQR_INSTALL_PACKAGES")
  
  if (install_required_packages == "1") {
    # List of Bioconductor packages the package depends on
    bioconductor_packages <- c("Biostrings", "GenomeInfoDb", "GenomicRanges", "IRanges", "plyranges", "S4Vectors", "VariantAnnotation")
    
    missing_packages <- bioconductor_packages[!sapply(bioconductor_packages, requireNamespace, quietly = TRUE)]
    
    if (length(missing_packages) > 0) {
      auto_install_options <- c("Yes, install Bioconductor packages", "No, do not install")
      
      user_choice <- menu("Would you like to automatically install required Bioconductor packages?",
                          choices = auto_install_options, title = "Package Installation")
      
      if (user_choice == 1) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        
        packageStartupMessage("Installing missing Bioconductor packages...")
        BiocManager::install(missing_packages, dependencies = TRUE)
        message("Bioconductor packages installed.")
      } else {
        warning("Some required Bioconductor packages are missing. Please install them manually.")
      }
    } else {
      packageStartupMessage("All required Bioconductor packages are already installed.")
    }
  }
}

# In your package installation code, set the environment variable to "1" to trigger package installation
# Sys.setenv(DUPSEQR_INSTALL_PACKAGES = "1")
# install.packages("YourPackage")

