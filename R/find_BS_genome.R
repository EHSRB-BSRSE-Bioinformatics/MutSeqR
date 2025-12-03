#' Find the appropriate BS genome for the specified organism and genome.
#' @description This function will browse available BSgenomes, indicating
#' which one should be installed for the specified organism and genome assembly
#' version. If you cannot specify both organism and genome, the function can
#' return a list of available genomes for a specified species.
#' @param organism the name of the organism for which to install the reference
#' genome. This can be the scientific name or a common name. For example Homo
#' Sapiens, H. sapiens, or human
#' @param genome The reference genome assembly version. Ex. hg18, mm10, rn6.
#' @param masked Logical value. Whether to search for the 'masked' BSgenome.
#' Default is FALSE.
#' @export
#' @importFrom BSgenome available.genomes installed.genomes
#' @return a BSgenome package name or a dataframe of possibilities
#' @examples
#' # Find the reference genome for Mouse, mm10 assembly:
#' mouse_mm10 <- find_BS_genome("mouse", "mm10")
#' # Find all possible mouse BS genomes:
#' mouse_all <- find_BS_genome("mouse")
find_BS_genome <- function(organism, genome, masked = FALSE) {

    # Map
    map_keys <- names(MutSeqR::BS_org_map)
    map_vals <- MutSeqR::BS_org_map
    all_canonicals <- rep(map_keys, times = lengths(map_vals))
    all_synonyms <- unlist(map_vals)
    all_canonicals <- c(all_canonicals, map_keys)
    all_synonyms   <- c(all_synonyms, map_keys)
    org_lookup <- setNames(all_canonicals, tolower(all_synonyms))
    input_clean <- gsub("\\.\\s", "", tolower(organism))
    # Perform Lookup
    # We check both the direct input and the cleaned input
    organism_name <- org_lookup[tolower(organism)]
    if (is.na(organism_name)) organism_name <- org_lookup[input_clean] # try clean

    if (is.na(organism_name)) {
        stop("Unrecognized organism name: '", organism, 
            "'. Please consult BSgenome::available.genomes for valid names.")
    }

    # Search Available Genomes
    available <- BSgenome::available.genomes(splitNameParts = TRUE)
    # Filter by Organism and Masked status
    possible_genomes <- available[
        available$organism == organism_name & available$masked == masked, 
    ]

    if (nrow(possible_genomes) == 0) {
        stop("No genomes found for organism '", organism_name, "' with masked = ", masked, ".")
    }

    # Handle Null Genome Argument
    if (missing(genome) || is.null(genome) || nchar(genome) == 0) {
        message(
        "Possible BS genomes for '", organism_name, "' (masked=", masked, "):\n",
        paste(possible_genomes$pkgname, collapse = "\n"),
        "\n\nPlease install one using BiocManager::install('pkgname')."
        )
        # Return the dataframe of possibilities
        return(possible_genomes[, c("pkgname", "organism", "genome", "masked")])
    }

    # Filter by genome
    selected_genome <- possible_genomes[possible_genomes$genome == genome, ]

    if (nrow(selected_genome) == 0) {
        stop(
            "No BS genome found for '", organism_name, "' assembly '", genome,
            "' (masked=", masked, ").\n", "Available assemblies: ",
            paste(unique(possible_genomes$genome), collapse = ", ")
        )
    }
    ref_genome_pkg <- selected_genome$pkgname[1] # first match if duplicated
  
    message("Selected reference genome: ", ref_genome_pkg)
    if (ref_genome_pkg %in% BSgenome::installed.genomes()) {
        message("Reference genome is already installed.")
    } else {
        message("Reference genome is NOT installed. Install using:\n",
            "  BiocManager::install('", ref_genome_pkg, "')")
    }
    message("Once installed, supply '", ref_genome_pkg, "' as the BS_genome parameter.")  
    return(ref_genome_pkg)
}
