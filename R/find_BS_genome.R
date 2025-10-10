#' Find the appropriate BS genome for the specified organism and genome.
#' @description This function will browse available BSgenomes, indicating
#' which one should be installed for the specified organism and genome assembly
#' version. If you cannot specify both organism and genome, the function can return
#' a list of available genomes for a specified species.
#' @param organism the name of the organism for which to install the reference genome.
#' This can be the scientific name or a common name. For example Homo Sapiens, H. sapiens, or human
#' @param genome The reference genome assembly version. Ex. hg18, mm10, rn6.
#' @param masked Logical value. Whether to search for the 'masked' BSgenome. Default is FALSE.
#' @export
#' @importFrom BSgenome available.genomes installed.genomes
#' @return a BSgenome object
#' @examples
#' # Find the reference genome for Mouse, mm10 assembly:
#' mouse_mm10 <- find_BS_genome("mouse", "mm10")
#' # Find all possible mouse BS genomesL
#' mouse_all <- find_BS_genome("mouse")
find_BS_genome <- function(organism, genome, masked = FALSE) {
  # Common name mapping
  name_map <- list(
    "Alyrata" = "arabidopsis lyrata",
    "Amellifera" = c("apis mellifera", "honey bee"),
    "Aofficinalis" = c("asparagus officinalis", "asparagus"),
    "Athaliana" = c("arabidopsis thaliana", "arabidopsis"),
    "Btaurus" = c("bos Taurus", "cow"),
    "Carietinum" = c("cicer arietinum", "chickpea"),
    "Celegans" = c("caenorhabditis elegans", "roundworm", "nematode", "worm"),
    "Cfamiliaris" = c("canis lupus familiaris", "dog"),
    "Cjacchus" = c("callithrix jacchus", "marmoset"),
    "CneoformansVarGrubiiKN99" = "cryptococcus neoformans var. grubii KN99",
    "Creinhardtii" = "chlamydomonas reinhardtii",
    "Dmelanogaster" = "drosophila melanogaster",
    "Drerio" = c("danio rerio", "zebrafish"),
    "Dvirilis" = c("drosophila virilis"),
    "Ecoli" = "escherichia coli",
    "Gaculeatus" = c("gasterosteus aculeatus", "stickleback", "three-spined stickleback"),
    "Ggallus" = c("gallus gallus", "chicken"),
    "Gmax" = c("glycine max", "soybean"),
    "Hsapiens" = c("homo sapiens", "homo sapiens sapiens", "human"),
    "Mdomestica" = c("monodelphis domestica", "opossum", "gray short-tailed opossum"),
    "Mfascicularis" = c("macaca fascicularis", "long-tailed macaque", "crab-eating macaque"),
    "Mfuro" = c("mustela putorius furo", "ferret"),
    "Mmulatta" = c("macaca mulatta", "rhesus macaque"),
    "Mmusculus" = c("mus musculus", "mouse", "house mouse"),
    "Osativa" = c("oryza sativa", "rice"),
    "Ppaniscus" = c("pan paniscus", "bonobo"),
    "Ptroglodytes" = c("pan troglodytes", "chimp", "chimpanzee"),
    "Rnorvegicus" = c("rattus norvegicus", "rat", "brown rat"),
    "Scerevisiae" = c("saccharomyces cerevisiae", "yeast", "brewer's yeast"),
    "Sscrofa" = c("sus scrofa", "pig", "wild boar"),
    "Tgondii" = "toxoplasma gondii",
    "Tguttata" = c("taeniopygia guttata", "zebra finch"),
    "Vvinifera" = c("vitis vinifera", "grape")
  )
  # Map the input name to the organism name in available.genomes
  organism <- gsub("\\.\\s", "", organism)
  convertToOrganismName <- function(name) {
    for (org_name in names(name_map)) {
      if (tolower(name) %in% c(tolower(org_name), tolower(name_map[[org_name]]))) {
        return(org_name)
      }
    }
    stop("Unrecognized organism name: ", name,
         ". Please consult BSgenome::available.genomes for valid names.")
  }
  organism_name <- convertToOrganismName(organism)
  # Search available genomes
  available_genomes <- BSgenome::available.genomes(splitNameParts = TRUE)
  
  possible_genomes <- available_genomes[
    available_genomes$organism == organism_name & 
      available_genomes$masked == masked, ]
  if (nrow(possible_genomes) == 0) {
      stop("No genomes found for specified organism.")
  } 
  # If genome argument is missing, return possibilities
  if (missing(genome) || is.null(genome) || nchar(genome)==0) {
    message("Possible BS genomes for organism = '", organism, "', masked = ", masked, ":")
    print(possible_genomes[, c("pkgname", "genome", "masked")])
    message("Please install one of the possible BS genomes using BiocManager::install('pkgname') and provide the pkgname to import_mut/vcf_data().")
    return(possible_genomes[, c("pkgname", "organism", "genome", "masked")])
  }
  # If genome specified, filter further for genome assembly
  selected_genome <- possible_genomes[
    possible_genomes$genome == genome, ]
  if (nrow(selected_genome) == 0) {
    stop("No BS genome found for the specified organism, assembly version and masked setting.\n",
      "Available assemblies for this organism (masked = ", masked, ") are:\n",
      paste(unique(possible_genomes$genome), collapse=", "))
  }

  ref_genome <- selected_genome$pkgname
  message("Selected reference genome: ", ref_genome)
  # Install the reference genome
  installed_BS_genomes <- BSgenome::installed.genomes()
  if (ref_genome %in% installed_BS_genomes) {
    message("Reference genome is already installed.")
  } else {
    message("Reference genome is not installed. Please install using BiocManager::install(", ref_genome, ")")
  }
  message("Once installed, supply: ", ref_genome, " as the BS_genome parameter in import_mut/vcf_data()")
  return(ref_genome)
}