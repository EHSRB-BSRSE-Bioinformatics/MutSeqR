#' Install the reference genome for the specified organism.
#' @description This function will use BSgenome to install the reference genome
#' for a specified organism and assembly version.
#' @param organism the name of the organism for which to install the reference genome.
#' This can be the scientific name or a common name. For example Homo Sapiens, H. sapiens, or human
#' @param genome The reference genome assembly version. Ex. hg18, mm10, rn6.
#' @param masked Logical value. Whether to search for the 'masked' BSgenome. Default is FALSE.
#' @importFrom BiocManager install
#' @export
#' @return a BSgenome object
install_ref_genome <- function(organism, genome, masked = FALSE) {
  if (!requireNamespace("BSgenome", quietly = TRUE)) {
    stop("Package BSgenome is required. Please install from Bioconductor.")
  }
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
  # Map the input name to the Organism name in available.genomes
  organism <- gsub("\\.\\s", "", organism)  # Clean up the input: collapse scientifc names
  convertToOrganismName <- function(name) {
    for (org_name in names(name_map)) {
      # Convert both user input and name_map keys to lowercase for case-insensitive comparison
      if (tolower(name) %in% c(tolower(org_name),
        tolower(name_map[[org_name]]))) {
        return(org_name)
      }
    }
    stop("Unrecognized organism name: ", name, "Please consult BSgenome::available.genomes for a full list of the available genomes and their associated organism names.")
  }
  organism_name <- convertToOrganismName(organism)

  # Choose the genome associated with the organism and genome assembly version
  available_genomes <- BSgenome::available.genomes(splitNameParts = TRUE)
  selected_genome <- available_genomes[available_genomes$organism == organism_name & available_genomes$genome == genome, ]  
  if (nrow(selected_genome) == 0) {
    stop("No reference genome found for the specified organism and assembly version. Please consult BSgenome::available.genomes for a full list of the available genomes and their associated organism names.")
  } else if (nrow(selected_genome) > 1) {
    if(masked) {
      selected_genome <- selected_genome[selected_genome$masked == TRUE, ]
    } else {
      selected_genome <- selected_genome[selected_genome$masked == FALSE, ]
    }
  } else {
    selected_genome <- selected_genome
  }

  ref_genome <- selected_genome$pkgname
  # Install the reference genome
  installed_BS_genomes <- BSgenome::installed.genomes()
  if (ref_genome %in% installed_BS_genomes) {
    message("Reference genome already installed.")
  } else {
    message("Installing reference genome: ", ref_genome, "from Bioconductor.")
    BiocManager::install(ref_genome)
  }
  message("Loading reference genome: ", ref_genome, ".")
  reference_genome <- suppressMessages(BSgenome::getBSgenome(ref_genome))

  return(reference_genome)
}