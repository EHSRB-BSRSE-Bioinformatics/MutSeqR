#' get target sequence from ensembl
#' 
#'@param species "human" or "mouse"
#'@param chr chromosome of target region "chr1"
#'@param start start location of target region. Numeric
#'@param end end location of target region. Numeric
#'@return plain text character set of target sequence
#'@examples
#' get_seq("human", "chr11", 108510787, 108513187)
#' @export
get_seq <- function(species, chr, start, end) {
  ext <- paste0("https://rest.ensembl.org/sequence/region/", species, "/", chr, ":", start, "..", end)
  r <-  httr::GET(paste(ext, sep = ""), httr::content_type("text/plain"))
  print(httr::content(r)) 
}


#ext <- "/sequence/region/human/chr11:108510787-108513187?coord_system_version=GRCh38"
#ext <- "/sequence/region/mouse/chr13:74030071-74032471?coord_system_version=GRCm38"

#https://rest.ensembl.org/sequence/region/human/chr11:108510787-108513187