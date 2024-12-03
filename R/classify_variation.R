#' classify_variation
#' @description Classify the variation type of a mutation based on its ref and
#' alt values.
#' @param ref The reference allele.
#' @param alt The alternate allele.
#' @return A character indicating the type of variation.
#'
classify_variation <- function(ref, alt) {
  no_variant_indicators <- c(".", "", "<NON_REF>")
  structural_indicators <- c("<DEL>", "<INS>", "<DUP>", "<INV>", "<FUS>",
                             "<CNV>", "<CNV:TR>", "<DUP:TANDEM>", "<DEL:ME>",
                             "<INS:ME>")
  iupac_indicators <- c("R", "K", "S", "Y", "M", "W", "B", "H", "N", "D", "V")

  # Case: No variant site
  # GVCF files sometimes list no_variant sites as <NON_REF> (GATK)
  # We will have to assume that anytime we see <NON_REF> alone, that
  # there is no variant, and if <NON_REF> is followed by an alt allele,
  # there is a variant.
  alt <- gsub("(^|,)<NON_REF>(,|$)", "", alt)
  alt <- gsub("^,|,$", "", alt)  # Trim leading/trailing commas

  if (alt %in% no_variant_indicators || alt == ref) {
    return("no_variant")
  }

  # Case: Structural variants
  if (alt %in% structural_indicators) {
    return("sv")
  }
  # Case: IUPAC ambiguity codes
  if (alt %in% iupac_indicators) {
    return("ambiguous")
  }
  # Case: SNV (Single Nucleotide Variant)
  if (nchar(ref) == 1 && nchar(alt) == 1 && ref != alt) {
    return("snv")
  }
  # Case: MNV (Multi-Nucleotide Variant)
  if (nchar(ref) > 1 && nchar(ref) == nchar(alt) && ref != alt) {
    return("mnv")
  }
   # Case: Insertion
  if (nchar(ref) < nchar(alt) && startsWith(alt, ref)) {
    return("insertion")
  }
  # Case: Deletion
  if (nchar(ref) > nchar(alt) && nchar(alt) == 1 && startsWith(ref, alt)) {
    return("deletion")
  }
  # Case: Complex; ref and alt have different lengths and different base compositions
  if (nchar(ref) != nchar(alt) && !grepl(paste0("^", ref), alt) && !grepl(paste0("^", alt), ref)) {
    return("complex")
  }

  # Otherwise, uncategorized
  return("uncategorized")
}
