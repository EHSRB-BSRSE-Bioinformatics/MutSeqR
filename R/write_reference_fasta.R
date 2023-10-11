#' Write FASTA file of reference sequences
#'
#' In some cases, you might want to generate an arbitrary multi-sequence FASTA
#' file from GRanges including the reference sequences. This function allows you
#' to do this.
#' @param ref_ranges A GRanges object including the sequences of the reference
#' regions included for the data.
#' @param fasta_out The path for the FASTA file output.
#' @returns Writes a FASTA reference file to be used in downstream applications.
#' @importFrom Biostrings writeXStringSet
#' @export
write_reference_fasta <- function(ref_ranges,
                                  fasta_out = "./reference_output.fasta") {
  ref <- ref_ranges$sequence
  names(ref) <- make.unique(as.vector(seqnames(ref_ranges)))
  Biostrings::writeXStringSet(ref,
                  fasta_out,
                  append=FALSE,
                  format="fasta")
}
