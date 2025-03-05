#' Write FASTA file of reference sequences.
#'
#' @details Generate an arbitrary multi-sequence FASTA file from GRanges
#' including the reference sequences.
#' @param regions_gr A GRanges object including the sequences of the reference
#' regions included for the data. This can be generated from the `get_seq`
#' function.
#' @param output_path The directory where the FASTA file should be written.
#' Default is NULL, which will write the file to the current working directory.
#' @returns Writes a FASTA reference file "reference_output.fasta".
#' If multiple ranges are included in the GRanges object, the sequences will be
#' written to a single FASTA file. Sequences names will be the seqnames
#' (contig) of the range.
#' @examples
#' \dontrun{
#' # Write FASTA files for the 20 genomic target sequences
#' # of TwinStrand's Mouse Mutagenesis Panel.
#' rg <- get_seq("TSpanel_mouse")
#' write_reference_fasta(rg, output_path = NULL)
#' }
#' @importFrom Biostrings writeXStringSet
#' @export
write_reference_fasta <- function(regions_gr,
                                  output_path = NULL) {
  if (is.null(output_path)) {
    output_dir <- file.path(here::here(), "reference_output.fasta")
  } else {
    output_dir <- file.path(output_path, "reference_output.fasta")
  }
  if (!dir.exists(output_path)) {
    dir.create(output_path)
    output_dir <- file.path(output_path, "reference_output.fasta")
  }
  ref <- regions_gr$sequence
  # Convert to DNAStringSet
  ref <- Biostrings::DNAStringSet(ref)
  names(ref) <- make.unique(as.vector(seqnames(regions_gr)))
  Biostrings::writeXStringSet(ref,
                              output_dir,
                              append = FALSE,
                              format = "fasta")
}
