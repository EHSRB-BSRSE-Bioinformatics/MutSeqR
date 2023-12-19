#' Write FASTA file of reference sequences
#'
#' In some cases, you might want to generate an arbitrary multi-sequence FASTA
#' file from GRanges including the reference sequences. This function allows you
#' to do this.
#' @param mutation_data A GRanges object imported from a .mut file
#' @param vcf_out The path for the VCF file output.
#' @importFrom VariantAnnotation makeVRangesFromGRanges writeVcf asVCF
#' @returns Writes a VCF file of mutations to be used in downstream applications.
#' @export
write_VCF_from_mut <- function(mutation_data,
                               vcf_out = "./mutation_output.vcf") {
 
   # Check if data is provided as data frame: if so, convert to GRanges.
  if (inherits(mutation_data, "data.frame")) { 
   # Turn into GRanges
      mutation_data <- makeGRangesFromDataFrame(
        df = mutation_data,
        keep.extra.columns = T,
        seqnames.field = "contig",
        start.field = "start",
        end.field = "end",
        starts.in.df.are.0based = FALSE
      ) }
  
  muts_for_vcf <-
    VariantAnnotation::makeVRangesFromGRanges(mutation_data,
                                              sampleNames.field = "sample",
                                              ref.field = "ref",
                                              alt.field = "alt",
                                              altDepth.field = "alt_depth",
                                              totalDepth.field = "total_depth",
                                              refDepth.field = "ref_depth",
                                              keep.extra.columns = T)
  VariantAnnotation::writeVcf(VariantAnnotation::asVCF(muts_for_vcf),
                              filename = vcf_out)
}