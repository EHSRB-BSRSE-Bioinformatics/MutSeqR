#' Write mutation_data to a VCF file
#'
#' Export your mutation_data to a VCF file for downstream applications.
#' @param mutation_data A data frame of a GRanges object containing your
#' mutation data. This can be the output of `import_mut_data`,
#' `import_vcf_data`, or `filter_mut.` Coordinates must be 1-based.
#' Required columns are "contig", "start", "end", "ref", "alt", "sample",
#' "alt_depth", "total_depth", and "ref_depth". Additional columns are allowed.
#' @param output_path The directory where the VCF file should be written.
#' Default is NULL, which will write the file to the current working directory.
#' @importFrom VariantAnnotation makeVRangesFromGRanges writeVcf asVCF
#' @returns Writes a VCF file of mutations "mutation_output.vcf".
#' @examples
#' \dontrun{
#' example_file <- system.file("extdata", "Example_files",
#'                             "example_mutation_data_filtered.rds",
#'                             package = "MutSeqR")
#' example_data <- readRDS(example_file)
#' write_vcf_from_mut(example_data)
#' }
#' @export

### Would be nice if some of the extra columns were put into an INFO field
### This would require a bit more work to get the data into the right format
### And I don't want to do it right now.
write_vcf_from_mut <- function(mutation_data,
                               output_path = NULL) {
  if (is.null(output_path)) {
    output_dir <- file.path(here::here(), "mutation_output.vcf")
  } else {
    output_dir <- file.path(output_path, "mutation_output.vcf")
  }
  if (!dir.exists(output_path)) {
    dir.create(output_path)
    output_dir <- file.path(output_path, "mutation_output.vcf")
  }

  # Convert to GRanges.
  if (inherits(mutation_data, "data.frame")) {
    mutation_data <- makeGRangesFromDataFrame(
      df = mutation_data,
      keep.extra.columns = TRUE,
      seqnames.field = "contig",
      start.field = "start",
      end.field = "end"
    )
  }
  if ("total_depth" %in% colnames(mutation_data)) {
    muts_for_vcf <-
    VariantAnnotation::makeVRangesFromGRanges(mutation_data,
                                              sampleNames.field = "sample",
                                              ref.field = "ref",
                                              alt.field = "alt",
                                              altDepth.field = "alt_depth",
                                              totalDepth.field = "total_depth",
                                              refDepth.field = "ref_depth",
                                              keep.extra.columns = TRUE)
  } else {
    muts_for_vcf <-
      VariantAnnotation::makeVRangesFromGRanges(mutation_data,
                                                sampleNames.field = "sample",
                                                ref.field = "ref",
                                                alt.field = "alt",
                                                altDepth.field = "alt_depth",
                                                keep.extra.columns = TRUE)
  }
  VariantAnnotation::writeVcf(VariantAnnotation::asVCF(muts_for_vcf),
                              filename = output_dir)
}