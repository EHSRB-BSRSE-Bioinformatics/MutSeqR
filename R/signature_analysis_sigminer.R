#' Run COSMIC signatures comparison
#'
#' After cleaning the mutation data input, runs several Alexandrov Lab tools for COSMIC signature analysis (assigns signatures to best explain the input data).
#' @param mutation_data A data frame imported from a .mut file
#' @param project_name The name of the project; used to get mutation data into the required .txt format for SigProfiler
#' @param project_genome A string describing the reference genome to use; e.g., GRCh38
#' @param group The column in the mutation data used to aggregate groups (e.g., sample ID, tissue, dose)
#' @param run_bootstrapping TRUE or FALSE. Default FALSE. Determines if the
#' sig_fit_bootstrap_batch() function should be run. This *should* be done, but
#' the process is slow, so it's best to confirm that the rest of the analysis is
#' working as expected first.
#' @param ...
#' @returns Creates a subfolder in the output directory with SigProfiler tools results.
#'  Suggests: sigminer
#' @importFrom here here
#' @importFrom dplyr select rename mutate
#' @importFrom utils write.table
#' @import reticulate
#' @import stringr
#' @export
signature_analysis_sigminer <- function(mutations = mutation_data,
                                        project_name = "Default",
                                        project_genome = "BSgenome.Mmusculus.UCSC.hg38",
                                        group = "sample",
                                        run_bootstrapping = F,
                                        ...) {
  
  sigminer_input <- as.data.frame(mutations) |>
    dplyr::filter(!variation_type %in% "no_variant") |>
    dplyr::select(all_of(group), variation_type, seqnames,
                  start, end, ref, alt) |>
    dplyr::rename( # TODO add a column for "gene", but use locus with TS data.. should be Hugo_Symbol for MAF compatibility
      "Tumor_Sample_Barcode" = group,
      "Chromosome" = "seqnames",
      "Start_Position" = "start",
      "End_Position" = "end",
      "Reference_Allele" = "ref",
      "Tumor_Seq_Allele2" = "alt"
      )
  
  output_path <- file.path(
    here::here(),
    "output",
    "SigMiner",
    group
  )
  
  if (!dir.exists(output_path)) {
    dir.create(file.path(output_path, "matrices"), recursive = T)
  }
  
  write.table(sigminer_input,
              file = file.path(output_path, "matrices", "mutations.sigminer.txt"),
              sep = "\t", row.names = F, quote = F
  )

  
  sigminer_results <- list()
  
  sigminer_results$signature_data <-
    sigminer::read_maf_minimal(sigminer_input)
  
  sigminer_results$tally_results <-
    sigminer::sig_tally(sigminer_results$signature_data,
                        ref_genome = project_genome,
                        ...)
  # add get_sig_similarity() to compare catalogs to signatures; TODO
  
  sigminer_results$fitting_results <-
    sigminer::sig_fit(t(sigminer_results$tally_results$nmf_matrix),
                      sig_db = "SBS", sig_index = "ALL")
  
  if (run_bootstrapping == T) {
  sigminer_results$bootstrap_results <-
    sigminer::sig_fit_bootstrap_batch(t(sigminer_results$tally_results$nmf_matrix),
                                      sig_db = "SBS", sig_index = "ALL")
  }
  return(sigminer_results)
  
}
