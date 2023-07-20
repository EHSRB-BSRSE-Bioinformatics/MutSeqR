#' Analyse transition/transversion ratios
#' 
#' Function to convert mutation data to a format that can be passed to GenVisR
#' @param mutations Mutation data, either as GRanges or data frame
#' @param group Experimental group (a column in the mutation data) by which
#' to aggregate samples
#' @param y Default expected proportions for each mutation subtype. See 
#' documentation for `GenVisR::TvTi`.
#' @param ... Additional arguments sent to GenVisR::TvTi (e.g., out = "data",
#' sample_order_input, sort, type = "Proportion" or "Frequency")
#' @importFrom GenVisR TvTi
#' @importFrom dplyr mutate select filter
#' @export
tvti_plot <- function(mutations = mutation_data,
                      group = "sample",
                      y = c(`A->C or T->G (TV)` = 0.066,
                            `A->G or T->C (TI)` = 0.217,
                            `A->T or T->A (TV)` = 0.065,
                            `G->A or C->T (TI)` = 0.4945,
                            `G->C or C->G (TV)` = 0.0645,
                            `G->T or C->A (TV)` = 0.093),
                      ...) {
  # If mutation data is GRanges, convert to data frame
  if (inherits(mutations, "GRanges")) { mutations <- as.data.frame(mutations) }
  snvs_genvisr <- mutations %>%
    dplyr::mutate(
      sample = .data[[group]],
      reference = ref,
      variant = alt) %>%
    dplyr::filter(variation_type == "snv") %>%
    dplyr::select(sample, reference, variant)
  GenVisR::TvTi(snvs_genvisr, fileType = "MGI",  progress = FALSE, ...)
 
  #tv_ti_ratio <- tv_ti_table$main %>% dplyr::mutate(Class = str_extract(string=trans_tranv, pattern=regex("\\(\\w++\\)")))
  #tv_ti_ratio %>% dplyr::group_by(Class) %>% dplyr::mutate(class_count = sum(Freq))
}
