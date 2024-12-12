#' Column names for mut tables
#'
#' A list of column specifications
#'
#' @format A list with potential variable column names
#' @export
op <- list()
op$column <- list()
op$column$alternate <- "alt"
op$column$alt_value <- "alt"
op$column$alt_depth <- "alt_depth"
op$column$alt_read_depth <- "alt_depth"
op$column$var_depth <- "alt_depth"
op$column$variant_depth <- "alt_depth"
op$column$vd <- "alt_depth"
op$column$sequence_context <- "context"
op$column$trinucleotide_context <- "context"
op$column$flanking_sequence <- "context"
op$column$contig <- "contig"
op$column$chr <- "contig"
op$column$chromosome <- "contig"
op$column$seqnames <- "contig"
op$column$depth <- "depth"
op$column$dp <- "depth"
op$column$end <- "end"
op$column$end_position <- "end"
op$column$is_snp <- "is_snp"
op$column$lower_ci <- "lower_ci"
op$column$mut_depth <- "mut_depth"
op$column$mut_freq <- "mut_freq"
op$column$n_calls <- "no_calls"
op$column$n_depth <- "no_calls"
op$column$no_depth <- "no_calls"
op$column$reference <- "ref"
op$column$ref_value <- "ref"
op$column$sample <- "sample"
op$column$sample_name <- "sample"
op$column$sample_id <- "sample"
op$column$start <- "start"
op$column$pos <- "start"
op$column$position <- "start"
op$column$mutation_subtype <- "subtype"
op$column$total_depth <- "total_depth"
op$column$informative_somatic_depth <- "total_depth"
op$column$upper_ci <- "upper_ci"
op$column$vaf <- "vaf"
op$column$type <- "variation_type"
op$column$mutation_type <- "variation_type"
op$column$variant_type <- "variation_type"
op$site$columns <- c("contig", "start")
op$mut_count_method <- "min"
op$processed_required_mut_cols <-
  c("alt_depth",
    "variation_type",
    "subtype",
    "normalized_subtype",
    "context_with_mutation",
    "normalized_context_with_mutation",
    "normalized_ref",
    "short_ref",
    "context",
    "normalized_context")
op$base_required_mut_cols <-
  c("contig",
    "start",
    "end",
    "sample",
    "ref",
    "alt")
op$default_vaf_cutoffs <- c(0.3, 0.7, 0.9)

#' Values accepted for mutation subtypes
#'
#' These values are used to enable user input to translate to columns in a
#' mut file
#'
#' @format A vector with corresponding values
#' @export
subtype_dict <- c(
  "none" = NA,
  "type" = "variation_type",
  "base_6" = "normalized_subtype",
  "base_12" = "subtype",
  "base_96" = "normalized_context_with_mutation",
  "base_192" = "context_with_mutation"
)

#' A list of mutation subtypes at different resolutions
#' @format A list with corresponding values
#' @export
subtype_list <- list(
  type = c("no_variant", "snv", "deletion", "insertion", "complex", "mnv",
           "sv", "ambiguous", "uncategorized"),
  base_6 = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
  base_12 = c("A>C", "A>G", "A>T", "C>A", "C>G", "C>T",
              "G>A", "G>C", "G>T", "T>A", "T>C", "T>G"),
  base_96 = c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "A[C>G]A", "A[C>G]C",
              "A[C>G]G", "A[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T",
              "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "A[T>C]A", "A[T>C]C",
              "A[T>C]G", "A[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T",
              "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "C[C>G]A", "C[C>G]C",
              "C[C>G]G", "C[C>G]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T",
              "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "C[T>C]A", "C[T>C]C",
              "C[T>C]G", "C[T>C]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T",
              "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "G[C>G]A", "G[C>G]C",
              "G[C>G]G", "G[C>G]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T",
              "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "G[T>C]A", "G[T>C]C",
              "G[T>C]G", "G[T>C]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T",
              "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "T[C>G]A", "T[C>G]C",
              "T[C>G]G", "T[C>G]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
              "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "T[T>C]A", "T[T>C]C",
              "T[T>C]G", "T[T>C]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T"),
  base_192 = c("A[A>C]A", "A[A>C]C", "A[A>C]G", "A[A>C]T", "A[A>G]A", "A[A>G]C",
               "A[A>G]G", "A[A>G]T", "A[A>T]A", "A[A>T]C", "A[A>T]G", "A[A>T]T", 
               "A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "A[C>G]A", "A[C>G]C",
               "A[C>G]G", "A[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T",
               "A[G>A]A", "A[G>A]C", "A[G>A]G", "A[G>A]T", "A[G>C]A", "A[G>C]C",
               "A[G>C]G", "A[G>C]T", "A[G>T]A", "A[G>T]C", "A[G>T]G", "A[G>T]T",
               "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "A[T>C]A", "A[T>C]C",
               "A[T>C]G", "A[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T",
               "C[A>C]A", "C[A>C]C", "C[A>C]G", "C[A>C]T", "C[A>G]A", "C[A>G]C",
               "C[A>G]G", "C[A>G]T", "C[A>T]A", "C[A>T]C", "C[A>T]G", "C[A>T]T",
               "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "C[C>G]A", "C[C>G]C",
               "C[C>G]G", "C[C>G]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T",
               "C[G>A]A", "C[G>A]C", "C[G>A]G", "C[G>A]T", "C[G>C]A", "C[G>C]C",
               "C[G>C]G", "C[G>C]T", "C[G>T]A", "C[G>T]C", "C[G>T]G", "C[G>T]T",
               "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "C[T>C]A", "C[T>C]C",
               "C[T>C]G", "C[T>C]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T",
               "G[A>C]A", "G[A>C]C", "G[A>C]G", "G[A>C]T", "G[A>G]A", "G[A>G]C",
               "G[A>G]G", "G[A>G]T", "G[A>T]A", "G[A>T]C", "G[A>T]G", "G[A>T]T",
               "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "G[C>G]A", "G[C>G]C",
               "G[C>G]G", "G[C>G]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T",
               "G[G>A]A", "G[G>A]C", "G[G>A]G", "G[G>A]T", "G[G>C]A", "G[G>C]C",
               "G[G>C]G", "G[G>C]T", "G[G>T]A", "G[G>T]C", "G[G>T]G", "G[G>T]T",
               "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "G[T>C]A", "G[T>C]C",
               "G[T>C]G", "G[T>C]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T",
               "T[A>C]A", "T[A>C]C", "T[A>C]G", "T[A>C]T", "T[A>G]A", "T[A>G]C",
               "T[A>G]G", "T[A>G]T", "T[A>T]A", "T[A>T]C", "T[A>T]G", "T[A>T]T",
               "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "T[C>G]A", "T[C>G]C",
               "T[C>G]G", "T[C>G]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
               "T[G>A]A", "T[G>A]C", "T[G>A]G", "T[G>A]T", "T[G>C]A", "T[G>C]C",
               "T[G>C]G", "T[G>C]T", "T[G>T]A", "T[G>T]C", "T[G>T]G", "T[G>T]T",
               "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "T[T>C]A", "T[T>C]C",
               "T[T>C]G", "T[T>C]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T")
)


#' Values used for denominators in frequency calculations
#'
#' These values are used to cross reference base substitution types to their
#' appropriate denominators for calculations. That is", "for example, the 6 base
#' substitution frequency should be subsetted based on the normalized_ref 
#' column which would contain only T or C (i.e., the pyrimidine context for
#' base substitutions).
#' @format A vector with corresponding values
#' @export
denominator_dict <- c(
  "none" = NA,
  "type" = NA,
  "base_6" = "normalized_ref",
  "base_12" = "short_ref",
  "base_96" = "normalized_context",
  "base_192" = "context"
)