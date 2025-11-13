#' Dictionary of column name synonyms for required mutation data columns
#'
#' A list of column specifications
#'
#' @format A list with potential variable column names
#' @examples
#' # examples of "alt" synonyms
#' op$column$alternate
#' op$column$alt_value
#' op$column$tumor_seq_allele2
#' @export
op <- list()
op$column <- list()
op$column$alternate <- "alt"
op$column$alt_value <- "alt"
op$column$tumor_seq_allele2 <- "alt" # CODEC MAF
op$column$alt_read_depth <- "alt_depth"
op$column$n_alt_count <- "alt_depth" # CODEC MAF
op$column$t_alt_count <- "alt_depth" # CODEC MAF
op$column$var_depth <- "alt_depth"
op$column$variant_depth <- "alt_depth"
op$column$vd <- "alt_depth"
op$column$sequence_context <- "context"
op$column$trinucleotide_context <- "context"
op$column$flanking_sequence <- "context"
op$column$chr <- "contig"
op$column$chromosome <- "contig"
op$column$seqnames <- "contig"
op$column$dp <- "depth"
op$column$end_position <- "end"
op$column$n_calls <- "no_calls"
op$column$n_depth <- "no_calls"
op$column$no_depth <- "no_calls"
op$column$reference <- "ref"
op$column$reference_allele <- "ref"
op$column$ref_value <- "ref"
op$column$sample_name <- "sample"
op$column$sample_id <- "sample"
op$column$pos <- "start"
op$column$position <- "start"
op$column$start_position <- "start"
op$column$total_depth <- "total_depth"
op$column$informative_somatic_depth <- "total_depth"
op$site$columns <- c("contig", "start")
op$mut_count_method <- "min"
op$processed_required_mut_cols <-
    c(
        "alt_depth",
        "variation_type",
        "subtype",
        "normalized_subtype",
        "context_with_mutation",
        "normalized_context_with_mutation",
        "normalized_ref",
        "short_ref",
        "context",
        "normalized_context",
        "filter_mut"
    )
op$base_required_mut_cols <-
    c(
        "contig",
        "start",
        "end",
        "sample",
        "ref",
        "alt"
    )
op$default_vaf_cutoffs <- c(0.3, 0.7, 0.9)

#' Subtype Resolution Dictionary
#'
#' Associates the subtype resolution names with their corresponding column
#' names in the mutation data.
#'
#' @format A vector with corresponding values
#' @examples
#' subtype_dict["base_96"]
#' subtype_dict["type"]
#' @export
subtype_dict <- c(
    "none" = NA,
    "type" = "variation_type",
    "base_6" = "normalized_subtype",
    "base_12" = "subtype",
    "base_96" = "normalized_context_with_mutation",
    "base_192" = "context_with_mutation"
)

#' A comprehensive list of mutation subtypes at different resolutions
#' @format A list with corresponding values
#' @examples
#' subtype_list[["type"]]
#' subtype_list[["base_6"]]
#' @export
subtype_list <- list(
    type = c(
        "ambiguous", "complex", "deletion", "insertion",
        "mnv", "no_variant", "sv", "snv", "uncategorized"
    ),
    base_6 = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"),
    base_12 = c(
        "A>C", "A>G", "A>T", "C>A", "C>G", "C>T",
        "G>A", "G>C", "G>T", "T>A", "T>C", "T>G"
    ),
    base_96 = c(
        "A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "A[C>G]A", "A[C>G]C",
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
        "T[T>C]G", "T[T>C]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T"
    ),
    base_192 = c(
        "A[A>C]A", "A[A>C]C", "A[A>C]G", "A[A>C]T", "A[A>G]A", "A[A>G]C",
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
        "T[T>C]G", "T[T>C]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T"
    )
)


#' Values used for denominators in frequency calculations
#'
#' These values are used to cross reference base substitution types to their
#' appropriate denominators for calculations. That is", "for example, the 6 base
#' substitution frequency should be subsetted based on the normalized_ref
#' column which would contain only T or C (i.e., the pyrimidine context for
#' base substitutions).
#' @format A vector with corresponding values
#' @examples
#' denominator_dict["base_96"]
#' denominator_dict["type"]
#' @export
denominator_dict <- c(
    "none" = NA,
    "type" = NA,
    "base_6" = "normalized_ref",
    "base_12" = "short_ref",
    "base_96" = "normalized_context",
    "base_192" = "context"
)
#' A list of reference contexts at different resolutions
#' @format A list with corresponding values
#' @examples
#' context_list[["base_6"]]
#' context_list[["base_12"]]
#' @export
context_list <- list(
    none = NA,
    type = NA,
    base_6 = c("C", "T"),
    base_12 = c("A", "C", "G", "T"),
    base_96 = c(
        "ACA", "ACC", "ACG", "ACT", "ATA", "ATC", "ATG", "ATT",
        "CCA", "CCC", "CCG", "CCT", "CTA", "CTC", "CTG", "CTT",
        "GCA", "GCC", "GCG", "GCT", "GTA", "GTC", "GTG", "GTT",
        "TCA", "TCC", "TCG", "TCT", "TTA", "TTC", "TTG", "TTT"
    ),
    base_192 = c(
        "AAA", "AAC", "AAG", "AAT",
        "ACA", "ACC", "ACG", "ACT",
        "AGA", "AGC", "AGG", "AGT",
        "ATA", "ATC", "ATG", "ATT",
        "CAA", "CAC", "CAG", "CAT",
        "CCA", "CCC", "CCG", "CCT",
        "CGA", "CGC", "CGG", "CGT",
        "CTA", "CTC", "CTG", "CTT",
        "GAA", "GAC", "GAG", "GAT",
        "GCA", "GCC", "GCG", "GCT",
        "GGA", "GGC", "GGG", "GGT",
        "GTA", "GTC", "GTG", "GTT",
        "TAA", "TAC", "TAG", "TAT",
        "TCA", "TCC", "TCG", "TCT",
        "TGA", "TGC", "TGG", "TGT",
        "TTA", "TTC", "TTG", "TTT"
    )
)

#' BS genome organism dictionary
#' @description A dictionary to map common organism names to their
#' corresponding BSgenome organism names.
#' @format A list with organism names as keys and corresponding BSgenome
#' organism names as values.
#' @export
#' @examples
#' BS_org_map["Hsapiens"]
BS_org_map <- list(
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
    "Gaculeatus" = c(
        "gasterosteus aculeatus",
        "stickleback",
        "three-spined stickleback"
    ),
    "Ggallus" = c("gallus gallus", "chicken"),
    "Gmax" = c("glycine max", "soybean"),
    "Hsapiens" = c("homo sapiens", "homo sapiens sapiens", "human"),
    "Mdomestica" = c(
        "monodelphis domestica",
        "opossum",
        "gray short-tailed opossum"
    ),
    "Mfascicularis" = c(
        "macaca fascicularis",
        "long-tailed macaque",
        "crab-eating macaque"
    ),
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
