#' Join Sample Metadata
#' @description This function imports the sample metadata and joins it with the
#' mutation data.
#' @param mutation_data A data frame containing mutation data.
#' @param sample_data The path to the file containing the sample metadata.
#' Alternatively, a data frame can be provided directly.
#' @param sd_sep The separator used in the sample metadata file.
#' Default is tab (`\t`).
#' @return A data frame that combines the mutation data with the sample
#' metadata.
#' @importFrom dplyr left_join
import_sample_data <- function(mutation_data, sample_data, sd_sep = "\t") {
    if (is.data.frame(sample_data)) {
        sd <- sample_data
        if (nrow(sd) == 0) {
            stop("The sample data frame you've provided is empty")
        }
    } else if (is.character(sample_data)) {
        sample_file <- file.path(sample_data)
        if (!file.exists(sample_file)) {
            stop("The sample data file path you've specified is invalid")
        }
        if (file.info(sample_file)$size == 0) {
            stop("You are trying to import an empty sample data file")
        }
        sd <- read.delim(sample_file, sep = sd_sep, header = TRUE)
        if (ncol(sd) <= 1) {
            stop(
                "Your imported sample data only has one column. You may want",
                " to set sd_sep to properly reflect the delimiter used for",
                " the data you are importing."
            )
        }
    } else {
        stop("sample_data must be a character string or a data frame")
    }
    # Join
    joined_data <- dplyr::left_join(mutation_data, sd, suffix = c("", ".sd"))
    message("Sample metadata successfully joined to mutation data\n")
    return(joined_data)
}

#' Join Regions Metadata
#' @description This function imports the regions metadata and joins it with
#' the mutation data.
#' @param mutation_granges A data frame containing mutation data.
#' @param regions The path to the file containing the regions metadata.
#' Alternatively, a data frame can be provided directly.
#' @param rg_sep The separator used in the regions metadata file.
#' Default is tab (`\t`).
#' @param is_0_based_rg A logical value indicating whether the regions file is
#' 0-based (TRUE) or 1-based (FALSE). Default is FALSE.
#' @param padding An integer value indicating the number of base pairs to pad
#' the regions on either side. Default is 0.
#' @return A GRanges object that combines the mutation data with the regions
#' metadata.
#' @importFrom plyranges join_overlap_left_within_directed mutate
#' @importFrom BiocGenerics start end
import_regions_metadata <- function(mutation_granges,
                                    regions,
                                    rg_sep,
                                    is_0_based_rg,
                                    padding) {
    # load regions file as GRanges
    regions_gr <- MutSeqR::load_regions_file(regions, rg_sep, is_0_based_rg)
    regions_gr$in_regions <- TRUE

    # Apply padding
    BiocGenerics::start(regions_gr) <- pmax(
        BiocGenerics::start(regions_gr) - padding, 1
    )
    BiocGenerics::end(regions_gr) <- BiocGenerics::end(regions_gr) + padding

    # Join mutation data and region data using overlap
    mutation_granges <- plyranges::join_overlap_left_within_directed(
        mutation_granges,
        regions_gr,
        suffix = c("", "_regions")
    )
    message("Regions metadata successfully joined to mutation data\n")
    # Count the rows that did not overlap
    mutation_granges <- mutation_granges %>%
        plyranges::mutate(in_regions = ifelse(is.na(in_regions), FALSE, TRUE))
    false_count <- sum(mutation_granges$in_regions == FALSE)
    if (false_count > 0) {
        warning(
            false_count,
            " rows were outside of the specified regions.",
            " To remove these rows, use the filter_mut() function\n"
        )
    }
    return(mutation_granges)
}
#' Populate Sequence context
#' @description This function populates the trinucleotide context for each
#' mutation in the mutation data.
#' @param mutation_granges A GRanges object containing mutation data.
#' @param BS_genome The name of the Bioconductor BSgenome package to use for
#' retrieving the reference genome sequence.
#' @param n An integer value indicating the number of base pairs to include
#' on either side of the mutation for context. Default is 1 (trinucleotide
#' context including the mutation).
#' @return A GRanges object with an additional column for the trinucleotide
#' context.
#' @importFrom Biostrings getSeq
#' @importFrom BSgenome installed.genomes getBSgenome
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom Seqinfo seqnames
#' @importFrom BiocGenerics start end strand
populate_sequence_context <- function(mutation_granges, BS_genome, n = 1) {
    if (is.null(BS_genome)) {
        stop(
            "The trinuceotide context is populated from BS genomes.",
            " Please install the appropriate BS genome and indicate the",
            " pkgname with the BS_genome parameter. If you are not sure",
            " which BS genome to use, please provide the species and",
            " reference genome to find_BS_genome()."
        )
    }
    installed_BS_genomes <- BSgenome::installed.genomes()
    if (!(BS_genome %in% installed_BS_genomes)) {
        stop(
            "The specified BS genome is not installed. Please install the",
            " appropriate BS genome using BiocManager::install('pkgname')",
            " where pkgname is the name of the BSgenome package. If you are",
            " not sure which BS genome to use, please provide the species and",
            " reference genome to find_BS_genome()."
        )
    }
    message("Loading reference genome: ", BS_genome, ".")
    ref_genome <- BSgenome::getBSgenome(BS_genome)
    extract_context <- function(mut_gr, bsgenome) {
        # Resize the mutation_granges to include the context
        expanded_ranges <- GenomicRanges::GRanges(
        seqnames = Seqinfo::seqnames(mut_gr),
        ranges = IRanges::IRanges(
            start = BiocGenerics::start(mut_gr) - n,
            end = BiocGenerics::start(mut_gr) + n
        ),
        strand = BiocGenerics::strand(mut_gr)
        )
        # Extract the sequences from the BSgenome
        sequences <- Biostrings::getSeq(bsgenome, expanded_ranges)
        return(sequences)
    }
    message("Retrieving context sequences from BSgenome")
    context <- extract_context(mutation_granges, ref_genome)
    mutation_granges$context <- context
    return(mutation_granges)
}

#' Characterize Variants
#' @description This function generates additional columns for the mutation
#' data, including a breakdown of the mutation subtypes at various resolutions.
#' @param mutation_data A data frame containing mutation data.
#' @return A data frame with additional columns for variant characterization.
#' @importFrom dplyr mutate rename case_when
#' @importFrom stringr str_sub str_count
characterize_variants <- function(mutation_data) {
    # RSIDS
    if ("id" %in% colnames(mutation_data)) {
        mutation_data <- mutation_data %>% dplyr::mutate(
            is_known = ifelse(!.data$id == ".", "TRUE", "FALSE")
        )
    }
    # variation_type
    if ("variation_type" %in% colnames(mutation_data)) {
        mutation_data <- dplyr::rename(mutation_data,
        original_variation_type = "variation_type"
        )
    }
    mutation_data$variation_type <- mapply(
        MutSeqR::classify_variation,
        mutation_data$ref,
        mutation_data$alt
    )

    # Define substitution dictionary to normalize to pyrimidine context
    sub_dict <- c(
        "G>T" = "C>A", "G>A" = "C>T", "G>C" = "C>G",
        "A>G" = "T>C", "A>C" = "T>G", "A>T" = "T>A"
    )
    # Calculate columns:
    # nchar_ref, nchar_alt, varlen, short_ref, normalized_ref, subtype,
    # normalized_subtype, normalized_context, context_with_mutation,
    # normalized_context_with_mutation, gc_content
    mutation_data <- mutation_data %>%
        dplyr::mutate(
        nchar_ref = nchar(ref),
        nchar_alt = ifelse(!(.data$variation_type %in% c(
            "no_variant",
            "sv",
            "ambiguous",
            "uncategorized"
        )),
        nchar(alt), NA
        ),
        varlen =
            ifelse(.data$variation_type %in%
                c("insertion", "deletion", "complex"),
                    .data$nchar_alt - .data$nchar_ref,
                        ifelse(.data$variation_type %in% c("snv", "mnv"),
                            .data$nchar_ref, NA
                        )
            ),
        short_ref = substr(.data$ref, 1, 1),
        normalized_ref = dplyr::case_when(
            substr(.data$ref, 1, 1) == "A" ~ "T",
            substr(.data$ref, 1, 1) == "G" ~ "C",
            substr(.data$ref, 1, 1) == "C" ~ "C",
            substr(.data$ref, 1, 1) == "T" ~ "T"
        ),
        subtype =
            ifelse(.data$variation_type == "snv",
            paste0(.data$ref, ">", .data$alt),
            .data$variation_type
            ),
        normalized_subtype = ifelse(.data$subtype %in% names(sub_dict),
            sub_dict[.data$subtype],
            .data$subtype
        ),
        normalized_context = ifelse(
            stringr::str_sub(.data$context, 2, 2) %in% c("G", "A"),
            mapply(
            function(x) reverseComplement(x, case = "upper"),
            .data$context
            ),
            .data$context
        ),
        context_with_mutation =
            ifelse(.data$variation_type == "snv",
            paste0(
                stringr::str_sub(.data$context, 1, 1),
                "[", .data$subtype, "]",
                stringr::str_sub(.data$context, 3, 3)
            ),
            .data$variation_type
            ),
        normalized_context_with_mutation =
            ifelse(.data$variation_type == "snv",
            paste0(
                stringr::str_sub(.data$normalized_context, 1, 1),
                "[", .data$normalized_subtype, "]",
                stringr::str_sub(.data$normalized_context, 3, 3)
            ),
            .data$variation_type
            ),
        gc_content = (stringr::str_count(string = .data$context, pattern = "G") +
            stringr::str_count(string = .data$context, pattern = "C"))
        / stringr::str_count(.data$context),
        filter_mut = FALSE
        )
    return(mutation_data)
}

#' A utility function that will return the reference context of a mutation
#' @param mut_string the mutation. Ex. T>C, `A[G>T]C`
#' @return the reference context of the mutation
get_ref_of_mut <- function(mut_string) {
  a <- str_extract(mut_string, ".*(?=\\s*>)")
  # Remove non-letter characters
  b <- str_replace_all(a, "[^a-zA-Z]", "")
  # Extract letter characters after square bracket
  c <- str_extract(mut_string, "\\](.*)") %>% str_replace_all("[^a-zA-Z]", "")
  if (is.na(c)) {
    return(b)
  } else {
    return(paste0(b, c))
  }
}

#' Get the reverse complement of a DNA or RNA sequence.
#'
#' @param x A character vector of DNA or RNA sequences.
#' @param content c("dna", "rna") The type of sequence to be reversed.
#' @param case c("lower", "upper", "as is") The case of the output sequence.
#' @details This file is part of the source code for
#' SPGS: an R package for identifying statistical patterns in genomic sequences.
#' Copyright (C) 2015  Universidad de Chile and INRIA-Chile
#
#' A copy of Version 2 of the GNU Public License is available in the
#' share/licenses/gpl-2 file in the R installation directory or from
#' http://www.R-project.org/Licenses/GPL-2.
#' reverseComplement.R
#' @return A character vector of the reverse complement sequences.

reverseComplement <- function(x,
                              content = c("dna", "rna"),
                              case = c("lower", "upper", "as is")) {
  # reverse character vector
  strreverse <- function(x) {
    if (!is.character(x)) {
      stop("x must be a character vector")
    }
    vapply(
      strsplit(x, ""),
      function(y) paste(rev(y), collapse = ""),
      character(1)
    )
  }
  # Check arguments
  if (!is.character(x)) x <- as.character(x) # coerse x to a character vector
  content <- match.arg(content)
  case <- match.arg(case)
  if (length(x) == 0 || (length(x) == 1 && nchar(x) == 0)) {
    return(x)
  } # bail if input is empty
  if (case == "lower") x <- tolower(x)
  if (case == "upper") x <- toupper(x)
  if (content == "dna") {
    src <- "acgturykmswbdhvnxACGTURYKMSWBDHVNX-"
    dest <- "tgcaayrmkswvhdbnxTGCAAyRMKSWVHDBNX-"
  } else {
    src <- "acgturykmswbdhvnxACGTURYKMSWBDHVNX-"
    dest <- "ugcaayrmkswvhdbnxUGCAAyRMKSWVHDBNX-"
  } # if
  if (max(nchar(x)) > 1) {
    return(chartr(src, dest, strreverse(x)))
  }
  # x is not a single string, so process it as a vector
  chartr(src, dest, rev(x))
} # function

#' classify_variation
#' @description Classify the variation type of a mutation based on its ref and
#' alt values.
#' @param ref The reference allele.
#' @param alt The alternate allele.
#' @return A character indicating the type of variation.
#' @export
#' @examples
#' df <- data.frame(
#'   ref = c("A", "CAGT", "GCC", "T", "ACG", "C", "G", "T", "A"),
#'   alt = c("R", "TGA", "G", "TC", "TAC", "C", "<DEL>", "G", "???")
#' )
#' df$variation_type <- mapply(classify_variation, df$ref, df$alt)
#' df
classify_variation <- function(ref, alt) {
    stopifnot(is.character(ref), is.character(alt))
    no_variant_indicators <- c(".", "", "<NON_REF>")
    structural_indicators <- c(
        "<DEL>", "<INS>", "<DUP>", "<INV>", "<FUS>",
        "<CNV>", "<CNV:TR>", "<DUP:TANDEM>", "<DEL:ME>",
        "<INS:ME>"
    )
    iupac_indicators <- c(
        "R", "K", "S", "Y", "M",
        "W", "B", "H", "N", "D", "V"
    )

    # Case: No variant site
    # GVCF files sometimes list no_variant sites as <NON_REF> (GATK)
    alt <- gsub("(^|,)<NON_REF>(,|$)", "", alt)
    alt <- gsub("^,|,$", "", alt) # Trim leading/trailing commas
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
    # Case: Complex; ref and alt diff lengths & diff base compositions
    if (nchar(ref) != nchar(alt) &&
        !grepl(paste0("^", ref), alt) &&
        !grepl(paste0("^", alt), ref)
    ) {
        return("complex")
    }
    # Otherwise, uncategorized
    return("uncategorized")
}
