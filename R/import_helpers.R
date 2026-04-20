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
import_sample_data <- function(sample_data, sd_sep = "\t") {
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
        sd <- read.delim(
            sample_file,
            sep = sd_sep,
            header = TRUE,
            check.names = FALSE
        )
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

    # DEFENSIVE CHECK: Ensure "sample" column exists in the metadata
    if (!"sample" %in% colnames(sd)) {
        if (any(tolower(colnames(sd)) == "sample")) {
            stop(
                "Error merging sample metadata: A column exactly named 'sample' is required, ",
                "but found a variation with different casing. Column names are case-sensitive."
            )
        } else {
            available_cols <- paste(head(colnames(sd), 10), collapse = ", ")
            stop(
                "Error merging sample metadata: The required column 'sample' was not found.\n",
                "Please ensure your sample metadata contains a column named 'sample'.\n",
                "Columns found: ",
                available_cols
            )
        }
    }

    return(sd)
}


#' Join imported mutation data with sample metadata
#' @param dat A data frame containing imported mutation data.
#' @param sample_df A data frame containing sample metadata.
#' @param source_label A short label used in user-facing error messages.
#' @param mismatch_hint Extra guidance appended to sample mismatch errors.
#' @param remove_sample_suffix Optional regex suffix to strip from sample names.
#' @importFrom dplyr left_join
join_import_sample_metadata <- function(
    dat,
    sample_df = NULL,
    source_label = "mutation data",
    mismatch_hint = "Please check for trailing suffixes or typos in your metadata file.",
    remove_sample_suffix = NULL
) {
    if (is.null(sample_df)) {
        return(dat)
    }

    if (!"sample" %in% colnames(dat)) {
        stop(
            "Error in mutation data: 'sample' column is missing prior to joining sample metadata."
        )
    }

    if (is.list(dat$sample)) {
        dat$sample <- vapply(
            dat$sample,
            function(x) paste(x, collapse = ","),
            character(1)
        )
    }

    dat$sample <- as.character(dat$sample)

    if (!is.null(remove_sample_suffix)) {
        dat$sample <- gsub(
            pattern = remove_sample_suffix,
            replacement = "",
            x = dat$sample
        )
    }

    sample_df$sample <- as.character(sample_df$sample)

    imported_samples <- unique(dat$sample)
    meta_samples <- unique(sample_df$sample)
    missing_in_meta <- setdiff(imported_samples, meta_samples)

    if (length(missing_in_meta) > 0) {
        stop(
            "Mismatch in sample names: Some samples in your ",
            source_label,
            " are MISSING from the metadata.\n",
            "Sample names must match EXACTLY. ",
            mismatch_hint,
            "\n\n",
            "Unmatched samples in ",
            source_label,
            ": ",
            paste(utils::head(missing_in_meta, 3), collapse = ", "),
            "\n",
            "Available samples in metadata: ",
            paste(utils::head(meta_samples, 3), collapse = ", "),
            "\n",
            call. = FALSE
        )
    }

    dat <- dplyr::left_join(
        dat,
        sample_df,
        by = "sample",
        suffix = c("", ".sd")
    )

    message("Sample metadata successfully joined to mutation data\n")
    dat
}


#' Standardize and validate imported mutation data
#' @param dat A data frame containing imported mutation data.
#' @param custom_column_names Optional custom column name mapping overrides.
#' @param sample_df Optional sample metadata as a data frame.
#' @param source_label A short label used in user-facing error messages.
#' @param mismatch_hint Extra guidance appended to sample mismatch errors.
#' @param remove_sample_suffix Optional regex suffix to strip from sample names.
#' @param required_columns Required columns that must be present.
#' @param allow_na_columns Required columns allowed to contain NA values.
#' @param BS_genome BSgenome package name, used when context needs populating.
prepare_imported_mutation_data <- function(
    dat,
    custom_column_names = NULL,
    sample_df = NULL,
    source_label = "mutation data",
    mismatch_hint = "Please check for trailing suffixes or typos in your metadata file.",
    remove_sample_suffix = NULL,
    required_columns = MutSeqR::op$base_required_mut_cols,
    allow_na_columns = character(0),
    BS_genome = NULL
) {
    if (!is.null(custom_column_names)) {
        cols <- modifyList(MutSeqR::op$column, custom_column_names)
        dat <- rename_columns(dat, cols)
    } else {
        dat <- rename_columns(dat)
    }

    dat <- join_import_sample_metadata(
        dat = dat,
        sample_df = sample_df,
        source_label = source_label,
        mismatch_hint = mismatch_hint,
        remove_sample_suffix = remove_sample_suffix
    )

    dat <- check_required_columns(dat, required_columns)
    context_exists <- "context" %in% colnames(dat)

    required_no_na <- setdiff(required_columns, allow_na_columns)
    na_columns_required <- required_no_na[
        vapply(dat[required_no_na], function(x) any(is.na(x)), logical(1))
    ]

    if (length(na_columns_required) > 0) {
        stop(
            "NA values were found within the following required column(s): ",
            paste(na_columns_required, collapse = ", "),
            ". Please confirm that your data is complete before proceeding."
        )
    }

    if (context_exists && any(is.na(dat$context))) {
        context_exists <- FALSE
    }

    if (!context_exists) {
        validate_BS_genome(BS_genome)
    }

    list(dat = dat, context_exists = context_exists)
}


#' Convert imported mutation data to GRanges and enrich it
#' @param dat A data frame containing imported mutation data.
#' @param context_exists Whether an existing context column can be trusted.
#' @param regions Optional regions input.
#' @param rg_sep Region separator.
#' @param is_0_based_rg Whether region coordinates are 0-based.
#' @param padding Region padding.
#' @param BS_genome BSgenome package name used to populate sequence context.
#' @param starts_in_df_are_0based Whether imported mutation starts are 0-based.
build_imported_mutation_ranges <- function(
    dat,
    context_exists,
    regions = NULL,
    rg_sep = "\t",
    is_0_based_rg = FALSE,
    padding = 0,
    BS_genome = NULL,
    starts_in_df_are_0based = FALSE
) {
    mut_ranges <- GenomicRanges::makeGRangesFromDataFrame(
        df = as.data.frame(dat),
        keep.extra.columns = TRUE,
        seqnames.field = "contig",
        start.field = "start",
        end.field = "end",
        starts.in.df.are.0based = starts_in_df_are_0based
    )

    if (!is.null(regions)) {
        mut_ranges <- import_regions_metadata(
            mutation_granges = mut_ranges,
            regions = regions,
            rg_sep = rg_sep,
            is_0_based_rg = is_0_based_rg,
            padding = padding
        )
    }

    if (!context_exists) {
        mut_ranges <- populate_sequence_context(
            mutation_granges = mut_ranges,
            BS_genome = BS_genome
        )
    }

    mut_ranges
}


#' Mark duplicated genomic positions within samples
#' @param dat A data frame containing mutation data.
#' @importFrom dplyr group_by mutate ungroup
mark_duplicate_mutation_rows <- function(dat) {
    dat <- dat %>%
        dplyr::group_by(.data$sample, .data$contig, .data$start) %>%
        dplyr::mutate(row_has_duplicate = dplyr::n() > 1) %>%
        dplyr::ungroup()

    if (sum(dat$row_has_duplicate) > 0) {
        warning(
            sum(dat$row_has_duplicate),
            " rows were found whose position was the same as that of at least one other row for the same sample."
        )

        if ("total_depth" %in% colnames(dat)) {
            warning(
                "The total_depth may be double-counted in some instances due to overlapping positions. Set the correct_depth parameter in calculate_mf() to correct the total_depth for these instances."
            )
        }
    }

    dat
}


#' Add VAF and reference depth columns when total depth is available
#' @param dat A data frame containing mutation data.
#' @importFrom dplyr mutate
add_vaf_columns <- function(dat) {
    if ("total_depth" %in% colnames(dat)) {
        dat <- dat %>%
            dplyr::mutate(
                vaf = .data$alt_depth / .data$total_depth,
                ref_depth = .data$total_depth - .data$alt_depth
            )
    }

    dat
}


#' Finalize imported mutation data after source-specific parsing
#' @param mut_ranges A GRanges object containing imported mutation data.
#' @param depth_resolver A function that adds or adjusts depth-related columns.
#' @param output_granges Whether to return a GRanges object.
finalize_imported_mutation_data <- function(
    mut_ranges,
    depth_resolver = identity,
    output_granges = FALSE
) {
    dat <- as.data.frame(mut_ranges) %>%
        dplyr::rename(contig = "seqnames")
    dat <- characterize_variants(dat)

    if (!"alt_depth" %in% colnames(dat)) {
        dat$alt_depth <- 1
    }

    dat <- depth_resolver(dat)
    dat <- mark_duplicate_mutation_rows(dat)
    dat <- add_vaf_columns(dat)

    if (output_granges) {
        return(GenomicRanges::makeGRangesFromDataFrame(
            df = dat,
            keep.extra.columns = TRUE,
            seqnames.field = "contig",
            start.field = "start",
            end.field = "end",
            starts.in.df.are.0based = FALSE
        ))
    }

    dat
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
#' @importFrom plyranges join_overlap_left_within_directed
#' @importFrom BiocGenerics start end
#' @importFrom S4Vectors mcols
import_regions_metadata <- function(
    mutation_granges,
    regions,
    rg_sep,
    is_0_based_rg,
    padding
) {
    # load regions file as GRanges
    regions_gr <- MutSeqR::load_regions_file(regions, rg_sep, is_0_based_rg)
    regions_gr$in_regions <- TRUE

    # Apply padding
    BiocGenerics::start(regions_gr) <- pmax(
        BiocGenerics::start(regions_gr) - padding,
        1
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
    S4Vectors::mcols(mutation_granges)$in_regions[is.na(
        S4Vectors::mcols(mutation_granges)$in_regions
    )] <- FALSE
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
        mutation_data <- mutation_data %>%
            dplyr::mutate(
                is_known = ifelse(!.data$id == ".", "TRUE", "FALSE")
            )
    }
    # variation_type
    if ("variation_type" %in% colnames(mutation_data)) {
        mutation_data <- dplyr::rename(
            mutation_data,
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
        "G>T" = "C>A",
        "G>A" = "C>T",
        "G>C" = "C>G",
        "A>G" = "T>C",
        "A>C" = "T>G",
        "A>T" = "T>A"
    )
    # Calculate columns:
    # nchar_ref, nchar_alt, varlen, short_ref, normalized_ref, subtype,
    # normalized_subtype, normalized_context, context_with_mutation,
    # normalized_context_with_mutation, gc_content
    mutation_data <- mutation_data %>%
        dplyr::mutate(
            nchar_ref = nchar(ref),
            nchar_alt = ifelse(
                !(.data$variation_type %in%
                    c(
                        "no_variant",
                        "sv",
                        "ambiguous",
                        "uncategorized"
                    )),
                nchar(alt),
                NA
            ),
            varlen = ifelse(
                .data$variation_type %in%
                    c("insertion", "deletion", "complex"),
                .data$nchar_alt - .data$nchar_ref,
                ifelse(
                    .data$variation_type %in% c("snv", "mnv"),
                    .data$nchar_ref,
                    NA
                )
            ),
            short_ref = substr(.data$ref, 1, 1),
            normalized_ref = dplyr::case_when(
                substr(.data$ref, 1, 1) == "A" ~ "T",
                substr(.data$ref, 1, 1) == "G" ~ "C",
                substr(.data$ref, 1, 1) == "C" ~ "C",
                substr(.data$ref, 1, 1) == "T" ~ "T"
            ),
            subtype = ifelse(
                .data$variation_type == "snv",
                paste0(.data$ref, ">", .data$alt),
                .data$variation_type
            ),
            normalized_subtype = ifelse(
                .data$subtype %in% names(sub_dict),
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
            context_with_mutation = ifelse(
                .data$variation_type == "snv",
                paste0(
                    stringr::str_sub(.data$context, 1, 1),
                    "[",
                    .data$subtype,
                    "]",
                    stringr::str_sub(.data$context, 3, 3)
                ),
                .data$variation_type
            ),
            normalized_context_with_mutation = ifelse(
                .data$variation_type == "snv",
                paste0(
                    stringr::str_sub(.data$normalized_context, 1, 1),
                    "[",
                    .data$normalized_subtype,
                    "]",
                    stringr::str_sub(.data$normalized_context, 3, 3)
                ),
                .data$variation_type
            ),
            gc_content = (stringr::str_count(
                string = .data$context,
                pattern = "G"
            ) +
                stringr::str_count(string = .data$context, pattern = "C")) /
                stringr::str_count(.data$context),
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

reverseComplement <- function(
    x,
    content = c("dna", "rna"),
    case = c("lower", "upper", "as is")
) {
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
    if (!is.character(x)) {
        x <- as.character(x)
    } # coerse x to a character vector
    content <- match.arg(content)
    case <- match.arg(case)
    if (length(x) == 0 || (length(x) == 1 && nchar(x) == 0)) {
        return(x)
    } # bail if input is empty
    if (case == "lower") {
        x <- tolower(x)
    }
    if (case == "upper") {
        x <- toupper(x)
    }
    if (content == "dna") {
        src <- "acgturykmswbdhvnxACGTURYKMSWBDHVNX-"
        dest <- "tgcaayrmkswvhdbnxTGCAAyRMKSWVHDBNX-"
    } else {
        src <- "acgturykmswbdhvnxACGTURYKMSWBDHVNX-"
        dest <- "ugcaayrmkswvhdbnxUGCAAyRMKSWVHDBNX-"
    }
    if (max(nchar(x)) > 1) {
        return(chartr(src, dest, strreverse(x)))
    }
    # x is not a single string, so process it as a vector
    chartr(src, dest, rev(x))
}

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
        "<DEL>",
        "<INS>",
        "<DUP>",
        "<INV>",
        "<FUS>",
        "<CNV>",
        "<CNV:TR>",
        "<DUP:TANDEM>",
        "<DEL:ME>",
        "<INS:ME>"
    )
    iupac_indicators <- c(
        "R",
        "K",
        "S",
        "Y",
        "M",
        "W",
        "B",
        "H",
        "N",
        "D",
        "V"
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
    if (
        nchar(ref) != nchar(alt) &&
            !grepl(paste0("^", ref), alt) &&
            !grepl(paste0("^", alt), ref)
    ) {
        return("complex")
    }
    # Otherwise, uncategorized
    return("uncategorized")
}

#' Map column names of mutation data to default column names.

#' A utility function that renames columns of mutation data to default columns
#' names.
#' @param data mutation data
#' @param column_map a list that maps synonymous column names to their default.
#' @returns the mutation data with column names changed to match default.
#' @examples
#' df <- data.frame(
#'   chromosome = c("chr1", "chr2", "chr3"),
#'   pos = c(100, 200, 300),
#'   end = c(100, 200, 300),
#'   sample_id = c("S1", "S2", "S3"),
#'   reference = c("G", "C", "T"),
#'   alternate = c("A", "T", "G")
#' )
#' renamed_data <- rename_columns(df, column_map = op$column)
#' @export

rename_columns <- function(data, column_map = op$column) {
    original_colnames <- colnames(data)

    # normalized names (clean regex)
    # remove leading X or dots, trailing dots,
    # and replace inner dots with underscores
    norm_names <- tolower(original_colnames)
    norm_names <- gsub("^((x\\.+)|(\\.+))?", "", norm_names) # Leading
    norm_names <- gsub("(\\.+)?$", "", norm_names) # Trailing
    norm_names <- gsub("\\.+", "_", norm_names) # Middle dots to _

    map_synonyms <- names(column_map)
    map_targets <- unlist(column_map)

    # Handle existing defaults (casing)
    is_target <- norm_names %in% map_targets
    if (any(is_target)) {
        target_indices <- match(norm_names[is_target], map_targets)
        original_colnames[is_target] <- map_targets[target_indices]
    }

    # Identify targets that are still missing from the data
    # Only rename synonyms if the target doesn't exist yet
    present_targets <- original_colnames[original_colnames %in% map_targets]
    targets_needed <- setdiff(unique(map_targets), present_targets)

    # Find synonyms for needed targets
    synonym <- map_synonyms %in% norm_names
    target_is_needed <- map_targets %in% targets_needed
    candidate_indices <- which(synonym & target_is_needed)

    if (length(candidate_indices > 0)) {
        selected_indices <- candidate_indices[
            !duplicated(map_targets[candidate_indices])
        ]
        final_synonyms <- map_synonyms[selected_indices]
        final_targets <- map_targets[selected_indices]
        col_indices <- match(final_synonyms, norm_names)
        invisible(Map(
            function(orig, new) {
                message(
                    "Expected '",
                    new,
                    "' but found '",
                    original_colnames[orig],
                    "', renaming it."
                )
            },
            col_indices,
            final_targets
        ))

        # Update names
        original_colnames[col_indices] <- final_targets
    }
    colnames(data) <- original_colnames
    return(data)
}

#' Check that all required columns are present before proceeding with the function
#'
#' A utility function that will check that all required columns are present.
#' @param data mutation data
#' @param required_columns a list of required column names.
#' @returns an error
#' @examples
#' df <- data.frame(
#'   contig = c("chr1", "chr2", "chr3"),
#'   start = c(100, 200, 300),
#'   end = c(100, 200, 300),
#'   sample = c("S1", "S2", "S3"),
#'   ref = c("G", "C", "T"),
#'   alt = c("A", "T", "G")
#' )
#' check_required_columns(df, required_columns = op$base_required_mut_cols)
#' @export

check_required_columns <- function(data, required_columns) {
    missing_columns <- setdiff(tolower(required_columns), tolower(names(data)))

    if (length(missing_columns) > 0) {
        missing_col_names <- paste(missing_columns, collapse = ", ")
        stop(
            "Some required columns are missing",
            "or their synonyms are not found: ",
            missing_col_names
        )
    } else {
        return(data)
    }
}

#' Retrieve the sample column from VCF files
#' @description Checks to find the sample name of the vcf in the INFO field or
#' in the FORMAT header. Can also handle sample name synonyms.
#' @param vcf The imported VCF
#' @importFrom VariantAnnotation info
#' @importFrom SummarizedExperiment colData
#' @returns The vcf with sample column name corrected
vcf_sample_fix <- function(vcf) {
    # Check INFO for Sample column (Incl synonyms)
    original_names <- names(VariantAnnotation::info(vcf))
    # Normalize names
    norm_names <- tolower(original_names)
    norm_names <- gsub("[ .]", "_", norm_names)
    # check for synonyms
    synonyms <- c("sample", "sample_name", "sample_id")
    match_idx <- match(synonyms, norm_names)
    found_idx <- match_idx[!is.na(match_idx)]
    if (length(found_idx) > 0) {
        # Rename the first match found
        names(VariantAnnotation::info(vcf))[found_idx[1]] <- "sample"
    } else if (!"sample" %in% norm_names) {
        # Fallback to colData rownames (VCF header sample name)
        # Must have 1 sample per file as per docs
        sample_name <- rownames(SummarizedExperiment::colData(vcf))
        VariantAnnotation::info(vcf)$sample <- sample_name
    }
    return(vcf)
}


#' Validate BSgenome Input
#'
#' @description
#' Internal utility function to validate the \code{BS_genome} argument prior to
#' sequence context extraction. Ensures that the provided genome is a valid
#' \pkg{BSgenome} package name and that it is installed locally.
#'
#' @param BS_genome A character string specifying the package name of a
#' \pkg{BSgenome} object (e.g., \code{"BSgenome.Hsapiens.UCSC.hg38"}), or
#' \code{NULL}.
#'
#' @details
#' This function performs three checks:
#' \enumerate{
#'   \item If \code{BS_genome} is \code{NULL}, an error is thrown indicating that
#'   a genome must be provided when sequence context is required.
#'   \item If \code{BS_genome} is not among the available \pkg{BSgenome}
#'   packages, an error is thrown.
#'   \item If \code{BS_genome} is valid but not installed locally, an error is
#'   thrown with instructions to install it via \code{BiocManager::install()}.
#' }
#'
#' This function is intended to be called only when sequence context needs to be
#' populated (i.e., when a \code{context} column is absent or incomplete).
#'
#' @return Invisibly returns \code{TRUE} if validation passes; otherwise, an
#' error is raised.
#'
#' @keywords internal
validate_BS_genome <- function(BS_genome) {
    if (is.null(BS_genome)) {
        stop(
            "The context column is missing, and no BS_genome was provided. ",
            "Please install and specify an appropriate BSgenome package."
        )
    }

    available_BS_genomes <- BSgenome::available.genomes(
        splitNameParts = TRUE
    )$pkgname

    if (!BS_genome %in% available_BS_genomes) {
        stop(
            "The specified BS genome ('",
            BS_genome,
            "') is not a recognized BSgenome package."
        )
    }

    if (!BS_genome %in% BSgenome::installed.genomes()) {
        stop(
            "The specified BS genome ('",
            BS_genome,
            "') is valid, but is not installed locally. ",
            "Please install it using BiocManager::install('",
            BS_genome,
            "') before proceeding."
        )
    }
}
