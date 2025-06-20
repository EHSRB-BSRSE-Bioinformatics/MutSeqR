#' Import a VCF file
#'
#' @description The function reads VCF file(s) and extracts the
#' data into a dataframe.
#' @param vcf_file The path to the .vcf (.gvcf, gzip, bgzip) to be
#' imported. If you specify a directory, the function will
#' attempt to read all files in the directory and combine them into
#' a single table. VCF files should follow the VCF specifications,
#' version 4.5. Multisample VCF files are not supported; VCF files
#' must contain one sample each. Required fields are listed below.
#' \itemize{
#' \item `FIXED FIELDS`
#' \item `CHROM`: The name of the reference sequence. Equivalent to `contig`.
#' \item `POS`: The 1-based start position of the feature. Equivalent to  `start`.
#' \item `REF`: The reference allele at this position.
#' \item `ALT`: The left-aligned, normalized, alternate allele at this position.
#' Multiple alt alleles called for a single position should be represented as
#' separate rows in the table.
#' \item `INFO FIELDS`
#' \item `END`: The half-open end position of the feature.
#' \item `sample`: An identifying field for your samples; either in the INFO
#' field or as the header to the FORMAT field.
#' }
#' The following FORMAT fields are not required, but are recommended for full
#' package functionality:
#' \itemize{
#' \item `AD`: The allelic depths for the reference and alternate allele
#' in the order listed. The sum of AD is equivalent to the `total_depth`
#' (read depth at this position excluding N-calls).
#'  \item `DP`: The read depth at this position (including N-calls).
#' Equivalent to `depth`. Note that in many VCF files, the DP field
#' is defined as `total_depth`. However, in most cases, the DP field
#' includes N-calls.
#'  \item `VD`: The read depth supporting the alternate allele. If
#' not included, the function will add this column, assuming a value of 1.
#' Equivalent to `alt_depth`.
#' }
#' We recommend that files include a record for every sequenced
#' position, regardless of whether a variant was called, along with the
#' `AD` for each record. This enables site-specific depth calculations
#' required for some downstream analyses. AD is used to calculate the `total_depth`
#' (the read depth excluding No-calls). If AD is not available, the `DP` field
#' will be used as the `total_depth`.
#' @param sample_data An optional file containing additional sample
#' metadata (dose, timepoint, etc.). This can be a data frame or a file path.
#' Metadata will be joined with the mutation data based on the sample column.
#' Required columns are `sample` and any additional columns you wish to
#' include.
#' @param sd_sep The delimiter for importing sample metadata tables.
#' Default is tab-delimited.
#' @param regions An optional file containing metadata of genomic regions.
#' Region metadata will be joined with mutation data and variants will be
#' checked for overlap with the regions. `regions` can be either a file path,
#' a data frame, or a GRanges object. File paths will be read using the rg_sep.
#' Users can also choose from the built-in TwinStrand's Mutagenesis Panels by
#' inputting "TSpanel_human",  "TSpanel_mouse", or "TSpanel_rat". Required
#' columns for the regions file are "contig", "start", and "end". For a GRanges
#' object, the required columns are "seqnames", "start", and "end". Default is
#' NULL.
#' @param rg_sep The delimiter for importing the custom_regions. The default is
#' tab-delimited "\\t".
#' @param is_0_based_rg A logical variable. Indicates whether the position
#' coordinates in `regions` are 0 based (TRUE) or 1 based (FALSE).
#' If TRUE, positions will be converted to 1-based (start + 1).
#' Need not be supplied for TSpanels. Default is TRUE.
#' @param padding Extend the range of your regions
#' in both directions by the given amount. Ex. Structural variants and
#' indels may start outside of the regions. Adjust the `padding` to
#' include these variants in your region's ranges.
#' @param genome The genome assembly version of the reference genome. This is
#' required if your data does not include a context column. The
#' function will install a BS genome for the given species/genome/masked
#' arguments to populate the context column.
#' Ex.Human GRCh38 = hg38 | Human GRCh37 = hg19 | Mouse GRCm38 = mm10 |
#' Mouse GRCm39 = mm39 | Rat RGSC 6.0 = rn6 | Rat mRatBN7.2 = rn7
#' @param species The species. Required if your data does not include a
#' context column. The function will install a BS genome for the given
#' species/genome/masked to populate the context column. The species can
#' be the common name of the species or the scientific name.
#' Ex. "human" or "Homo sapiens".
#' @param masked_BS_genome A logical value. Required when using a BS genome
#' to poulate the context column. Whether to use the masked version of the
#' BS genome (TRUE) or not (FALSE). Default is FALSE.
#' @param output_granges `TRUE` or `FALSE`; whether you want the mutation
#' data to output as a GRanges object. Default output is as a dataframe.
#' @returns A table where each row is a mutation, and columns indicate the
#' location, type, and other data. If `output_granges` is set to TRUE, the
#' mutation data will be returned as a GRanges object, otherwise mutation
#' data is returned as a dataframe.
#'
#' Output Column Definitions:
#' \itemize{
#' \item `short_ref`: The reference base at the start position.
#' \item `normalized_ref`: The short_ref in C/T-base notation for
#' this position (e.g. A -> T, G -> C).
#' \item `context` The trinucleotide context at this position. Consists
#' of the reference base and the two flanking bases (e.g. TAC).
#' \item `normalized_context`: The trinucleotide context in C/T base
#' notation for this position (e.g. TAG -> CTA).
#'  \item `variation_type` The type of variant (snv, mnv, insertion,
#' deletion, complex, sv, no_variant, ambiguous, uncategorized).
#' \item `subtype` The substitution type for the snv variant (12-base spectrum;
#' e.g. A>C).
#' \item `normalized_subtype` The C/T-based substitution type for the snv
#' variant (6-base spectrum; e.g. A>C -> T>G).
#' \item `context_with_mutation`: The substitution type for the snv variant
#' including the two flanking nucleotides (192-trinucleotide spectrum;
#' e.g. `T[A>C]G`)
#' \item `normalized_context_with_mutation`: The C/T-based substitution
#' type for the snv variant including the two flanking nucleotides
#' (96-base spectrum e.g. `T[A>C]G` -> `C[T>G]A`).
#' \item `nchar_ref`: The length (in bp) of the reference allele.
#' \item `nchar_alt`: The length (in bp) of the alternate allele.
#' \item `varlen`: The length (in bp) of the variant.
#' \item `ref_depth`: The depth of the reference allele. Calculated as
#' `total_depth` - `alt_depth`, if applicable.
#' \item `vaf` : The variant allele fraction. Calculated as
#' `alt_depth`/`total_depth`.
#' \item `gc_content`: % GC of the trinucleotide context at this position.
#' \item `is_known`: TRUE or FALSE. Flags known variants (ID != ".").
#' \item `row_has_duplicate`: TRUE or FALSE. Flags rows whose position is
#' the same as that of at least one other row for the same sample.
#' }
#' @examples
#' # Example: Import a single bg-zipped vcf file. This library was sequenced
#' # with Duplex Sequencing using the TwinStrand Mouse Mutagenesis Panel which
#' # consists of 20 2.4kb targets = 48kb of sequence.
#' example_file <- system.file("extdata", "Example_files",
#'                             "example_import_vcf_data_cleaned.vcf.bgz",
#'                             package = "MutSeqR")
#' # We will create an example metadata table for this data.
#' sample_meta <- data.frame(sample = "dna00996.1",
#'                           dose = "50",
#'                           dose_group = "High")
#' # Import the data
#' imported_example_data <- import_vcf_data(vcf_file = example_file,
#'                                          sample_data = sample_meta,
#'                                          regions = "TSpanel_mouse",
#'                                          genome = "mm10",
#'                                          species = "mouse",
#'                                          masked_BS_genome = FALSE)
#' @importFrom  VariantAnnotation alt info geno readVcf ref rbind
#' @importFrom dplyr filter group_by left_join mutate rename select summarize ungroup
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stringr str_sub str_count
#' @importFrom SummarizedExperiment colData
#' @importFrom plyranges join_overlap_left
#' @importFrom Biostrings getSeq
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom BiocGenerics strand start end
#' @importFrom GenomeInfoDb seqnames
#' @export
#'
import_vcf_data <- function(vcf_file,
                            sample_data = NULL,
                            sd_sep = "\t",
                            regions = NULL,
                            rg_sep = "\t",
                            is_0_based_rg = FALSE,
                            padding = 0,
                            genome = NULL,
                            species = NULL,
                            masked_BS_genome = FALSE,
                            output_granges = FALSE) {
  vcf_file <- file.path(vcf_file)

  # Check if a sample identifier is already present in the INFO field
  check_and_rename_sample <- function(vcf) {
    # Define possible variations of sample identifier names
    possible_sample_names <- c("sample", "sample_name", "sample_id")
    # Initialize the sample_info column
    sample_info <- NULL
    # Search for variations of sample identifier names
    for (sample_name_var in possible_sample_names) {
      sample_name_var <- tolower(sample_name_var)  # Make the comparison case-insensitive
      if (sample_name_var %in% tolower(names(VariantAnnotation::info(vcf)))) {
        # If found, rename it to "sample" and break the loop
        names(VariantAnnotation::info(vcf))[tolower(names(VariantAnnotation::info(vcf))) == sample_name_var] <- "sample"
        break
      }
    }
    # If "sample" still doesn't exist, create it using colData
    if (!"sample" %in% colnames(VariantAnnotation::info(vcf))) {
      sample_info <- rownames(SummarizedExperiment::colData(vcf))
      VariantAnnotation::info(vcf)$sample <- sample_info
    }
    return(vcf)
  }

  # Read and bind vcfs from folder
  if (file.info(vcf_file)$isdir == TRUE) {
    vcf_files <- list.files(path = vcf_file,
                            pattern = "\\.g?vcf(\\.bgz|\\.gz)?$",
                            full.names = TRUE)
    # FIX: add check for empty file list.
    # Initialize an empty VCF object to store the combined data
    vcf <- NULL
    # Read and combine VCF files
    for (file in vcf_files) {
      vcf_list <- VariantAnnotation::readVcf(file)
      # Rename or create the "sample" column in the INFO field
      vcf_list <- suppressWarnings(check_and_rename_sample(vcf_list))
      # Ensure consistent column names
      rownames(SummarizedExperiment::colData(vcf_list)) <- "sample_info"
      # Combine the VCF data
      if (is.null(vcf)) {
        vcf <- vcf_list
      } else {
        vcf <- suppressWarnings(VariantAnnotation::rbind(vcf, vcf_list))
      }
    }
  } else {
    # Read a single vcf file
    vcf <- VariantAnnotation::readVcf(vcf_file)
    # Rename or create the "sample" column in the INFO field
    vcf <- suppressWarnings(check_and_rename_sample(vcf))
  }
  # Extract and Clean alt column
  ## May want to use the expand function to unlist ALT column of a CollapsedVCF object to one row per ALT value.
  alt <- VariantAnnotation::alt(vcf)
 # alt_values_clean <- lapply(alt, function(x) x[x != "<NON_REF>"])
  alt <- IRanges::CharacterList(alt)

  # Extract mutation data into a dataframe
  dat <- data.frame(
    contig = SummarizedExperiment::seqnames(vcf),
    start = SummarizedExperiment::start(vcf),
    ref = VariantAnnotation::ref(vcf),
    alt = alt
  )
  # Retain all INFO fields
  info <- as.data.frame(VariantAnnotation::info(vcf))
  # Extract GENO fields depending on the type of data
  geno <- VariantAnnotation::geno(vcf)
  geno_df <- data.frame(row.names = seq_len(nrow(geno[[1]])))
  for (field_name in names(geno)) {
    field <- geno[[field_name]]
    if (is.list(field)) { # Ex. AD
      max_length <- max(sapply(field, length))
      expanded_field <- do.call(rbind, lapply(field, function(x) {
        c(x, rep(NA, max_length - length(x)))
      }))
      colnames(expanded_field) <- paste(field_name, seq_len(max_length), sep = "_")
      geno_df <- cbind(geno_df, expanded_field)
    } else if (is.matrix(field)) { # Ex. GT, DP, VD
      geno_df[[field_name]] <- as.vector(field)
    } else if (is.array(field) && length(dim(field)) == 3) { # Ex. RD, ALD
      # Collapse the array over the 2nd and 3rd dimensions
      collapsed_field <- apply(field, c(1), function(x) as.vector(x))
      collapsed_field <- as.data.frame(t(collapsed_field))
      ncols <- ncol(collapsed_field)
      colnames(collapsed_field) <- paste(field_name, seq_len(ncols), sep = "_")
      geno_df <- cbind(geno_df, collapsed_field)
    } else {
      geno_df[[field_name]] <- field
    }
  }
  # Ensure info and geno do not have the same columns
  common_cols <- intersect(colnames(info), colnames(geno_df))
  info <- info[, !(colnames(info) %in% common_cols), drop = FALSE]

  # Combine data frames
  dat <- cbind(dat, geno_df, info)
  row.names(dat) <- NULL

  # Join with sample metadata if provided
  if (!is.null(sample_data)) {
    if (is.data.frame(sample_data)) {
      sampledata <- sample_data
      if (nrow(sampledata) == 0) {
        stop("Error: The sample data frame you've provided is empty")
      }
    } else if (is.character(sample_data)) {
      sample_file <- file.path(sample_data)
      if (!file.exists(sample_file)) {
        stop("Error: The sample data file path you've specified is invalid")
      }
      if (file.info(sample_file)$size == 0) {
        stop("Error: You are trying to import an empty sample data file")
      }
      sampledata <- read.delim(file.path(sample_data),
                               sep = sd_sep,
                               header = TRUE)
      if (ncol(sampledata) <= 1) {
        stop("Your imported sample data only has one column.
                            You may want to set sd_sep to properly reflect
                            the delimiter used for the data you are importing.")
      }
    } else {
      stop("Error: sample_data must be a character string or a data frame")
    }
    # Join
    dat <- dplyr::left_join(dat, sampledata, suffix = c("", ".sampledata"))
  }

  # Rename columns to default names
  dat <- rename_columns(dat)
  # Check for all required columns before proceeding
  dat <- MutSeqR::check_required_columns(dat, op$base_required_mut_cols)
  context_exists <- "context" %in% colnames(dat)

  # Check for NA values in required columns.
  # Except for the alt column, which can have NA values.
  required_columns <- setdiff(op$base_required_mut_cols, "alt")
  columns_with_na <- colnames(dat)[apply(dat, 2, function(x) any(is.na(x)))]
  na_columns_required <- intersect(columns_with_na,
                                   required_columns)
  if (length(na_columns_required) > 0) {
    stop(paste0("Error: NA values were found within the following required
                column(s): ", paste(na_columns_required, collapse = ", "),
                ".
                Please confirm that your data is complete before proceeding."))
  }
  # Check for NA values in the context column. If so, will populate it.
  if (context_exists) {
    if ("context" %in% columns_with_na) {
      context_exists <- FALSE
    }
  }
  # Join Regions
  # Turn mutation data into GRanges
  mut_ranges <- GenomicRanges::makeGRangesFromDataFrame(
    df = as.data.frame(dat),
    keep.extra.columns = TRUE,
    seqnames.field = "contig",
    start.field = "start",
    end.field = "end"
  )

  if (!is.null(regions)) {

    # load regions file
    regions_gr <- MutSeqR::load_regions_file(regions,
                                             rg_sep,
                                             is_0_based_rg)
    regions_gr$in_regions <- TRUE
    # Apply range buffer
    BiocGenerics::start(regions_gr) <- pmax(BiocGenerics::start(regions_gr) - padding, 1)
    BiocGenerics::end(regions_gr) <- BiocGenerics::end(regions_gr) + padding

    # Join mutation data and region data using overlap
    mut_ranges <- plyranges::join_overlap_left_within_directed(mut_ranges,
                                                               regions_gr,
                                                               suffix = c("",
                                                                          "_regions"))

    mut_ranges <- mut_ranges %>%
      plyranges::mutate(in_regions = ifelse(is.na(in_regions), FALSE, TRUE))

    false_count <- sum(mut_ranges$in_regions == FALSE)
    if (false_count > 0) {
      warning("Warning: ", false_count, " rows were outside of the specified regions. To remove these rows, use the filter_mut() function")
    }
  }
  # Create a context column, if needed
  if (!context_exists) {
    if (is.null(genome) || is.null(species)) {
      stop("Error: We need to populate the context column for your data. Please provide a genome and species so that we can retrieve the sequences.")
    }
    ref_genome <- install_ref_genome(organism = species,
                                     genome = genome,
                                     masked = masked_BS_genome)

    extract_context <- function(mut_gr,
                                bsgenome) {
      # Resize the mut_ranges to include the context
      expanded_ranges <- GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(mut_gr),
                                                ranges = IRanges::IRanges(
                                                  start = BiocGenerics::start(mut_gr) - 1,
                                                  end = BiocGenerics::start(mut_gr) + 1
                                                ),
                                                strand = BiocGenerics::strand(mut_gr))
      # Extract the sequences from the BSgenome
      sequences <- Biostrings::getSeq(bsgenome, expanded_ranges)
      # Return the sequences
      return(sequences)
    }
    message("Retrieving context sequences from the reference genome")
    context <- extract_context(mut_ranges, ref_genome)
    mut_ranges$context <- context
  }
  dat <- as.data.frame(mut_ranges) %>%
    dplyr::rename(contig = "seqnames")

  # Create is_known based on ID col, if present
  if ("id" %in% colnames(dat)) {
    dat <- dat %>% dplyr::mutate(is_known = ifelse(!.data$id == ".", "Y", "N"))
  }
  # Create variation_type
  if (!"variation_type" %in% colnames(dat)) {
    dat$variation_type <- mapply(MutSeqR::classify_variation, dat$ref, dat$alt)
  } else {
    dat <- dplyr::rename(dat, original_variation_type = "variation_type")
    dat$variation_type <- mapply(MutSeqR::classify_variation, dat$ref, dat$alt)
  }

  # Define substitution dictionary to normalize to pyrimidine context
  sub_dict <- c(
    "G>T" = "C>A", "G>A" = "C>T", "G>C" = "C>G",
    "A>G" = "T>C", "A>C" = "T>G", "A>T" = "T>A"
  )
  # Calculate columns:
  # nchar_ref, nchar_alt, varlen, short_ref, normalized_ref, subtype,
  # normalized_subtype, normalized_context, context_with_mutation,
  # normalized_context_with_mutation, gc_content
  dat <- dat %>%
    dplyr::mutate(
      nchar_ref = nchar(ref),
      nchar_alt = ifelse(!(.data$variation_type %in% c("no_variant",
                                                       "sv",
                                                       "ambiguous",
                                                       "uncategorized")),
                         nchar(alt), NA),
      varlen =
        ifelse(.data$variation_type %in% c("insertion", "deletion", "complex"),
          .data$nchar_alt - .data$nchar_ref,
          ifelse(.data$variation_type %in% c("snv", "mnv"), .data$nchar_ref,
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
      subtype =
        ifelse(.data$variation_type == "snv",
          paste0(.data$ref, ">", .data$alt),
          .data$variation_type
        ),
      normalized_subtype = ifelse(.data$subtype %in% names(sub_dict),
                                  sub_dict[.data$subtype],
                                  .data$subtype),
      normalized_context = ifelse(
        stringr::str_sub(.data$context, 2, 2) %in% c("G", "A"),
        mapply(function(x) MutSeqR::reverseComplement(x, case = "upper"),
               .data$context),
        .data$context),
      context_with_mutation =
        ifelse(.data$variation_type == "snv",
               paste0(stringr::str_sub(.data$context, 1, 1),
                      "[", .data$subtype, "]",
                      stringr::str_sub(.data$context, 3, 3)),
               .data$variation_type),
      normalized_context_with_mutation =
        ifelse(.data$variation_type == "snv",
               paste0(stringr::str_sub(.data$normalized_context, 1, 1),
                      "[", .data$normalized_subtype, "]",
                      stringr::str_sub(.data$normalized_context, 3, 3)),
               .data$variation_type),
      gc_content = (stringr::str_count(string = .data$context, pattern = "G") +
                    stringr::str_count(string = .data$context, pattern = "C"))
      / stringr::str_count(.data$context),
      filter_mut = FALSE
    )

  # Depth
  # Add alt_depth column, if it doesn't exist
  if (!"alt_depth" %in% colnames(dat)) {
    dat$alt_depth <- 1
  }
  # Create a total_depth column, if able
  # Create total_depth and no_calls columns based on set parameter depth_calc.
  # Requires AD field in FORMAT of vcf. If this field is missing, we use depth instead of total_depth
  total_depth_exists <- "total_depth" %in% colnames(dat)
  depth_exists <- "depth" %in% colnames(dat)
  no_calls_exists <- "no_calls" %in% colnames(dat)
  ad_columns <- grep("^AD_", colnames(dat), value = TRUE)

  if (!total_depth_exists) {
    if (no_calls_exists && depth_exists) {
      dat <- dat %>%
        dplyr::mutate(total_depth = .data$depth - .data$no_calls)
    } else if (length(ad_columns) > 0) { # create total_depth from AD
      dat$total_depth <- rowSums(dat[, ad_columns], na.rm = TRUE)
    } else { # use the DP field
      if (depth_exists) {
        dat <- dat %>%
          dplyr::mutate(
            total_depth = .data$depth)
        warning("Could not find total_depth column and cannot calculate. The 'total_depth' will be set to DP. You can review the diffference in the README")
      } else {
        warning("Could not find an appropriate depth column. Some package functionality may be limited.\n")
      }
    }
  }

  # Check for duplicated rows
  dat <- dat %>%
    dplyr::group_by(.data$sample, .data$contig, .data$start) %>%
    dplyr::mutate(row_has_duplicate = dplyr::n() > 1) %>%
    dplyr::ungroup()

  if (sum(dat$row_has_duplicate) > 0) {
    warning(sum(dat$row_has_duplicate), " rows were found whose position was the same as that of at least one other row for the same sample.")

    # Warn about the depth for the duplicated rows
    if ("total_depth" %in% colnames(dat)) {
      warning("The total_depth may be double-counted in some instances due to overlapping positions. Set the correct_depth parameter in calculate_mf() to correct the total_depth for these instances.")
    }
  }

  # Make VAF and ref_depth columns, if depth exists
  if ("total_depth" %in% colnames(dat)) {
    dat <- dat %>%
      dplyr::mutate(vaf = .data$alt_depth / .data$total_depth,
                    ref_depth = .data$total_depth - .data$alt_depth)
  }

  # Add empty filter column
  if (!"filter_mut" %in% colnames(dat)) {
    dat$filter_mut <- FALSE
  }

  if (output_granges) {
    gr <-  GenomicRanges::makeGRangesFromDataFrame(
      df = dat,
      keep.extra.columns = TRUE,
      seqnames.field = "contig",
      start.field = "start",
      end.field = "end",
      starts.in.df.are.0based =  FALSE)
    return(gr)
  } else {
    return(dat)
  }
}
