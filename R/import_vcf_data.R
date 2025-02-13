#' Import a vcf file
#'
#' @description The function reads the genomic vcf file(s) and extracts the
#' data into a dataframe. The function also reads in sample metadata
#' if provided and joins it with the mutation data. An interval list of genomic
#' regions can be provided to filter out variants that occur outside of the
#' defined regions' ranges. The function will use the reference genome to
#' extract the trinucleotide context of every position in the mutation data.
#' The function can output the mutation data as a dataframe or a granges
#' object.
#' @param vcf_file The path to the genomic .vcf or .vcf.gz file(s)  to be imported. If
#' you specify a folder, the function will attempt to read all files in the
#' folder and combine them into on dataset. Multisample vcf files are not
#' supported; vcf files must contain one sample each. Required fields are
#' listed below
#'  - FIXED FIELDS:
#'  - `CHROM`: The reference sequence name.Equivalent to `contig`
#'  - `POS`: 0-based start position of the feature in contig.
#'  - `REF`: The reference allele at this position
#'  - `ALT`: The left-aligned, normalized, alternate allele at this position.
#' - INFO FIELDS
#'  - `END`: The half-open end position of the feature in contig.
#'  - `sample`: An identifying field for your samples; either in the INFO
#' field or as the header to the FORMAT field.
#' - SUGGESTED FIELDS:
#' - FORMAT `AD`: The allelic depths for the reference and alternate alleles in the
#' order listed.
#'  - FORMAT `DP`: The total read depth at this position (including N-calls).
#' Equivalent to `depth`.
#'  - FORMAT `VD`: Variant Depth. Equivalent to `alt_depth`.
#'  - INFO `SVTYPE`: Structural variant types; INV DUP DEL INS FUS.
#'  - INFO `SVLEN`: Length of the structural variant in base pairs.
#' @param sample_data An optional file containing additional sample
#' metadata (dose, timepoint, etc.). This can be a data frame or a file path.
#' @param sd_sep The delimiter for importing sample metadata tables.
#' Default is tab-delimited
#' @param regions "TSpanel_human", "TSpanel_mouse", "TSpanel_rat" ,
#' "custom" or "none". The 'TSpanel_' argument refers to the TS
#' Mutagenesis panel of the specified species, or to a custom regions
#' interval file. If set to 'custom', please provide the file path in
#' custom_regions and the genome assembly version of the reference
#' genome using the 'genome' parameter. If you are not using a targeted
#' approach, set regions to none, and supply the species and genome
#' assembly of the reference genome using the 'species' and 'genome'
#' parameters respectively.
#' @param custom_regions "filepath". If regions is set to custom,
#'  provide the file path for the file containing regions metadata.
#'  Required columns are "contig", "start", and "end"
#' @param rg_sep The delimiter for importing the custom_regions.
#' Default is tab-delimited
#' @param genome The genome assembly of the reference genome.
#' Ex.Human GRCh38 = hg38 | Human GRCh37 = hg19 | Mouse GRCm38 = mm10 |
#' Mouse GRCm39 = mm39 | Rat RGSC 6.0 = rn6 | Rat mRatBN7.2 = rn7
#' @param species The species of the reference genome. Required if
#' regions is set to none. The value can be the common name of the species
#' or the scientific name. Ex. "human" or "Homo sapiens".
#' @param range_buffer An integer >= 0 .Required if using a targetted
#' approach.  Use the range-buffer to extend the range outside
#' of a region within which a variant can occur. The default is 0 nucleotides
#' outside of region ranges. Ex. Structural variants and indels may start
#' outside of the regions. Adjust the range_buffer to include these variants
#' in the target sequencies.
#' @param output_granges `TRUE` or `FALSE`; whether you want the mutation
#' data to output as a GRanges object. Default output is as a dataframe.
#' @returns A data frame or a GRanges object where each row is a mutation,
#' and columns indicate the location, type, and other data.
#' @importFrom  VariantAnnotation alt info geno readVcf ref rbind
#' @importFrom dplyr filter group_by left_join mutate rename select summarize ungroup
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stringr str_sub str_count
#' @importFrom SummarizedExperiment colData
#' @importFrom plyranges join_overlap_left
#' @importFrom Biostrings getSeq
#' @export
#'
import_vcf_data <- function(
    vcf_file,
    sample_data = NULL,
    sd_sep = "\t",
    regions = c("TSpanel_human", "TSpanel_mouse", "TSpanel_rat", "custom", "none"),
    custom_regions = NULL,
    rg_sep = "\t",
    range_buffer = 0,
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
    vcf_files <- list.files(path = vcf_file, pattern = "\\.(vcf|gvcf)$", full.names = TRUE)
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
  alt_values_clean <- lapply(alt, function(x) x[x != "<NON_REF>"])
  alt <- IRanges::CharacterList(alt_values_clean)

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

  if (regions != "none") {

    # load regions file
      regions_df <- MutSeqR::load_regions_file(regions, custom_regions, rg_sep)
      regions_df$in_regions <- TRUE

    # Apply range buffer
    regions_df <- regions_df %>%
      dplyr::mutate(start = .data$start - range_buffer,
                    end = .data$end + range_buffer)

    # adjust start position to be 1-based for TSpanels
    if (regions %in% c("TSpanel_mouse", "TSpanel_human", "TSpanel_rat")) {
      is_0_based_rg <- TRUE
    }

    # Turn region data into GRanges
    region_ranges <- GenomicRanges::makeGRangesFromDataFrame(
      df = regions_df,
      keep.extra.columns = TRUE,
      seqnames.field = "contig",
      start.field = "start",
      end.field = "end",
      starts.in.df.are.0based = is_0_based_rg
    )

    # Join mutation data and region data using overlap
    mut_ranges <- plyranges::join_overlap_left_within_directed(mut_ranges,
                                                               region_ranges,
                                                               suffix = c("",
                                                                          "_regions"))

    mut_ranges <- mut_ranges %>%
      plyranges::mutate(in_regions = ifelse(is.na(in_regions), FALSE, TRUE))

    false_count <- sum(mut_ranges$in_regions == FALSE)
    if (false_count > 0) {
      warning("Warning: ", false_count, " rows were outside of the specified regions.\n
        To remove these rows, use the filter_mut() function")
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

    extract_context <- function(mutations,
                                bsgenome,
                                upstream = 1,
                                downstream = 1) {
      # Resize the mut_ranges to include the context
      expanded_ranges <- GenomicRanges::resize(x = mutations,
                                               width = upstream + downstream + 1,
                                               fix = "center")
      # Extract the sequences from the BSgenome
      sequences <- Biostrings::getSeq(bsgenome, expanded_ranges)
      # Return the sequences
      return(sequences)
    }
    context <- extract_context(mut_ranges, ref_genome)
    mut_ranges$context <- context
    dat <- as.data.frame(mut_ranges) %>%
      dplyr::rename(contig = "seqnames")
  }

  # Create is_known based on ID col, if present
  if ("id" %in% colnames(dat)) {
    dat <- dat %>% dplyr::mutate(is_known = ifelse(!.data$id == ".", "Y", "N"))
  }
  # Create variation_type
  if (!"variation_type" %in% colnames(dat)) {
    dat$variation_type <- mapply(classify_variation, dat$ref, dat$alt)
  } else {
    dat <- dplyr::rename(dat, original_variation_type = "variation_type")
    dat$variation_type <- mapply(classify_variation, dat$ref, dat$alt)
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
      / stringr::str_count(.data$context)
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
  ad_columns <- grep("^ad_", colnames(dat), value = TRUE)

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
            total_depth = .data$depth,
            vaf = .data$alt_depth / .data$total_depth)
        warning("Could not find total_depth column.\n
          Could not calculate total_depth\n
          No Allelic Depth (AD) field.\n
          The 'total_depth' will be set to 'depth' (DP)\n
          You can review the definitions of each column in the README")
      } else {
        warning("Could not find an appropriate depth column.\n
            Some package functionality may be limited.\n")
      }
    }
  }

  # Check for duplicated rows
  dat <- dat %>%
    dplyr::group_by(.data$sample, .data$contig, .data$start) %>%
    dplyr::mutate(row_has_duplicate = n() > 1) %>%
    dplyr::ungroup()

  if (sum(dat$row_has_duplicate) > 0) {
    warning(sum(dat$row_has_duplicate), " rows were found whose
    position was the same as that of at least one other row for the same
    sample.")

    # Warn about the depth for the duplicated rows
    if ("total_depth" %in% colnames(dat)) {
      warning("The total_depth may be double-counted in some instances due to
      overlapping positions. Use the filter_mut() function to correct the
      total_depth for these instances.")
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
      seqnames.field = "seqnames",
      start.field = "start",
      end.field = "end",
      starts.in.df.are.0based =  FALSE)
      return(gr)
  } else {
    return(dat)
  }
}
