#' Import tabular mutation data
#'
#' Imports tabular mutation file into the local R environment.
#' @param mut_file The mutation data file(s) to be imported.
#' This can be either a data frame object or a filepath
#' to a file or directory. If you specify a directory, the function will
#' attempt to read all files in the directory and combine them into
#' a single data frame. Mutation data should consist of a row for each
#' variant. Required columns are listed in details.
#' @param mut_sep The delimiter for importing the mutation file.
#' Default is tab-delimited.
#' @param is_0_based_mut A logical variable. Indicates whether the
#' position coordinates in the mutation data are 0 based (TRUE) or
#' 1 based (FALSE). If TRUE, positions will be converted to 1-based.
#' @param sample_data An optional file containing additional sample
#' metadata (dose, timepoint, etc.). This can be a data frame or a file path.
#' Metadata will be joined with the mutation data based on the sample column.
#' Required columns are `sample` and any additional columns you wish to
#' include.
#' @param sd_sep The delimiter for importing sample data.
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
#' @param padding An integer >= 0. Extend the range of your regions
#' in both directions by the given amount. Ex. Structural variants and
#' indels may start outside of the regions. Adjust the padding to
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
#' @param custom_column_names A list of names to specify the meaning of column
#'  headers. Since column names can vary with data, this might be necessary to
#'  digest the mutation data properly. Typical defaults are set, but can
#'  be substituted in the form of `list(my_custom_contig_name = "contig",
#'  my_custom_sample_column_name = "sample")`. You can change one or more of
#' these. Set column synonyms are defined in MutSeqR::op$column and will
#' automatically be changed to their default value.
#' @param output_granges A logical variable; whether you want the mutation
#' data to output as a GRanges object. Default output (FALSE) is as a dataframe.
#' @details Required columns for mut files are:
#' \itemize{
#'      \item `contig`: The name of the reference sequence.
#'      \item `start`: The start position of the feature.
#'      \item `end`: The half-open end position of the feature.
#'      \item `sample`: The sample name.
#'      \item `ref`: The reference allele at this position
#'      \item `alt`: The left-aligned, normalized, alternate allele at this
#' position. Multiple alt alleles called for a single position should be
#' represented as separate rows in the table.
#' }
#' The following columns are not required, but are recommended for full
#' package functionality:
#' \itemize{
#'   \item `alt_depth`: The read depth supporting the alternate allele. If
#' not included, the function will add this column, assuming a value of 1.
#'    \item `total_depth`: The total read depth at this position, excluding
#' no-calls (N calls). If not present, the function will attempt to calculate
#' the `total_depth` as `depth` - `no_calls`. If no_calls is not present, the
#' function will use `depth` as the `total_depth.`
#'    \item `depth`: The total read depth at this position, including no-calls.
#'    \item `no_calls`: The number of no-calls (N-calls) at this position.
#' }
#' We recommend that files include a record for every sequenced
#' position, regardless of whether a variant was called, along with the
#' `total_depth` for each record. This enables site-specific depth calculations
#' required for some downstream analyses.
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
#' \item `filter_mut` : A logical value, initially set to FALSE that indicates
#' to calculte_mf() if the variant should be excluded from mutation counts.
#' See the filter_mut function for more detail.
#' }
#' @examples
#' # Example: Import a single mutation file. This library was sequenced with
#' # Duplex Sequencing using the TwinStrand Mouse Mutagenesis Panel which
#' # consists of 20 2.4kb targets = 48kb of sequence.
#' example_file <- system.file("extdata", "Example_files",
#'                             "example_import_mut_data.rds",
#'                             package = "MutSeqR")
#' example_data <- readRDS(example_file)
#' # We will create an example metadata table for this data.
#' sample_meta <- data.frame(sample = "dna00996.1",
#'                           dose = "50",
#'                           dose_group = "High")
#' # Import the data
#' imported_example_data <- import_mut_data(mut_file = example_data,
#'                                          sample_data = sample_meta,
#'                                          regions = "TSpanel_mouse",
#'                                          genome = "mm10",
#'                                          species = "mouse",
#'                                          masked_BS_genome = FALSE)
#' @importFrom dplyr bind_rows mutate left_join case_when
#' @importFrom magrittr %>%
#' @importFrom stringr str_sub str_count
#' @importFrom plyranges join_overlap_left
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom utils read.delim read.table
#' @importFrom rlang .data
#' @importFrom BiocGenerics strand start end
#' @importFrom IRanges IRanges
#' @importFrom Biostrings getSeq
#' @importFrom GenomeInfoDb seqnames
#' @export
import_mut_data <- function(mut_file,
                            mut_sep = "\t",
                            is_0_based_mut = TRUE,
                            sample_data = NULL,
                            sd_sep = "\t",
                            regions = NULL,
                            rg_sep = "\t",
                            is_0_based_rg = TRUE,
                            padding = 0,
                            genome = NULL,
                            species = NULL,
                            masked_BS_genome = FALSE,
                            custom_column_names = NULL,
                            output_granges = FALSE) {                             

  if (!is.numeric(padding) || padding < 0) {
    stop("Error: The range buffer must be a non-negative number")
  }
  if (!is.logical(is_0_based_mut) || !is.logical(is_0_based_rg)) {
    stop("Error: is_0_based must be a logical variable")
  }

  if (!is.null(custom_column_names)) {
    if (!is.list(custom_column_names)) {
      stop("Error: custom_column_names must be a list")
    }
  }
  if (!is.logical(output_granges)) {
    stop("Error: output_granges must be a logical variable")
  }

  # Import the mut files: data frame or file path
  if (is.data.frame(mut_file)) {
    dat <- mut_file
    if (nrow(dat) == 0) {
      stop("Error: The data frame you've provided is empty")
    }
  } else if (is.character(mut_file)) {
    mut_file <- file.path(mut_file)
    # Validate file/folder input
    if (!file.exists(mut_file)) {
      stop("Error: The file path you've specified is invalid")
    }
    file_info <- file.info(mut_file)
    if (file_info$isdir == TRUE) {
      # Handle the case where mut_file is a directory
      mut_files <- list.files(path = mut_file, full.names = TRUE, no.. = TRUE)

      if (length(mut_files) == 0) {
        stop("Error: The folder you've specified is empty")
      }
      # Warning if any of the files in folder are empty
      files_info_all <- file.info(mut_files)

      empty_indices <- is.na(files_info_all$size) | files_info_all$size == 0
      empty_list <- basename(mut_files[empty_indices])

      empty_list_str <- paste(empty_list, collapse = ", ")

      if (length(empty_list) == length(mut_files)) {
        stop("Error: All the files in the specified directory are empty")
      }
      if (length(empty_list) != 0) {
        warning(paste("Warning: The following files in the specified directory are empty and will not be imported: ", empty_list_str))
      }

      # Remove empty files from mut_files
      mut_files <- mut_files[!empty_indices]

      # Read in the files and bind them together
      dat <- lapply(mut_files, function(file) {
        read.table(file,
          header = TRUE, sep = mut_sep,
          fileEncoding = "UTF-8-BOM"
        )
      }) %>% dplyr::bind_rows()
    } else {
      # Handle the case where mut_file exists and is a file
      if (file_info$size == 0 || is.na(file_info$size)) {
        stop("Error: You are trying to import an empty file")
      }
      dat <- read.table(mut_file,
        header = TRUE, sep = mut_sep,
        fileEncoding = "UTF-8-BOM"
      )
    }
    if (ncol(dat) <= 1) {
      stop("Your imported data only has one column.
           You may want to set mut_sep to properly reflect
           the delimiter used for the data you are importing.")
    }
  } else {
    stop("Error: mut_file must be a character string or a data frame")
  }
  ## Sample Data File
  # Validate and join sample data file if provided
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

  # Rename columns to default.
  # Add custom column names to default list
  if (!is.null(custom_column_names)) {
    cols <- modifyList(MutSeqR::op$column, custom_column_names)
    dat <- rename_columns(dat, cols)
  } else {
    dat <- rename_columns(dat)
  }
  # Check that all required columns are present
  dat <- check_required_columns(dat, op$base_required_mut_cols)
  context_exists <- "context" %in% colnames(dat)

  # Check for NA values in required columns.
  columns_with_na <- colnames(dat)[apply(dat, 2, function(x) any(is.na(x)))]
  na_columns_required <- intersect(columns_with_na,
                                   MutSeqR::op$base_required_mut_cols)
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
    end.field = "end",
    starts.in.df.are.0based = is_0_based_mut
  )

  if (!is.null(regions)) {
    # load regions file
    regions_gr <- MutSeqR::load_regions_file(regions,
                                             rg_sep,
                                             is_0_based_rg)
    regions_gr$in_regions <- TRUE

    # Apply padding
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
      warning("Warning: ", false_count, " rows were outside of the specified regions. To remove these rows, use the filter_mut() function\n")
    }
  }
  # Create a context column, if needed: BSGenome
  if (!context_exists) {
    if (is.null(genome) || is.null(species)) {
      stop("Error: We need to calculate the context column for your data. Please provide a genome and species so that we can retrieve the appropriate BS genome.")
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
      return(sequences)
    }
    message("Retrieving context sequences from BSgenome")
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

  # Set Depth column as total_depth or depth
  total_depth_exists <- "total_depth" %in% colnames(dat)
  depth_exists <- "depth" %in% colnames(dat)
  no_calls_exists <- "no_calls" %in% colnames(dat)

  if (!total_depth_exists && no_calls_exists && depth_exists) {
    dat <- dat %>%
      dplyr::mutate(total_depth = .data$depth - .data$no_calls)
  }
  if (!total_depth_exists && !no_calls_exists && depth_exists) {
    dat <- dplyr::rename(dat, total_depth = "depth")
    warning("Could not find total_depth column and cannot calculate. Will use depth column as total_depth. Renamed 'depth' to 'total_depth'. Review the differences in the README. \n")
  }
  if (!total_depth_exists && !depth_exists) {
    warning("Could not find an appropriate depth column. Some package functionality may be limited.\n")
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

  if (output_granges) {
    gr <- GenomicRanges::makeGRangesFromDataFrame(
      df = dat,
      keep.extra.columns = TRUE,
      seqnames.field = "contig",
      start.field = "start",
      end.field = "end",
      starts.in.df.are.0based = FALSE
    )
    return(gr)
  } else {
    return(dat)
  }
}
