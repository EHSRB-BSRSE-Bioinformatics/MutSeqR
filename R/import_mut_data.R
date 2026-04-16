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
#' @param BS_genome The pkgname of a BS genome. A BS genome must be installed
#' prior to import to populate the context column (trinucleotide context for
#' each position). Only required if data does not already include a context
#' column. Please install the appropriate BS genome using
#' BiocManager::install("pkgname") where pkgname is the name of the BSgenome
#' package. The pkgname can be found using the find_BS_genome() function, which
#' requires the species and assembly version. Ex."BSgenome.Hsapiens.UCSC.hg38"
#' | "BSgenome.Hsapiens.UCSC.hg19" | "BSgenome.Mmusculus.UCSC.mm10" |
#' "BSgenome.Mmusculus.UCSC.mm39" | "BSgenome.Rnorvegicus.UCSC.rn6".
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
#' # Mutation data is just for example purposes. It does not reflect real data
#' file <- system.file("extdata", "Example_files",
#'                    "simple_mut_import.txt", package = "MutSeqR")
#' # Import the data
#' imported_example_data <- import_mut_data(mut_file = file)
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
#' @importFrom Seqinfo seqnames
#' @importFrom BSgenome getBSgenome installed.genomes
#' @export
import_mut_data <- function(
  mut_file,
  mut_sep = "\t",
  is_0_based_mut = TRUE,
  sample_data = NULL,
  sd_sep = "\t",
  regions = NULL,
  rg_sep = "\t",
  is_0_based_rg = TRUE,
  padding = 0,
  BS_genome = NULL,
  custom_column_names = NULL,
  output_granges = FALSE
) {
  stopifnot(
    "mut_file is required" = !missing(mut_file),
    "mut_file must be a character indicating a filepath or a data frame" = is.character(
      mut_file
    ) ||
      is.data.frame(mut_file),
    "mut_sep must be a character string" = is.character(mut_sep),
    "is_0_based_mut must be a logical variable" = is.logical(is_0_based_mut),
    "sample_data must be NULL, a character indicating a filepath, or a data frame" = is.null(
      sample_data
    ) ||
      is.character(sample_data) ||
      is.data.frame(sample_data),
    "sd_sep must be a character string" = is.character(sd_sep),
    "regions must be NULL, a character indicating a filepath, a data frame, or a GRanges object" = is.null(
      regions
    ) ||
      is.character(regions) ||
      is.data.frame(regions) ||
      methods::is(regions, "GRanges"),
    "rg_sep must be a character string" = is.character(rg_sep),
    "is_0_based_rg must be a logical variable" = is.logical(is_0_based_rg),
    "padding must be a non-negative integer" = is.numeric(padding) &&
      padding >= 0 &&
      (padding %% 1 == 0),
    "BS_genome must be NULL or a character string" = is.null(BS_genome) ||
      is.character(BS_genome),
    "custom_column_names must be NULL or a list" = is.null(
      custom_column_names
    ) ||
      is.list(custom_column_names),
    "output_granges must be a logical variable" = is.logical(output_granges)
  )
  BS_genome <- match.arg(
    BS_genome,
    choices = c(
      NULL,
      BSgenome::available.genomes(splitNameParts = TRUE)$pkgname
    )
  )

  # Check if BS_genome is actually installed locally
  if (!is.null(BS_genome)) {
    if (!(BS_genome %in% BSgenome::installed.genomes())) {
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

  # Load and validate sample metadata before heavy lifting
  sample_df <- NULL
  if (!is.null(sample_data)) {
    sample_df <- import_sample_data(sample_data, sd_sep)
  }

  # Import the mut files: data frame or file path
  if (is.data.frame(mut_file)) {
    dat <- mut_file
    if (nrow(dat) == 0) {
      stop("The data frame you've provided is empty")
    }
  } else if (is.character(mut_file)) {
    mut_file <- file.path(mut_file)
    # Validate file/folder input
    if (!file.exists(mut_file)) {
      stop("The file path you've specified is invalid")
    }
    file_info <- file.info(mut_file)
    if (file_info$isdir == TRUE) {
      # Handle the case where mut_file is a directory
      mut_files <- list.files(path = mut_file, full.names = TRUE, no.. = TRUE)

      if (length(mut_files) == 0) {
        stop("The folder you've specified is empty")
      }
      # Warning if any of the files in folder are empty
      files_info_all <- file.info(mut_files)

      empty_indices <- is.na(files_info_all$size) | files_info_all$size == 0
      empty_list <- basename(mut_files[empty_indices])

      empty_list_str <- paste(empty_list, collapse = ", ")

      if (length(empty_list) == length(mut_files)) {
        stop("All the files in the specified directory are empty")
      }
      if (length(empty_list) != 0) {
        warning(
          "The following files in the specified directory are empty and will not be imported: ",
          empty_list_str
        )
      }

      # Remove empty files from mut_files
      mut_files <- mut_files[!empty_indices]

      # Read in the files and bind them together
      dat <- lapply(mut_files, function(file) {
        read.table(
          file,
          header = TRUE,
          sep = mut_sep,
          fileEncoding = "UTF-8-BOM"
        )
      }) %>%
        dplyr::bind_rows()
    } else {
      # Handle the case where mut_file exists and is a file
      if (file_info$size == 0 || is.na(file_info$size)) {
        stop("You are trying to import an empty file")
      }
      dat <- read.table(
        mut_file,
        header = TRUE,
        sep = mut_sep,
        fileEncoding = "UTF-8-BOM"
      )
    }
    if (ncol(dat) <= 1) {
      stop(
        "Your imported data only has one column.
           You may want to set mut_sep to properly reflect
           the delimiter used for the data you are importing."
      )
    }
  }

  # Rename columns to default (including custom names)
  if (!is.null(custom_column_names)) {
    cols <- modifyList(MutSeqR::op$column, custom_column_names)
    dat <- rename_columns(dat, cols)
  } else {
    dat <- rename_columns(dat)
  }

  ## Join with sample metadata if provided
  if (!is.null(sample_df)) {
    if (!"sample" %in% colnames(dat)) {
      stop(
        "Error in mutation data: 'sample' column is missing prior to joining sample metadata."
      )
    }
    dat <- dplyr::left_join(
      dat,
      sample_df,
      by = "sample",
      suffix = c("", ".sd")
    )
    message("Sample metadata successfully joined to mutation data\n")
  }

  # Check that all required columns are present
  dat <- check_required_columns(dat, op$base_required_mut_cols)
  context_exists <- "context" %in% colnames(dat)

  # Check for NA values in required columns.
  columns_with_na <- colnames(dat)[apply(dat, 2, function(x) any(is.na(x)))]
  na_columns_required <- intersect(
    columns_with_na,
    MutSeqR::op$base_required_mut_cols
  )
  if (length(na_columns_required) > 0) {
    stop(
      "NA values were found within the following required column(s): ",
      paste(na_columns_required, collapse = ", "),
      ". Please confirm that your data is complete before proceeding."
    )
  }
  # Check for NA values in the context column. If so, will populate it.
  if (context_exists) {
    if ("context" %in% columns_with_na) {
      context_exists <- FALSE
    }
  }

  # Turn mutation data into GRanges
  mut_ranges <- GenomicRanges::makeGRangesFromDataFrame(
    df = as.data.frame(dat),
    keep.extra.columns = TRUE,
    seqnames.field = "contig",
    start.field = "start",
    end.field = "end",
    starts.in.df.are.0based = is_0_based_mut
  )
  # Join Regions Metadata
  if (!is.null(regions)) {
    mut_ranges <- import_regions_metadata(
      mutation_granges = mut_ranges,
      regions = regions,
      rg_sep = rg_sep,
      is_0_based_rg = is_0_based_rg,
      padding = padding
    )
  }
  # Populate Context (if not present)
  if (!context_exists) {
    mut_ranges <- populate_sequence_context(
      mutation_granges = mut_ranges,
      BS_genome = BS_genome
    )
  }

  # Characterize variants
  dat <- as.data.frame(mut_ranges) %>%
    dplyr::rename(contig = "seqnames")
  dat <- characterize_variants(dat)

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
    warning(
      "Could not find total_depth column and cannot calculate. Will use depth column as total_depth. Renamed 'depth' to 'total_depth'. Review the differences in the README. \n"
    )
  }
  if (!total_depth_exists && !depth_exists) {
    warning(
      "Could not find an appropriate depth column. Some package functionality may be limited.\n"
    )
  }

  # Check for duplicated rows
  dat <- dat %>%
    dplyr::group_by(.data$sample, .data$contig, .data$start) %>%
    dplyr::mutate(row_has_duplicate = dplyr::n() > 1) %>%
    dplyr::ungroup()

  if (sum(dat$row_has_duplicate) > 0) {
    warning(
      sum(dat$row_has_duplicate),
      " rows were found whose position was the same as that of at least one other row for the same sample."
    )

    # Warn about the depth for the duplicated rows
    if ("total_depth" %in% colnames(dat)) {
      warning(
        "The total_depth may be double-counted in some instances due to overlapping positions. Set the correct_depth parameter in calculate_mf() to correct the total_depth for these instances."
      )
    }
  }

  # Make VAF and ref_depth columns, if depth exists
  if ("total_depth" %in% colnames(dat)) {
    dat <- dat %>%
      dplyr::mutate(
        vaf = .data$alt_depth / .data$total_depth,
        ref_depth = .data$total_depth - .data$alt_depth
      )
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
