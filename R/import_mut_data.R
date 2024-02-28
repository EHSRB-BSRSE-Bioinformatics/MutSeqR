#' Import a .mut file
#'
#' Imports a .mut file into the local R environment.
#' @param mut_file "filepath". The .mut file containing mutation
#' data to be imported. If you specify a folder, the function will
#' attempt to read all files in the folder and combine them into
#' a single data frame.Required columns are listed below.
#' Synonymous names for these columns are accepted.
#' \itemize{
#'      \item `contig`: The reference sequence name.
#'      \item `start`: 0-based start position of the feature in contig.
#'      \item `end`: half-open end position of the feature in contig.
#'      \item `sample`: The sample name.
#'      \item `ref`: The reference allele at this position
#'      \item `alt`: The left-aligned, normalized, alternate allele at this
#' position.
#'      \item `alt_depth`: The read depth supporting the alternate allele.
#'      \item depth col: The total read depth at this position. This column can
#' be `total_depth` (excluding N-calls) or `depth`(including N-calls; if
#' `total_depth` is not available).
#'      \item `variation_type`: The category to which this variant is assigned.
#'      \item `context`: The local reference trinucleotide context at this
#' position (e.g. ATC - not necessarily the transcript codon)
#' }
#' @param mut_sep The delimiter for importing the .mut file.
#' Default is tab-delimited.
#' @param rsids A logical variable; whether or not the .mut file
#' contains rsID information (existing SNPs).
#' @param sample_data_file "filepath". An optional file containing
#' additional sample metadata (dose, tissue, timepoint, etc.).
#' @param sd_sep The delimiter for importing sample metadata table.
#' Default is tab-delimited.
#' @param vaf_cutoff The function will add `is_germline` column
#' that identifies ostensibly germline variants using a cutoff
#' for variant allele fraction (VAF). There is no default value
#' provided, but generally a value of 0.1 (i.e., 10%) is a good
#' starting point. Setting this will flag variants that are
#' present at a frequency greater than this value at a given site.
#' @param regions Values are `c("human", "mouse", "custom")`.
#' Indicates the target panel used for Duplex Sequencing.
#' The argument refers to the TS Mutagenesis panel of the
#' specified species, or to a custom panel. If "custom",
#' provide the file path of your regions file in
#' `custom_regions_file`.
#' @param custom_regions_file "filepath". If `regions` is set to "custom",
#' provide the file path for the file containing regions metadata.
#' Required columns are `contig`, `start`, and `end`.
#' @param rg_sep The delimiter for importing the `custom_regions_file`.
#' Default is tab-delimited.
#' @param is_0_based A logical variable. Indicates whether the
#' `custom_regions_file` target region coordinates are 0 based (TRUE)
#' or 1 based (FALSE). If TRUE, ranges will be converted to 1-based.
#' @param range_buffer An integer >= 0. Variants that occur outside of
#' the defined regions' ranges will be filtered out. Use the range-buffer
#' to extend the range within which a variant can occur. The default is
#' 1 nucleotide outside of region ranges. Ex. Structural variants and
#' indels may start outside of the regions. Adjust the range_buffer to
#' include these variants. Rows that were filtered out of the mutation
#' data will be returned in a separate data frame.
#' @param depth_calc Values are `c("take_del", "take_mean")`. In the instance
#' when there are two or more calls at the same location within a sample, and
#' the depths differ, this parameter chooses the method for resolving the
#' difference. This occurs when a deletion is called in the data. It will be
#' called alongside a no_variant. "take_mean" calculates the depth column by
#' taking the mean of all depths in the group. "take_del" calculates the depth
#' column by choosing only the depth of the deletion in the group, or if no
#' deletion is present, the complex variant. If there is no deletion or
#' complex variant, then it takes the mean of the depths within the group.
#' Default is "take_del". depth_col = `total_depth` or `depth`.
#' @param custom_column_names A list of names to specify the meaning of column
#'  headers. Since column names can vary with data, this might be necessary to
#'  digest the mutation data table properly. Typical defaults are set, but can
#'  be substituted in the form of `list(total_depth = "my_custom_depth_name",
#'  sample = "my_custom_sample_column_name")`. For a comprehensive list, see
#'  examples. You can change one or more of these.
#' @param output_granges A logical variable; whether you want the mutation
#' data to output as a GRanges object. Default output (FALSE) is as a dataframe.
#' @returns A table where each row is a mutation, and columns indicate the
#' location, type, and other data. If `output_granges` is set to TRUE, the
#' mutation data will be returned as a GRanges object, otherwise mutation
#' data is returned as a dataframe. If the mutation data contains rows that
#' are filtered out of the specified regions, results will be returned as a
#' list. The first element will be the mutation data called `mut_dat` and
#' the second element will be the rows that were filtered from the data called
#' `rows_outside_regions`. If no rows are filtered out, the results will be
#' returned as a single dataframe or GRanges object.
#'
#'
#' Output Column Definitions:
#' \itemize{
#'      \item `nchar_ref`: The length (in bp) of the reference allele.
#'      \item `nchar_alt`: The length (in bp) of the alternate allele.
#'      \item `varlen`: The length (in bp) of the variant.
#'      \item `total_depth`: The total read depth at this position, excluding
#' N-calls.
#'      \item `vaf`: The variant allele fraction. Calculated as
#' `alt_depth`/`depth_col` where `depth_col` can be `total_depth` or `depth`.
#'      \item `is_germline`: TRUE or FALSE. Flags ostensible germline mutations
#' (`vaf` > `vaf_cutoff`).
#'      \item `ref_depth`: The total read depth at the position calling for the
#' reference allele. Calculated as `depth_col` - `alt_depth` where
#' `depth_col` can be `total_depth`or `depth`.
#'      \item `subtype`: The substitution type for the snv variant
#' (12-base spectrum; e.g. A>C)
#'      \item `short_ref`: The reference base at this position.
#'      \item `normalized_subtype`: The C/T-based substitution type for the
#' snv variant (6-base spectrum; e.g. A>C -> T>G).
#'      \item `normalized_ref`: The reference base in C/T-base notation for
#' this position (e.g. A -> T).
#'      \item `context_with_mutation`: The substitution type fo the snv variant
#' including the two flanking nucleotides (192-trinucleotide spectrum;
#' e.g. `T[A>C]G`)
#'      \item `normalized_context_with_mutation`: The C/T-based substitution
#' type for the snv variant including the two flanking nucleotide
#' (96-base spectrum e.g. `T[A>C]G` -> `C[T>G]A`)
#'      \item `normalized_context`: The trinucleotide context in C/T base
#' notation for this position (e.g. TAG -> CTA).
#'      \item `gc_content`: % GC of the trinucleotide context at this position.
#' }
#' @importFrom dplyr bind_rows mutate left_join case_when
#' @importFrom magrittr %>%
#' @importFrom stringr str_sub str_count
#' @importFrom plyranges join_overlap_left
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom utils read.delim read.table
#' @importFrom rlang .data
#' @export

import_mut_data <- function(mut_file,
                            mut_sep = "\t",
                            rsids = FALSE,
                            sample_data_file = NULL,
                            sd_sep = "\t",
                            vaf_cutoff,
                            range_buffer = 0,
                            regions = c("human", "mouse", "custom"),
                            custom_regions_file = NULL,
                            rg_sep = "\t",
                            is_0_based = TRUE,
                            depth_calc = "take_del",
                            custom_column_names = NULL,
                            output_granges = FALSE) {
  mut_file <- file.path(mut_file)

  # Validate file/folder input
  if (file.exists(mut_file)) {
    file_info <- file.info(mut_file)

    if (file_info$isdir == TRUE) {
      # Handle the case where mut_file exists and is a directory
      mut_files <- list.files(path = mut_file, full.names = TRUE, no.. = TRUE)

      if (length(mut_files) == 0) {
        stop("Error: The folder you've specified is empty")
      }

      # Warning/error if any of the files in folder are empty
      files_info_all <- file.info(mut_files)

      empty_indices <- is.na(files_info_all$size) | files_info_all$size == 0
      empty_list <- basename(mut_files[empty_indices])

      empty_list_str <- paste(empty_list, collapse = ", ")

      if (length(empty_list) == length(mut_files)) {
        stop("Error: All the files in the specified directory are empty")
      }
      if (length(empty_list) != 0) {
        warning(paste("Warning: The following files in the specified
                        directory are empty:", empty_list_str))
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
      # Handle the case where mut_file exists and is not a directory (a file)
      if (file_info$size == 0 || is.na(file_info$size)) {
        stop("Error: You are trying to import an empty file")
      }

      dat <- read.table(mut_file,
        header = TRUE, sep = mut_sep,
        fileEncoding = "UTF-8-BOM"
      )
    }
  } else {
    # Handle the case where mut_file does not exist
    stop("Error: The file path you've specified is invalid")
  }

  if (ncol(dat) <= 1) {
    stop("Your imported data only has one column.
                           You may want to set mut_sep to properly reflect
                           the delimiter used for the data you are importing.")
  }
  if (rsids == TRUE) {
    if (!"id" %in% colnames(dat)) {
      stop("Error: you have set rsids to TRUE,
      but there is no id column in the mut file!")
    }
    # If we have rs IDs, add a column indicating whether the mutation
    # is a known SNP

    dat <- dat %>% dplyr::mutate(is_known = ifelse(!.data$id == ".", "Y", "N"))
  }

  # Rename columns to default .
  # Add custom column names to default list
  if (!is.null(custom_column_names)) {
    cols <- modifyList(MutSeqR::op$column, custom_column_names)
    dat <- rename_columns(dat, cols)
  } else {
    dat <- rename_columns(dat)
  }
  # Check that all required columns are present
  dat <- check_required_columns(dat, c(op$base_required_mut_cols, "context"))

  # Check for NA values
  # If there are NA values in required columns. stop
  columns_with_na <- colnames(dat)[apply(dat, 2, function(x) any(is.na(x)))]
  # Check if there are any NA columns in base_required_mut_cols
  na_columns_required <- intersect(columns_with_na,
                                    c(MutSeqR::op$base_required_mut_cols,
                                      "context"))

  if (length(na_columns_required) > 0) {
    stop(paste(
      "Function stopped: NA values were found within the following required
      column(s): ",
      paste(na_columns_required, collapse = ", "),
      ". Please confirm that your data is complete before proceeding."
    ))
  } else {
    print("Data checked for NA values in required columns: PASS")
  }

  # Read in sample data if it's provided
  #  #Trim and lowercase column headings
  # read.delim check.names = TRUE adds an X
  if (!is.null(sample_data_file)) {
    colnames(dat) <- tolower(gsub("\\.+", "", # deals with middle periods
      gsub(
        "(\\.+)?$", "", # deals with trailing periods
        gsub(
          "^((X\\.+)|(\\.+))?", "", # deals with beginning X. and periods
          colnames(dat)
        )
      ),
      perl = TRUE
    ))

    sampledata <- read.delim(file.path(sample_data_file),
      sep = sd_sep,
      header = TRUE
    )
    # Join
    dat <- dplyr::left_join(dat, sampledata, suffix = c("", ".sampledata"))
  }
  # Clean Up Data
  # Modify variation_type, create varlen and subtype columns.
  # NOTE: cases may occur where there is an indel and ref or alt are 1 bp,
  # but there is a "substitution" as well. Ex GC -> T. Unclear how to define
  # this. Example seen in Jonatan data.
  dat <- dat %>%
    dplyr::mutate(
      nchar_ref = nchar(ref),
      nchar_alt = ifelse(.data$variation_type != "symbolic" |
                           .data$variation_type != "sv", nchar(alt), NA),
      variation_type = tolower(dat$variation_type),
      variation_type =
        ifelse(.data$variation_type == "sv", "symbolic",
          ifelse(.data$variation_type == "indel" &
              .data$nchar_ref > .data$nchar_alt &
              .data$nchar_alt == 1, "deletion",
            ifelse(.data$variation_type == "indel" &
                .data$nchar_ref < .data$nchar_alt &
                .data$nchar_ref == 1, "insertion",
              ifelse(.data$variation_type == "indel" &
                  .data$nchar_ref != .data$nchar_alt &
                  .data$nchar_alt > 1 & .data$nchar_ref > 1, "complex",
                .data$variation_type
              )
            )
          )
        ),
      varlen =
        ifelse(.data$variation_type %in% c("insertion", "deletion", "complex"),
          .data$nchar_alt - .data$nchar_ref,
          ifelse(.data$variation_type %in% c("snv", "mnv"), .data$nchar_ref,
            NA
          )
        ),
      subtype =
        ifelse(.data$variation_type == "snv",
          paste0(.data$ref, ">", .data$alt),
          "."
        )
    )

  # Clean up depth column
  # Set Depth column as total_depth or depth
  has_total_depth <- "total_depth" %in% colnames(dat)
  has_depth <- "depth" %in% colnames(dat)
  has_no_calls <- "no_calls" %in% colnames(dat)

  if (has_total_depth) {
    depth_col <- "total_depth"
  }
  if (!has_total_depth && has_no_calls && has_depth) {
    dat <- dat %>%
      mutate(total_depth = .data$depth - .data$no_calls)
    depth_col <- "total_depth"
  }
  if (!has_total_depth && !has_no_calls && has_depth) {
    depth_col <- "depth"
    warning(" 'total_depth' column was not found and cannot be
            calculated without the 'no_calls' column.'Depth_col'
            will be set as 'depth'. Please review the definitions
            of each column ~here~ (README) before proceeding")
  }
  if (!has_total_depth && !has_depth && !has_no_calls) {
    stop("Required columns are missing or could not be determined:
          depth column ('depth' and 'no_calls' OR 'total_depth')")
  }
  if (!has_total_depth && !has_depth && has_no_calls) {
    stop("Required columns are missing or could not be determined:
          depth column ('depth' OR 'total_depth')")
  }

  # Calculate depth_col for the duplicated rows
  if (depth_calc == "take_del") {
    # Group by sample, contig, and start
    df_grouped <- dat %>%
      dplyr::group_by(.data$sample, .data$contig, .data$start)

    # Define a function to prioritize variation types
    prioritize_depth <- function(variation_type, total_depth) {
      if ("deletion" %in% variation_type) {
        # Choose the total_depth of the first deletion
        total_depth[variation_type == "deletion"][1]
      } else if ("complex" %in% variation_type) {
        # Choose the total_depth of the first complex
        total_depth[variation_type == "complex"][1]
      } else {
        # Choose the total_depth of the first row if no deletion or complex
        total_depth[1]
      }
    }
    # Apply the function to each group
    dat <- df_grouped %>%
      dplyr::mutate(new_depth_col =
                      prioritize_depth(variation_type = .data$variation_type,
                                       .data[[depth_col]])) %>%
      dplyr::ungroup() %>%
      dplyr::select(-{{ depth_col }}) %>%
      dplyr::rename(!!depth_col := "new_depth_col")
  } else if (depth_calc == "take_mean") {
    df_grouped <- dat %>%
      dplyr::group_by(.data$sample, .data$contig, .data$start) %>%
      dplyr::filter(n() > 1) %>%
      dplyr::summarize(
        new_depth_col = round(mean(.data[[depth_col]], na.rm = TRUE))
      ) %>%
      dplyr::ungroup()

    dat <- dat %>%
      dplyr::left_join(df_grouped, by = c("sample", "contig", "start"))

    dat <- dat %>%
      dplyr::mutate(final_depth_col = ifelse(is.na(.data$new_depth_col),
                                              .data[[depth_col]],
                                              .data$new_depth_col)) %>%
      dplyr::select(-{{ depth_col }}, -"new_depth_col") %>%
      dplyr::rename(!!depth_col := "final_depth_col")
  } else {
    stop("Invalid depth_calc input. Please choose 'take_mean' or 'take_del'.")
  }

  # Create VAF, is_germline, and ref_depth columns
  dat <- dat %>%
    dplyr::mutate(VAF = .data$alt_depth / .data[[depth_col]]) %>%
    dplyr::mutate(
      is_germline = ifelse(.data$VAF < vaf_cutoff, FALSE, TRUE),
      ref_depth = .data[[depth_col]] - .data$alt_depth
    )

  # Define substitution dictionary to normalize to pyrimidine context
  sub_dict <- c(
    "G>T" = "C>A", "G>A" = "C>T", "G>C" = "C>G",
    "A>G" = "T>C", "A>C" = "T>G", "A>T" = "T>A"
  )

  # Get reverse complement of sequence context where mutation is
  # listed on purine context
  # Change all purine substitutions to pyrimidine substitutions
  # Make new column with COSMIC-style 96 base context
  # TO DO - describe better
  dat <- dat %>%
    dplyr::mutate(
      context_with_mutation =
        ifelse(.data$subtype != ".",
          paste0(
            stringr::str_sub(.data$context, 1, 1),
            "[", .data$subtype, "]",
            stringr::str_sub(.data$context, 3, 3)
          ),
          .data$variation_type
        ),
      normalized_context = ifelse(
        stringr::str_sub(.data$context, 2, 2) %in% c("G", "A"),
        mapply(function(x) MutSeqR::reverseComplement(x, case = "upper"),
               .data$context),
        .data$context
      ),
      normalized_subtype = ifelse(
        .data$subtype %in% names(sub_dict),
        sub_dict[.data$subtype],
        .data$subtype
      ),
      short_ref = substr(.data$ref, 1, 1),
      normalized_ref = dplyr::case_when(
        substr(.data$ref, 1, 1) == "A" ~ "T",
        substr(.data$ref, 1, 1) == "G" ~ "C",
        substr(.data$ref, 1, 1) == "C" ~ "C",
        substr(.data$ref, 1, 1) == "T" ~ "T"
      )
    ) %>%
    dplyr::mutate(
      normalized_context_with_mutation =
        ifelse(.data$subtype != ".",
          paste0(
            stringr::str_sub(.data$normalized_context, 1, 1),
            "[", .data$normalized_subtype, "]",
            stringr::str_sub(.data$normalized_context, 3, 3)
          ),
          .data$variation_type
        ),
      gc_content = (stringr::str_count(string = .data$context, pattern = "G") +
                    stringr::str_count(string = .data$context, pattern = "C"))
      / stringr::str_count(.data$context)
    ) %>%
    dplyr::mutate(
      normalized_subtype = ifelse(
        .data$normalized_subtype == ".",
        .data$variation_type,
        .data$normalized_subtype
      ),
      subtype = ifelse(
        .data$subtype == ".",
        .data$variation_type,
        .data$subtype
      )
    )

  # Turn into GRanges
  mut_ranges <- GenomicRanges::makeGRangesFromDataFrame(
    df = as.data.frame(dat),
    keep.extra.columns = TRUE,
    seqnames.field = "contig",
    start.field = "start",
    end.field = "end",
    starts.in.df.are.0based = TRUE
  )

  regions_df <- load_regions_file(regions, custom_regions_file, rg_sep)

  regions_df <- regions_df %>%
    dplyr::mutate(regions_start_buffered =
                    as.numeric(paste0(.data$start)) - range_buffer)

  if (regions %in% c("mouse", "human", "rat")) {
    is_0_based <- TRUE
  }

  if (is_0_based) {
    regions_df$regions_start_buffered <- regions_df$regions_start_buffered + 1
  }

  region_ranges <- GenomicRanges::makeGRangesFromDataFrame(
    df = regions_df,
    keep.extra.columns = TRUE,
    seqnames.field = "contig",
    start.field = "start",
    end.field = "end",
    starts.in.df.are.0based = is_0_based
  )

  ranges_joined <- plyranges::join_overlap_left(mut_ranges,
    region_ranges,
    maxgap = range_buffer,
    suffix = c("_mut", "_regions")
  )
  dat <- as.data.frame(ranges_joined)
  dat <- dat %>%
    dplyr::mutate(bp_outside_rg = .data$regions_start_buffered - .data$start)
  ranges_outside_regions <- dat %>%
    dplyr::filter(is.na(.data$bp_outside_rg) | .data$bp_outside_rg > 0) %>%
    dplyr::select("sample", "seqnames", "start", "end", "ref", "alt")

  # Display the ranges that were filtered out of the data
  if (nrow(ranges_outside_regions) > 0) {
    print(paste(
      nrow(ranges_outside_regions),
      "rows of data were outside specified regions and were filtered out.
      Mutation data and filtered rows will be returned in seperate
      data frames."
    ))
  }

  dat <- dat %>%
    dplyr::filter(!is.na(.data$bp_outside_rg) & .data$bp_outside_rg <= 0) %>%
    plyranges::select(-"regions_start_buffered", -"bp_outside_rg")

  if (output_granges) {
    gr <- GenomicRanges::makeGRangesFromDataFrame(
      df = dat,
      keep.extra.columns = TRUE,
      seqnames.field = "seqnames",
      start.field = "start",
      end.field = "end",
      starts.in.df.are.0based = FALSE
    )
    if (nrow(ranges_outside_regions) > 0) {
      ls <- list(
        mut_dat = gr,
        rows_outside_regions = ranges_outside_regions
      )
      return(ls)
    } else {
      return(gr)
    }
  } else {
    dat <- dat %>%
      dplyr::rename(contig = "seqnames")

    if (nrow(ranges_outside_regions) > 0) {
      ls <- list(
        mut_dat = dat,
        rows_outside_regions = ranges_outside_regions
      )
      return(ls)
    } else {
      return(dat)
    }
  }
}
