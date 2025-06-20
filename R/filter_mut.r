#' Filter your mutation data
#' @description This function creates a `filter_mut`` column that will be read
#' by the \code{calculate_mf} function and other downstream functions.
#' Variants with `filter_mut == TRUE`` will be excluded from group mutation
#' counts. This function may also remove records upon on user specification.
#' Running this function again on the same data will not overide the previous
#' filters. To reset previous filters, set the filter_mut column values to
#' FALSE.
#' @param mutation_data Your mutation data. This can be a data frame or a
#' GRanges object.
#' @param vaf_cutoff Filter out ostensibly germline variants using a cutoff for
#' variant allele fraction (VAF). Any variant with a \code{vaf} larger than
#' the cutoff will be filtered. The default is 1 (no filtering). It is
#' recommended to use a value of 0.01 (i.e. 1%) as a conservative approach
#' to retain only somatic variants.
#' @param snv_in_germ_mnv Filter out snv variants that overlap with
#' germline mnv variants within the same samples. mnv variants will be
#' considered germline if their vaf > vaf_cutoff. Default is FALSE.
#' @param rm_abnormal_vaf A logical value. If TRUE, rows in
#' \code{mutation_data} with a variant allele fraction (VAF) between 0.05 and
#' 0.45 or between 0.55 and 0.95 will be removed. We expect variants to have a
#' VAF ~0. 0.5, or 1, reflecting rare somatic mutations, heterozygous germline
#' mutations, and homozygous germline mutations, respectively. Default is
#' FALSE.
#' @param custom_filter_col The name of the column in mutation_data to apply a
#' custom filter to. This column will be checked for specific values, as defined
#' by \code{custom_filter_val}. If any row in this column contains one of the
#' specified values, that row will either be flagged in the
#' \code{filter_mut column} or, if specified by \code{custom_filter_rm},
#' removed from mutation_data.
#' @param custom_filter_val A set of values used to filter rows in
#' \code{mutation_data} based on \code{custom_filter_col}. If a row in
#' \code{custom_filter_col} matches any value in \code{custom_filter_val},
#' it will either be set to TRUE in the \code{filter_mut} column or removed,
#' depending on \code{custom_filter_rm}.
#' @param custom_filter_rm A logical value. If TRUE, rows in custom_filter_col
#' that match any value in custom_filter_val will be removed from the
#' mutation_data. If FALSE, \code{filter_mut} will be set to TRUE for those
#' rows.
#' @param regions Remove rows that are within/outside of specified regions.
#' `regions` can be either a file path, a data frame, or a GRanges object
#' containing the genomic ranges by which to filter. File paths will be read
#' using the rg_sep. Users can also choose from the built-in TwinStrand's
#' Mutagenesis Panels by inputting "TSpanel_human",  "TSpanel_mouse", or
#' "TSpanel_rat". Required columns for the regions file are "contig", "start",
#' and "end". In a GRanges object, the required columns are "seqnames",
#' "start", and "end".
#' @param regions_filter Specifies how the provided \code{regions} should be
#' applied to \code{mutation_data}. Acceptable values are "remove_within" or
#' "keep_within". If set to "remove_within", records that fall within the
#' specified regions wil be removed from mutation_data. If set to
#' "keep_within", only records within the specified regions will be kept in
#' mutation_data, and all other records will be removed.
#' @param allow_half_overlap A logical value. If TRUE, records that start or
#' end in your \code{regions}, but extend outside of them in either direction
#' will be included in the filter. If FALSE, only records that start and end
#' within the \code{regions} will be included in the filter. Default is FALSE.
#' @param rg_sep The delimiter for importing the custom_regions. The default is
#' tab-delimited "\\t".
#' @param is_0_based_rg A logical variable. Indicates whether the position
#' coordinates in `regions` are 0 based (TRUE) or 1 based (FALSE).
#' If TRUE, positions will be converted to 1-based (start + 1).
#' Need not be supplied for TSpanels. Default is TRUE.
#' @param rm_filtered_mut_from_depth A logical value. If TRUE, the function will
#' subtract the \code{alt_depth} of records that were flagged by the
#' \code{filter_mut} column from their \code{total_depth}. This will treat
#' flagged variants as No-calls. This will not apply to variants flagged as
#' germline by the \code{vaf_cutoff}. However, if the germline variant
#' has additional filters applied, then the subtraction will still occur.
#' If FALSE, the \code{alt_depth} will be retained in the
#' \code{total_depth} for all variants.  Default is FALSE.
#' @param return_filtered_rows A logical value. If TRUE, the function will
#' return both the filtered mutation data and the records that were
#' removed/flagged in a seperate data frame. The two dataframes will be
#' returned inside a list, with names \code{mutation_data} and
#' \code{filtered_rows}. Default is FALSE.
#' @examples
#' # Load example data
#' example_file <- system.file("extdata", "Example_files",
#'                             "example_mutation_data.rds",
#'                             package = "MutSeqR")
#' example_data <- readRDS(example_file)
#' # Filter the data
#' # Basic Usage: Filter out putative germline variants
#' filter_example_1 <- filter_mut(mutation_data = example_data,
#'                                vaf_cutoff = 0.01)
#' # Remove rows outside of the TwinStand Mouse Mutagenesis Panel regions
#' filter_example_2 <- filter_mut(mutation_data = example_data,
#'                                vaf_cutoff = 0.01,
#'                                regions = "TSpanel_mouse",
#'                                regions_filter = "keep_within")
#' # Apply a custom filter to flag rows with "EndRepairFillInArtifact"
#' # in the column 'filter'
#' filter_example_3 <- filter_mut(mutation_data = example_data,
#'                                vaf_cutoff = 0.01,
#'                                regions = "TSpanel_mouse",
#'                                regions_filter = "keep_within",
#'                                custom_filter_col = "filter",
#'                                custom_filter_val = "EndRepairFillInArtifact",
#'                                custom_filter_rm = FALSE)
#' # Flag snv variants that overlap with germline mnv variants.
#' # Subtract the alt_depth of these variants from their total_depth
#' # (treat them as No-calls).
#' # Return all the flagged/removed rows in a seperate data frame
#' filter_example_4 <- filter_mut(mutation_data = example_data,
#'                                vaf_cutoff = 0.01,
#'                                regions = "TSpanel_mouse",
#'                                regions_filter = "keep_within",
#'                                custom_filter_col = "filter",
#'                                custom_filter_val = "EndRepairFillInArtifact",
#'                                custom_filter_rm = FALSE,
#'                                snv_in_germ_mnv = TRUE,
#'                                rm_filtered_mut_from_depth = TRUE,
#'                                return_filtered_rows = TRUE)
#' # Flagging germline mutations...
#' # Found 612 germline mutations.
#' # Flagging SNVs overlapping with germline MNVs...
#' # Found 20 SNVs overlapping with germline MNVs.
#' # Applying custom filter...
#' # Flagged 2021 rows with values in <filter> column that matched EndRepairFillInArtifact
#' # Applying region filter...
#' # Removed 22 rows based on regions.
#' # Correcting depth...
#' # 909 rows had their total_depth corrected.
#' # Removing filtered mutations from the total_depth...
#' # Filtering complete.
#' # Returning a list: mutation_data and filtered_rows.
#' filtered_rows <- filter_example_4$filtered_rows
#' filtered_example_mutation_data <- filter_example_4$mutation_data
#' @importFrom dplyr group_by mutate ungroup select filter starts_with
#' n_distinct first case_when if_else
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom S4Vectors queryHits mcols
#' @importFrom plyranges join_overlap_left_directed join_overlap_left_within_directed
#' @export
filter_mut <- function(mutation_data,
                       vaf_cutoff = 1,
                       snv_in_germ_mnv = FALSE,
                       # SNV GERM INDELS by 10bp
                       rm_abnormal_vaf = FALSE,
                       custom_filter_col = NULL,
                       custom_filter_val = NULL,
                       custom_filter_rm = FALSE,
                       regions = NULL,
                       regions_filter,
                       allow_half_overlap = FALSE,
                       rg_sep = "\t",
                       is_0_based_rg = TRUE,
                       rm_filtered_mut_from_depth = FALSE,
                       return_filtered_rows = FALSE) {
# clonality_cutoff NOT CURRENTLY IMPLEMENTED! Up for consideration.
# This value determines the fraction of reads that
# is considered a constitutional variant. If a mutation is present at a
# fraction higher than this value, the reference base will be swapped,
# and the alt_depth recalculated. 0.3 (30%) would be a sane default?

  if (return_filtered_rows) {
    rm_rows <- data.frame()
  }
  if (inherits(mutation_data, "GRanges")) {
    mutation_data <- as.data.frame(mutation_data)
    mutation_data <- dplyr::rename(mutation_data, contig = seqnames)
  }

  # Create a new Filter Column
  if (!("filter_mut" %in% colnames(mutation_data))) {
    mutation_data$filter_mut <- FALSE
  }
  if (!is.logical(mutation_data$filter_mut)) {
    mutation_data$filter_mut <- as.logical(mutation_data$filter_mut)
    if (any(is.na(mutation_data$filter_mut))) {
      stop("NAs or non-logical values were found in the filter_mut column")
    }
  }
  if (!("filter_reason" %in% colnames(mutation_data))) {
    mutation_data$filter_reason <- ""
  }

  # Quality check for depth TO ADD:
  ## read depth is not abnormally high or low (corrected for regional GC content)
  ## total_depth is within 1.5 fold of local mean
  ## total depth ratio is >0.95 (5% of all reads were excluded as no calls)

  if (vaf_cutoff < 0 || vaf_cutoff > 1) {
    stop("Error: The VAF cutoff must be between 0 and 1")
  }

  ######## VAF Filter #########################################################
  if (vaf_cutoff < 1) {
    if (!("vaf" %in% colnames(mutation_data))) {
      stop("\nError: You have set a vaf_cutoff but there is no 'vaf' column in
      your mutation_data. vaf = alt_depth/total_depth or vaf = alt_depth/depth,
      if the total_depth column is not available.")
    }
    message("Flagging germline mutations...")
    mutation_data <- mutation_data %>%
      dplyr::mutate(
        filter_mut = ifelse(.data$vaf > vaf_cutoff, TRUE, .data$filter_mut),
        is_germline = ifelse(.data$vaf > vaf_cutoff, TRUE, FALSE),
        filter_reason = ifelse(.data$vaf > vaf_cutoff,
          ifelse(.data$filter_reason == "", "germline",
            paste0(.data$filter_reason, "|germline")
          ), .data$filter_reason
        )
      )
    vaf_filtered_count <- sum(as.numeric(mutation_data$is_germline == TRUE))
    message("Found ", vaf_filtered_count, " germline mutations.")

    ###### snv_in_germ_mnv Filter #############################################
    if (snv_in_germ_mnv) {
      message("Flagging SNVs overlapping with germline MNVs...")
      mnv_ranges <- mutation_data %>%
        dplyr::filter(.data$variation_type == "mnv" & .data$is_germline == TRUE)

      hits_indices <- integer(0)

      # Only proceed if there are any germline MNVs to check against
      if (nrow(mnv_ranges) > 0) {
        for (sample_name in unique(mutation_data$sample)) {
          snv_subset <- mutation_data %>%
            dplyr::filter(.data$sample == sample_name, .data$variation_type == "snv")
          germ_mnv_subset <- mnv_ranges %>%
            dplyr::filter(.data$sample == sample_name)

          # If there are no SNVs or no germline MNVs for this sample, skip to the next
          if (nrow(snv_subset) == 0 || nrow(germ_mnv_subset) == 0) {
            next
          }
          snv_gr <- GenomicRanges::makeGRangesFromDataFrame(snv_subset,
                                                            seqnames.field = "contig",
                                                            keep.extra.columns = TRUE)
          germ_mnv_gr <- GenomicRanges::makeGRangesFromDataFrame(germ_mnv_subset,
                                                                  seqnames.field = "contig",
                                                                  keep.extra.columns = TRUE)
          overlaps <- GenomicRanges::findOverlaps(query = snv_gr, subject = germ_mnv_gr)
          hits <- S4Vectors::queryHits(overlaps)

          # Ensure 'hits' is not empty before using it as an index
          if (length(hits) > 0) {
                original_indices <- which(mutation_data$sample == sample_name & mutation_data$variation_type == "snv")
                hits_indices <- c(hits_indices, original_indices[hits])
          }
        }
      } # end of if(nrow(mnv_ranges) > 0)
      mutation_data$snv_in_germ_mnv <- FALSE
      # Check if hits_indices has any values before trying to index
      if (length(hits_indices) > 0) {
        mutation_data$snv_in_germ_mnv[hits_indices] <- TRUE
      }

      mutation_data <- mutation_data %>%
        dplyr::mutate(filter_mut = ifelse(.data$snv_in_germ_mnv == TRUE,
                                          TRUE, .data$filter_mut),
                      filter_reason = ifelse(.data$snv_in_germ_mnv == TRUE,
                                             ifelse(.data$filter_reason == "", "snv_in_germ_mnv",
                                                    paste0(.data$filter_reason, "|snv_in_germ_mnv")),
                                            .data$filter_reason))
      snv_in_germ_mnv_count <- sum(mutation_data$snv_in_germ_mnv == TRUE)
      message("Found ", snv_in_germ_mnv_count, " SNVs overlapping with germline MNVs.")
    }
  }
  ######## rm_abnormal_vaf Filter #############################################
  if (rm_abnormal_vaf) {
    if (!("vaf" %in% colnames(mutation_data))) {
      stop("\nError: You have set rm_abnormal_vaf to TRUE but there is no 'vaf'
      column in your mutation_data. vaf = alt_depth/total_depth or vaf =
      alt_depth/depth, if the total_depth column is not available.")
    }
    message("Removing rows with abnormal VAF...")
    original_row_count <- nrow(mutation_data)

    if (return_filtered_rows) {
      rm_abnormal_vaf <- mutation_data %>%
        dplyr::filter((.data$vaf > 0.05 & .data$vaf < 0.45) | (.data$vaf > 0.55 & .data$vaf < 0.95)) %>%
        dplyr::mutate(filter_reason = ifelse(.data$filter_reason == "",
                                             "abnormal_vaf",
                                             paste0(filter_reason, "|abnormal_vaf")))

      rm_rows <- rbind(rm_rows, rm_abnormal_vaf)
    }
    mutation_data <- mutation_data %>%
      dplyr::filter(!(.data$vaf > 0.05 & .data$vaf < 0.45) & !(.data$vaf > 0.55 & .data$vaf < 0.95))
    abnormal_vaf_count <- original_row_count - nrow(mutation_data)
    message("Removed ", abnormal_vaf_count, " rows with abnormal VAF.")
  }

  ######## Custom Filter ######################################################
  if (!is.null(custom_filter_col)) {
    message("Applying custom filter...")
    if (is.null(custom_filter_val)) {
      stop("Error: You provided a custom filter column but did not specify the
      filter value(s). Please provide the value(s) within the custom filter
      column that should be used to apply the filter")
    }
    if (!(custom_filter_col) %in% colnames(mutation_data)) {
      stop(paste("Error: could not find", custom_filter_col, "in mutation_data"))
    }
    pattern <- paste(custom_filter_val, collapse = "|")
    custom_filtered_rows <- grepl(pattern, mutation_data[[custom_filter_col]])
    custom_filtered_count <- sum(custom_filtered_rows)
    if (custom_filter_rm) {
      if (return_filtered_rows) {
        rm_custom <- mutation_data[custom_filtered_rows, ]
        matching_filter_values <- custom_filter_val[sapply(custom_filter_val, function(val) grepl(val, rm_custom[[custom_filter_col]]))]
        rm_custom <- rm_custom %>%
          dplyr::mutate(filter_reason =
            ifelse(.data$filter_reason == "",
              .data[[custom_filter_col]],
              paste0(.data$filter_reason, "|", .data[[custom_filter_col]])
            )
          )
        rm_rows <- rbind(rm_rows, rm_custom)
      }
      mutation_data <- mutation_data[!custom_filtered_rows, ]
      message("Removed ", custom_filtered_count, " rows with values in <", custom_filter_col, "> that contained ",
             custom_filter_val, " from mutation_data")
    } else {
      mutation_data$filter_mut[custom_filtered_rows] <- TRUE

      mutation_data <- mutation_data %>%
        dplyr::mutate(filter_reason = ifelse(custom_filtered_rows,
          ifelse(.data$filter_reason == "",
            .data[[custom_filter_col]],
            paste0(.data$filter_reason, "|", .data[[custom_filter_col]])
          ),
          .data$filter_reason
        ))
      message("Flagged ", custom_filtered_count, " rows with values in <", custom_filter_col, "> column that matched ", custom_filter_val)
    }
  }
  ######## Regions Filter #####################################################
  if (!is.null(regions)) {
    message("Applying region filter...")
    regions_gr <- MutSeqR::load_regions_file(regions, rg_sep, is_0_based_rg)
    regions_gr$in_regions <- TRUE
    colnames(S4Vectors::mcols(regions_gr)) <- paste0("TO_REMOVE_", colnames(S4Vectors::mcols(regions_gr)))

    mut_ranges <- GenomicRanges::makeGRangesFromDataFrame(
      df = as.data.frame(mutation_data),
      keep.extra.columns = TRUE,
      seqnames.field = "contig",
      start.field = "start",
      end.field = "end",
      starts.in.df.are.0based = FALSE
    )
    if (allow_half_overlap) {
      ranges_joined <- plyranges::join_overlap_left_directed(mut_ranges,
                                                             regions_gr,
                                                             suffix = c("", "_regions"))
    } else {
      ranges_joined <- plyranges::join_overlap_left_within_directed(mut_ranges, regions_gr, suffix = c("", "_regions"))
    }
    mutation_data <- as.data.frame(ranges_joined)
    mutation_data$TO_REMOVE_in_regions[is.na(mutation_data$TO_REMOVE_in_regions)] <- FALSE
    mutation_data <- mutation_data %>%
      dplyr::rename(contig = seqnames)
    original_row_count <- nrow(mutation_data)
    if (regions_filter == "remove_within") {
      if (return_filtered_rows) {
        rm_regions <- mutation_data %>%
          dplyr::filter(.data$TO_REMOVE_in_regions == TRUE) %>%
          dplyr::mutate(filter_reason = ifelse(.data$filter_reason == "", "regions",
                                               paste0(.data$filter_reason, "|regions"))) %>%
          dplyr::select(-dplyr::starts_with("TO_REMOVE_"))
        rm_rows <- rbind(rm_rows, rm_regions)
      }
      mutation_data <- mutation_data %>%
        dplyr::filter(.data$TO_REMOVE_in_regions == FALSE)
    } else if (regions_filter == "keep_within") {
      if (return_filtered_rows) {
        rm_regions <- mutation_data %>%
          dplyr::filter(.data$TO_REMOVE_in_regions == FALSE) %>%
          dplyr::mutate(filter_reason = ifelse(.data$filter_reason == "", "regions",
                                               paste0(.data$filter_reason, "|regions"))) %>%
          dplyr::select(-dplyr::starts_with("TO_REMOVE_"))
        rm_rows <- rbind(rm_rows, rm_regions)
      }

      mutation_data <- mutation_data %>%
        dplyr::filter(.data$TO_REMOVE_in_regions == TRUE)
    } else {
      stop("regions_filter must be either 'remove_within' or 'keep_within'.")
    }
    mutation_data <- mutation_data %>%
      dplyr::select(-dplyr::starts_with("TO_REMOVE_"))
    region_filtered_count <- original_row_count - nrow(mutation_data)
    message("Removed ", region_filtered_count, " rows based on regions.")
  }

  if (rm_filtered_mut_from_depth) {
    message("Removing filtered mutations from the total_depth...")
    mutation_data <- mutation_data %>%
      dplyr::mutate(total_depth =
        dplyr::if_else(.data$filter_mut &
                       .data$filter_reason != "germline" &
                       .data$total_depth != 0,
                       .data$total_depth - .data$alt_depth,
                       .data$total_depth)
      )
  }
  message("Filtering complete.")
  if (return_filtered_rows) {
    filtered_muts <- mutation_data %>%
      dplyr::filter(filter_mut == TRUE)
    filter_rows_return <- rbind(rm_rows, filtered_muts)
    message("Returning a list: mutation_data and filtered_rows.")
    return(list(mutation_data = mutation_data, filtered_rows = filter_rows_return))
  } else {
    return(mutation_data)
  }
}
