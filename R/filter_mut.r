#' Filter your mutation data
#' @description This function creates a filter_mut column that will be read by
#' the \code{calculate_mut_freq} function. Variants with filter == TRUE will
#' not be included in final mutation counts. This function may also remove
#' records of given loci from the mutation data based on user specification.
#' Running this function again on the same data will not overide the previous
#' filters. To reset previous filters, set the filter_mut column values to
#' FALSE.
#' @param mutation_data Your mutation data.
#' @param correct_depth A logical value. If TRUE, the function will correct the
#' \code{total_depth} column in \code{mutation_data} in order to prevent
#' double-counting the \code{total_depth} values for the same genomic position.
#' For rows with the same sample contig, and start values, the \code{total_depth}
#' will be retained for only one row. All other rows in the group will have their
#' \code{total_depth} set to 0. The default is FALSE
#' @param correct_depth_to_indel A logical value. If TRUE, during depth
#' correction, should there be different \code{total_depth} values within a
#' group of rows with the same sample, contig, and start values, the
#' \code{total_depth} value for the row with the highest priority
#' \code{variation_type} will be retained, while the other rows will have their
#' \code{total_depth} set to 0. \code{variation_type} priority order is:
#' deletion, complex, insertion, snv, mnv, sv, uncategorised, no_variant.
#' If FALSE, the \code{total_depth} value for the first row in the group will
#' be retained, while the other rows will have their \code{total_depth} set to
#' 0. The default is TRUE.
#'@param vaf_cutoff Filter out ostensibly germline variants using a cutoff for
#' variant allele fraction (VAF). Any variant with a \code{vaf} larger than
#' the cutoff will be filtered. The default is 1 (no filtering). It is
#' recommended to use a value of 0.01 (i.e. 1%) to retain only somatic
#' variants.
#' @param snv_in_germ_mnv Filter out snv variants that overlap with
#' germline mnv variants within the same samples. mnv variants will be
#' considered germline if their vaf > vaf_cutoff. Default is FALSE.
#' @param custom_filter_col The name of the column in mutation_data to apply a
#' custom filter to. This column will be checked for specific values, as defined
#' by \code{custom_filter_val}. If any row in this column contains one of the
#' specified values, that row will either be flagged in the
#' \code{filter_mut column} or, if specified by \code{custom_filter_rm},
#' removed from mutation_data.
#' @param rm_abnormal_vaf A logical value. If TRUE, rows in
#' \code{mutation_data} with a variant allele fraction (VAF) between 0.05 and
#' 0.45 or between 0.55 and 0.95 will be removed. We expect variants to have a
#' VAF ~0. 0.5, or 1, reflecting rare somatic mutations, heterozygous germline
#' mutations, and homozygous germline mutations, respectively. Default is
#' FALSE.
#' @param custom_filter_val A set of values used to filter rows in
#' \code{mutation_data} based on \code{custom_filter_col}. If a row in
#' \code{custom_filter_col} matches any value in \code{custom_filter_val},
#' it will either be set to TRUE in the \code{filter_mut} column or removed,
#' depending on \code{custom_filter_rm}.
#' @param custom_filter_rm A logical value. If TRUE, rows in custom_filter_col
#' that match any value in custom_filter_val will be removed from the
#' mutation_data. If FALSE, filter_mut will be set to TRUE for those rows.
#' @param regions Remove rows that are within or outside of specified regions.
#' Provide either a data frame or a file path of the specified intervals. Your
#' regions must contain "contig", "start", and "end".  Use the
#' \code{regions_filter} parameter to specify whether rows within the regions
#' should be kept versus removed. To use one of the TSpanels, set the regions
#' parameter as follows: \code{regions = load_regions_file("TSpanel_mouse")}.
#' Change the species as needed for human/rat.
#' @param regions_filter Specifies how the provided \code{regions} should be
#' applied to \code{mutation_data}. Acceptable values are "remove_within" or
#' "keep_within". If set to "remove_within", any rows that fall within the
#' specified regions wil be removed from mutation_data. If set to
#' "keep_within", only the rows within the specified regions will be kept in
#' mutation_data, and all other rows will be removed.
#' @param allow_half_overlap A logical value. If TRUE, rows that start or end
#' in your \code{regions}, but extend outside of them in either direction will
#' be included in the filter. If FALSE, only rows that start and end within the
#' \code{regions} will be included in the filter. Default is FALSE.
#' @param rg_sep The delimiter for importing the \code{regions} file, if
#' applicable. Default is tab-delimited.
#' @param is_0_based_rg A logical value indicating whether the genomic intervals
#' you provided in \code{regions} are 0_based. Set this to TRUE if you are
#' using one of the TSpanels.
#' @param rm_filtered_mut_from_depth A logical value. If TRUE, the function will
#' subtract the \code{alt_depth} of rows that were flagged by the
#' \code{filter_mut} column from their \code{total_depth}. This will treat
#' flagged variants as N-calls. This will not apply to variants flagged as
#' germline by the \code{vaf_cutoff}. However, if the germline variant
#' has additional filters applied, then the subtraction will still occur.
#' If FALSE, the \code{alt_depth} will be retained in the
#' \code{total_depth} for all variants.  Default is FALSE.
#' @param return_filtered_rows A logical value. If TRUE, the function will
#' return both the filtered mutation data and the rows that were
#' removed/flagged in a seperate data frame. The two dataframes will be
#' returned inside a list, with names \code{mutation_data} and
#' \code{filtered_rows}. Default is FALSE.
#' @details
#' Depth correction is important for preventing double-counting of reads in
#' mutation data when summing the total_depth across samples or other groups.
#' Generally, when several mutations have been detected at the same genomic
#' position, within a sample, the total_depth value will be the same for all of
#' them. However, in some datasets, whenever a deletion is detected, the data
#' may contain an additional row with the same genomic position calling a
#' "no_variant". The total_depth will differ between the deletion and the
#' no_variant. In these cases, correct_depth_to_indel == TRUE will ensure that
#' the total_depth value for the deletion is retained, while the total_depth
#' value for the no_variant is removed.
#' @importFrom dplyr group_by mutate ungroup select filter starts_with
#' n_distinct first case_when if_else
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom plyranges join_overlap_left_directed
#' join_overlap_left_within_directed
#' @export
filter_mut <- function(mutation_data,
                       correct_depth = FALSE,
                       correct_depth_to_indel = TRUE,
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
  # import mut will add in_regions column that can be used to filter. Make sure this doesn't interfere with regions filter
  if (return_filtered_rows) {
    rm_rows <- data.frame()
  }

  # Create a new Filter Column
  if (!("filter_mut" %in% colnames(mutation_data))) {
    mutation_data$filter_mut <- FALSE
  }
  if (!("filter_reason" %in% colnames(mutation_data))) {
    mutation_data$filter_reason <- ""
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
      dplyr::mutate(filter_mut = ifelse(.data$vaf > vaf_cutoff,
                                        TRUE, .data$filter_mut),
                    is_germline = ifelse(.data$vaf > vaf_cutoff, TRUE,
                                         FALSE),
                    filter_reason = ifelse(.data$vaf > vaf_cutoff,
                                            ifelse(.data$filter_reason == "", "germline",
                                                   paste0(.data$filter_reason, "|germline")),
                                                   .data$filter_reason))
    vaf_filtered_count <- sum(as.numeric(mutation_data$is_germline == TRUE))
    message("Found ", vaf_filtered_count, " germline mutations.")

  ######## snv_in_germ_mnv Filter #############################################
    if (snv_in_germ_mnv) {
      message("Flagging SNVs overlapping with germline MNVs...")
      mnv_ranges <- mutation_data %>%
        dplyr::filter(.data$variation_type == "mnv" & .data$is_germline == TRUE)

      hits_indices <- integer(0)

      for (sample_name in unique(mutation_data$sample)) {
        snv_subset <- mutation_data %>%
          dplyr::filter(.data$sample == sample_name, .data$variation_type == "snv")
        germ_mnv_subset <- mnv_ranges %>%
          dplyr::filter(.data$sample == sample_name)
        snv_gr <- GenomicRanges::makeGRangesFromDataFrame(snv_subset,
                                                          seqnames.field = "contig",
                                                          keep.extra.columns = TRUE)
        germ_mnv_gr <- GenomicRanges::makeGRangesFromDataFrame(germ_mnv_subset,
                                                               seqnames.field = "contig",
                                                               keep.extra.columns = TRUE)
        overlaps <- GenomicRanges::findOverlaps(query = snv_gr, subject = germ_mnv_gr)
        hits <- S4Vectors::queryHits(overlaps)
        hits_indices <- c(hits_indices, which(mutation_data$sample == sample_name & mutation_data$variation_type == "snv")[hits])
      }
      mutation_data$snv_mnv_overlaps <- FALSE
      mutation_data$snv_mnv_overlaps[hits_indices] <- TRUE
      mutation_data <- mutation_data %>%
        dplyr::mutate(filter_mut = ifelse(.data$snv_mnv_overlaps == TRUE,
                                          TRUE, .data$filter_mut),
                      filter_reason = ifelse(.data$snv_mnv_overlaps == TRUE,
                                             ifelse(.data$filter_reason == "", "snv_in_germ_mnv",
                                                    paste0(.data$filter_reason, "|snv_in_germ_mnv")),
                                            .data$filter_reason))
      snv_in_germ_mnv_count <- sum(mutation_data$snv_mnv_overlaps == TRUE)
      message("Found ", snv_in_germ_mnv_count, " SNVs overlapping with germline MNVs.")
    }
  } else {
    mutation_data$is_germline <- FALSE
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
        dplyr::filter((vaf > 0.05 & vaf < 0.45) | (vaf > 0.55 & vaf < 0.95)) %>%
        dplyr::mutate(filter_reason = ifelse(filter_reason == "",
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
    regions_df <- MutSeqR::load_regions_file("custom_interval", regions, rg_sep)
    regions_df$in_regions <- TRUE
    colnames(regions_df) <- paste0("TO_REMOVE_", colnames(regions_df))

    region_ranges <- GenomicRanges::makeGRangesFromDataFrame(
      df = regions_df,
      keep.extra.columns = TRUE,
      seqnames.field = "TO_REMOVE_contig",
      start.field = "TO_REMOVE_start",
      end.field = "TO_REMOVE_end",
      starts.in.df.are.0based = is_0_based_rg
    )
    mut_ranges <- GenomicRanges::makeGRangesFromDataFrame(
      df = as.data.frame(mutation_data),
      keep.extra.columns = TRUE,
      seqnames.field = "contig",
      start.field = "start",
      end.field = "end",
      starts.in.df.are.0based = FALSE ### FIX: double check that mutation data will always be 1-based.
    )
    if (allow_half_overlap) {
      ranges_joined <- plyranges::join_overlap_left_directed(mut_ranges,
                                                             region_ranges,
                                                             suffix = c("", "_regions"))
    } else {
      ranges_joined <- plyranges::join_overlap_left_within_directed(mut_ranges, region_ranges, suffix = c("", "_regions"))
    }
    mutation_data <- as.data.frame(ranges_joined)
    mutation_data$TO_REMOVE_in_regions[is.na(mutation_data$TO_REMOVE_in_regions)] <- FALSE
    mutation_data <- dplyr::rename(mutation_data, contig = seqnames)
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

  ######## Depth Correction ###################################################
  if (correct_depth) {
    if (!("total_depth" %in% colnames(mutation_data))) {
      stop("Error: You have set correct_depth to TRUE but there is no
      'total_depth' column in your mutation_data.")
    }
    message("Correcting depth...")
    if (correct_depth_to_indel) {

      # priority list for variation types
      variation_priority <- c("deletion", "complex", "insertion", "snv", "mnv",
                              "sv", "uncategorized", "no_variant")

      # Step 1: Identify the highest priority type for each unique group
      priority_data <- mutation_data %>%
        dplyr::group_by(.data$sample, .data$contig, .data$start) %>%
        dplyr::reframe(
          all_depths_same = dplyr::n_distinct(.data$total_depth) == 1,
          highest_priority_type = variation_type[which.min(match(variation_type, variation_priority))])
      # Step 2: Join back to the original data for efficient updates
      mutation_data <- mutation_data %>%
        dplyr::left_join(priority_data, by = c("sample", "contig", "start")) %>%
        dplyr::group_by(sample, contig, start) %>%
        dplyr::mutate(first_row = dplyr::row_number() == 1) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(total_depth = dplyr::case_when(
                        all_depths_same ~ dplyr::if_else(first_row, total_depth, 0),
                        variation_type == highest_priority_type ~ total_depth,
                        TRUE ~ 0)) %>%
        dplyr::select(-"all_depths_same",
                      -"highest_priority_type",
                      -"first_row")
    } else {
      mutation_data <- mutation_data %>%
        dplyr::group_by(.data$sample, .data$contig, .data$start) %>%
        dplyr::mutate(
          total_depth = dplyr::if_else(row_number() == 1, .data$total_depth, 0)) %>%
        dplyr::ungroup()
    }
    corrected_depth_count <- sum(mutation_data$total_depth == 0)
    message(corrected_depth_count, " rows had their total_depth corrected.")
  }
  if (rm_filtered_mut_from_depth) {
    ####TO DO Rethink this: what if a germline mutation is problematic? Can that happen? If VAF = 1, is_germline is not created = error
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
    ## If filter_rows_return is empty, return just mutation data, no list.
    filtered_muts <- mutation_data %>%
      dplyr::filter(filter_mut == TRUE)
    filter_rows_return <- rbind(rm_rows, filtered_muts)
    message("Returning a list: mutation_data and filtered_rows.")
    return(list(mutation_data = mutation_data, filtered_rows = filter_rows_return))
  } else {
    return(mutation_data)
  }
}
