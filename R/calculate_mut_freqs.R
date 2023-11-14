#' Calculate mutation frequency
#'
#' Calculates the mutation frequency for arbitrary groupings and creates a new
#' dataframe with the results. Mutation frequency is # mutations / total bases,
#' but this can be subset in different ways: e.g., by mutation context. In this
#' case, it is necessary to change the denominator of total bases to reflect
#' the sequencing depth at the proper reference bases under consideration.
#' 
#' Additionally, the operation is run by default using clonally expanded as well
#' as unique mutations.
#' @param data The data frame to be processed containing mutation data.
#' @param cols_to_group A vector of grouping variables: this should be the groups
#' of interest that you want to calculate a frequency for. For instance, getting
#' the frequency by `sample`. Other options might include `locus`, or,
#' `c("sample","locus")`. Must be a column in the mutation data table.
#' @param subtype_resolution At what resolution should the frequencies be
#' calculated? This is an important consideration because the reference bases
#' sequencing depth will differ depending on which mutation subtypes are being
#' considered. Options are "none", 6base", "12base", "96base", and "192base".
#' Could add 24 base option eventually. By selecting "none", this will use the
#' total depth across all the selected metadata columns.
#' @param vaf_cutoff Exclude rows from frequency calculations
#' with a variant allele fraction (VAF) greater than this value. Default is 0.1.
#' @param summary TRUE or FALSE, whether to return a summary table (i.e., where
#' only relevant columns for frequencies and groupings are returned). Setting 
#' this to false returns all columns in the original data, which might make
#' plotting more difficult, but may provide additional flexibility to power
#' users.
#' @param  clonality_cutoff NOT CURRENTLY IMPLEMENTED! Up for consideration.
#' This value determines the fraction of reads that
#' is considered a constitutional variant. If a mutation is present at a 
#' fraction higher than this value, the reference base will be swapped,
#' and the alt_depth recalculated. 0.3 (30%) would be a sane default?
#' @param variant_types Include these variant types. A vector of one or more
#'  "snv", "indel", "sv", "mnv", "no_variant"
#' @param custom_column_names A list of names to specify the meaning of column
#'  headers. Since column names can vary with data, this might be necessary to
#'  digest the mutation data table properly. Typical defaults are set, but can
#'  be substituted in the form of `list(total_depth = "my_custom_depth_name", 
#'  sample = "my_custom_sample_column_name")`. For a comprehensive list, see 
#'  examples. You can change one or more of these. The most likely column to
#'  need changing is the "chromosome name" (chr), which by default is "seqnames"
#'  but could be "contig", "chr", or others.  (TODO - MAKE EXAMPLES)
#' @returns A data frame with the mutation frequency calculated.
#' @importFrom dplyr across all_of filter group_by mutate n row_number select distinct ungroup  
#' @importFrom magrittr %>%
#' @importFrom rlang := .data
#' @importFrom utils modifyList
#' @export
calculate_mut_freq <- function(data,
                               cols_to_group = "sample",
                               subtype_resolution = "6base",
                               vaf_cutoff = 0.1,
                               clonality_cutoff = 0.3,
                               summary = TRUE,
                               variant_types = c("snv","indel","mnv","sv"),
                               custom_column_names = list(chr = "seqnames")) {

  # Determine columns to use for different parts
  default_columns <- DupSeqR::op$column
  cols <- modifyList(default_columns, custom_column_names)
  
  freq_col_prefix <- paste0(cols_to_group, collapse = "_")
  
  if (!subtype_resolution %in% names(DupSeqR::subtype_dict)) {
    stop(paste0("Error: you need to set subtype_resolution to one of: ",
                paste(names(DupSeqR::subtype_dict), collapse = ", ")))
  }

  numerator_groups <- c(cols_to_group, DupSeqR::subtype_dict[[subtype_resolution]])
  numerator_groups <- numerator_groups[!is.na(numerator_groups)]
  denominator_groups <- c(cols_to_group, DupSeqR::denominator_dict[[subtype_resolution]])
  denominator_groups <- denominator_groups[!is.na(denominator_groups)]
  # Check if data is provided as GRanges: if so, convert to data frame.
  if (inherits(data, "GRanges")) { data <- as.data.frame(data) }
  if (!inherits(data, "data.frame")) { warning("You should probably use a 
                                             data frame as input here.")}

  
  #############################################################
  ###  Dealing with Total  Depth and  Duplicated  Rows ####
  #  ISSUES: Group_by cols_to_group -  if this  is  set  to "dose", then  this 
  # code will not  work
  # 1. Take the deletion total_depth
  # 2. Take the mean of the total_depths (rounded)
  #  read_vaf will have already made  sure that  all   duplicates have the  same total_depth
  # import_mut_file allows for different total_depths amongst the duplicates. 
    # SO  we will need  to implement that code into the import_mut_file function. 
    # Calculate mutation frequencies
  mut_freq_table <- data %>%
    # Identify duplicate entries prior to depth calculation
    dplyr::group_by(dplyr::across(dplyr::all_of(c(cols$sample, cols$chr, cols$start)))) %>%
    mutate(num_dups = dplyr::n(), 
           dup_id = dplyr::row_number()) %>% 
    dplyr::ungroup() %>%
    mutate(is_duplicated = .data$num_dups > 1) %>%
    mutate(depth_undupes = ifelse(.data$dup_id>1, 0, .data[[cols$total_depth]])) %>%
    #Identify values with the same start, sample, and depth
    dplyr::group_by(dplyr::across(dplyr::all_of(c(cols_to_group, cols$chr, cols$start, cols$total_depth)))) %>%
    mutate(is_depth_duplicated = dplyr::n() > 1) %>%
    dplyr::ungroup() %>%
    # For sites that share start positions within a group, de-duplicate depth.
    # For duplicated depths, set all instance but one to zero; for
    # indels/svs/mnvs, take the 'no_variant' depth instead.
    mutate(depth_final = ifelse(.data$is_duplicated == TRUE & .data$is_depth_duplicated == FALSE,
                                ifelse(.data$variation_type == "no_variant", .data[[cols$total_depth]], 0),
                                ifelse(.data$is_depth_duplicated == TRUE, .data$depth_undupes, .data$total_depth)))
   
    
    #####################################################################
    
    
     dplyr::group_by(dplyr::across(dplyr::all_of(c(numerator_groups)))) %>%
    mutate(!!paste0(freq_col_prefix, "_sum_clonal") :=
        sum(.data$alt_depth[!.data$variation_type == "no_variant" & .data$VAF < vaf_cutoff])) %>%
    mutate(!!paste0(freq_col_prefix, "_sum_unique") :=
             length(.data$alt_depth[!.data$variation_type == "no_variant" & .data$VAF < vaf_cutoff])) %>%
    # Calculate denominator (same for clonal and unique mutations)
    dplyr::group_by(dplyr::across(dplyr::all_of(c(denominator_groups)))) %>%
    mutate(!!paste0(freq_col_prefix, "_depth") := sum(.data$depth_final)) %>%
    #mutate(!!paste0(freq_col_prefix, "_depth") := sum(.data$total_depth[!is_duplicate])) %>%
    # Calculate frequencies
    mutate(!!paste0(freq_col_prefix, "_MF_clonal") :=
             .data[[paste0(freq_col_prefix, "_sum_clonal")]] /
             .data[[paste0(freq_col_prefix, "_depth")]]) %>%
    mutate(!!paste0(freq_col_prefix, "_MF_unique") :=
             .data[[paste0(freq_col_prefix, "_sum_unique")]] /
             .data[[paste0(freq_col_prefix, "_depth")]]) %>%
    dplyr::ungroup()

  # Define columns of interest for summary table
  summary_cols <- c(
    numerator_groups,
    paste0(freq_col_prefix, "_sum_unique"),
    paste0(freq_col_prefix, "_sum_clonal"),
    paste0(freq_col_prefix, "_depth"),
    paste0(freq_col_prefix, "_MF_clonal"),
    paste0(freq_col_prefix, "_MF_unique")
  )

  # Make summary table of frequencies
  # This is also where subtype proportions are calculated
  summary_table <- mut_freq_table %>%
    dplyr::filter(.data$variation_type %in% variant_types) %>%
    dplyr::select({{ summary_cols }}) %>%
    dplyr::distinct() %>%
    mutate(freq_clonal =
             .data[[paste0(freq_col_prefix, "_sum_clonal")]] /
             sum(.data[[paste0(freq_col_prefix, "_sum_clonal")]]) /
             .data[[paste0(freq_col_prefix, "_depth")]] ) %>%
    mutate(prop_clonal = .data$freq_clonal / sum(.data$freq_clonal)) %>%
    mutate(freq_unique =
             .data[[paste0(freq_col_prefix, "_sum_unique")]] /
             sum(.data[[paste0(freq_col_prefix, "_sum_unique")]]) /
             .data[[paste0(freq_col_prefix, "_depth")]] ) %>%
    mutate(prop_unique = .data$freq_unique / sum(.data$freq_unique))

  if (!summary) {
    return(mut_freq_table)
  } else if (summary) {
    return(summary_table)
  } else {
    stop("Error: summary must be TRUE or FALSE.")
  }
}

# For testing
#cols_to_group <- c("sample")
#cols_to_group <- c("sample","description")
#cols_to_group <- c("sample","variation_type")
#cols_to_group <- c("sample", "CpG_site")

# To do... locate and enumerate recurrent mutations?
