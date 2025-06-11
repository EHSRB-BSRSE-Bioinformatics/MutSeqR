#' Calculate mutation frequency
#'
#' Calculates the mutation frequency for arbitrary groupings and creates a new
#' dataframe with the results. Mutation frequency is # mutations / total bases,
#' but this can be subset in different ways: e.g., by mutation context. In this
#' case, it is necessary to change the denominator of total bases to reflect
#' the sequencing depth at the proper reference bases under consideration.
#' Additionally, by default, the operation is run by default using both the minimum 
#' and maximum independent methods for counting mutations.
#' @param mutation_data The data frame to be processed containing mutation data.
#' Required columns are listed below.
#' Synonymous names for these columns are accepted.
#' \itemize{
#'      \item `contig`: The reference sequence name.
#'      \item `start`: 1-based start position of the feature.
#'      \item `alt_depth`: The read depth supporting the alternate allele.
#'      \item `variation_type`: The category to which this variant is assigned.
#'      \item subtype_col: The column containing the mutation subtype. This
#' column depends on the `subtype_resolution` parameter.
#'     \item reference context: The column containing the referene base(s) for the
#' mutation. This column depends on the `subtype_resolution` parameter.
#'    \item cols to group: all columns across which you want to calculate
#' the mutation frequency. Ex. `c("tissue", "dose")`. These columns should be listed
#' in cols_to_group.
#' }
#' It is also required to include the total_depth column if you are calculating
#' depth from the mutation data. If you are using precalculated depth data, the
#' total_depth column is not required.
#' 
#' @param cols_to_group A vector of grouping variables: this should be the
#' groups of interest that you want to calculate a frequency for.
#' For instance, getting the frequency by `sample`. Other options might
#' include `dose`, `locus`, or, `c("sample","locus")`. All listed variables
#' must be a column in the mutation_data.
#' @param subtype_resolution The resolution of the mutation subtypes. Options are
#'  \itemize{
#'        \item "none" calculates mutation frequencies across all selected
#' grouping columns.
#'         \item "type" calculates mutation frequencies across all selected
#' grouping columns for each `variation_type` seperately; snv, mnv, deletion,
#' insertion, complex, sv, ambiguous, uncategorized.
#'          \item "base_6" calculates mutation frequencies across all selected
#' grouping columns for each variation_type with snv mutations separated by
#' `normalized_subtype`; C>A, C>G, C>T, T>A, T>C, T>G. The reference context is
#' `normalized_ref`.
#'          \item "base_12" calculates mutation frequencies across all
#' selected grouping columns for each variation_type with snv mutations
#' separated by `subtype`; A>C, A>G, A>T, C>A, C>G, C>T, G>A, G>C, G>T,
#' T>A, T>C, T>G. The reference context is `short_ref`.
#'           \item "base_96" calculates mutation frequencies across all
#' selected grouping columns for each variation_type with snv mutations
#' separated by `normalized_context_with_mutation`, i.e. the 96-base
#' trinucleotide context. Ex. A\\[C>T\\]A. The reference context is
#' `normalized_context`.
#'           \item "base_192" calculates mutation frequencies across all
#' selected grouping columns for each variation_type with snv mutations
#' separated by `context_with_mutation`, i.e. the 192-base trinucleotide
#' context. Ex A\\[G>A\\]A. The reference context is `context`.
#'  }
#' @param variant_types Use this parameter to choose which variation types
#' to include in the mutation counts. Provide a character vector of the
#' variation types that you want to include. Alternatively, provide a
#' character vector of the variation types that you want to exclude preceded
#' by "-". Options are: "snv", "complex", "deletion", "insertion", "mnv", "sv",
#' "ambiguous", "uncategorized". Ex. inclusion: "snv", exclusion: "-snv".
#' Default includes all variants. For `calculate_depth = TRUE`: Regardless of
#' whether or not a variant is included in the mutation counts, the total_depth
#' for that position will be counted.
#' @param correct_depth A logical value, only used if `calculate_depth = TRUE`.
#' If TRUE (default), an internal correction is applied to prevent double-counting
#' sequencing depth at genomic sites with multiple mutations. Set to FALSE only if
#' your data has already been processed with `correct_depth()`.
#' @param correct_depth_by_indel_priority A logical value, passed to the internal
#' depth correction. If TRUE, depth is retained for the mutation with the highest
#' indel/variant priority. If FALSE (default), depth is retained for the first
#' mutation encountered at a site.
#' @param calculate_depth A logical variable, whether to calculate the
#' per-group total_depth from the mutation data. If set to TRUE, the mutation
#' data must contain a total_depth value for every sequenced base (including
#' variants AND no-variant calls). If set to FALSE, pre-calculated per-group
#' total_depth values may be supplied at the desired subtype_resolution
#' using the precalc_depth_data parameter. Alternatively, if no per-group
#' total_depth is available, per-group mutation counts will be calculated,
#' but mutation frequency will not. In such cases, mutation subtype proportions
#' will not be normalized to the total_depth.
#' @param precalc_depth_data A data frame or a file path to a text file
#' containing pre-calculated per-group total_depth values. This data frame
#' should contain the columns for the desired grouping variable(s)
#' and the reference context at the desired subtype resolution (if applicable).
#' The precalculated total_depth column(s) should be called one of
#' `group_depth` and `subtype_depth`. `group_depth` is used for subtype
#' resolutions of "none", "type", and all non-snv mutations in "base_6",
#' "base_12", "base_96", and "base_192". `subtype_depth` is used for snv
#' mutations in "base_6", "base_12", "base_96", and "base_192". You can
#' access a list of context values for each subtype resolution using
#' `MutSeqR::context_list$your_subtype_resolution`.
#' @param d_sep The delimiter used in the precalc_depth_data, if applicable.
#' Default is tab-delimited.
#' @param summary A logical variable, whether to return a summary table
#' (i.e., where only relevant columns for frequencies and groupings are
#' returned). Setting this to false returns all columns in the original
#' mutation_data, which might make plotting more difficult, but may provide
#' additional flexibility to power users.
#' @param retain_metadata_cols a character vector that contains the names of the
#' metadata columns that you would like to retain in the summary table.
#' This may be useful for plotting your summary data. Ex. retain the "dose"
#' column when summarising by "sample".
#' @returns A data frame with the mutation frequency calculated. If summary
#' is set to TRUE, the data frame will be a summary table with the mutation
#' frequency calculated for each group. If summary is set to FALSE, the
#' mutation frequency will be appended to each row of the original
#' mutation_data.
#'  \itemize{
#'      \item `sum_min`: The sum of all mutations within the group, calculated
#' using the "min" method for mutation counting. All identical mutations
#' within a samples are assumed to be the result of clonal expansion and are
#' thus only counted once.
#'      \item `sum_max`: The sum of all mutations within the group, calculated
#' using the "max" method for mutaiton counting. All identical mutations
#' within a sample are assumed to be idenpendant mutational evens and are
#' included in the mutation frequency calculation.
#'       \item `mf_min`: The mutation frequency calculated using the "min"
#' method for mutation counting. mf_min = sum_min / depth.
#'        \item `mf_max`: The mutation frequency calculated using the "max"
#' method for mutation counting. mf_max = sum_max / depth.
#'     \item `proportion_min`: The proportion of each mutation
#' subtype within the group, normalized to the depth. Calculated
#' using the "min" method. This is only calculated if `subtype_resolution`
#' is not "none". If no depth is calculated or provided, proportion is
#' calculated without normalization to the depth.
#'     \item `proportion_max`: The proportion of each mutation
#' subtype within the group, normalized to its read depth. Calculated
#' using the "max" method. This is only calculated if `subtype_resolution`
#' is not "none". If no depth is calculated or provided, proportion is
#' calculated without normalization to the depth.
#' }
#' @examples
#' # Load example data
#' example_file <- system.file("extdata", "Example_files",
#'                             "example_mutation_data_filtered.rds",
#'                             package = "MutSeqR")
#' example_data <- readRDS(example_file)
#'
#' # Example 1 Calculate mutation frequency by sample.
#' # Calculate depth from the mutation data
#' mf_example <- calculate_mf(mutation_data = example_data,
#'                            cols_to_group = "sample")
#' # Example 2: Calculate the trinucleotide mutation proportions for each dose
#' mf_96_example <- calculate_mf(mutation_data = example_data,
#'                               cols_to_group = "dose",
#'                               subtype_resolution = "base_96",
#'                               variant_types  = "snv")
#' # Example 3: Calculate the mean mutation frequency for each 6 base subtype
#' # per dose
#' # calculate_mf does not calculate mean mutation frequency for
#' # groups; this function only sums mutations across groups. Thus, if you are
#' # interested in calculating the mean of a group, this must be done
#' # separately.
#' # First, calculate 6 base MF per sample. Retain the dose column.
#' mf_6_example <- calculate_mf(mutation_data = example_data,
#'                              cols_to_group = "sample",
#'                              subtype_resolution = "base_6",
#'                              retain_metadata_cols = "dose")
#' # Note: NA values in retain_metadata_cols. When we create a summary table
#' # that includes mutation subtypes, there may occasionally be NA values in
#' # the metadata columns. This is because the mutation data does not contain
#' # any mutations (filtered or not) for that particular subtype within the
#' # given group. For example, our example_data does not contain any ambiguous
#' # or uncategorized mutations, so the dose column is NA for all those
#' # mutations in the summary table. This will not affect downstream analyses.
#'
#' # Calculate the mean MF for each 6 base subtype per dose
#' mf_6_mean_example <- mf_6_example %>%
#'  dplyr::group_by(dose, normalized_subtype) %>%
#'  dplyr::summarise(mean_mf_min = mean(mf_min),
#'                   se_mf_min = sd(mf_min) / sqrt(dplyr::n()),
#'                   mean_mf_max = mean(mf_max),
#'                   se_mf_max = sd(mf_max) / sqrt(dplyr::n()))
#'
#' # Example 4: Calculate MF using precalculated depth data
#' sample_depth_example <- data.frame(sample = c("dna00973.1", "dna00974.1",
#'                                               "dna00975.1", "dna00976.1",
#'                                               "dna00977.1", "dna00978.1",
#'                                               "dna00979.1", "dna00980.1",
#'                                               "dna00981.1", "dna00982.1",
#'                                               "dna00983.1", "dna00984.1",
#'                                               "dna00985.1", "dna00986.1",
#'                                               "dna00987.1", "dna00988.1",
#'                                               "dna00989.1", "dna00990.1",
#'                                               "dna00991.1", "dna00992.1",
#'                                               "dna00993.1", "dna00994.1",
#'                                               "dna00995.1", "dna00996.1"),
#'                                    group_depth = c(565395266, 755574283,
#'                                                    639909215, 675090988,
#'                                                    598104021, 611295330,
#'                                                    648531765, 713240735,
#'                                                    669734626, 684951248,
#'                                                    716913381, 692323218,
#'                                                    297661400, 172863681,
#'                                                    672259724, 740901132,
#'                                                    558051386, 733727643,
#'                                                    703349287, 884821671,
#'                                                    743311822, 799605045,
#'                                                    677693752, 701163532))
#' mf_example_precalc <- calculate_mf(mutation_data = example_data,
#'                                    cols_to_group = "sample",
#'                                    calculate_depth = FALSE,
#'                                    precalc_depth_data = sample_depth_example)
#'
#' # Example 5: Calculate MF using precalculated depth data for 6 base
#' # mutation subtypes per sample.
#' # The base_6 resolution uses reference context 'normalized_ref'; C or T.
#' # Our precalc_depth_data needs group_depth (depth per sample) and the
#' # subtype_depth (depth per sample AND per normalized_ref)
#' # We will create the example precalc_depth data for the base_6 resolution
#' # from Example 3 results for simplicity.
#' sample_subtype_depth_example <- mf_6_example %>%
#'  dplyr::select(sample, normalized_ref, group_depth, subtype_depth) %>%
#'  unique() %>%
#' dplyr::filter(normalized_ref != "N")
#' mf_6_example_precalc <- calculate_mf(mutation_data = example_data,
#'                                      cols_to_group = "sample",
#'                                      subtype_resolution = "base_6",
#'                                      calculate_depth = FALSE,
#'                                      precalc_depth_data = sample_subtype_depth_example)
#' @importFrom dplyr across all_of filter group_by mutate n row_number
#' select distinct ungroup
#' @importFrom magrittr %>%
#' @importFrom rlang := .data
#' @importFrom utils modifyList
#' @importFrom stats na.omit
#' @export

calculate_mf <- function(mutation_data,
                         cols_to_group = "sample",
                         subtype_resolution = "none",
                         variant_types = c("snv",
                                           "deletion",
                                           "insertion",
                                           "complex",
                                           "mnv",
                                           "sv",
                                           "ambiguous",
                                           "uncategorized"),
                         calculate_depth = TRUE,
                         correct_depth = TRUE,
                         correct_depth_by_indel_priority = FALSE,
                         precalc_depth_data = NULL,
                         d_sep = "\t",
                         summary = TRUE,
                         retain_metadata_cols = NULL) {

  # Variant list
  all_variant_types <- setdiff(MutSeqR::subtype_list$type, "no_variant")
  filter_variants <- function(selected_types, all_variant_types) {
    # Convert to character vector if it's a single string input
    if (is.character(selected_types) && length(selected_types) == 1) {
      selected_types <- unlist(strsplit(selected_types, ","))
    }
    # Check if any input starts with "-"
    exclusions <- grepl("^-", selected_types)
    if (any(exclusions)) {
      # Remove "-" and exclude those variants
      excluded_variants <- sub("^-", "", selected_types[exclusions])
      selected_variants <- setdiff(all_variant_types, excluded_variants)
    } else {
      # Only include the specified variants
      selected_variants <- intersect(all_variant_types, selected_types)
    }
    return(selected_variants)
  }

  variant_types <- filter_variants(variant_types, all_variant_types)

  # Validate Parameters
  # Check if data is provided as GRanges: if so, convert to data frame.
  if (inherits(mutation_data, "GRanges")) {
    mutation_data <- as.data.frame(mutation_data)
  }
  if (!inherits(mutation_data, "data.frame")) {
    warning("You should use a data frame as input here.")
  }
  if (!subtype_resolution %in% names(MutSeqR::subtype_dict)) {
    stop(paste0(
      "Error: you need to set subtype_resolution to one of:
      none, type, base_6, base_12, base_96, base_192")
    )
  }
  if (any(!variant_types %in% MutSeqR::subtype_list$type)) {
    stop(paste0(
      "Error: you need to set variant_types to one or more of: ",
      paste(MutSeqR::subtype_list$type, collapse = ", "),
      ". Variation_types outside of this list will not be included in the mutation frequency calculation."
    ))
  }
  if (!is.logical(summary)) {
    stop("summary must be a logical variable.")
  }
  if (!is.null(retain_metadata_cols) && !is.character(retain_metadata_cols)) {
    stop("retain_metadata_cols must be a character vector.")
  }

  if (calculate_depth && correct_depth) {
    if (!"total_depth" %in% colnames(mutation_data)) {
      stop("Error: `correct_depth` is TRUE but 'total_depth' column not found in mutation_data.")
    }

    message("Performing internal depth correction to prevent double-counting...")
    dt <- data.table::as.data.table(mutation_data)

    if (correct_depth_by_indel_priority) {
      variation_priority <- c("deletion", "complex", "insertion", "snv", "mnv", "sv", "uncategorized", "no_variant")
      dt[, priority_order := factor(variation_type, levels = variation_priority, ordered = TRUE)]
      dt[, total_depth := {
        group_order <- order(priority_order, na.last = TRUE)
        corrected_depths <- rep(0, .N)
        corrected_depths[group_order[1]] <- total_depth[group_order[1]]
        corrected_depths
      }, by = .(sample, contig, start)]
      dt[, priority_order := NULL]
    } else {
      dt[, total_depth := c(total_depth[1], rep(0, .N - 1)), by = .(sample, contig, start)]
    }
    
    # Overwrite the input data frame with the corrected version
    mutation_data <- as.data.frame(dt)
    message("Internal depth correction complete.")
  }

  # Rename columns in mutation_data to default
  mutation_data <- MutSeqR::rename_columns(mutation_data)
  # Check for all required columns
  required_columns <- c(
   # "sample",
    "alt_depth",
    "variation_type",
    "filter_mut",
    cols_to_group
  )
  if (calculate_depth) {
    required_columns <- c(required_columns, "total_depth")
  }
  if(subtype_resolution != "none") {
    required_columns <- c(required_columns,
                          MutSeqR::subtype_dict[[subtype_resolution]])
    if (subtype_resolution != "type") {
      required_columns <- c(required_columns,
                            MutSeqR::denominator_dict[[subtype_resolution]])
    }
  }

  if (!is.null(retain_metadata_cols)) {
    required_columns <- c(required_columns, retain_metadata_cols)
  }
  mutation_data <- MutSeqR::check_required_columns(mutation_data,
                                                   required_columns)

  if (subtype_resolution %in% c("base_6", "base_12", "base_96", "base_192")
      && !("snv" %in% variant_types)) {
    warning("Please include 'snv' in parameter 'variant_types' to calculate
            single-nucleotide variant subtype frequencies.")
  }

  numerator_groups <- c(cols_to_group,
                        MutSeqR::subtype_dict[[subtype_resolution]])
  numerator_groups <- numerator_groups[!is.na(numerator_groups)]
  denominator_groups <- c(cols_to_group,
                          MutSeqR::denominator_dict[[subtype_resolution]])
  denominator_groups <- denominator_groups[!is.na(denominator_groups)]

  # Calculate mutation counts groups
  mut_freq_table <- mutation_data %>%
    dplyr::mutate(alt_depth_min = ifelse(.data$alt_depth == 0, 0, 1)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(numerator_groups)))) %>%
    dplyr::mutate(sum_max =
                    sum(.data$alt_depth[.data$variation_type %in% variant_types
                                        & .data$filter_mut == FALSE])) %>%
    dplyr::mutate(sum_min =
                    sum(.data$alt_depth_min[.data$variation_type %in% variant_types
                                            & .data$filter_mut == FALSE])) %>%
    dplyr::ungroup() %>%
    dplyr::select(-"alt_depth_min")

  # Calculate denominator (same for max and min mutations)
  # snv depth: depth across groups and snv subtype resolution
  # group depth: depth across groups.
  if (calculate_depth) {
    mut_freq_table <- mut_freq_table %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(cols_to_group))) %>%
      dplyr::mutate(group_depth =
                      sum(.data$total_depth)) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(denominator_groups)))) %>%
      dplyr::mutate(subtype_depth =
                      sum(.data$total_depth)) %>%
      dplyr::ungroup()
    depth_exists <- TRUE
  } else {
    if (!is.null(precalc_depth_data)) {
      if (is.data.frame(precalc_depth_data)) {
        depth_df <- precalc_depth_data
      } else if (is.character(precalc_depth_data)) {
        depth_file <- file.path(precalc_depth_data)
        if (!file.exists(depth_file)) {
          stop("The precalc_depth_data does not exist.")
        }
        if (file.info(depth_file)$size == 0) {
          stop("Error: You are trying to import an empty precalc_depth_data")
        }
        depth_df <- read.delim(file.path(depth_file),
                               sep = d_sep,
                               header = TRUE)
        if (ncol(depth_df) <= 1) {
          stop("Your imported precalc only has one column.
              You may want to set d_sep to properly reflect
              the delimiter used for the data you are importing.")
        }
      } else {
        stop("precalc_depth_data must be NULL, a data frame, or a file path.")
      }
      # check for required columns in depth_df
      # If they are just using snvs, maybe don't require group_depth
      required_columns <- c(denominator_groups, "group_depth")
      if (subtype_resolution %in% c("base_6", "base_12", "base_96", "base_192")) {
        required_columns <- c(required_columns, "subtype_depth")
      } else {
        # add a subtype_depth column if it doesn't exist (for types)
        depth_df$subtype_depth <- depth_df$group_depth
      }
      missing_columns <- setdiff(required_columns, colnames(depth_df))
      if (length(missing_columns) > 0) {
        stop("Missing columns in precalc_depth_data: ",
             paste(missing_columns, collapse = ", "), "\n")
      }
      # Merge depth_df with mut_freq_table
      mut_freq_table <- dplyr::left_join(mut_freq_table, depth_df,
                                         by = c(denominator_groups))
      depth_exists <- TRUE
    } else {
      warning("No depth data provided. Mutation frequencies will not be calculated. Mutation subtype proportions will not be normalized to the total_depth.")
      depth_exists <- FALSE
    }
  }

  # Calculate frequencies
  if (depth_exists) {
    if (subtype_resolution == "none") {
      mut_freq_table <- mut_freq_table %>%
        dplyr::mutate(mf_max = .data$sum_max / .data$group_depth,
                      mf_min = .data$sum_min / .data$group_depth) %>%
        dplyr::ungroup()
    } else {
      mut_freq_table <- mut_freq_table %>%
        dplyr::mutate(mf_max = ifelse(!!rlang::sym(MutSeqR::subtype_dict[[subtype_resolution]])  %in% MutSeqR::subtype_list$type,
                                      .data$sum_max / .data$group_depth,
                                      .data$sum_max / .data$subtype_depth),
                      mf_min = ifelse(!!rlang::sym(MutSeqR::subtype_dict[[subtype_resolution]]) %in% MutSeqR::subtype_list$type,
                                      .data$sum_min / .data$group_depth,
                                      .data$sum_min / .data$subtype_depth)) %>%
        dplyr::ungroup()
    }
    summary_cols <- c(numerator_groups,
                      "sum_min",
                      "sum_max",
                      "mf_min",
                      "mf_max"
                    )
  } else {
    summary_cols <- c(numerator_groups,
                      "sum_min",
                      "sum_max"
                    )
  }

  # Create a summary table
  group_list <- lapply(cols_to_group, function(col) unique(mut_freq_table[[col]]))
  names(group_list) <- cols_to_group
  group_df <- do.call(expand.grid, group_list)
  non_snvs <- setdiff(MutSeqR::subtype_list$type, "snv")

  # Extract the depth values from mutation data, if we calculated them
  if (calculate_depth) {
    depth_df <- mut_freq_table %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(denominator_groups))) %>%
      dplyr::summarise(group_depth = dplyr::first(.data$group_depth),
                       subtype_depth = dplyr::first(.data$subtype_depth),
                       .groups = "drop")

    # Make sure we have all a complete list of groups and context rows for the summary table
    # Do not depend on the data to have all the possible context rows; we will grab from set list in context_list
    # Missing data will be 0s eventually.
    if (subtype_resolution %in% c("base_6", "base_12", "base_96", "base_192")) {
      context_rows_snv <- MutSeqR::context_list[[subtype_resolution]]
      context_rows_snv <- tidyr::expand_grid(group_df, !!paste(MutSeqR::denominator_dict[[subtype_resolution]]) := context_rows_snv)

      depth_df <- dplyr::left_join(context_rows_snv, depth_df, by = denominator_groups)

      # re-extract the group depth seperately, to make sure we have a complete list of all groups
      # Do not depend on extracting all group x context combinations from the data as we may end up missing some groups
      depth_df <- dplyr::select(depth_df, -"group_depth")
      group_depth_df <- mut_freq_table %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(cols_to_group))) %>%
        dplyr::summarise(group_depth = dplyr::first(.data$group_depth),
                         .groups = "drop")
      depth_df <- dplyr::left_join(depth_df, group_depth_df, by = cols_to_group)
      # make a list of depths for non-snv subtypes to add to the depth_df (reference = N)      
      if (any(non_snvs %in% variant_types)) {
        group_depth_df[[MutSeqR::denominator_dict[[subtype_resolution]]]] <- "N"
        group_depth_df$subtype_depth <- group_depth_df$group_depth
        depth_df <- rbind(group_depth_df, depth_df)
      }
    } else { # none or type
      depth_df$subtype_depth <- depth_df$group_depth
    }
  }

  # Make Summary Rows for Mutations
  subset_type <- list(MutSeqR::subtype_list$type[MutSeqR::subtype_list$type %in% variant_types])
  if (subtype_resolution != "none") {
    names(subset_type) <- MutSeqR::subtype_dict[[subtype_resolution]]
    summary_rows <- do.call(expand.grid, c(subset_type, group_list)) # full list of non-snv subtypes
    if (subtype_resolution %in% c("base_6", "base_12", "base_96", "base_192")) {
      snv_subtype <- list(MutSeqR::subtype_list[[subtype_resolution]])
      names(snv_subtype) <- MutSeqR::subtype_dict[[subtype_resolution]]
      summary_rows_snv <- do.call(expand.grid, c(snv_subtype, group_list)) # full list of snv subtypes
      summary_rows <- rbind(summary_rows, summary_rows_snv)
      summary_rows <- summary_rows %>% # add the context column
        dplyr::rowwise() %>%
        dplyr::mutate(!!paste(denominator_dict[[subtype_resolution]]) :=
                      get_ref_of_mut(get(subtype_dict[[subtype_resolution]])))
      summary_rows <- summary_rows %>%
        dplyr::mutate(!!paste0(MutSeqR::denominator_dict[[subtype_resolution]])
                     := if_else(is.na(!!sym(paste0(MutSeqR::denominator_dict[[subtype_resolution]]))), 
                                "N",
                                !!sym(paste0(MutSeqR::denominator_dict[[subtype_resolution]]))))
    }
  } else {
    summary_rows <- group_df
  }

  # grab the metadata columns
  if (!is.null(retain_metadata_cols)) {
    metadata <- mut_freq_table %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(numerator_groups))) %>%
      dplyr::summarise(across(all_of(retain_metadata_cols), ~ first(.)), .groups = "drop")
    summary_rows <- dplyr::left_join(summary_rows, metadata)
  }

  # Grab the data
  summary_data <- mut_freq_table %>%
    dplyr::filter(.data$variation_type %in% unlist(subset_type)) %>%
    dplyr::select(
      {{ summary_cols }}) %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(c(numerator_groups))),
                    .keep_all = TRUE)

  # Merge summary rows and data cols.
  # this makes sure that every group and subtype has a row in the summary table, regardless of if there is a variant in the data.
  summary_table <- merge(summary_rows, summary_data, all = TRUE)

  if (subtype_resolution %in% c("base_6", "base_12", "base_96", "base_192")) {
    # Remove snv rows since they are now represented by their subtypes
    summary_table <- filter(summary_table, get(subtype_dict[[subtype_resolution]]) != "snv")
  }

  if (depth_exists) {
    if (!is.null(precalc_depth_data) && !is.na(MutSeqR::denominator_dict[[subtype_resolution]])) {
      ref_values <- unique(depth_df[[MutSeqR::denominator_dict[[subtype_resolution]]]])
      if (!"N" %in% ref_values) {
        n_depth <- depth_df %>%
          dplyr::select(-MutSeqR::denominator_dict[[subtype_resolution]], -"subtype_depth") %>%
          unique()
        n_depth[[MutSeqR::denominator_dict[[subtype_resolution]]]] <- "N"
        n_depth$subtype_depth <- n_depth$group_depth
        depth_df <- rbind(n_depth, depth_df)
      }
    }
    summary_table <- dplyr::left_join(summary_table, depth_df)
  }

  # Replace NAs in data with 0s (sum and MF cols only)
  summary_table <- summary_table %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(summary_cols), ~ tidyr::replace_na(., 0)))


  # Calculate the proportions of each subtype
  if (subtype_resolution != "none") {
    # Get total sum of mutations across groups
    proportions <- summary_table %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(cols_to_group)))) %>%
      dplyr::mutate(
        total_group_mut_sum_min =
          sum(.data$sum_min),
        total_group_mut_sum_max =
          sum(.data$sum_max)
      ) %>%
      dplyr::ungroup()

    if (depth_exists) {
      # freq = mut_sum / total_mut_sum / depth
      # subtype_depth is = group_depth for non-snv mutations
      proportions <- proportions %>%
        dplyr::mutate(freq_min = .data$sum_min / .data$total_group_mut_sum_min / .data$subtype_depth,
                      freq_max = .data$sum_max / .data$total_group_mut_sum_max / .data$subtype_depth
        )

      # proportion = freq / total_freq  
      summary_table <- proportions %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(cols_to_group)))) %>%
        dplyr::mutate(
          total_freq_min = sum(.data$freq_min),
          total_freq_max = sum(.data$freq_max)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          proportion_min = .data$freq_min / .data$total_freq_min,
          proportion_max = .data$freq_max / .data$total_freq_max
        ) %>% # Remove extra columns
        dplyr::select(
          -"total_group_mut_sum_min", -"total_freq_min",
          -"total_group_mut_sum_max", -"total_freq_max",
          -"freq_min", -"freq_max"
        )
    } else {
      summary_table <- proportions %>%
        dplyr::mutate(proportion_min = .data$sum_min / .data$total_group_mut_sum_min,
                      proportion_max = .data$sum_max / .data$total_group_mut_sum_max) %>%
        dplyr::select(-"total_group_mut_sum_min", -"total_group_mut_sum_max")
    }
    # replace NAs with 0s (when sum of group = 0, dividing by 0 = NaN)
    summary_table <- summary_table %>%
      dplyr::mutate(dplyr::across(c(proportion_min, proportion_max), ~ tidyr::replace_na(., 0)))
  }
  # Remove the subtype_depth column if we are not using it
  if (depth_exists && subtype_resolution %in% c("none", "type")) {
    summary_table <- summary_table %>%
      dplyr::select(-"subtype_depth")
  }

  if (!summary) {
    return(mut_freq_table)
  } else {
    return(summary_table)
  }
}


#' Correct total depth for unique genomic positions
#'
#' @description
#' This function prevents double-counting of sequencing depth when multiple
#' mutations are reported at the same genomic position within the same sample.
#' For each group of rows with the same `sample`, `contig`, and `start`,
#' it retains the `total_depth` value for only one row and sets the rest to 0.
#' This is a necessary pre-processing step before calculating mutation frequency
#' to ensure the denominator (total bases sequenced) is not inflated.
#'
#' @details
#' Running this operation can be computationally intensive on large datasets.
#' By isolating it in a dedicated function, you can run it once and use the
#' resulting data frame for multiple downstream analyses, such as calls to
#' `calculate_mf` with different grouping variables.
#'
#' @param mutation_data Your mutation data. Must be a data.frame.
#' @param by_indel_priority A logical value. If TRUE, the `total_depth` value for
#' the row with the highest priority `variation_type` is retained. Priority is:
#' deletion, complex, insertion, snv, mnv, sv, uncategorised, no_variant.
#' If FALSE (default), the `total_depth` of the first row in the group is retained.
#'
#' @return A data.frame with the `total_depth` column corrected.
#' @examples
#' # Load example data
#' example_file <- system.file("extdata", "Example_files",
#'                             "example_mutation_data.rds",
#'                             package = "MutSeqR")
#' example_data <- readRDS(example_file)
#'
#' # Correct the depth
#' corrected_data <- correct_depth(mutation_data = example_data, by_indel_priority = TRUE)
#' # Now `corrected_data` can be used for filtering and frequency calculation
#'
#' @importFrom data.table as.data.table :=
#' @export
correct_depth <- function(mutation_data, by_indel_priority = FALSE) {

  if (!"total_depth" %in% colnames(mutation_data)) {
    stop("Error: 'total_depth' column not found in mutation_data.")
  }

  # Use data.table for high performance
  dt <- data.table::as.data.table(mutation_data)

  message("Correcting depth for unique genomic positions...")

  if (by_indel_priority) {
    # Define priority list for variation types
    variation_priority <- c("deletion", "complex", "insertion", "snv", "mnv",
                            "sv", "uncategorized", "no_variant")

    # Create a factor for ordering
    dt[, priority_order := factor(variation_type, levels = variation_priority, ordered = TRUE)]

    # Within each group, set total_depth to 0 for all but the highest priority row
    # The order() function on the priority_order column ensures the highest priority is first
    dt[, total_depth := {
      group_order <- order(priority_order, na.last = TRUE)
      corrected_depths <- rep(0, .N)
      corrected_depths[group_order[1]] <- total_depth[group_order[1]]
      corrected_depths
    }, by = .(sample, contig, start)]

    # Remove the temporary priority column
    dt[, priority_order := NULL]

  } else {
    # Simpler case: keep the first row's depth and zero out the rest
    # This is extremely fast in data.table
    dt[, total_depth := c(total_depth[1], rep(0, .N - 1)), by = .(sample, contig, start)]
  }

  corrected_rows <- nrow(dt[total_depth == 0]) - nrow(mutation_data[mutation_data$total_depth == 0, ])
  message(corrected_rows, " rows had their total_depth set to 0 to prevent double-counting.")
  
  # Return a data.frame for consistency with the rest of the package
  return(as.data.frame(dt))
}