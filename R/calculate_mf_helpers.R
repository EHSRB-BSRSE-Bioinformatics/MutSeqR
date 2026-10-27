#' Read precalculated depth input
#' @param depth_data file path or data frame
#' @param d_sep seperator for reading depth data from file.
#' @keywords internal
#' @noRd
read_precalc_depth <- function(depth_data, d_sep) {
    if (is.data.frame(depth_data)) {
        return(depth_data)
    }
    if (!is.character(depth_data) || length(depth_data) != 1L ||
        is.na(depth_data) || !nzchar(depth_data)) {
        stop("precalc_depth_data must be a data frame or a single file path.",
             call. = FALSE)
    }
    if (!file.exists(depth_data)) {
        stop("The precalc_depth_data does not exist.", call. = FALSE)
    }
    if (file.info(depth_data)$size == 0) {
        stop("You are trying to import an empty precalc_depth_data",
          call. = FALSE
        )
    }
    depth_df <- utils::read.delim(
        depth_data,
        sep = d_sep,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    if (ncol(depth_df) <= 1L) {
        stop("Your imported precalc only has one column. You may want to set ",
             "d_sep to properly reflect the delimiter used for the data.",
             call. = FALSE)
    }
  return(depth_df)
}

#' Infer and map the grouping columns for precalculated depth.
#' @description Adds the grouping columns to the precal depth data.
#' @keywords internal
#' @noRd
map_precalc_depth_groups <- function(depth_df, mutation_data, cols_to_group) {
    if (all(cols_to_group %in% names(depth_df))) {
        return(depth_df)
    }
    if (!"sample" %in% names(depth_df)) {
        stop(
            "precalc_depth_data must contain all cols_to_group columns or a ",
            "'sample' column for mapping: ", paste(cols_to_group, collapse = ", "),
            call. = FALSE
        )
    }

    sample_groups <- dplyr::distinct(
        mutation_data[unique(c("sample", cols_to_group))]
    )
    sample_groups$sample <- as.character(sample_groups$sample)
    depth_df$sample <- as.character(depth_df$sample)

    group_count <- nrow(dplyr::distinct(mutation_data[cols_to_group]))
    sample_count <- nrow(dplyr::distinct(mutation_data["sample"]))
    if (group_count > sample_count) {
        stop(
            "Cannot map sample-level precalc_depth_data: cols_to_group has ",
            "a higher resolution than sample. Supply depth data at the ",
            "requested groups.",
            call. = FALSE
        )
    }

    ambiguous_samples <- sample_groups %>%
        dplyr::group_by(.data$sample) %>%
        dplyr::summarise(n_groups = dplyr::n(), .groups = "drop") %>%
        dplyr::filter(.data$n_groups > 1L)
    if (nrow(ambiguous_samples)) {
        stop(
            "Cannot map sample-level precalc_depth_data: each sample must ",
            "identify exactly one combination of cols_to_group (",
            paste(cols_to_group, collapse = ", "), "). Ambiguous sample(s): ",
            paste(utils::head(ambiguous_samples$sample, 10L), collapse = ", "),
            call. = FALSE
        )
    }

    absent_samples <- setdiff(unique(depth_df$sample), sample_groups$sample)
    if (length(absent_samples)) {
        warning(
            "samples in precalc_depth_data are ",
            "absent from mutation_data: ",
            paste(utils::head(absent_samples, 10L), collapse = ", "),
            call. = FALSE
        )
    }
    missing_samples <- setdiff(sample_groups$sample, unique(depth_df$sample))
    if (length(missing_samples)) {
        stop(
            "Cannot map sample-level precalc_depth_data: sample(s) in ",
            "mutation_data are absent from precalc_depth_data: ",
            paste(utils::head(missing_samples, 10L), collapse = ", "),
            call. = FALSE
        )
    }

    supplied_groups <- setdiff(intersect(cols_to_group, names(depth_df)), "sample")
    mapped <- sample_groups[match(depth_df$sample, sample_groups$sample), ,
                            drop = FALSE]
    for (column in supplied_groups) {
        supplied <- depth_df[[column]]
        expected <- mapped[[column]]
        if (is.factor(supplied)) supplied <- as.character(supplied)
        if (is.factor(expected)) expected <- as.character(expected)
        equal <- (is.na(supplied) & is.na(expected)) |
            (!is.na(supplied) & !is.na(expected) & supplied == expected)
        if (any(!equal)) {
            stop(
                "Cannot map sample-level precalc_depth_data: grouping column '",
                column, "' disagrees with mutation_data for sample(s): ",
                paste(utils::head(unique(depth_df$sample[!equal]), 10L),
                      collapse = ", "),
                call. = FALSE
            )
        }
    }

    depth_df <- depth_df[
        setdiff(names(depth_df), setdiff(supplied_groups, "sample"))
    ]
    depth_df <- dplyr::left_join(depth_df, sample_groups, by = "sample")
  return(depth_df)
}

#' Prepare precalculated depth at the requested group and subtype resolutions
#' @keywords internal
#' @noRd
prepare_precalc_depth <- function(depth_df, mutation_data, cols_to_group,
                                  subtype_resolution, d_sep) {
    depth_df <- map_precalc_depth_groups(
        depth_df, mutation_data, cols_to_group
    )

    input_resolutions <- c("base_192", "base_96", "base_12", "base_6")
    context_columns <- unname(
        MutSeqR::denominator_dict[input_resolutions]
    )
    context_columns <- context_columns[!is.na(context_columns)]
    present_context_columns <- intersect(context_columns, names(depth_df))
    if (length(present_context_columns) > 1L) {
        stop(
            "precalc_depth_data must contain exactly one context column when ",
            "using automatic resolution conversion.",
            call. = FALSE
        )
    }

    if (length(present_context_columns) == 1L) {
        context_column <- present_context_columns[[1L]]
        if (!"subtype_depth" %in% names(depth_df)) {
            stop("Context-resolved precalc_depth_data must contain ",
                 "'subtype_depth'.", call. = FALSE)
        }
        if (!is.numeric(depth_df$subtype_depth) ||
            any(!is.finite(depth_df$subtype_depth)) ||
            any(depth_df$subtype_depth < 0)) {
            stop("subtype_depth must contain finite, non-negative numeric values.",
                 call. = FALSE)
        }
        depth_df <- depth_df %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(
                c(cols_to_group, context_column)
            ))) %>%
            dplyr::summarise(
                subtype_depth = sum(.data$subtype_depth),
                .groups = "drop"
            )
        converted_depth <- suppressWarnings(
            MutSeqR::write_depth_data(
                depth_df,
                d_sep = d_sep,
                group_cols = cols_to_group
            )
        )
        depth_key <- switch(
            subtype_resolution,
            none = "global",
            type = "global",
            base_6 = "base6",
            base_12 = "base12",
            base_96 = "base96",
            base_192 = "base192"
        )
        depth_df <- converted_depth[[depth_key]]
        if (is.null(depth_df)) {
            requested_context_column <-
                MutSeqR::denominator_dict[[subtype_resolution]]
            stop(
                "precalc_depth_data at a finer context resolution cannot ",
                "generate depths for subtype_resolution '", subtype_resolution,
                "'. Supply data at 'base_192' or another resolution fine ",
                "enough to derive '", requested_context_column, "'.",
                call. = FALSE
            )
        }
    } else {
        if (!"group_depth" %in% names(depth_df)) {
            stop("precalc_depth_data is missing 'group_depth'.", call. = FALSE)
        }
        if (!is.numeric(depth_df$group_depth) ||
            any(!is.finite(depth_df$group_depth)) ||
            any(depth_df$group_depth < 0)) {
            stop("group_depth must contain finite, non-negative numeric values.",
                 call. = FALSE)
        }
        depth_df <- depth_df %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(cols_to_group))) %>%
            dplyr::summarise(
                group_depth = sum(.data$group_depth),
                .groups = "drop"
            )
        if (any(!is.finite(depth_df$group_depth))) {
            stop("Summed group_depth is not finite for one or more groups.",
                 call. = FALSE)
        }
    }

    expected_groups <- dplyr::distinct(mutation_data[cols_to_group]) %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
    depth_groups <- dplyr::distinct(depth_df[cols_to_group]) %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
    missing_groups <- dplyr::anti_join(
        expected_groups, depth_groups, by = cols_to_group
    )
    if (nrow(missing_groups)) {
        missing_labels <- apply(as.data.frame(missing_groups), 1L, paste,
                                collapse = ", ")
        stop(
            "precalc_depth_data has no depth value for cols_to_group group(s): ",
            paste(utils::head(missing_labels, 10L), collapse = "; "),
            call. = FALSE
        )
    }
  return(depth_df)
}