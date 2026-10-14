#' Read precalculated depth input
#' @keywords internal
#' @noRd
read_precalc_depth <- function(depth_data, d_sep) {
    if (is.data.frame(depth_data)) {
        depth_df <- depth_data
    } else if (is.character(depth_data)) {
        depth_file <- file.path(depth_data)
        if (!file.exists(depth_file)) {
            stop("The precalc_depth_data does not exist.")
        }
        if (file.info(depth_file)$size == 0) {
            stop(
                "You are trying to import an empty precalc_depth_data"
            )
        }
        depth_df <- read.delim(file.path(depth_file),
            sep = d_sep,
            header = TRUE
        )
        if (ncol(depth_df) <= 1) {
            stop("Your imported precalc only has one column.
                You may want to set d_sep to properly reflect
                the delimiter used for the data you are importing.")
        }
    }
    depth_df
}

#' Prepare requested depth resolution and grouping
#' @keywords internal
#' @noRd
prepare_precalc_depth <- function(depth_df, mutation_data, cols_to_group,
                                   subtype_resolution, d_sep,
                                   depth_group_cols = NULL) {
    # Determine the resolution of the input depth data
    input_resolutions <- c(
        "base_192", "base_96", "base_12", "base_6"
    )
    input_context_columns <- unname(
        MutSeqR::denominator_dict[input_resolutions]
    )
    input_context_columns <- input_context_columns[
        !is.na(input_context_columns)
    ]
    present_context_columns <- intersect(
        input_context_columns, names(depth_df)
    )
    explicit_groups <- !is.null(depth_group_cols)
    sample_input <- explicit_groups && identical(depth_group_cols, "sample")
    supplied_depth <- depth_df
    if (explicit_groups) {
        missing_keys <- setdiff(depth_group_cols, names(depth_df))
        if (length(missing_keys)) {
            stop("precalc_depth_data is missing depth_group_cols: ",
                 paste(missing_keys, collapse = ", "), call. = FALSE)
        }
        if (!sample_input && !setequal(depth_group_cols, cols_to_group)) {
            stop("Automatic pooling only supports depth_group_cols = 'sample'. ",
                 "Summarize the depth table to cols_to_group first: ",
                 paste(cols_to_group, collapse = ", "), call. = FALSE)
        }
    }
    if (length(present_context_columns) == 0L) {
        if (sample_input && subtype_resolution %in% c("none", "type")) {
            depth_df <- pool_sample_depth(
                depth_df, depth_df, mutation_data, cols_to_group,
                validate_global = TRUE, supplied_depth = supplied_depth,
                warn_missing = TRUE
            )
        }
        return(depth_df)
    }
    if (length(present_context_columns) > 0L) {
        if (length(present_context_columns) != 1L) {
            stop(
                "precalc_depth_data must contain exactly one context ",
                "column when using automatic resolution conversion.",
                call. = FALSE
            )
        }
        input_context_column <- present_context_columns[[1L]]
        input_resolution <- input_resolutions[
            match(
                input_context_column,
                unname(MutSeqR::denominator_dict[input_resolutions])
            )
        ]
        if (explicit_groups) {
            input_group_cols <- depth_group_cols
        } else if (all(cols_to_group %in% names(depth_df))) {
            input_group_cols <- cols_to_group
        } else if ("sample" %in% names(depth_df)) {
            input_group_cols <- "sample"
        } else {
            stop(
                "Context-resolved precalc_depth_data must contain ",
                "either 'sample' or all cols_to_group columns: ",
                paste(cols_to_group, collapse = ", "),
                call. = FALSE
            )
        }
        if (!explicit_groups && setequal(input_group_cols, cols_to_group) &&
            anyDuplicated(depth_df[c(input_group_cols, input_context_column)])) {
            stop("Your precalc_depth_data is at a finer grouping than required ",
                 "(contains multiple rows per combination of: ",
                 paste(c(cols_to_group, input_context_column), collapse = ", "),
                 "). Adding metadata columns does not aggregate sample depths. ",
                 "For sample-level depths, set depth_group_cols = 'sample'; ",
                 "otherwise summarize the depth first.", call. = FALSE)
        }
        converted_depth <- suppressWarnings(
            MutSeqR::write_depth_data(
                depth_df,
                d_sep = d_sep,
                group_cols = input_group_cols
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
                "precalc_depth_data at '", input_resolution,
                "' resolution cannot generate depths for subtype_resolution '",
                subtype_resolution, "'. Supply data at ",
                "'base_192' or another resolution fine enough to ",
                "derive '", requested_context_column, "'.",
                call. = FALSE
            )
        }

        if (sample_input || !all(cols_to_group %in% names(depth_df))) {
            context_column <- if (depth_key == "global") NULL else
                MutSeqR::denominator_dict[[subtype_resolution]]
            depth_df <- pool_sample_depth(
                depth_df, converted_depth$global, mutation_data,
                cols_to_group, context_column,
                validate_global = sample_input,
                supplied_depth = if (sample_input) supplied_depth else NULL,
                warn_missing = sample_input
            )
        }
    }
    return(depth_df)
}

#' Pool sample depths using a unique sample-to-group mapping
#' @keywords internal
#' @noRd
pool_sample_depth <- function(depth_df, global_depth, mutation_data,
                              cols_to_group, context_column = NULL,
                              validate_global = FALSE, supplied_depth = NULL,
                              warn_missing = FALSE) {
    if (validate_global) {
        duplicate_samples <- duplicated(global_depth$sample) |
            duplicated(global_depth$sample, fromLast = TRUE)
        if (any(duplicate_samples)) {
            stop("Cannot aggregate per-sample precalc_depth_data: supply exactly ",
                 "one depth row per sample; summarize other finer depth tables first. ",
                 "Duplicate sample(s): ",
                 paste(utils::head(unique(global_depth$sample[duplicate_samples]), 10L),
                       collapse = ", "),
                 call. = FALSE)
        }
        values <- global_depth$group_depth
        if (is.null(values) || !is.numeric(values) ||
            any(!is.finite(values)) || any(values < 0)) {
            stop("Cannot aggregate per-sample precalc_depth_data: group_depth ",
                 "must contain finite, non-negative numeric values.", call. = FALSE)
        }
    }
    group_map_columns <- unique(c("sample", cols_to_group))
    missing_group_columns <- setdiff(
        group_map_columns, names(mutation_data)
    )
    if (length(missing_group_columns)) {
        stop(
            "Cannot aggregate per-sample precalc_depth_data ",
            "to cols_to_group; mutation_data is missing: ",
            paste(missing_group_columns, collapse = ", "),
            call. = FALSE
        )
    }
    sample_groups <- mutation_data %>%
        dplyr::select(dplyr::all_of(group_map_columns)) %>%
        dplyr::distinct()
    ambiguous_samples <- sample_groups %>%
        dplyr::group_by(.data$sample) %>%
        dplyr::summarise(n_groups = dplyr::n(),
                         .groups = "drop") %>%
        dplyr::filter(.data$n_groups > 1L)
    if (nrow(ambiguous_samples)) {
        stop(
            "Cannot aggregate per-sample precalc_depth_data: ",
            "each sample must map to exactly one combination ",
            "of cols_to_group (", paste(cols_to_group, collapse = ", "),
            "). Supply depths at those groups. Ambiguous sample(s): ",
            paste(utils::head(ambiguous_samples$sample, 10L),
                  collapse = ", "),
            call. = FALSE
        )
    }
    unmapped_samples <- setdiff(
        as.character(global_depth$sample),
        as.character(sample_groups$sample)
    )
    if (length(unmapped_samples)) {
        stop(
            "Cannot aggregate per-sample precalc_depth_data: ",
            "sample(s) are absent from mutation_data: ",
            paste(utils::head(unmapped_samples, 10L),
                  collapse = ", "),
            call. = FALSE
        )
    }

    if (validate_global) {
        if (is.null(supplied_depth)) supplied_depth <- depth_df
        supplied_groups <- setdiff(
            intersect(cols_to_group, names(supplied_depth)), "sample"
        )
        if (length(supplied_groups)) {
            mapped <- sample_groups[
                match(supplied_depth$sample, sample_groups$sample), , drop = FALSE
            ]
            conflicts <- rep(FALSE, nrow(supplied_depth))
            for (column in supplied_groups) {
                supplied <- supplied_depth[[column]]
                expected <- mapped[[column]]
                # Factor level sets need not be identical across the two inputs.
                if (is.factor(supplied)) supplied <- as.character(supplied)
                if (is.factor(expected)) expected <- as.character(expected)
                equal <- (is.na(supplied) & is.na(expected)) |
                    (!is.na(supplied) & !is.na(expected) & supplied == expected)
                conflicts <- conflicts | !equal
            }
            if (any(conflicts)) {
                stop("Cannot aggregate per-sample precalc_depth_data: grouping ",
                     "columns disagree with mutation_data for sample(s): ",
                     paste(utils::head(supplied_depth$sample[conflicts], 10L), collapse = ", "),
                     ". Requested groups: ", paste(cols_to_group, collapse = ", "),
                     call. = FALSE)
            }
        }
        # Discard redundant metadata before joining the authoritative mapping.
        global_depth <- global_depth[c("sample", "group_depth")]
        if (!is.null(context_column)) {
            depth_df <- depth_df[c("sample", context_column, "subtype_depth")]
        }
    }
    if (warn_missing) {
        missing_samples <- setdiff(as.character(sample_groups$sample),
                                   as.character(global_depth$sample))
        if (length(missing_samples)) {
            affected <- sample_groups[
                as.character(sample_groups$sample) %in% missing_samples, , drop = FALSE
            ]
            labels <- apply(as.data.frame(affected[cols_to_group]), 1L,
                            function(row) paste(row, collapse = "/"))
            warning("No precalculated depth supplied for mutation sample(s): ",
                    paste(utils::head(missing_samples, 10L), collapse = ", "),
                    ". Affected groups (", paste(cols_to_group, collapse = ", "),
                    "): ", paste(utils::head(unique(labels), 10L), collapse = ", "),
                    ". Supplied depths alone contribute to pooled denominators; ",
                    "groups without supplied depth retain missing-depth behavior.",
                    call. = FALSE)
        }
    }
    global_depth <- dplyr::left_join(
        global_depth,
        sample_groups,
        by = "sample"
    ) %>%
        dplyr::group_by(
            dplyr::across(dplyr::all_of(cols_to_group))
        ) %>%
        dplyr::summarise(
            group_depth = sum(.data$group_depth),
            .groups = "drop"
        )
    if (validate_global && any(!is.finite(global_depth$group_depth))) {
        stop("Cannot aggregate per-sample precalc_depth_data: pooled group_depth ",
             "is not finite. Check depths for groups: ",
             paste(cols_to_group, collapse = ", "), call. = FALSE)
    }
    if (is.null(context_column)) {
        depth_df <- global_depth
    } else {
        depth_df <- dplyr::left_join(
            depth_df, sample_groups, by = "sample"
        ) %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(
                c(cols_to_group, context_column)
            ))) %>%
            dplyr::summarise(
                subtype_depth = sum(.data$subtype_depth),
                .groups = "drop"
            ) %>%
            dplyr::left_join(global_depth, by = cols_to_group)
    }
    depth_df
}
