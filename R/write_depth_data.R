#' Prepare depth tables from a single resolution-specific input
#'
#' @description Reads and validates a precalculated depth data frame
#' or delimited file, detects its subtype resolution from the
#' denominator context column, and derives all lower-resolution depth
#' tables that can be recovered from that input.
#' @param depth_data A data frame or path to a delimited file. It must contain
#'   the columns named in `group_cols`, `subtype_depth`, and exactly one context
#'   column from `denominator_dict` for `base_6`, `base_12`, `base_96`, or
#'   `base_192`.
#' @param d_sep A single-character delimiter used to read `depth_data` when it
#'   is a file path.
#' @param group_cols Character vector of columns that identify independent
#'   groups in `depth_data`. Defaults to `"sample"`; use other columns when
#'   depth data is already aggregated to those groups.
#'
#' @details
#' The input must contain the grouping columns, `subtype_depth`, and the
#' appropriate context column. Use `MutSeqR::denominator_dict` to see the
#' context column name associated with each subtype resolution. Each group must
#' have exactly one row for every context in `MutSeqR::context_list` at the
#' given resolution; contexts with zero depth must be included explicitly.
#' Duplicate group-context pairs, unknown contexts, missing contexts, or
#' non-finite per-group depth totals produce an error.
#'
#' *Base192* input produces Base96, Base12, Base6, and global depth tables.
#' *Base96* input produces Base6 and global tables, but cannot produce
#' Base12 because its normalized contexts have already combined
#' reverse-complement strands.
#' *Base12* input produces Base6 and global tables.
#' *Base6* input produces global depth.
#'
#' @return A named list of the depth tables at available subtype resolutions:
#' "base192", "base96", "base12", "base6", "global" (none). Each subtype
#' table contains the `group_cols`, its `denominator_dict` context column,
#' `subtype_depth`, and `group_depth`; the global table contains `group_cols`
#' and `group_depth`.
#' @importFrom dplyr across anti_join distinct group_by left_join summarise
#' @importFrom rlang .data
#' @importFrom tidyr crossing
#' @importFrom utils head read.delim
#' @export
write_depth_data <- function(depth_data, d_sep = "\t", group_cols = "sample") {
    if (is.character(depth_data)) {
        if (length(depth_data) != 1L || is.na(depth_data) ||
            !nzchar(depth_data)) {
            stop("depth_data must be a single file path.", call. = FALSE)
        }
        if (!file.exists(depth_data)) {
            stop("The depth_data file does not exist.", call. = FALSE)
        }
        if (file.info(depth_data)$size == 0) {
            stop("The depth_data file is empty.", call. = FALSE)
        }
        depth_data <- utils::read.delim(
            depth_data,
            sep = d_sep,
            header = TRUE,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    } else if (!is.data.frame(depth_data)) {
        stop("depth_data must be a data frame or a file path.", call. = FALSE)
    }

    if (nrow(depth_data) == 0L) {
        stop("depth_data must contain at least one row.", call. = FALSE)
    }
    if (anyDuplicated(names(depth_data))) {
        stop("depth_data must not contain duplicate column names.",
             call. = FALSE)
    }
    if (anyNA(depth_data)) {
        stop("depth_data must not contain any NA values.", call. = FALSE)
    }

    if (!is.character(group_cols) || length(group_cols) == 0L ||
        anyNA(group_cols) || any(!nzchar(group_cols)) ||
        anyDuplicated(group_cols)) {
        stop("group_cols must be a non-empty vector of unique column names.",
             call. = FALSE)
    }
    required <- c(group_cols, "subtype_depth")
    missing <- setdiff(required, names(depth_data))
    if (length(missing)) {
        stop("depth_data is missing required columns: ",
             paste(missing, collapse = ", "), call. = FALSE)
    }
    resolutions <- c("base_192", "base_96", "base_12", "base_6")
    context_columns <- unname(MutSeqR::denominator_dict[resolutions])
    context_columns <- context_columns[!is.na(context_columns)]
    present_context_columns <- intersect(context_columns, names(depth_data))
    if (length(present_context_columns) != 1L) {
        stop(
            "depth_data must contain exactly one resolution context column: ",
            paste(context_columns, collapse = ", "),
            call. = FALSE
        )
    }
    if (any(group_cols %in% c(present_context_columns, "group_depth", "subtype_depth"))) {
        stop("group_cols must not include a context or depth column.",
             call. = FALSE)
    }

    context_column <- present_context_columns[[1L]]
    resolution <- resolutions[
        match(context_column, unname(MutSeqR::denominator_dict[resolutions]))
    ]
    context_values <- as.character(depth_data[[context_column]])
    group_values <- depth_data[group_cols]
    subtype_depth <- depth_data$subtype_depth
    if (!is.numeric(subtype_depth) || any(!is.finite(subtype_depth)) ||
        any(subtype_depth < 0)) {
        stop("subtype_depth must contain finite, non-negative numeric values.",
             call. = FALSE)
    }
    empty_groups <- vapply(group_values, function(values) {
        is.character(values) && any(!nzchar(trimws(values)))
    }, logical(1))
    if (any(empty_groups)) {
        stop("Grouping values must not be empty.", call. = FALSE)
    }

    contexts <- as.character(MutSeqR::context_list[[resolution]])
    if (any(!context_values %in% contexts)) {
        unknown <- unique(context_values[!context_values %in% contexts])
        stop("depth_data contains unknown ", context_column, " values: ",
             paste(unknown, collapse = ", "), call. = FALSE)
    }
    keys <- group_values
    keys$context <- context_values
    if (anyDuplicated(keys)) {
        stop("depth_data must contain exactly one row per group and context; duplicate combinations found.",
             call. = FALSE)
    }
    groups <- dplyr::distinct(group_values)
    expected_keys <- tidyr::crossing(groups, context = contexts)
    missing_keys <- dplyr::anti_join(
        expected_keys, keys, by = c(group_cols, "context")
    )
    if (nrow(missing_keys)) {
        missing_groups <- do.call(
            paste,
            c(lapply(missing_keys[group_cols], as.character), sep = ", ")
        )
        missing_labels <- paste0(
            missing_groups, ": ", missing_keys$context
        )
        shown_labels <- utils::head(missing_labels, 10L)
        remaining <- length(missing_labels) - length(shown_labels)
        if (remaining > 0L) {
            shown_labels <- c(shown_labels,
                              paste0("and ", remaining, " more"))
        }
        stop(
            "depth_data is missing context values (group: context): ",
            paste(shown_labels, collapse = ", "),
            call. = FALSE
        )
    }

    source <- data.frame(
        group_values,
        context = context_values,
        subtype_depth = subtype_depth,
        check.names = FALSE
    )
    global <- source %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
        dplyr::summarise(group_depth = sum(.data$subtype_depth),
                         .groups = "drop")
    invalid_groups <- global[!is.finite(global$group_depth), group_cols,
                             drop = FALSE]
    if (nrow(invalid_groups)) {
        invalid_labels <- do.call(
            paste,
            c(lapply(invalid_groups, as.character), sep = ", ")
        )
        shown_groups <- utils::head(invalid_labels, 10L)
        remaining <- length(invalid_labels) - length(shown_groups)
        if (remaining > 0L) {
            shown_groups <- c(shown_groups,
                              paste0("and ", remaining, " more"))
        }
        stop(
            "Could not derive a finite group_depth: summing subtype_depth ",
            "overflowed for group(s) ", paste(shown_groups, collapse = "; "),
            ". Check these groups' subtype_depth values and retry.",
            call. = FALSE
        )
    }

    build_table <- function(context, output_resolution) {
        context_name <- MutSeqR::denominator_dict[[output_resolution]]
        table <- data.frame(
            source[group_cols],
            output_context = context,
            subtype_depth = source$subtype_depth,
            check.names = FALSE
        ) %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(
                c(group_cols, "output_context")
            ))) %>%
            dplyr::summarise(subtype_depth = sum(.data$subtype_depth),
                             .groups = "drop")
        names(table)[names(table) == "output_context"] <- context_name
        dplyr::left_join(table, global, by = group_cols)
    }
    pyrimidine <- function(context) {
        purine <- context %in% c("A", "G")
        context[purine] <- vapply(
            context[purine],
            function(base) reverseComplement(base, case = "upper"),
            character(1)
        )
        context
    }

    prepared <- list()
    if (resolution == "base_192") {
        prepared$base192 <- build_table(source$context, "base_192")
        normalized_context <- source$context
        purine_center <- substr(normalized_context, 2L, 2L) %in% c("A", "G")
        normalized_context[purine_center] <- reverseComplement(
            normalized_context[purine_center], case = "upper"
        )
        prepared$base96 <- build_table(normalized_context, "base_96")
        prepared$base12 <- build_table(substr(source$context, 2L, 2L),
                                        "base_12")
        prepared$base6 <- build_table(pyrimidine(substr(source$context, 2L, 2L)),
                                      "base_6")
    } else if (resolution == "base_96") {
        prepared$base96 <- build_table(source$context, "base_96")
        prepared$base6 <- build_table(substr(source$context, 2L, 2L), "base_6")
        warning(
            "Base12 depth cannot be derived from normalized Base96 contexts; ",
            "purine and pyrimidine strand counts are already combined.",
            call. = FALSE
        )
    } else if (resolution == "base_12") {
        prepared$base12 <- build_table(source$context, "base_12")
        prepared$base6 <- build_table(pyrimidine(source$context), "base_6")
    } else {
        prepared$base6 <- build_table(source$context, "base_6")
    }
    prepared$global <- global
    prepared
}

#' Prepare precalculated depth tables for the summary report
#'
#' @description Validates context-resolved precalculated depth inputs and
#' derives the global, base_6, and base_96 tables used by the summary report.
#' Explicit inputs at a requested resolution take precedence over derived
#' tables. When an input is omitted, the finest available context-resolved
#' table is used to derive it.
#'
#' @param base192 Optional base_192 depth data. Must be a data frame accepted
#'   by [write_depth_data()].
#' @param base96 Optional base_96 depth data. Must be a data frame accepted
#'   by [write_depth_data()].
#' @param global Optional global depth data with `sample` and `group_depth`
#'   columns.
#' @param base6 Optional base_6 depth data. Must be a data frame accepted by
#'   [write_depth_data()].
#'
#' @return A list with `global`, `base6`, `base96`, and `base192` elements.
#'   The context-resolved tables include `group_depth` in addition to their
#'   context and `subtype_depth` columns. Elements are `NULL` when neither an
#'   explicit input nor a compatible higher-resolution input is available.
#'
#' @keywords internal
prepare_report_depth_data <- function(base192 = NULL, base96 = NULL,
                                      global = NULL, base6 = NULL) {
    prepare_context_depth <- function(depth_data) {
        if (is.null(depth_data)) {
            return(NULL)
        }
        suppressWarnings(MutSeqR::write_depth_data(depth_data))
    }
    first_available <- function(...) {
        values <- list(...)
        available <- which(!vapply(values, is.null, logical(1)))
        if (!length(available)) {
            return(NULL)
        }
        values[[available[[1L]]]]
    }

    base192_depth <- prepare_context_depth(base192)
    base96_depth <- prepare_context_depth(base96)
    base6_depth <- prepare_context_depth(base6)

    list(
        global = first_available(
            global,
            base192_depth$global,
            base96_depth$global,
            base6_depth$global
        ),
        base6 = first_available(
            base6_depth$base6,
            base192_depth$base6,
            base96_depth$base6
        ),
        base96 = first_available(
            base96_depth$base96,
            base192_depth$base96
        ),
        base192 = if (is.null(base192_depth)) NULL else base192_depth$base192
    )
}
