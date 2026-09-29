#' Prepare depth tables from a single resolution-specific input
#'
#' @description Reads and validates a precalculated depth data frame
#' or delimited file, detects its subtype resolution from the
#' denominator context column, and derives all lower-resolution depth
#' tables that can be recovered from that input.
#' @param depth_data A data frame or path to a delimited file. It must contain
#'   `sample`, `subtype_depth`, and exactly one context column from
#'   `denominator_dict` for `base_6`, `base_12`, `base_96`, or `base_192`.
#' @param d_sep A single-character delimiter used to read `depth_data` when it
#'   is a file path.
#'
#' @details
#' The input must contain columns "sample", "subtype_depth" and the appropriate
#' context column. Use `MutSeqR::denominator_dict` to see the context column
#' name associated with each subtype resolution. Each sample must have exactly
#' one row for every context in `MutSeqR::context_list` at the given
#' resolution; contexts with zero depth must be included explicitly. Duplicate
#' sample-context pairs, unknown contexts, missing contexts, or non-finite
#' per-sample depth totals produce an error.
#'
#' *Base192* input produces Base96, Base12, Base6, and global depth tables.
#' *Base96* input produces Base6 and global tables, but cannot produce
#' Base12 because its normalized contexts have already combined
#' reverse-complement strands.
#' *Base12* input produces Base6 and global tables.
#' *Base6* input produces global depth.
#'
#' @return A named list of the depth tables at available subtype resolutions:
#' "base192", "base96", "base12", "base6", "global" (none). Each subtype table
#' contains `sample`, its `denominator_dict` context column, `subtype_depth`,
#' and `group_depth`; the global table contains `sample` and `group_depth`.
#' @importFrom dplyr anti_join group_by left_join summarise
#' @importFrom rlang .data
#' @importFrom tidyr expand_grid
#' @importFrom utils head read.delim
#' @export
write_depth_data <- function(depth_data, d_sep) {
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

    required <- c("sample", "subtype_depth")
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

    context_column <- present_context_columns[[1L]]
    resolution <- resolutions[
        match(context_column, unname(MutSeqR::denominator_dict[resolutions]))
    ]
    context_values <- as.character(depth_data[[context_column]])
    samples <- as.character(depth_data$sample)
    subtype_depth <- depth_data$subtype_depth
    if (!is.numeric(subtype_depth) || any(!is.finite(subtype_depth)) ||
        any(subtype_depth < 0)) {
        stop("subtype_depth must contain finite, non-negative numeric values.",
             call. = FALSE)
    }
    if (any(!nzchar(trimws(samples)))) {
        stop("sample values must not be empty.", call. = FALSE)
    }

    contexts <- as.character(MutSeqR::context_list[[resolution]])
    if (any(!context_values %in% contexts)) {
        unknown <- unique(context_values[!context_values %in% contexts])
        stop("depth_data contains unknown ", context_column, " values: ",
             paste(unknown, collapse = ", "), call. = FALSE)
    }
    keys <- data.frame(sample = samples, context = context_values)
    if (anyDuplicated(keys)) {
        stop("depth_data must contain exactly one row per sample and context; duplicate combinations found.",
             call. = FALSE)
    }
    expected_keys <- tidyr::expand_grid(
        sample = unique(samples),
        context = contexts
    )
    missing_keys <- dplyr::anti_join(expected_keys, keys,
                                     by = c("sample", "context"))
    if (nrow(missing_keys)) {
        missing_labels <- paste0(
            missing_keys$sample, ": ", missing_keys$context
        )
        shown_labels <- utils::head(missing_labels, 10L)
        remaining <- length(missing_labels) - length(shown_labels)
        if (remaining > 0L) {
            shown_labels <- c(shown_labels,
                              paste0("and ", remaining, " more"))
        }
        stop(
            "depth_data is missing context values (sample: context): ",
            paste(shown_labels, collapse = ", "),
            call. = FALSE
        )
    }

    source <- data.frame(
        sample = samples,
        context = context_values,
        subtype_depth = subtype_depth,
        stringsAsFactors = FALSE
    )
    global <- source %>%
        dplyr::group_by(.data$sample) %>%
        dplyr::summarise(group_depth = sum(.data$subtype_depth),
                         .groups = "drop")
    invalid_samples <- as.character(
        global$sample[!is.finite(global$group_depth)]
    )
    if (length(invalid_samples)) {
        shown_samples <- utils::head(invalid_samples, 10L)
        remaining <- length(invalid_samples) - length(shown_samples)
        if (remaining > 0L) {
            shown_samples <- c(shown_samples,
                               paste0("and ", remaining, " more"))
        }
        stop(
            "Could not derive a finite group_depth: summing subtype_depth ",
            "overflowed for sample(s) ", paste(shown_samples, collapse = ", "),
            ". Check these samples' subtype_depth values and retry.",
            call. = FALSE
        )
    }

    build_table <- function(context, output_resolution) {
        context_name <- MutSeqR::denominator_dict[[output_resolution]]
        table <- data.frame(
            sample = source$sample,
            output_context = context,
            subtype_depth = source$subtype_depth,
            stringsAsFactors = FALSE
        ) %>%
            dplyr::group_by(.data$sample, .data$output_context) %>%
            dplyr::summarise(subtype_depth = sum(.data$subtype_depth),
                             .groups = "drop")
        names(table)[names(table) == "output_context"] <- context_name
        dplyr::left_join(table, global, by = "sample")
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