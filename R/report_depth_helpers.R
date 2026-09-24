# Prepare report depth tables at all supported mutation-context resolutions.
#
#' @keywords internal
prepare_report_depth_data <- function(base96 = NULL,
                                      global = NULL,
                                      base6 = NULL) {
    if (is.null(base96)) {
        return(list(global = global, base6 = base6, base96 = NULL))
    }

    required <- c("sample", "normalized_context", "subtype_depth")
    missing <- setdiff(required, names(base96))
    if (length(missing)) {
        stop("Base96 depth data is missing required columns: ",
             paste(missing, collapse = ", "), call. = FALSE)
    }
    if (nrow(base96) == 0L) {
        stop("Base96 depth data must contain at least one sample.",
             call. = FALSE)
    }
    if (!is.numeric(base96$subtype_depth) ||
        any(!is.finite(base96$subtype_depth)) ||
        any(base96$subtype_depth < 0)) {
        stop("Base96 subtype_depth must contain finite, non-negative numeric values.",
             call. = FALSE)
    }
    if (anyNA(base96$sample) || anyNA(base96$normalized_context)) {
        stop("Base96 sample and normalized_context values cannot be missing.",
             call. = FALSE)
    }

    contexts <- MutSeqR::context_list$base_96
    unknown <- setdiff(unique(as.character(base96$normalized_context)), contexts)
    if (length(unknown)) {
        stop("Base96 depth data contains unknown normalized_context values: ",
             paste(unknown, collapse = ", "), call. = FALSE)
    }
    base96$sample <- as.character(base96$sample)
    base96$normalized_context <- as.character(base96$normalized_context)
    duplicate <- duplicated(base96[c("sample", "normalized_context")])
    if (any(duplicate)) {
        stop("Base96 depth data must contain exactly one row per sample and normalized_context; duplicate combinations found.",
             call. = FALSE)
    }
    missing_context <- base96 |>
        dplyr::distinct(sample) |>
        tidyr::crossing(normalized_context = contexts) |>
        dplyr::anti_join(base96, by = c("sample", "normalized_context"))
    if (nrow(missing_context)) {
        stop("Base96 depth data must contain all 32 normalized_context values for every sample; missing combinations found.",
             call. = FALSE)
    }

    if ("group_depth" %in% names(base96) &&
        (!is.numeric(base96$group_depth) ||
         any(!is.finite(base96$group_depth)) ||
         any(base96$group_depth < 0))) {
        stop("Base96 group_depth must contain finite, non-negative numeric values.",
             call. = FALSE)
    }

    global_derived <- base96 |>
        dplyr::group_by(sample) |>
        dplyr::summarise(group_depth = sum(subtype_depth), .groups = "drop")

    if (is.null(global)) {
        global <- global_derived
    }

    if (is.null(base6)) {
        base6 <- base96 |>
            dplyr::mutate(normalized_ref = substr(normalized_context, 2, 2)) |>
            dplyr::group_by(sample, normalized_ref) |>
            dplyr::summarise(subtype_depth = sum(subtype_depth), .groups = "drop") |>
            dplyr::left_join(global, by = "sample")
    } else if (!"group_depth" %in% names(base6)) {
        base6 <- dplyr::left_join(base6, global, by = "sample")
    }

    # Older base96 files can omit group_depth. Keep a provided denominator
    # unchanged for compatibility with files that already carry one.
    if (!"group_depth" %in% names(base96)) {
        base96 <- dplyr::left_join(base96, global, by = "sample")
    }

    list(global = global, base6 = base6, base96 = base96)
}
