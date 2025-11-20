#' Cleveland Plot
#' @description Make a Cleveland plot for the PROAST results.
#' @param results PROAST results object.
#' @param covariate_col Covariate column name.
#' @param output_path Output path for the plot. If the output_path doesn't
#' exist, it will be created. If NULL, the plots will not be exported.
#' @return A single ggplot object with facets for each response.
#' @importFrom dplyr arrange filter mutate pull rename group_by summarise ungroup if_else
#' @importFrom ggplot2 ggplot aes geom_errorbar geom_point scale_size_continuous theme element_blank element_line xlab ylab ggtitle facet_wrap ggsave
#' @importFrom rlang sym
cleveland_plot <- function(results,
                           covariate_col = NULL,
                           output_path = NULL) {

    results_df <- results[[2]] %>%
        dplyr::mutate(
            CED = as.numeric(.data$CED),
            CEDL = as.numeric(.data$CEDL),
            CEDU = as.numeric(.data$CEDU),
            weights = dplyr::if_else(.data$Selected.Model == "Model averaging", 1, .data$weights)
        ) %>%
        dplyr::rename(Model = "Selected.Model")

    if (!is.null(covariate_col)) {
        results_df <- results_df %>%
            dplyr::mutate(
                Model = paste(.data$Model, .data[[covariate_col]])
            )
    }

    model_order_global <- results_df %>%
        dplyr::group_by(.data$Model) %>%
        dplyr::summarise(mean_weight = mean(.data$weights, na.rm = TRUE)) %>%
        dplyr::arrange(.data$mean_weight) %>%
        dplyr::pull(.data$Model) %>%
        unique()

    results_df <- results_df %>%
        dplyr::mutate(Model = factor(.data$Model, levels = model_order_global))

    p_faceted <- ggplot(results_df, aes(x = .data$CED, y = .data$Model)) +
        geom_errorbar(aes(xmin = .data$CEDL, xmax = .data$CEDU),
            color = "gray",
            width = 0.1
        ) +
        geom_point(aes(size = .data$weights), color = "red") +
        ggplot2::scale_size_continuous(guide = "none", range = c(0.5, 3)) +
        ggplot2::theme(
            panel.background = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(),
            panel.grid = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_line(),
            axis.ticks.y = ggplot2::element_line(),
            strip.background = element_rect(fill = "lightgray", color = NA),
            strip.text = element_text(face = "bold")
        ) +
        ggplot2::xlab("BMD Estimate") +
        ggplot2::ylab("Model") +
        ggplot2::ggtitle("BMD by Selected Model (Sorted by Weights) per Response") +
        ggplot2::facet_wrap(~Response, scales = "free_x", ncol = 2)

    if (!is.null(output_path)) {
        # Create directory if it doesn't exist
        if (!dir.exists(output_path)) {
            dir.create(output_path, recursive = TRUE)
        }
        file_name <- file.path(output_path, "PROAST_AllResponses_ClevelandPlot.svg")
        ggsave(
            filename = file_name,
            plot = p_faceted,
            device = "svg",
            width = 12,
            height = 8 + (length(unique(results_df$Response)) * 0.5),
            units = "in",
            create.dir = FALSE
        )
    }

    return(p_faceted)
}