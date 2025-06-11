#' plot_ci
#' @description Plot confidence intervals
#' @param data A data frame with the results of the BMD analysis.
#' Data must contain columns "Response", "BMD", "BMDL", and "BMDU".
#' BMD values can be NA.
#' @param order Indicates how the responses should be ordered. Options are
#' "none" (default), "asc" for ascending BMD values, "desc" for descending
#' BMD values, or a custom order.
#' @param custom_order A character vector with the custom order of the
#' Responses.
#' @param nudge A numeric value to nudge the text labels away from points.
#' Default is 0.3.
#' @param log_scale A logical value indicating if the x-axis should be in
#' log10 scale. Default is false.
#' @param x_lab A character string with the x-axis label. Default is "BMD" or
#' "log10(BMD)" if log_scale is TRUE.
#' @param y_lab A character string with the y-axis label. Default is "Response".
#' @param title A character string with the plot title. Default is "BMD with
#' 90% Confidence Intervals".
#' @return a ggplot object
#' @examples
#' # Plot results from PROAST and ToxicR
#' dat <- data.frame(Response = c("PROAST MF Min", "PROAST MF Max", "ToxicR MF Min", "ToxicR MF Max"),
#'                   BMD = c(NA, NA, 9.641894, 8.100164),
#'                   BMDL = c(7.38, 2.98, 8.032936,5.463013),
#'                   BMDU = c(10.9, 7.68, 10.97636, 10.04638))
#' plot <- plot_ci(dat)
#' @export
#' @importFrom dplyr arrange pull mutate group_by ungroup across where desc
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line geom_point theme_minimal theme
#' scale_color_manual scale_x_continuous geom_text labs ggplot2::element_text
#' element_rect element_blank element_line
plot_ci <- function(data,
                    order = "none",
                    custom_order = NULL,
                    nudge = 0.3,
                    log_scale = FALSE,
                    x_lab = NULL,
                    y_lab = NULL,
                    title = NULL) {

  if (order == "asc") {
    response_order <- data %>%
      dplyr::arrange(.data$BMD) %>%
      dplyr::pull(.data$Response)
    data$Response <- factor(data$Response, levels = response_order)
  } else if (order == "desc") {
    response_order <- data %>%
      dplyr::arrange(dplyr::desc(.data$BMD)) %>%
      dplyr::pull(.data$Response)
    data$Response <- factor(data$Response, levels = response_order)
  } else if (order == "custom") {
    if (is.null(custom_order)) {
      stop("order = 'custom' was selected but no custom_order was provided.")
    }
    data$Response <- factor(data$Response, levels = custom_order)
  }
  data <- data %>%
    dplyr::mutate(BMD = as.numeric(.data$BMD),
                  BMDL = as.numeric(.data$BMDL),
                  BMDU = as.numeric(.data$BMDU))
  if (log_scale) {
    data <- data %>%
      dplyr::mutate(BMD = log10(.data$BMD),
                    BMDL = log10(.data$BMDL),
                    BMDU = log10(.data$BMDU))
    results_bmd_df_plot <- data %>%
      dplyr::group_by(.data$Response) %>%
      dplyr::mutate(max = .data$BMDU) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, 2))) %>%
      tidyr::pivot_longer(cols = c("BMD", "BMDL", "BMDU"))
    if (is.null(x_lab)) {x_lab <- "log10(BMD)"}
  } else {
    results_bmd_df_plot <- data %>%
      dplyr::group_by(.data$Response) %>%
      dplyr::mutate(max = .data$BMDU) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, 1))) %>%
      tidyr::pivot_longer(cols = c("BMD", "BMDL", "BMDU"))
    if (is.null(x_lab)) {x_lab <- "BMD"}
  }
  nudge_value <- nudge

  if (is.null(y_lab)) {y_lab <- "Response"}
  if (is.null(title)) {title <- "BMD with 90% Confidence Intervals"}


  g <- ggplot(results_bmd_df_plot,
              ggplot2::aes(x = results_bmd_df_plot$value,
                           y = results_bmd_df_plot$Response,
                           color = results_bmd_df_plot$name)) +
    ggplot2::geom_line(ggplot2::aes(group = results_bmd_df_plot$Response),
                       color = "#b8b8b8", linewidth = 3.5, na.rm = TRUE) +
    ggplot2::geom_point(size = 3, na.rm = TRUE) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.y = ggplot2::element_text(color = "black"),
                   axis.text.x = ggplot2::element_text(color = "#000000"),
                   panel.border = ggplot2::element_rect(colour = "black",
                                                        fill = NA,
                                                        size = 1),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::scale_color_manual(values = c("black", "#BF2F24", "#436685")) +
    ggplot2::scale_x_continuous() +
    ggplot2::geom_text(ggplot2::aes(label = results_bmd_df_plot$value,
                                    color = results_bmd_df_plot$name),
                       size = 3.25,
                       nudge_x = dplyr::if_else(results_bmd_df_plot$value == results_bmd_df_plot$max,
                                                nudge_value, -nudge_value),
                       hjust = dplyr::if_else(results_bmd_df_plot$value == results_bmd_df_plot$max,
                                              0, 1), na.rm = TRUE) +
    ggplot2::labs(x = x_lab, y = y_lab,
                  title = title,
                  color = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0,
                                                       vjust = 0.5,
                                                       hjust = 0.5),
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(axis.ticks = ggplot2::element_line(color = "black",
                                                      linewidth = 0.5))

  return(g)
}