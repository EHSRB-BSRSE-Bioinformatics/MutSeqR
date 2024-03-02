#library(ggplot2)
#library(tidyverse)

#' Create a heatmap plot
#'
#' This function creates a heatmap plot using the provided data file.
#'
#' @param data The path to the data file.
#' @param group_var The variable to group by.
#' @param mut_proportion_scale The scale option for the mutation proportion. Options are passed to viridis::scale_fill_viridis_c. One of # inferno, magma, plasma, viridis, cividis, turbo, mako, or rocket. We highly reccomend the default for its ability to disciminate hard to see patterns. (Default: "turbo")
#' @param max Maximum value used for plotting the relative contributions.
#' Contributions that are higher will have the maximum colour. (Default: 0.2)
#' @param rescale_data Logical value indicating whether to rescale the mutation proportions to increase the dynamic range of colors shown on the plot. (Default: TRUE)
#' @param condensed More condensed plotting format. Default = F.
#' @import ggplot2
#' @importFrom dplyr %>% group_by summarise mutate
#' @importFrom scales rescale
#'
#' @return A ggplot object representing the heatmap plot.
#'
#' @examples
#' create_heatmap(mutation_data, dose, "inferno")
spectra_comparison_heatmap <- function(data,
                           group_var,
                           mut_proportion_scale = "turbo",
                           max = 0.2,
                           rescale_data = TRUE,
                           condensed = FALSE) { # nolint: line_length_linter.

    ### For production replace this with an object that already exists,
    # i.e., data frame with mutation proportions
    df <- read.table(
        file = data,
        header = TRUE,
        sep = "\t"
    )
    df <- df %>%
        separate(SampleName, c("Lab", "Sample", "Replicate", "Downsampling"), sep = "_", remove = FALSE) %>%
        mutate(Downsampling = ifelse(is.na(Downsampling), "None", Downsampling))
    df <- df[df$Downsampling == "None", ]
    ########################################################################################

    context_size <- stringr::str_length(df$Context[1])
    if (context_size == 1) {
        axis_size <- 10
    } else if (context_size == 2) {
        axis_size <- 8
    } else if (context_size == 3) {
        axis_size <- 4
    } else {
        axis_size <- 3
    }

    # Change plotting parameters based on whether plot should be condensed.
    if (condensed == TRUE) {
        spacing <- 0
    } else {
        spacing <- 0.5
    }

    # Get counts of mutations per subtyope for facet labels
    # TODO refactor (replace MutDepth and Subtype with right column names)
    mut_counts <- df %>%
        dplyr::group_by(Subtype) %>%
        dplyr::summarise(nrmuts = sum(MutDepth), .groups = "drop_last")
    # Count number muts per sample_group
    mut_counts_groups <- df %>%
        dplyr::group_by(!!ensym(group_var)) %>%
        dplyr::summarise(nrmuts = sum(MutDepth), .groups = "drop_last")

    # Create x facet labels
    facet_labs_x <- stringr::str_c(mut_counts$Subtype, " (n = ", mut_counts$nrmuts, ")")
    names(facet_labs_x) <- mut_counts$Subtype
    # Create y facet labels
    facet_labs_y <- stringr::str_c(mut_counts_groups[[group_var]], " (n = ", mut_counts_groups$nrmuts, ")")
    names(facet_labs_y) <- mut_counts_groups[[group_var]]

    # If user specifies a mutation proportion max, then if value is higher than max,
    # change it to max (i.e., cut off the values at max)
    if (max < 1 && !rescale_data) {
    message(paste0("Cutting off at maximum mutation proportion value of ", max))
    df2 <- df %>%
        dplyr::mutate(ProportionPlot = ifelse(Proportion > max, max, Proportion))
    } else if (max == 1 && rescale_data) {
    # If user specifies scaling to max value (the default), rescale the values to 0-1
        message(paste0("Rescaling bewteen 0 and 1"))
        df2 <- df %>%
            dplyr::mutate(ProportionPlot = scales::rescale(Proportion, to = c(0, 1)))
    } else if (rescale_data && max < 1) {
    # If user specifies both scaling and cutting off at max, then do both
        message(paste0("Rescaling and cutting off at maximum mutation proportion value of ", max))
        df2 <- df %>%
            dplyr::mutate(ProportionPlot = scales::rescale(Proportion, to = c(0, 1))) %>%
            dplyr::mutate(ProportionPlot = ifelse(Proportion > max, max, Proportion))
    } else {
        message(paste0("No scaling or maximum mutation proportion value applied"))
        df2$ProportionPlot <- df$Proportion
    }

    fig <- ggplot(df2, aes(
        x = Context,
        y = SampleName,
        fill = ProportionPlot
    )) +
        geom_raster() +
        scale_fill_viridis_c(
            name = "Relative proportion", limits = c(0, max),
            option = mut_proportion_scale,
            na.value = "white"
        ) +
        facet_grid(Lab ~ Subtype, labeller = labeller(Lab = facet_labs_y, Subtype = facet_labs_x), scales = "free") +
        labs(x = "Trinucleotide context", y = "Lab") +
        theme_minimal() +
        theme(
            axis.text.y = element_text(size = axis_size),
            axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = axis_size, family = "mono"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.spacing.x = unit(spacing, "lines"),
            panel.spacing.y = unit(spacing, "lines")
        )

    return(fig)
}

# Example usage:
# spectra_comparison_heatmap(data = "data/trinucleotide_spectra_data_down_sampled.txt",
#                                       condensed = FALSE,
#                                       mut_proportion_scale = "mako",
#                                       max = 1,
#                                       group_var = "Lab",
#                                       rescale_data = FALSE)