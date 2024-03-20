#' Create a heatmap plot
#'
#' This function creates a heatmap plot using the provided data file.
#'
#' @param mf_data The data file
#' @param group_var The variable to group by.
#' @param mf_type The type of mutation frequency to plot. Options are "unique" or "clonal". (Default: "unique")
#' @param mut_proportion_scale The scale option for the mutation proportion. Options are passed to viridis::scale_fill_viridis_c.
#'  One of # inferno, magma, plasma, viridis, cividis, turbo, mako, or rocket. We highly reccomend the default for its ability to
#' disciminate hard to see patterns. (Default: "turbo")
#' @param max Maximum value used for plotting the relative contributions.
#' Contributions that are higher will have the maximum colour. (Default: 0.2)
#' @param rescale_data Logical value indicating whether to rescale the mutation proportions to increase the dynamic range of colors shown on the plot. (Default: TRUE)
#' @param condensed More condensed plotting format. Default = F.
#' @import ggplot2
#' @importFrom dplyr %>% group_by summarise mutate
#' @importFrom scales rescale
#' @return A ggplot object representing the heatmap plot.
#' @export
#' @examples
#' create_heatmap(mutation_data, dose, "inferno")
#' 
# Some options for data format:
    # mf_data: single column with mutation subtype, context, and proportion. May have variant types mixed in. 
    # A data set with context and subtype seperated into different columns, simple spectrum  
spectra_comparison_heatmap <- function(mf_data,
                                       group_var = "dose",
                                       mf_type = "unique",
                                       mut_proportion_scale = "turbo",
                                       max = 0.2,
                                       rescale_data = FALSE,
                                       condensed = FALSE) {

mf_data <- MutSeqR::rename_columns(mf_data)
# Check for sample colum
sample_col <- c("sample", "sample_id", "sample_name")
if (!"sample" %in% colnames(mf_data)) {
  stop("The dataframe does not contain a 'sample' column.")
}

# Check for proportion column
proportion_col <- c(paste0("proportion_", mf_type), "proportion", "prop")
found_prop_col <- proportion_col[proportion_col %in% tolower(colnames(mf_data))] 
if(length(found_prop_col) == 1) {
   mf_data <- dplyr::rename(mf_data, proportion = dplyr::all_of(found_prop_col))
} else if (length(found_prop_col) > 1) {
  stop(paste("More than one possible proportion column name found in mf_data: ",
             paste(found_prop_col, collapse = ", "),
             " Please remove columns that are not the proportion column to be plotted."))
} else if (length(found_prop_col) == 0) {
 stop(paste0("The dataframe does not contain a proportion column. Please add a column that contains the mutation proportion to be plotted or rename column to 'proportion'."))
}

# Check for mutation count column
count_cols <- c("mut_count", "alt_depth", "mut_depth")
sum_column_pattern <- paste0(".*(_sum_", mf_type, ")$")
count_col1 <- any(count_cols %in% tolower(colnames(mf_data)))
count_col2 <- any(grepl(sum_column_pattern, names(mf_data)))
if (count_col1 & !count_col2) {
  found_count_col <- count_cols[count_cols %in% tolower(colnames(mf_data))]
  mf_data <- dplyr::rename(mf_data, sum_column = dplyr::all_of(found_count_col))
} else if (count_col2 & !count_col1) {
  found_count_col <- names(mf_data)[grepl(sum_column_pattern, names(mf_data))]
  mf_data <- dplyr::rename(mf_data, sum_column = dplyr::all_of(found_count_col))
} else if (count_col1 && count_col2) {
  warning(paste0("Multiple count columns were found in the data. Using the column ending in '_sum_", mf_type, " for group counts."))
  found_count_col <- names(mf_data)[grepl(sum_column_pattern, names(mf_data))]
  mf_data <- dplyr::rename(mf_data, sum_column = dplyr::all_of(found_count_col))
} else if (!count_col1 && !count_col2) {
  warning("No count column was found in the data. The plot will not display group counts.")
}
# Check for subtype column
subtype_column_names <- c("normalized_context_with_mutation",
                          "context_with_mutation",
                          "normalized_subtype",
                          "subtype",
                          "variation_type",
                          "mutation_subtype")
found_subtype_cols <- subtype_column_names[subtype_column_names %in% colnames(mf_data)] 
if (length(found_subtype_cols) == 1) {
  # Rename the subtype column to "subtype"
  mf_data <- dplyr::rename(mf_data, subtype = dplyr::all_of(found_subtype_cols))
} else if (length(found_subtype_cols) > 1) {
 stop(paste("More than one possible subtype column name found in mf_data: ",
             paste(found_subtype_cols, collapse = ", "),
             " Please remove columns that are not the subtype column to be plotted.")) 
  } else if (length(found_subtype_cols) == 0) {
    stop(paste("No subtype column name found in mf_data from the following options: ", 
               paste(subtype_column_names, collapse = ", "),
               " Please add a column that contains the mutation subtype to be plotted or rename column to one of the listed options.")) 
  }

  # Context Column synonyms
  context_column_names <- c("normalized_context",
                            "context",
                            "mutation_context",
                            "trinucleotide_context",
                            "flanking_sequence")
found_context_cols <- context_column_names[context_column_names %in% colnames(mf_data)] 
if (length(found_context_cols) == 1) {
    mf_data$x_variable <- mf_data[[found_context_cols]]
    mf_data <- mf_data %>% 
      dplyr::mutate(x_variable = ifelse(subtype %in% MutSeqR::subtype_list$type, 
                                        subtype, 
                                        x_variable))
    plot_context = TRUE
} else if (length(found_context_cols) == 0) { # Create context column from mutation column if possible.
    if (any(grep(".*\\[([A-Z]?>[A-Z])\\].*", mf_data$subtype))) {
     mf_data <- mf_data %>% 
      dplyr::mutate(context = paste0(substr(mf_data$subtype, 1, 1),
                                     substr(mf_data$subtype, 3, 3),
                                     substr(mf_data$subtype, 7, 7)))
         mf_data$x_variable <- mf_data$context
         mf_data <- mf_data %>% 
          dplyr::mutate(x_variable = ifelse(subtype %in% MutSeqR::subtype_list$type, 
                                            subtype, 
                                            x_variable))
         plot_context = TRUE
        } else {
            print("No context column found in mf_data, plotting by subtype.") 
             mf_data$x_variable <- mf_data$subtype 
             plot_context = FALSE  
        }
  } else if (length(found_context_cols) > 1) {
 stop(paste("More than one possible context column name found in mf_data: ",
             paste(found_context_cols, collapse = ", "),
             " Please remove columns that are not the context column to be plotted.")) 
  }
  
  context_size <- max(stringr::str_length(mf_data$x_variable))
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

# Facet x
if(plot_context) {
 # Make facet labels for subtypes. 
  pattern <- "[A-Z]?>[A-Z]"
  mf_data$subtype_labels <- stringr::str_extract(mf_data$subtype, pattern)
  mf_data <- mf_data %>%
      mutate(subtype_labels = ifelse(mf_data$subtype %in% subtype_list$type, 
                                "other", 
                                subtype_labels))
    # Count number muts per subtype
  if ("sum_column" %in% colnames(mf_data)) {
    mut_counts <- mf_data %>%
          dplyr::group_by(subtype_labels) %>%
          dplyr::summarise(nrmuts = sum(sum_column), .groups = "drop_last") 
      facet_labs_x <- stringr::str_c(mut_counts$subtype_labels, " (n = ", mut_counts$nrmuts, ")")
      names(facet_labs_x) <- mut_counts$subtype_labels  
  } else {
    facet_labs_x <-unique(mf_data$subtype_labels)
    names(facet_labs_x) <- unique(mut_counts$subtype_labels)
  }
  facet_x_order <- c(MutSeqR::subtype_list$base_12, "other")
  mf_data <- mf_data %>%
      mutate(subtype_labels = factor(subtype_labels, levels = facet_x_order))
  x_label <- "context"                           
} else {
x_label <- "Subtype"
x_order <- c(MutSeqR::subtype_list$base_192, MutSeqR::subtype_list$base_12, MutSeqR::subtype_list$type)
mf_data <- mf_data %>%
  mutate(x_variable = factor(x_variable, levels = x_order))
}
 
# Facet y
if (!is.null(group_var)) {   
 if ("sum_column" %in% colnames(mf_data)) {
  # Count number muts per sample_group
  mut_counts_groups <- mf_data %>%
    dplyr::group_by(!!ensym(group_var)) %>%
    dplyr::summarise(nrmuts = sum(sum_column), .groups = "drop_last") 

  facet_labs_y <- stringr::str_c(mut_counts_groups[[group_var]], " (n = ", mut_counts_groups$nrmuts, ")")
  names(facet_labs_y) <- mut_counts_groups[[group_var]]
 } else {
  facet_labs_y <- unique(mf_data[[group_var]])
  names(facet_labs_y) <- unique(mf_data[[group_var]])
 }
  y_label <- paste(group_var)
  mf_data <- mf_data %>%
    rename(Group = !!ensym(group_var)) 
  mf_data$Group <- as.factor(mf_data$Group)
} else {
  y_label <- "Sample"
}

    # If user specifies a mutation proportion max, then if value is higher than max,
    # change it to max (i.e., cut off the values at max)
    if (max < 1 && !rescale_data) {
    message(paste0("Cutting off at maximum mutation proportion value of ", max))
    df <- mf_data %>%
        dplyr::mutate(ProportionPlot = ifelse(proportion > max, max, proportion))
    } else if (max == 1 && rescale_data) {
    # If user specifies scaling to max value (the default), rescale the values to 0-1
        message(paste0("Rescaling bewteen 0 and 1"))
        df <- mf_data %>%
            dplyr::mutate(ProportionPlot = scales::rescale(proportion, to = c(0, 1)))
    } else if (rescale_data && max < 1) {
    # If user specifies both scaling and cutting off at max, then do both
        message(paste0("Rescaling and cutting off at maximum mutation proportion value of ", max))
        df <- mf_data %>%
            dplyr::mutate(ProportionPlot = scales::rescale(proportion, to = c(0, 1))) %>%
            dplyr::mutate(ProportionPlot = ifelse(proportion > max, max, proportion))
    } else {
        message(paste0("No scaling or maximum mutation proportion value applied"))
        df <- mf_data
        df$ProportionPlot <- mf_data$proportion
    }

 # General figure, no facetting   
fig <- ggplot(df, aes(
            x = x_variable,
            y = sample,
            fill = ProportionPlot)) +
            geom_raster() +
            scale_fill_viridis_c(
                name = "Relative proportion", limits = c(0, max),
                option = mut_proportion_scale,
                na.value = "white") +
            theme_minimal() +
            labs(x = x_label, y = y_label) +
            theme(
                axis.text.y = element_text(size = axis_size),
                axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = axis_size, family = "mono"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.spacing.x = unit(spacing, "lines"),
                panel.spacing.y = unit(spacing, "lines")
            )
# facet x : plot_context == TRUE
# facet y : !is.null(group_var)
# facet x and y : plot_context == TRUE & !is.null(group_var)
# facet none : plot_context == FALSE & is.null(group_var)
if(plot_context & is.null(group_var)) {
  figfx <- fig +
    facet_grid(cols = vars(subtype_labels), scales = "free_x",
               labeller = labeller(subtype_labels = facet_labs_x))
  return(figfx)
} else if (plot_context == FALSE & !is.null(group_var)) {
  figfy <- fig +
    facet_grid(rows = vars(Group), scales = "free_y",
               labeller = labeller(Group = facet_labs_y))
  return(figfy)
} else if (plot_context & !is.null(group_var)) {
    figfxy <- fig +
      facet_grid(Group ~ subtype_labels, scales = "free",
      labeller = labeller(Group = facet_labs_y, subtype_labels = facet_labs_x))
    return(figfxy)
} else {
    return(fig)
}

}
