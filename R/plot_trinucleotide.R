#' Plot the trinucleotide spectrum
#' @description Creates barplots of the trinucleotide spectrum for all levels of
#' a given group based on the mutation data. All plots are exported. 
#' @param mutation_data A data frame containing mutation data. This can be obtained
#' using the 'import_mut_data' or 'read_vcf' functions.
#' @param response A character string specifying the type of response to plot.
#' Must be one of 'frequency', 'proportion', or 'sum'.
#' @param mf_type A character string specifying the mutation count method to
#' plot. Must be one of 'min' or 'max'. Default is 'min'.
#' @param group_col A character string specifying the column(s) in 'mutation_data'
#' to group the data by. Default is 'sample'. The sum, proportion, or frequency
#' will be calculated and a plot will be generated for all unique levels of this
#' group. You can specify more than one column to group by.
#' @param max_y A character string specifying the max response value for the y-axis.
#' Must be one of 'individual' or 'group'.'individual' will adjust the maximum y-axis
#' value for each level of the group independently of the others. 'group' will set the
#' maximum y-axis value based on the entire dataset such that all plots will have the
#' same scale. Default is 'group'.
#' @param output_path A character string specifying the path to save the output plot.
#' Default is NULL. This will create an output directory in the current working
#' directory.
#' @param output_type A character string specifying the type of output file.
#' Options are  'jpeg', 'pdf', 'png', 'svg', or 'tiff'. Default is 'svg'.
#' @param sum_totals A logical value specifying whether to sum the total mutations.
#' @importFrom dplyr arrange group_by mutate summarise
#' @importFrom stringr str_extract str_c

#' @details The function calculates the mutation frequency and plots the trinucleotide
#' spectrum for all levels of a given group based on the mutation data.
#' The function calculates the mutation frequency using the 'calculate_mut_freq'
#' function with "cols_to_group" set to 'group_col' and "subtype_resolution" set
#' to 'base_96'. For a given group, mutation counts and total informative duplex
#' bases are summed across all samples. Mutation frequency is calculated by
#' dividing the total mutation counts by the total number of duplex bases.
#' For a given mutation subtype, proportion is calculated as the proportion
#' of total mutation counts normalized to the total number of duplex bases
#' for a given group and subtype.
#' ^ should just explain this in calculate mutation freq and refer to that function.
#' @export

plot_trinucleotide <- function(mutation_data,
                               response = c("frequency", "proportion", "sum"),
                               mf_type = "min",
                               group_col = "dose",
                               max_y = c("individual", "group"),
                               sum_totals = TRUE,
                               output_path = NULL,
                               output_type = "svg") {
  # Output directory
  if (is.null(output_path)) {
    output_dir <- file.path(here::here(), "output")
  } else {
    output_dir <- file.path(output_path)
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # Calculate mutation frequency
  mf_96 <- MutSeqR::calculate_mut_freq(mutation_data = mutation_data,
                                       cols_to_group = group_col,
                                       subtype_resolution = "base_96",
                                       variant_types = "snv")

  group_col_prefix <- paste(group_col, collapse = "_")

  # Desginate the response column
  if (response == "proportion") {
    response_col <- paste0("proportion_", mf_type)
  } else if (response == "frequency") {
    response_col <- paste0("mf_", mf_type)
  } else if (response == "sum") {
    response_col <- paste0("sum_", mf_type)
  } else {
    stop("response must be one of 'frequency', 'proportion', or 'sum'")
  }
  # concatenate multiple group columns into one
  if (length(group_col) > 1) {
    mf_96$group <- apply(mf_96[group_col], 1, paste, collapse = "_")
  } else {
    mf_96$group <- mf_96[[group_col]]
  }
  # Select and rename required columns
  data <- mf_96 %>%
    dplyr::select("group",
                  "normalized_context_with_mutation",
                  "normalized_context",
                  dplyr::all_of(paste0("sum_", mf_type)),
                  dplyr::all_of(response_col)) %>%
    dplyr::rename(context = "normalized_context",
      subtype = "normalized_context_with_mutation",
      response = dplyr::all_of(response_col),
      sum = dplyr::all_of(paste0("sum_", mf_type))
    ) %>%
    dplyr::mutate(mutation = stringr::str_extract(.data$subtype, "(?<=\\[)[^\\]]+(?=\\])")) %>%
    dplyr::arrange(.data$mutation, .data$context) %>%
    dplyr::mutate(subtype = factor(.data$subtype, levels = unique(.data$subtype)))

  if (response == "sum") {
    data$response <- data$sum
  }
  group_levels <- unique(data$group)

  # Determine y_max for all levels of the group
  if (max_y == "group") {
    if (response == "proportion") {
      y_max <- ceiling(max(data$response) * 10) / 10
      y_lab <- "Proportion of Mutations"
    } else if (response == "frequency") {
      max_frequency <- max(data$response)
      max_frequency_string <- format(max_frequency, scientific = TRUE)
      split_frequency <- strsplit(max_frequency_string, "e", fixed = TRUE)[[1]]
      coefficient <- as.numeric(split_frequency[1])
      rounded_coefficient <- ceiling(coefficient * 10) / 10
      round_str <- paste(rounded_coefficient, "e", split_frequency[2], sep = "")
      y_max <- as.numeric(round_str)
      y_lab <- "Frequency of Mutations"
    } else if (response == "sum") {
      y_max <- ceiling(max(data$response) / 5) * 5
      y_lab <- "Sum of Mutations"
    }
  }
  # Loop over all levels in group
  for (i in seq_along(group_levels)) {
    plot_data <- data[data$group == group_levels[i], ]

    # Facet labels
    labels <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    # Not implemented yet
    if (sum_totals) {
      mut_counts <- plot_data %>%
        dplyr::group_by(.data$mutation) %>%
        dplyr::summarise(nrmuts = sum(sum), .groups = "drop_last")
      labels <- stringr::str_c(mut_counts$mutation, " (n = ", mut_counts$nrmuts, ")")
      names(labels) <- mut_counts$mutation
    }

    plotcolours <- c(grDevices::rgb(5, 195, 239, maxColorValue = 255),
      grDevices::rgb(0, 0, 0, maxColorValue = 255),
      grDevices::rgb(230, 47, 41, maxColorValue = 255),
      grDevices::rgb(208, 207, 207, maxColorValue = 255),
      grDevices::rgb(169, 212, 108, maxColorValue = 255),
      grDevices::rgb(238, 205, 204, maxColorValue = 255)
    )

    rearr.colours <- c(rep(plotcolours[1], 16),
      rep(plotcolours[2], 16),
      rep(plotcolours[3], 16),
      rep(plotcolours[4], 16),
      rep(plotcolours[5], 16),
      rep(plotcolours[6], 16)
    )

    xlabels <- plot_data$subtype
    # Define your title
    title <- group_levels[i]
    # Define your axis text size
    cex.axistext <- 0.4

    # Determine y_max for each level of the group seperately.
    if (max_y == "individual") {
      if (response == "proportion") {
        y_max <- ceiling(max(plot_data$response) * 10) / 10
        y_lab <- "Proportion of Mutations"
      } else if (response == "frequency") {
        max_frequency <- max(plot_data$response)
        max_frequency_string <- format(max_frequency, scientific = TRUE)
        split_frequency <- strsplit(max_frequency_string, "e", fixed = TRUE)[[1]]
        coefficient <- as.numeric(split_frequency[1])
        rounded_coefficient <- ceiling(coefficient * 10) / 10
        round_str <- paste(rounded_coefficient, "e", split_frequency[2], sep = "")
        y_max <- as.numeric(round_str)
        y_lab <- "Frequency of Mutations"
      } else if (response == "sum") {
        y_max <- ceiling(max(plot_data$response) / 5) * 5
        y_lab <- "Sum of Mutations"
      }
    }

    # Define ylim and the y-axis_labels baseed on y_max
    ylim <- c(0, y_max)
    y_axis_labels <- seq(0, y_max, by = y_max / 5)

    # Define the filename
    filename <- paste0(output_dir, "/plot_", group_levels[i], "_trinucleotide_", response)
    if (output_type == "jpeg") {
      filename <- paste0(filename, ".jpeg")
      grDevices::jpeg(filename, width = 7, height = 3.5, units = "in", res = 300)
    } else if (output_type == "pdf") {
      filename <- paste0(filename, ".pdf")
      grDevices::pdf(filename, width = 7, height = 3.5)
    } else if (output_type == "png") {
      filename <- paste0(filename, ".png")
      grDevices::png(filename, width = 7, height = 3.5, units = "in", res = 300)
    } else if (output_type == "tiff") {
      filename <- paste0(filename, ".tiff")
      grDevices::tiff(filename, width = 7, height = 3.5, units = "in", res = 300)
    } else if (output_type == "svg") {
      filename <- paste0(filename, ".svg")
      grDevices::svg(filename, width = 7, height = 3.5)
    }

    # Create the plot
    bp <- graphics::barplot(plot_data$response, # height of bars
      main = title,
      names.arg = NULL,
      xlab = "Mutation Type",
      ylab = y_lab,
      ylim = ylim,
      axes = FALSE, # don't draw axes yet
      col = NA, # no color for bars yet
      beside = TRUE,
      las = 2, # make labels perpendicular to axis
      border = NA,
      space = 2,
      cex.main = cex.axistext * 2
      )

    # Add horizontal gridlines
    xlim <- graphics::par("usr")[1:2]
    gridlines <- y_axis_labels
    for (j in gridlines) {
      graphics::lines(x = xlim, y = c(j, j), col = "#f3eeeea4", lty = "solid")
    }

    # Now draw the bars over the gridlines
    bp <- graphics::barplot(plot_data$response,
                  col = rearr.colours,
                  beside = TRUE,
                  border = NA,
                  space = 2,
                  add = TRUE,
                  axes = FALSE)

    # Add x-labels manually so they are closer to the bars
    graphics::text(x = bp,
         y = par("usr")[3] * 0.90,
         labels = xlabels,
         srt = 90, adj = 1,
         xpd = TRUE,
         cex = cex.axistext,
         family = "mono")
    # Add x-axis line
    graphics::lines(x = c(min(bp[,1]) - 1.2, max(bp[length(bp), ]) + 0.3),
          y = c(par("usr")[3], par("usr")[3]),
          lwd = 1,
          col = "gray")

    # Draw custom y-axis
    graphics::axis(side = 2,
         at = y_axis_labels,
         col = "gray",
         col.axis = "black",
         pos = min(bp[,1]) - 1.2,
         cex.axis = cex.axistext,
         las = 1)

    # Add subtype labels and coloured rectangles
    graphics::par(xpd = TRUE) # allow plotting outside the plot region
    usr <- graphics::par("usr") # plot region dimensions
    rect_top_relative <- 1.05
    rect_top <- usr[4] * rect_top_relative
    text_y <- rect_top * 1.02
    graphics::rect(xleft = bp[1], ybottom = usr[4], xright = bp[16], ytop = rect_top, col = plotcolours[1], border = NA)
    graphics::rect(xleft = bp[17], ybottom = usr[4], xright = bp[32], ytop = rect_top, col = plotcolours[2], border = NA)
    graphics::rect(xleft = bp[33], ybottom = usr[4], xright = bp[48], ytop = rect_top, col = plotcolours[3], border = NA)
    graphics::rect(xleft = bp[49], ybottom = usr[4], xright = bp[64], ytop = rect_top, col = plotcolours[4], border = NA)
    graphics::rect(xleft = bp[65], ybottom = usr[4], xright = bp[80], ytop = rect_top, col = plotcolours[5], border = NA)
    graphics::rect(xleft = bp[81], ybottom = usr[4], xright = bp[96], ytop = rect_top, col = plotcolours[6], border = NA)

    # Add text to the rectangles
    graphics::text(bp[8], text_y, labels[1], cex = cex.axistext*1.2, col = "black", font = 2)
    graphics::text(bp[24], text_y, labels[2], cex = cex.axistext*1.2, col = "black", font = 2)
    graphics::text(bp[40], text_y, labels[3], cex = cex.axistext*1.2, col = "black", font = 2)
    graphics::text(bp[56], text_y, labels[4], cex = cex.axistext*1.2, col = "black", font = 2)
    graphics::text(bp[72], text_y, labels[5], cex = cex.axistext*1.2, col = "black", font = 2)
    graphics::text(bp[88], text_y, labels[6], cex = cex.axistext*1.2, col = "black", font = 2)

    grDevices::dev.off()
  }
}