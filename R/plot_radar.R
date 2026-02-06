#' Create a radar plot
#' @description Create a radar plot
#' @param mf_data A data frame with the data to plot
#' @param response_col The column with the response values
#' @param label_col The column with the labels for the radar plot.
#' @param facet_col The column with the group to facet the radar plots.
#' @param indiv_y A logical indicating whether to use individual y-axis scales
#' for each plot.
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr select filter pull
#' @examples
#' # Plot the MFmin for each variation type, including the 6 SNV subtypes
#' # Facet the plots by dose group
#' mf_ex <- readRDS(system.file("extdata", "Example_files", "mf_data_6.rds",
#'      package = "MutSeqR")
#' )
#' plot <- plot_radar(
#'      mf_data = mf_ex,
#'      response_col = "mf_min",
#'      label_col = "normalized_subtype",
#'      facet_col = "dose_group",
#'      indiv_y = TRUE
#' )
#' @return A radar plot
#' @export
plot_radar <- function(mf_data,
                       response_col,
                       label_col,
                       facet_col,
                       indiv_y = TRUE) {
 # --- 1. Validation ---
  stopifnot(
    !missing(mf_data) && is.data.frame(mf_data),
    !missing(response_col), !missing(label_col)
  )
  if (!requireNamespace("fmsb", quietly = TRUE)) stop("Package fmsb required.")

  # Validate columns exist
  cols_to_check <- c(response_col, label_col)
  # Only add facet_col to check if it's not NULL
  if (!is.null(facet_col)) cols_to_check <- c(cols_to_check, facet_col)

  if (!all(cols_to_check %in% colnames(mf_data))) {
     stop("Columns not found: ", paste(setdiff(cols_to_check, colnames(mf_data)), collapse=", "))
  }

  # --- 2. Data Preparation ---
  # Handle the optional facet column
  if (is.null(facet_col)) {
    # If no facet provided, create a dummy column to treat the whole dataset as one group
    dummy_facet_name <- "All_Data"
    mf_data[[dummy_facet_name]] <- "All Data"
    facet_col <- dummy_facet_name # Update variable to use the dummy
    # Turn off individual y scaling since there is only 1 plot
    indiv_y <- FALSE 
  }

  # Ensure consistent ordering if label_col is a factor
  if (is.factor(mf_data[[label_col]])) {
    label_levels <- levels(mf_data[[label_col]])
  } else {
    label_levels <- sort(unique(mf_data[[label_col]]))
  }

  # Pivot to Wide Format (Rows=Facets, Cols=Labels)
  plot_data <- mf_data %>%
    dplyr::select(dplyr::all_of(c(facet_col, label_col, response_col))) %>%
    tidyr::pivot_wider(names_from = dplyr::all_of(label_col), 
                       values_from = dplyr::all_of(response_col),
                       values_fill = 0)

  # Reorder columns to match label levels
  # Check intersection to avoid errors if some levels are missing in data
  valid_levels <- intersect(label_levels, colnames(plot_data))
  plot_data <- plot_data[, c(facet_col, valid_levels)]

  # Ensure facet is a factor
  if (!is.factor(plot_data[[facet_col]])) {
      plot_data[[facet_col]] <- factor(plot_data[[facet_col]])
  }
  facet_levels <- levels(plot_data[[facet_col]])

  # --- 3. Scale Calculation ---
  # Calculate max values (excluding the facet column)
  data_matrix <- as.matrix(plot_data[, -1])

  if (!indiv_y) {
    global_max <- max(data_matrix, na.rm = TRUE) * 1.1
  } else {
    row_maxes <- apply(data_matrix, 1, max, na.rm = TRUE) * 1.1
    global_max <- NULL
  }

  # --- 4. Plot Layout ---
  n_plots <- length(facet_levels)

  # Only set up grid layout if we have multiple plots
  if (n_plots > 1) {
    n_cols <- 2
    n_rows <- ceiling(n_plots / n_cols)
    on.exit(graphics::layout(1)) # Reset on exit
    graphics::layout(matrix(seq_len(n_rows * n_cols), nrow = n_rows, ncol = n_cols, byrow = TRUE))
  }

  # --- 5. Plotting Loop ---
  lapply(seq_along(facet_levels), function(i) {

    facet_name <- facet_levels[i]
    current_vals <- data_matrix[i, ]

    # Determine Max Scale
    current_max <- if (indiv_y) row_maxes[i] else global_max

    # Prepare Dataframe for fmsb
    df_radar <- rbind(
      rep(current_max, length(current_vals)), # Max
      rep(0, length(current_vals)),           # Min
      current_vals                            # Data
    )
    df_radar <- as.data.frame(df_radar)
    colnames(df_radar) <- colnames(data_matrix)

    axis_labels <- sprintf("%.1e", seq(0, current_max, length.out = 5))

    # Determine Title
    # If dummy facet, use generic title, otherwise use Facet Name
    if (facet_col == "All_Data") {
      plot_title <- "Radar Plot"
    } else {
      plot_title <- paste(facet_col, facet_name)
    }

    fmsb::radarchart(
      df = df_radar,
      axistype = 1,
      caxislabels = axis_labels,
      seg = 4,
     # caxiscol = "grey",
      vlabels = names(df_radar),
      axislabcol = "#797777",
      vlcex = 1.2,
      title = plot_title,
      pcol = "black",
      plwd = 2,
      plty = 1,
      cglcol = "grey",
      cglty = 1,
      cglwd = 0.7
    )
  })
  
  invisible(NULL)
}