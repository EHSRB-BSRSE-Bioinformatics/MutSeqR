#' Write results to Excel tables
#'
#' @description Writes data to an Excel file.
#' @param data A data frame, a list of data frames, or model_mf output.
#' @param output_path Directory to write to. Defaults to current working
#' directory.
#' @param workbook_name Filename (without extension).
#' @param model_results Logical. Set to TRUE if data is output from model_mf.
#' @export
#' @returns A saved Excel workbook.
#' @examples
#' # Example data consists of 24 mouse bone marrow DNA samples imported
#' # using import_mut_data(), filtered with filter_mut, and summarized
#' # using calculate_mf().
#' outputpath <- tempdir()
#' mf_example <- readRDS(system.file("extdata/Example_files/mf_data_global.rds",
#'   package = "MutSeqR"
#' ))
#' mf_example2 <- readRDS(system.file("extdata/Example_files/mf_data_6.rds",
#'   package = "MutSeqR"
#' ))
#' mf_example3 <- readRDS(
#'  system.file("extdata/Example_files/mf_data_6_sample.rds",
#'      package = "MutSeqR"
#'  )
#' )
#' list <- list(mf_example, mf_example2, mf_example3)
#' names(list) <- c("Global MF", "Base 6 Spectra", "Base 6 Sample Spectra")
#'
#' # save a single data frame to an Excel file
#' write_excel(
#'  mf_example,
#'  output_path = outputpath,
#'  workbook_name = "test_single"
#' )
#'
#' # save a list of data frames to an Excel file
#' write_excel(list, output_path = outputpath, workbook_name = "test_list")
write_excel <- function(data,
                        output_path = NULL,
                        workbook_name,
                        model_results = FALSE) {

    if (!requireNamespace("openxlsx", quietly = TRUE)) {
        stop("openxlsx not installed")
    }
  if (is.null(output_path)) {
    output_dir <- file.path(here::here())
  } else {
    output_dir <- file.path(output_path)
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # Data Normalization
  data_list <- list()

  if (model_results) {
    # Handle Row Names for point estimates
    pe <- data$point_estimates
    pe$rows <- rownames(pe)
    data_list[["point_estimates"]] <- pe
    data_list[["model_data"]] <- data$model_data
    if ("pairwise_comparisons" %in% names(data)) {
      pc <- data$pairwise_comparisons
      pc$rows <- rownames(pc)
      data_list[["pairwise_comparisons"]] <- pc
    }

  } else if (is.data.frame(data)) {
    # Wrap single DF in a list
    data_list[[workbook_name]] <- data

  } else if (is.list(data)) {
    # Use existing list
    data_list <- data
    # Ensure names exist
    if (is.null(names(data_list))) {
      names(data_list) <- paste0("Sheet", seq_along(data_list))
    }
  } else {
    stop("Data must be a data frame, a list, or a model_mf object.")
  }

  # Style Setup
  hs1 <- openxlsx::createStyle(
    textDecoration = "Bold",
    border = "Bottom",
    fontColour = "black"
  )

  # Set options, rest on exit
  old_opts <- options(
    "openxlsx.borderColour" = "#4F80BD",
    "openxlsx.borderStyle" = "thin",
    "openxlsx.maxWidth" = 50
  )
  on.exit(options(old_opts), add = TRUE)

  # Workbook Writing
  wb <- openxlsx::createWorkbook()

  # ff(x) to write one sheet
  write_sheet <- function(sheet_name, df) {
    valid_name <- substr(sheet_name, 1, 31) # sheet names must not exceed 31

    openxlsx::addWorksheet(wb, valid_name)
    openxlsx::freezePane(
        wb, sheet = valid_name,
        firstRow = TRUE, firstActiveCol = 1
    )

    openxlsx::writeDataTable(
      wb,
      sheet = valid_name,
      x = as.data.frame(df), # Ensure it's a DF (handles tibbles etc)
      colNames = TRUE,
      rowNames = FALSE,
      tableStyle = "none",
      headerStyle = hs1,
      keepNA = TRUE,
      na.string = "NA"
    )

    openxlsx::setColWidths(
        wb, sheet = valid_name,
        cols = 1:ncol(df), widths = "auto"
    )
  }

  # Map the writing function over the data list and its names.
  invisible(Map(write_sheet, names(data_list), data_list))

  # Save
  fname <- file.path(output_path, paste0(workbook_name, ".xlsx"))
  openxlsx::saveWorkbook(wb, fname, overwrite = TRUE)

  message("Saved: ", fname)
}