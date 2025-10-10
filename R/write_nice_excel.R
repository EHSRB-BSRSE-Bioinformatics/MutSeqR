#' Write Excel tables
#'
#' @description Writes data to an Excel file.
#' @param data A data frame or a list of data frames. If a data frame,
#' it will be written to a single sheet in the Excel workbook. If a list, each
#' data frame will be written to a separate sheet in the Excel workbook.
#' Data may also be the output to model_mf, in which case set
#' `model_results = TRUE`.
#' @param output_path The directory where the Excel file should be written.
#' Default is NULL, which will write the file to the current working directory.
#' @param workbook_name The file name for the Excel file.
#' @param model_results A logical value indicating whether the data is the
#' output of model_mf. Default is FALSE. If TRUE, the function will grab the
#' model_data, point_estimates, and pairwise_comparisons data frames from the
#' model_mf output and write them to separate sheets in the Excel workbook.
#' @returns A saved Excel workbook.
#' @examples
#' # Example data consists of 24 mouse bone marrow DNA samples imported
#' # using import_mut_data(), filtered with filter_mut, and summarized
#' # using calculate_mf().
#' outputpath <- tempdir()
#' mf_example <- readRDS(system.file("extdata/Example_files/mf_data_global.rds",
#'                       package = "MutSeqR"))
#' mf_example2 <- readRDS(system.file("extdata/Example_files/mf_data_6.rds",
#'                        package = "MutSeqR"))
#' mf_example3 <- readRDS(system.file("extdata/Example_files/mf_data_6_sample.rds",
#' list <- list(mf_example, mf_example2, mf_example3)
#' names(list) <- c("Global MF", "Base 6 Spectra", "Base 6 Sample Spectra")
#'
#' # save a single data frame to an Excel file
#' write_excel(mf_example, output_path = outputpath, workbook_name = "test_single")
#'
#' # save a list of data frames to an Excel file
#' write_excel(list, output_path = outputpath, workbook_name = "test_list")
#' @export
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
  if (model_results) {
    data$point_estimates$rows <- rownames(data$point_estimates)
    if ("pairwise_comparisons" %in% names(data)) {
      data$pairwise_comparisons$rows <- rownames(data$pairwise_comparisons)
    }
    model_dfs <- list(
      model_data = data$model_data,
      point_estimates = data$point_estimates
    )
    if ("pairwise_comparisons" %in% names(data)) {
      model_dfs$pairwise_comparisons <- data$pairwise_comparisons
    }
    data <- model_dfs
  }
  if (is.data.frame(data)) {
    hs1 <- openxlsx::createStyle(
      textDecoration = "Bold",
      border = "Bottom",
      fontColour = "black"
    )
    options("openxlsx.borderColour" = "#4F80BD")
    options("openxlsx.borderStyle" = "thin")
    options("openxlsx.maxWidth" = 50)
    wb1 <- openxlsx::createWorkbook()
    dataToWrite <- as.data.frame(data)
    openxlsx::addWorksheet(wb1, workbook_name)
    openxlsx::freezePane(wb1, sheet = 1, firstRow = TRUE, firstActiveCol = 1)
    openxlsx::writeDataTable(wb1,
      sheet = 1,
      x = dataToWrite,
      colNames = TRUE,
      rowNames = FALSE,
      tableStyle = "none",
      headerStyle = hs1,
      keepNA = TRUE,
      na.string = "NA"
    )
    openxlsx::setColWidths(wb1, sheet = 1, cols = seq_len(ncol(dataToWrite)), widths = "auto")
    fname <- file.path(output_dir, paste0(workbook_name, ".xlsx"))
    openxlsx::saveWorkbook(wb1, fname, overwrite = TRUE)
  } else if (is.list(data)) {
    if (is.null(names(data)) || NA %in% names(data)) {
      names(data) <- seq_along(data)
    }
    hs1 <- openxlsx::createStyle(
      textDecoration = "Bold",
      border = "Bottom",
      fontColour = "black"
    )
    options("openxlsx.borderColour" = "#4F80BD")
    options("openxlsx.borderStyle" = "thin")
    options("openxlsx.maxWidth" = 50)
    wb1 <- openxlsx::createWorkbook()
    for (i in seq_along(data)) {
      dataToWrite <- as.data.frame(data[i])
      openxlsx::addWorksheet(wb1, names(data[i]))
      openxlsx::freezePane(wb1, sheet = i, firstRow = TRUE, firstActiveCol = 1)
      openxlsx::writeDataTable(wb1,
                               sheet = i,
                               x = dataToWrite,
                               colNames = TRUE,
                               rowNames = FALSE,
                               tableStyle = "none",
                               headerStyle = hs1,
                               keepNA = TRUE,
                               na.string = "NA")
      openxlsx::setColWidths(wb1, sheet = i, cols = seq_len(ncol(dataToWrite)),
                             widths = "auto")
    }
    fname <- file.path(output_dir, paste0(workbook_name, ".xlsx"))
    openxlsx::saveWorkbook(wb1, fname, overwrite = TRUE)
  } else {
    stop("Data must be a list or data frame")
  }
}
