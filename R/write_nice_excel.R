write_excel_from_list <- function(list_of_tables, output_path, workbook_name) {
  if (!require(openxlsx)) {
    stop("openxlsx not installed")
  }
  stopifnot(is.list(list_of_tables))
  hs1 <- createStyle(textDecoration = "Bold",
                     border = "Bottom",
                     fontColour = "black")
  hs2 <- createStyle(textDecoration = "Bold",
                     border = c("top", "bottom", "left", "right"),
                     fontColour = "black",
                     fgFill = "#C5D9F1")
  options("openxlsx.borderColour" = "#4F80BD")
  options("openxlsx.borderStyle" = "thin")
  options("openxlsx.maxWidth" = 50)
  wb1 <- createWorkbook()
  for (i in seq_along(list_of_tables)) {
    print(i)
    dataToWrite <- as.data.frame(list_of_tables[i])
    addWorksheet(wb1, names(list_of_tables[i]))
    freezePane(wb1, sheet = i, firstRow = TRUE, firstActiveCol = 1)
    writeDataTable(wb1,
                   sheet = i,
                   x = dataToWrite,
                   colNames = TRUE,
                   rowNames = F,
                   tableStyle = "none",
                   headerStyle = hs1,
                   keepNA = T,
                   na.string = "NA")
    setColWidths(wb1, sheet = i, cols = 1:ncol(dataToWrite), widths = "auto")
  }
  fname <- file.path(output_path, paste0(workbook_name,".xlsx"))
  saveWorkbook(wb1, fname, overwrite = TRUE)
}

