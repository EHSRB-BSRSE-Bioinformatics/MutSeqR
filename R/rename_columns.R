#' Map column names of mutation data to default column names.

#' A utility function that renames columns of mutation data to default columns names.
#' @param data mutation data
#' @param column_map a list that maps synonymous column names to their default. 
#' @returns the mutation data with column names changed to match default.
#' @export

rename_columns <- function(data,
                           column_map = op$column) {
  # Remove dots and replace with underscores
  colnames(data) <- tolower(gsub("\\.+", "_", #deals with middle periods
                                gsub("(\\.+)?$", "", #deals with trailing periods
                                     gsub("^((X\\.+)|(\\.+))?", "", #deals with beginning X. and periods
                                          colnames(data))),
                                perl = TRUE))

  # Map column names using the provided column mapping (case-insensitive)
  for (col in colnames(data)) {
    if (col %in% names(column_map)) {
      default_name <- column_map[[col]]
      
      if (col != default_name) {
        cat("Expected '", default_name, "' but found '", col, "', matching columns in input data\n")
        colnames(data)[colnames(data) == col] <- default_name
      }
    }
  }
  
  return(data)
}

#' Check that all required columns are present before proceeding with the function
#' 
#' A utility function that will check that all required columns are present. 
#' @param data mutation data
#' @param required_columns a list of required column names.
#' @returns an error
#' @export

check_required_columns <- function(data,
                           required_columns) {
missing_columns <- setdiff(tolower(required_columns), tolower(names(data)))

if (length(missing_columns) > 0) {
  missing_col_names <- paste(missing_columns, collapse = ", ")
  stop(paste("Some required columns are missing or their synonyms are not found: ", missing_col_names))
 } else {
    return(data)
 }
}


