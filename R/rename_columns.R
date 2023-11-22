#' Map column names of mutation data to default column names.

#' A utility function that renames columns of mutation data to default columns names.
#' @param data mutation data
#' @param column_map a list that maps synonymous column names to their default. 
#' @returns the mutation data with column names changed to match default.
#' @export

rename_columns <- function(data,
                           column_map = op$column) {
  # Map column names using the provided column mapping (case-insensitive)
  mapped_columns <- names(data)
  for (col in names(data)) {
    # Convert column name and mapping key to lowercase for case-insensitive comparison
    col_lower <- tolower(col)
    if (col_lower %in% tolower(names(column_map))) {
      mapped_col_name <- column_map[[tolower(col)]]
      if (col_lower != tolower(mapped_col_name)) {
        cat("Expected '", mapped_col_name, "' but found '", col, "', matching columns in input data\n")
      }
      mapped_columns[mapped_columns == col] <- mapped_col_name
    }
  }
  names(data) <- mapped_columns
  
  return(data)}


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


