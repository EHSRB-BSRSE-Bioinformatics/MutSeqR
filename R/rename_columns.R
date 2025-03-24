#' Map column names of mutation data to default column names.

#' A utility function that renames columns of mutation data to default columns names.
#' @param data mutation data
#' @param column_map a list that maps synonymous column names to their default.
#' @returns the mutation data with column names changed to match default.
#' @export

rename_columns <- function(data, column_map = op$column) {
  # store original names
  original_colnames <- colnames(data)
  # normalized column names for matching
  normalized_colnames  <- tolower(gsub("\\.+", "_",
                                 gsub("(\\.+)?$", "",
                                      gsub("^((X\\.+)|(\\.+))?", "",
                                           original_colnames)),
                                 perl = TRUE))
  # Identify existing default column names in the data
  existing_defaults <- unique(unlist(column_map))  # Extract all default names
  present_defaults <- existing_defaults[existing_defaults %in% normalized_colnames]

  # Create a mapping of normalized names to original names
  norm_to_orig <- setNames(original_colnames, normalized_colnames)
  
  # Initialize new column names (default to original names)
  new_colnames <- original_colnames
  
  # Track which default names have already been assigned
  assigned_defaults <- setNames(rep(FALSE, length(existing_defaults)), existing_defaults)
  
  # Apply renaming rules
  for (synonym in names(column_map)) {
    default_name <- column_map[[synonym]]

    # Only rename if the default name is NOT already in the data
    if (!(default_name %in% present_defaults) && synonym %in% normalized_colnames) {
      original_name <- norm_to_orig[[synonym]]  # Get original column name

      # Only rename the first synonym encountered
      if (!assigned_defaults[[default_name]]) {
        cat("Expected '", default_name, "' but found '", original_name, "', renaming it.\n")
        new_colnames[new_colnames == original_name] <- default_name
        assigned_defaults[[default_name]] <- TRUE  # Mark this default as assigned
      }
    }
  }

  # Ensure all existing default columns are in lowercase format
  for (default_col in present_defaults) {
    original_name <- norm_to_orig[[default_col]]
    if (!is.null(original_name)) {
      new_colnames[new_colnames == original_name] <- default_col
    }
  }

  # Apply final column names
  colnames(data) <- new_colnames
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
