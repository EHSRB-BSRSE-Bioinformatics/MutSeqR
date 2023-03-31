migrate_mut <- function(mut_table, op = DupSeqR::op) {
  mut_table <- data.table::as.data.table(mut_table)
  
  for (required_column in names(op)) {
    if (!op[[required_column]] %in% colnames(mut_table)) {
      matching_column_indices <- which(colnames(mut_table) %in% op[[required_column]])
      stopifnot(
        length(matching_column_indices) == 1,
        paste0(
          "Expected column '", op[[required_column]], "', but found ",
          length(matching_column_indices), " matching columns in input data."
        )
      )
      colnames(mut_table)[matching_column_indices] <- op[[required_column]]
    }
  }
  
  return(mut_table)
}

# migrate_mut <- function(mut_table, processed = FALSE) {
#   required_cols <- if (processed) 
#     c(op$base_required_mut_cols, processed_required_mut_cols)
#   else 
#     op$base_required_mut_cols
#   
#   mut_table <- as.data.table(mut_table)
#   
#   for (col in names(required_cols)) {
#     if (!col %in% colnames(mut_table)) {
#       alt_col <- required_cols[[col]]
#       stopifnot(length(grep(alt_col, colnames(mut_table))) == 1)
#       colnames(mut_table)[grep(alt_col, colnames(mut_table))] <- col
#     }
#   }
#   return(mut_table)
# }



# TS original code
# migrate_mut <- function (mut_table, processed = FALSE)
# {
#   required_cols <- if (processed)
#     c(op$base_required_mut_cols, processed_required_mut_cols)
#   else
#     op$base_required_mut_cols
#   mut_table <- data.table::as.data.table(mut_table)
#   for (required_column in names(required_cols)) {
#     if (!op[[required_column]] %in% colnames(mut_table)) {
#       matching_column_indices <- which(colnames(mut_table) %in%
#                                          required_cols[[required_column]])
#       assertthat::assert_that(
#         length(matching_column_indices) ==
#           1,
#         msg = paste0(
#           "Found ",
#           length(matching_column_indices),
#           " columns matching ",
#           required_column,
#           ":",
#           op[[required_column]],
#           ". Instead of 1."
#         )
#       )
#       colnames(mut_table)[matching_column_indices] <-
#         op[[required_column]]
#     }
#   }
#   mut_table
# }