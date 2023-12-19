#' A utility function that will return the reference context of a mutation
#' @param mut_string the mutation. Ex. T>C, A[G>T]C
#' @export
get_ref_of_mut <- function(mut_string) {
  a <- str_extract(mut_string, ".*(?=\\s*>)")
  # Remove non-letter characters
  b <- str_replace_all(a, "[^a-zA-Z]", "")
  # Extract letter characters after square bracket
  c <- str_extract(mut_string, "\\](.*)") %>% str_replace_all("[^a-zA-Z]", "")
  if(is.na(c)) {
    return(b)
  } else {
    return(paste0(b,c))
  }
}
