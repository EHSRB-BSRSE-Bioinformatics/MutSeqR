#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
#' @examples 
#' df <- data.frame(x = 1:5, y = rnorm(5))
#' df %>% dplyr::mutate(z = x + y)
#' df %>% head(3) %>% summary()
NULL
