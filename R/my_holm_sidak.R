#' Correct p-values for multiple comparisons
#' 
#' @param P p-value
#' @returns adjusted p-values
#' @export

# Holm-Sidak Correction for multiple comparison
my.holm.sidak <- function(P){
  m <- length(P)
  if(m > 1){
    Psort <- matrix(c(P, 1:m), 2, m, byrow = TRUE)
    Psort <- Psort[, order(Psort[1, ])]
    for (i in 1:m) {
      adjust <- m + 1 - i
      Psort[1, i] <- pmin(1, (1 - ((1 - Psort[1, i])^adjust)))
    }
    Psort <- Psort[, order(Psort[2, ])]
    P.adjust <- Psort[1, ]
    P.adjust
  } else{
    #No adjustment as there is only one comparison
    P.adjust <- P
    P.adjust
  }}