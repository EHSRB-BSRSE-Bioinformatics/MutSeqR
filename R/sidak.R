#' Correct p-values for multiple comparisons
#' 
#' @param vecP vector of p-values
#' @returns adjusted p-values
#' @details This function corrects a vector of probabilities for multiple testing
#' using the Bonferroni (1935) and Sidak (1967) corrections.
#' References: Bonferroni (1935), Sidak (1967), Wright (1992).
#' Bonferroni, C. E. 1935. Il calcolo delle assicurazioni su gruppi di teste.
#' Pp. 13-60 in: Studi in onore del Professore Salvatore Ortu Carboni. Roma.
#' Sidak, Z. 1967. Rectangular confidence regions for the means of multivariate
#' normal distributions. Journal of the American Statistical Association 62:626-633.
#' Wright, S. P. 1992. Adjusted P-values for simultaneous inference.
#' Biometrics 48: 1005-1013.
#' Pierre Legendre, May 2007
#' @export

sidak <- function(vecP) {
  k <- length(vecP)
  vecPB <- 0
  vecPS <- 0
  for(i in 1:k) {
    bonf = vecP[i] * k
      if(bonf > 1) bonf = 1
      vecPB <- c(vecPB, bonf)
      vecPS <- c(vecPS, (1- (1 - vecP[i])^k))
    }
    return(list(OriginalP = vecP, BonfP = vecPB[-1], SidakP = vecPS[-1]))
}