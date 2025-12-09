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
#' @examples
#' p_values <- c(0.01, 0.04, 0.03, 0.08, 0.05)
#' adjusted_p <- sidak(p_values)
#' adjusted_p$SidakP
sidak <- function(vecP) {
    k <- length(vecP)
    BonfP <- pmin(vecP * k, 1)
    SidakP <- 1 - (1 - vecP)^k
    return(list(OriginalP = vecP, BonfP = BonfP, SidakP = SidakP))
}
