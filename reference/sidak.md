# Correct p-values for multiple comparisons

Correct p-values for multiple comparisons

## Usage

``` r
sidak(vecP)
```

## Arguments

- vecP:

  vector of p-values

## Value

adjusted p-values

## Details

This function corrects a vector of probabilities for multiple testing
using the Bonferroni (1935) and Sidak (1967) corrections. References:
Bonferroni (1935), Sidak (1967), Wright (1992). Bonferroni, C. E. 1935.
Il calcolo delle assicurazioni su gruppi di teste. Pp. 13-60 in: Studi
in onore del Professore Salvatore Ortu Carboni. Roma. Sidak, Z. 1967.
Rectangular confidence regions for the means of multivariate normal
distributions. Journal of the American Statistical Association
62:626-633. Wright, S. P. 1992. Adjusted P-values for simultaneous
inference. Biometrics 48: 1005-1013. Pierre Legendre, May 2007

## Examples

``` r
p_values <- c(0.01, 0.04, 0.03, 0.08, 0.05)
adjusted_p <- sidak(p_values)
adjusted_p$SidakP
#> [1] 0.04900995 0.18462730 0.14126597 0.34091848 0.22621906
```
