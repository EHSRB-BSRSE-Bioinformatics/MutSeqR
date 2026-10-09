# rdm_fast

Rapidly simulates count matrices drawn from a Dirichlet-Multinomial
distribution. This function is highly optimized for parametric
bootstrapping, utilizing a Gamma-Multinomial approximation for speed.

## Usage

``` r
rdm_fast(depths, p, theta)
```

## Arguments

- depths:

  A numeric vector representing the total sequencing depth (total
  mutation counts) for each sample to be simulated.

- p:

  A numeric vector representing the expected probability of each
  subtype.

- theta:

  Numeric scalar. The overdispersion parameter (concentration
  parameter).

## Value

A numeric matrix of simulated mutation counts, with subtypes as rows and
simulated samples as columns.
