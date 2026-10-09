# dm_loglik_ridge_stable

Computes the robust, ridge-penalized log-likelihood for a
Dirichlet-Multinomial distribution. It safely handles massive counts
utilizing lgamma components and matrix operations. The ridge penalty
shrinks probability vectors toward a uniform distribution to stabilize
optimization for sparse or heavily overdispersed data.

## Usage

``` r
dm_loglik_ridge_stable(x, p, theta, lambda = 1)
```

## Arguments

- x:

  A numeric matrix of mutational counts (subtypes as rows, samples as
  columns).

- p:

  A numeric vector representing the expected probability of each
  subtype.

- theta:

  Numeric scalar. The overdispersion parameter (concentration
  parameter). Larger values of theta indicate lower overdispersion
  (approaching multinomial).

- lambda:

  Numeric scalar. The strength of the ridge penalty. Default is 1.0. Set
  to 0 for unpenalized (raw) likelihood calculation.

## Value

A numeric scalar representing the computed log-likelihood value,
adjusted by the penalty term.
