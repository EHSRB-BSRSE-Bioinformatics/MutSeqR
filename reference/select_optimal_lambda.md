# select_optimal_lambda

Performs leave-one-out cross-validation (LOOCV) to empirically select
the optimal ridge penalty (lambda) for a given count matrix. Models are
fit on \$N-1\$ samples using a specified lambda, and unpenalized
likelihood is evaluated on the held-out sample.

## Usage

``` r
select_optimal_lambda(count_matrix, lambdas = c(0.1, 0.5, 1, 5, 10, 50, 100))
```

## Arguments

- count_matrix:

  A numeric matrix of mutation counts with subtypes as rows and samples
  as columns.

- lambdas:

  A numeric vector of lambda penalty values to evaluate. Default is
  c(0.1, 0.5, 1, 5, 10, 50, 100).

## Value

A list containing two items:

- `best_lambda`: Numeric scalar. The lambda value from the provided
  vector that yielded the highest overall cross-validation
  log-likelihood.

- `scores`: A named numeric vector detailing the cumulative
  log-likelihood score achieved by each tested lambda.
