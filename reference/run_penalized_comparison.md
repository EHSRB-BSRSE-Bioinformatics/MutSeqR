# run_penalized_comparison

Core statistical engine for the RP-DMM test. Evaluates whether two count
matrices represent significantly different mutational spectra. It
utilizes a "hybrid logic": robust parameters (thetas) are estimated
utilizing a ridge penalty (QA phase), while the test statistic (LRT) is
evaluated using the unpenalized likelihood. Significance is evaluated
via parametric bootstrap.

## Usage

``` r
run_penalized_comparison(g1, g2, lambda = 1, n_boot = 500)
```

## Arguments

- g1:

  A numeric count matrix for Group 1 (subtypes as rows, samples as
  columns).

- g2:

  A numeric count matrix for Group 2 (subtypes as rows, samples as
  columns).

- lambda:

  Numeric. The ridge penalty strength used during parameter estimation.

- n_boot:

  Integer. The number of parametric bootstrap iterations to perform.

## Value

A list containing 5 items:

- `p_value`: Numeric. The bootstrapped p-value.

- `lrt`: Numeric. The observed unpenalized Likelihood Ratio Test
  statistic.

- `theta_g1`: Numeric. The optimized robust theta (overdispersion) for
  Group 1.

- `theta_g2`: Numeric. The optimized robust theta for Group 2.

- `theta_shared`: Numeric. The optimized robust theta under the Null
  hypothesis.
