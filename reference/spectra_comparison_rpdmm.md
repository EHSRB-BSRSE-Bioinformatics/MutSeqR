# spectra_comparison_rpdmm

Ridge-penalized Dirichlet-Multinomial test for comparing mutation
spectra between two groups. Robust for large mutation counts.

## Usage

``` r
spectra_comparison_rpdmm(
  mf_data,
  exp_variable,
  contrasts,
  cont_sep = "\t",
  mf_type = "min",
  lambda = 1,
  n_boot = 500
)
```

## Arguments

- mf_data:

  A data frame containing the MF data. This is the output from
  calculate_mf(). MF data should be calculate at the *sample level* and
  at the desired subtype resolution. Required columns are sample, the
  exp_variable column(s), the subtype column, and sum_min or sum_max.

- exp_variable:

  The column names of the experimental variable(s) to be compared.

- contrasts:

  A filepath (character) OR a data frame specifying the comparisons to
  be made. Must consist of exactly two columns. The level in the first
  column will be compared to the level in the second column for each
  row. If using multiple exp_variables, separate levels with a colon
  (e.g., "Drug:High").

- cont_sep:

  Character. The delimiter used to import the contrasts table if a
  filepath is provided. Default is tab.

- mf_type:

  Character. The type of mutation frequency count to use. Choices are
  "min" or "max". Default is "min" (recommended).

- lambda:

  Numeric. The strength of the ridge penalty. Default is 1.0. Higher
  values increase shrinkage towards a uniform distribution, stabilizing
  estimation in sparse datasets.

- n_boot:

  Integer. The number of bootstrap iterations to perform for p-value
  calculation. Default is 500.

## Value

A data frame containing one row per specified contrast. Columns include
the group names, comparison string, bootstrap p-value, observed
Likelihood Ratio Test (LRT) statistic, and robust overdispersion metrics
(theta) for group 1, group 2, and the shared null model.

## Details

Experimental: this function is under construction and its interface may
change.
