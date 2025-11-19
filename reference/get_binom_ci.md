# Add binomial confidence intervals to mutation frequencies.

Uses the binomial distribution to create confidence intervals for
mutation frequencies calculated from a single point estimate.
Calculating binomial confidence intervals for mutation frequencies is
not part of MutSeqR's recommended workflow, but is provided here for
users who wish to use it.

## Usage

``` r
get_binom_ci(
  mf_data,
  sum_col = "sum_min",
  depth_col = "group_depth",
  conf_level = 0.95,
  method = "wilson"
)
```

## Arguments

- mf_data:

  The data frame containing the mutation frequencies per sample.
  Obtained as an output from `calculate_mf`.

- sum_col:

  Column name that specifies the mutation count (e.g., sum_min)

- depth_col:

  Column name that specifies the sequencing depth (e.g., total_depth)

- conf_level:

  Confidence interval to calculate, default 95% (0.95)

- method:

  The method used by binom::binom.confint to calculate intervals.
  Default is "wilson" (recommended).

## Value

A mf data frame with added columns indicating the confidence intervals.

## Examples

``` r
if (requireNamespace("MutSeqRData", quietly = TRUE)) {
  # Example data consists of 24 mouse bone marrow DNA samples imported
  # using import_mut_data() and filtered with filter_mut as in Example 4.
  # Sequenced on TS Mouse Mutagenesis Panel. Example data is
  # retrieved from MutSeqRData, an ExperimentHub data package.
  library(ExperimentHub)
  eh <- ExperimentHub()
  example_data <- eh[["EH9861"]]

  mf <- calculate_mf(example_data)
  confint <- get_binom_ci(
    mf_data = mf,
    sum_col = "sum_min",
    depth_col = "group_depth"
  )
}
```
