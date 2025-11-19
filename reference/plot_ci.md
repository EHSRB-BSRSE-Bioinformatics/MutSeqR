# plot_ci

Plot confidence intervals

## Usage

``` r
plot_ci(
  data,
  order = "none",
  custom_order = NULL,
  nudge = 0.3,
  log_scale = FALSE,
  x_lab = NULL,
  y_lab = NULL,
  title = NULL
)
```

## Arguments

- data:

  A data frame with the results of the BMD analysis. Data must contain
  columns "Response", "BMD", "BMDL", and "BMDU". BMD values can be NA.

- order:

  Indicates how the responses should be ordered. Options are "none"
  (default), "asc" for ascending BMD values, "desc" for descending BMD
  values, or a custom order.

- custom_order:

  A character vector with the custom order of the Responses.

- nudge:

  A numeric value to nudge the text labels away from points. Default is
  0.3.

- log_scale:

  A logical value indicating if the x-axis should be in log10 scale.
  Default is false.

- x_lab:

  A character string with the x-axis label. Default is "BMD" or
  "log10(BMD)" if log_scale is TRUE.

- y_lab:

  A character string with the y-axis label. Default is "Response".

- title:

  A character string with the plot title. Default is "BMD with 90%
  Confidence Intervals".

## Value

a ggplot object

## Examples

``` r
if (requireNamespace("MutSeqRData", quietly = TRUE)) {
  # Plot results from PROAST
  dat <- data.frame(
    Response = c("PROAST MF Min", "PROAST MF Max"),
    BMD = c(NA, NA),
    BMDL = c(7.38, 2.98),
    BMDU = c(10.9, 7.68)
  )
  plot <- plot_ci(dat)
}
```
