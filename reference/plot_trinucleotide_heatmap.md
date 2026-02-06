# Create a heatmap plot of mutation subtype proportions.

This function creates a heatmap plot of subtype proportions for a given
grouping variable. The groups may be facetted by a second variable.
Mutation sums for each facet group and normalized subtype are calculated
and displayed.

## Usage

``` r
plot_trinucleotide_heatmap(
  mf_data,
  group_col = "sample",
  facet_col = "dose",
  mf_type = "min",
  mut_proportion_scale = "turbo",
  max = 0.2,
  rescale_data = FALSE,
  condensed = FALSE
)
```

## Arguments

- mf_data:

  A data frame containing the mutation frequency data at the desired
  base resolution. This is obtained using the 'calculate_mf' with
  subtype_resolution set to the desired resolution. cols_to_group should
  be the same as 'group_col'.

- group_col:

  The variable to group by.

- facet_col:

  The variable to facet by.

- mf_type:

  The type of mutation frequency to plot. Options are "min" or "max".
  (Default: "min")

- mut_proportion_scale:

  The scale option for the mutation proportion. Options are passed to
  viridis::scale_fill_viridis_c. One of \# inferno, magma, plasma,
  viridis, cividis, turbo, mako, or rocket. We highly reccomend the
  default for its ability to disciminate hard to see patterns. (Default:
  "turbo")

- max:

  Maximum value used for plotting the proportions. Proportions that are
  higher will have the maximum colour. (Default: 0.2)

- rescale_data:

  Logical value indicating whether to rescale the mutation proportions
  to increase the dynamic range of colors shown on the plot. (Default:
  TRUE)

- condensed:

  More condensed plotting format. Default = FALSE.

## Value

A ggplot object representing the heatmap plot.

## Examples

``` r
mf_96 <- readRDS(system.file("extdata/Example_files/mf_data_96_sample.rds",
package = "MutSeqR"))
# define dose_group order
mf_96$dose_group <- factor(mf_96$dose_group,
     levels = c("Control", "Low","Medium", "High")
)
plot <- plot_trinucleotide_heatmap(mf_96,
  group_col = "sample",
  facet_col = "dose_group"
)
#> Cutting off at maximum mutation proportion value of 0.2
```
