# Create a radar plot

Create a radar plot

## Usage

``` r
plot_radar(mf_data, response_col, label_col, facet_col, indiv_y = TRUE)
```

## Arguments

- mf_data:

  A data frame with the data to plot

- response_col:

  The column with the response values

- label_col:

  The column with the labels for the radar plot.

- facet_col:

  The column with the group to facet the radar plots.

- indiv_y:

  A logical indicating whether to use individual y-axis scales for each
  plot.

## Value

A radar plot

## Examples

``` r
# Plot the MFmin for each variation type, including the 6 SNV subtypes
# Facet the plots by dose group
mf_ex <- readRDS(system.file("extdata", "Example_files", "mf_data_6.rds",
     package = "MutSeqR")
)
plot <- plot_radar(
     mf_data = mf_ex,
     response_col = "mf_min",
     label_col = "normalized_subtype",
     facet_col = "dose_group",
     indiv_y = TRUE
)
```
