# Plot the Mutation Frequency

This function creates a plot of the mutation frequency.

## Usage

``` r
plot_mf(
  mf_data,
  group_col,
  plot_type = "bar",
  mf_type = "min",
  fill_col = NULL,
  custom_palette = NULL,
  group_order = "none",
  group_order_input = NULL,
  labels = "count",
  scale_y_axis = "linear",
  x_lab = NULL,
  y_lab = NULL,
  title = NULL,
  rotate_labels = FALSE,
  label_size = 3
)
```

## Arguments

- mf_data:

  A data frame containing the mutation frequency data. This is obtained
  from the calculate_mf function with SUMMARY = TRUE.

- group_col:

  The name of the column containing the sample/group names for the
  x-axis.

- plot_type:

  The type of plot to create. Options are "bar" or "point".

- mf_type:

  The type of mutation frequency to plot. Options are "min", "max",
  "both", or "stacked". If "both", the min and max mutation frequencies
  are plotted side by side. "stacked" can be chosen for bar plot_type
  only. If "stacked", the difference between the min and max MF is
  stacked on top of the min MF such that the total height of both bars
  represent the max MF.

- fill_col:

  The name of the column containing the fill variable.

- custom_palette:

  A character vector of colour codes to use for the plot. If NULL, a
  default palette is used

- group_order:

  The order of the samples/groups along the x-axis. ' Options include:

  - `none`: No ordering is performed. Default.

  - `smart`: Samples are ordered based on the sample names.

  - `arranged`: Samples are ordered based on one or more factor
    column(s) in mf_data. Factor column names are passed to the function
    using the `group_order_input`.

  - `custom`: Samples are ordered based on a custom vector of sample
    names. The custom vector is passed to the function using the
    `group_order_input`.

- group_order_input:

  The order of the samples/groups if group_order is "custom". The column
  name by which to arrange samples/groups if group_order is "arranged"

- labels:

  The data labels to display on the plot. Either "count", "MF", or
  "none". Count labels display the number of mutations, MF labels
  display the mutation frequency.

- scale_y_axis:

  The scale of the y axis. Either "linear" or "log".

- x_lab:

  The label for the x axis.

- y_lab:

  The label for the y axis.

- title:

  The title of the plot.

- rotate_labels:

  A logical value aplied when labels is not "none". Indicates whether
  the labels should be rotated 90 degrees. Default is FALSE.

- label_size:

  A numeric value that adjusts the size of the labels. Default is 3.

## Value

A ggplot object

## Examples

``` r
# Example data  consists of 24 mouse bone marrow
# samples exposed to three doses of BaP alongside vehicle controls.
# Libraries were sequenced with Duplex Sequencing using
# the TwinStrand Mouse Mutagenesis Panel which consists of 20 2.4kb
# targets = 48kb of sequence. Example data can be retrieved from
# MutSeqRData, an ExperimentHub data package:
## library(ExperimentHub)
## eh <- ExperimentHub()
## query(eh, "MutSeqRData")
# Mutation frequency data was precalculated using
## mf_data_global <- calculate_mf(mutation_data = eh[["EH9861"]],
##   cols_to_group = "sample",
##   retain_metadata_cols = c("dose_group", "dose"))
# Load Example MF data
mf_example <- readRDS(system.file("extdata/Example_files/mf_data_global.rds",
  package = "MutSeqR"
))
# Specify the order of the dose groups along the x-axis
mf_example$dose_group <- factor(mf_example$dose_group,
  levels = c(
    "Control", "Low",
    "Medium", "High"
  )
)
# Plot the min MF per sample as a bar plot with count labels
plot <- plot_mf(
  mf_data = mf_example,
  group_col = "sample",
  plot_type = "bar",
  mf_type = "min",
  fill_col = "dose_group",
  group_order = "arranged",
  group_order_input = "dose_group",
  labels = "count",
  title = "Mutation Frequency per Sample"
)
```
