# Plot the Mean Mutatation Frequency

This function calculates the mean mutation frequency across samples for
given groups and plots the results.

## Usage

``` r
plot_mean_mf(
  mf_data,
  group_col = "dose",
  fill_col = NULL,
  mf_type = "both",
  plot_type = "line",
  plot_error_bars = TRUE,
  plot_indiv_vals = TRUE,
  group_order = "none",
  group_order_input = NULL,
  add_labels = "mean_count",
  scale_y_axis = "linear",
  x_lab = NULL,
  y_lab = NULL,
  plot_title = NULL,
  custom_palette = NULL,
  plot_legend = TRUE,
  rotate_labels = FALSE,
  label_size = 3
)
```

## Arguments

- mf_data:

  A data frame containing the mutation frequency data. This is obtained
  from the calculate_mf function with SUMMARY = TRUE.

- group_col:

  The column(s) in mf_data by which to calculate the mean. When
  supplying more than one column, the values of all group columns will
  be concatenated into a single value by which to calculate the mean.
  Values will be displayed along the x-axis. Ex. "dose" or c("dose",
  "tissue").

- fill_col:

  An optional column name in the data used to define the fill aesthetic
  in the plot. If fill_col has multiple levels within each group_col
  level, the mean will be calculated for each level of fill_col
  (recommend plot_type = "line" for this use case). Default is NULL.

- mf_type:

  The type of mutation frequency to plot. Options are "min", "max",
  "both", or "stacked". If "both", the min and max mutation frequencies
  are plotted side by side. "stacked" can be chosen for bar plot_type
  only. If "stacked", the difference between the min and max MF is
  stacked on top of the min MF such that the total height of both bars
  represent the max MF. Default is "min".

- plot_type:

  The type of plot to create. Options are "bar" or "line". Default is
  "bar".

- plot_error_bars:

  Whether to plot the error bars. Default is TRUE. Error bars are
  standard error of the mean.

- plot_indiv_vals:

  Whether to plot the individual values as data points. Default is
  FALSE.

- group_order:

  The order of the groups along the x-axis. ' Options include:

  - `none`: No ordering is performed. Default.

  - `smart`: Groups are ordered based on the sample names.

  - `arranged`: Groups are ordered based on one or more factor column(s)
    in mf_data. Factor column names are passed to the function using the
    `group_order_input`.

  - `custom`: Groups are ordered based on a custom vector of group
    names. The custom vector is passed to the function using the
    `group_order_input`.

- group_order_input:

  The order of the groups if group_order is "custom". The column name by
  which to arrange groups if group_order is "arranged". If "custom", and
  using more than one group_col, values are concatenated in the order
  listed, separated by a "\_".

- add_labels:

  The data labels to display on the plot. Either "indiv_count",
  "indiv_MF", "mean_count", "mean_MF", or "none". Count labels display
  the number of mutations, MF labels display the mutation frequency.
  Mean plots the mean value. Indiv plots the labels for individual data
  points (only if plot_indiv_vals = TRUE). Default is "none".

- scale_y_axis:

  The scale of the y axis. Either "linear" or "log". Default is
  "linear".

- x_lab:

  The x-axis label. Default is the value of group_col.

- y_lab:

  The y-axis label. Default is "Mutation Frequency (mutations/bp)".

- plot_title:

  The title of the plot. Default is "Mean Mutation Frequency".

- custom_palette:

  A custom color palette to use for the plot. Input a character vector
  of colours. Input a named character vector to specify olours to
  specific groups. Fill labels will be constructed by the following
  components

  1.  "Mean/Individual" if plot_indiv_vals = TRUE, fill labels will
      specify Mean/Individual values.

  2.  "min/max" if mf_type = "both" or "stacked", fill labels will
      specify min/max values.

  3.  fill_col value. Name colours to match the fill labels. Default is
      NULL. If no custom_palette, a rainbow palette is generated.
      Min/Max values and Mean/Individual values will be the same colour,
      different shades.

- plot_legend:

  Logical. Whether to show the fill (and color) legend. Default is TRUE.

- rotate_labels:

  A logical value indicating whether data labels should be rotated 90
  degrees. Default is FALSE.

- label_size:

  A numeric value that controls the size of the data labels.

## Value

a ggplot object

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
# Plot the mean min MF per dose group as a bar plot with error bars
plot <- plot_mean_mf(
  mf_data = mf_example,
  group_col = "dose_group",
  mf_type = "min",
  plot_type = "line",
  fill_col = "dose_group",
  plot_error_bars = TRUE,
  plot_indiv_vals = TRUE,
  add_labels = "none"
)
```
