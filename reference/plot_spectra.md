# Plot spectra

Given mf data, construct a plot displaying the mutation subtypes
observed in a cohort.

## Usage

``` r
plot_spectra(
  mf_data,
  group_col = "sample",
  subtype_resolution = "base_6",
  response = "proportion",
  mf_type = "min",
  group_order = "none",
  group_order_input = NULL,
  dist = "cosine",
  cluster_method = "ward.D",
  custom_palette = NULL,
  x_lab = NULL,
  y_lab = NULL,
  rotate_xlabs = FALSE
)
```

## Arguments

- mf_data:

  A data frame containing the mutation frequency data at the desired
  subtype resolution. This is obtained using the 'calculate_mf' function
  with subtype_resolution set to the desired resolution. Data must
  include a column containing the group_col, a column containing the
  mutation subtypes, a column containing the desired response variable
  (mf, proportion, sum) for the desired mf_type (min or max), and if
  applicable, a column containing the variable by which to order the
  samples/groups.

- group_col:

  The name of the column(s) in the mf data that contains the
  sample/group names. This will generally be the same values used for
  the cols_to_group argument in the calculate_mf function. However, you
  may also use groups that are at a higher level of the aggregation in
  mf_data.

- subtype_resolution:

  The subtype resolution of the mf data. Options are `base_6`,
  `base_12`, `base_96`, `base_192`, or `type`. Default is `base_6`.

- response:

  The desired response variable to be plotted. Options are mf,
  proportion, or sum. Default is `proportion`. Your mf_data must contain
  columns with the name of your desired response: `mf_min`, `mf_max`,
  `proportion_min`, `proportion_max`, `sum_min`, and `sum_max`.

- mf_type:

  The mutation counting method to use. Options are min or max. Default
  is `min`.

- group_order:

  The method for ordering the samples within the plot. Options include:

  - `none`: No ordering is performed. Default.

  - `smart`: Groups are automatically ordered based on the group names
    (alphabetical, numerical)

  - `arranged`: Groups are ordered based on one or more factor column(s)
    in mf_data. Column names are passed to the function using the
    `group_order_input`.

  - `custom`: Groups are ordered based on a custom vector of group
    names. The custom vector is passed to the function using the
    `group_order_input`.

  - `clustered`: Groups are ordered based on hierarchical clustering.
    The dissimilarity matrix can be specified using the `dist` argument.
    The agglomeration method can be specified using the `cluster_method`
    argument.

- group_order_input:

  A character vector specifying details for the group order method. If
  `group_order` is `arranged`, `group_order_input` should contain the
  column name(s) to be used for ordering the samples. If `group_order`
  is `custom`, `group_order_input` should contain the custom vector of
  group names.

- dist:

  The dissimilarity matrix for hierarchical clustering. Options are
  `cosine`, `euclidean`, `maximum`, `manhattan`, `canberra`, `binary` or
  `minkowski`. The default is `cosine`. See
  [dist](https://rdrr.io/r/stats/dist.html) for details.

- cluster_method:

  The agglomeration method for hierarchical clustering. Options are
  `ward.D`, `ward.D2`, `single`, `complete`, `average` (= UPGMA),
  `mcquitty` (= WPGMA), `median` (= WPGMC) or `centroid` (= UPGMC). The
  default is `Ward.D`. See [hclust](https://rdrr.io/r/stats/hclust.html)
  for details.

- custom_palette:

  A named vector of colors to be used for the mutation subtypes. The
  names of the vector should correspond to the mutation subtypes in the
  data. Alternatively, you can specify a color palette from the
  RColorBrewer package. See
  [`brewer.pal`](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html)
  for palette options. You may visualize the palettes at the ColorBrewer
  website: <https://colorbrewer2.org/>. Default is `NULL`.

- x_lab:

  The label for the x-axis. Default is the value of `group_col`.

- y_lab:

  The label for the y-axis. Default is the value of `response_col`.

- rotate_xlabs:

  A logical value indicating whether the x-axis labels should be rotated
  90 degrees. Default is FALSE.

## Value

A ggplot object representing the mutation spectra plot.

## Examples

``` r
# Example data consists of 24 mouse bone marrow DNA samples imported
# using import_mut_data() and filtered with filter_mut. Filtered
# mutation data is available in the MutSeqRData ExperimentHub package:
# eh <- ExperimentHub::ExperimentHub()
# Example 1: Visualized the 6-base mutation proportions per dose group.
# Data was summarized per dose_group using:
# calculate_mf(mutation_data = eh[["EH9861"]],
#              cols_to_group = "dose_group",
#              subtype_resolution = "base_6")
# Load the example data
mf_example <- readRDS(system.file("extdata", "Example_files", "mf_data_6.rds",
  package = "MutSeqR"
))
# Convert dose_group to a factor with the desired order.
mf_example$dose_group <- factor(mf_example$dose_group,
  levels = c("Control", "Low", "Medium", "High")
)
# Plot the mutation spectra
plot <- plot_spectra(
  mf_data = mf_example,
  group_col = "dose_group",
  subtype_resolution = "base_6",
  response = "proportion",
  group_order = "arranged",
  group_order_input = "dose_group"
)

# Example 2: plot the proportion of 6-based mutation subtypes
# for each sample, ordered by hierarchical clustering:
# Data was summarized per dose_group using:
# calculate_mf(mutation_data = eh[["EH9861"]],
#              cols_to_group = "sample",
#              subtype_resolution = "base_6")
# Load the example data
mf_example2 <- readRDS(system.file("extdata", "Example_files", "mf_data_6_sample.rds",
  package = "MutSeqR"
))
plot <- plot_spectra(
  mf_data = mf_example2,
  group_col = "sample",
  subtype_resolution = "base_6",
  response = "proportion",
  group_order = "clustered"
)
```
