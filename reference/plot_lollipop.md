# Plot recurrent mutations in a lollipop plot using ggplot2

Plot recurrent mutations in a lollipop plot using ggplot2

## Usage

``` r
plot_lollipop(
  mutation_data,
  min_recurrence = 2,
  group_col = "contig",
  subtype_resolution = "base_6",
  custom_palette = NULL
)
```

## Arguments

- mutation_data:

  A data frame containing mutation data.

- min_recurrence:

  An integer specifying the minimum number of times a mutation must be
  observed at the same position to be plotted. Default is 2.

- group_col:

  A character vector specifying the column name(s) to group
  mutation_data by. Default is "contig".

- subtype_resolution:

  The subtype resolution by which to group and colour the mutations.
  Options are "none", "type", "base_6", "base_12", "base_96", and
  "base_192".

- custom_palette:

  A named vector of colors to be used for the mutation subtypes. The
  names of the vector should correspond to the mutation subtypes in the
  data. Alternatively, you can specify a color palette from the
  RColorBrewer package. See
  [`brewer.pal`](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html)
  for palette options. You may visualize the palettes at the ColorBrewer
  website: <https://colorbrewer2.org/>. Default is `NULL`.

## Value

A list of ggplot objects.

## Examples

``` r
# For this example, we will use a subset of the example mutation data.
# The subset contains mutations from target chr1 in samples from the high
# dose group (50mg).
 example_data <- readRDS(system.file("extdata", "Example_files",
     "variants_subset_d50_chr1.rds", package = "MutSeqR"))
 # We will plot mutations that recoccur in at least two samples, grouped
 # by the "label" column, which signifies the target region (chr1).
 # Mutations will be grouped and coloured by their base 6 subtype (default) 
 plot <- plot_lollipop(
     mutation_data = example_data,
     min_recurrence = 2,
     group_col = "label"
 )
```
