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
if (requireNamespace("MutSeqRData", quietly = TRUE)) {
  # Plot the mean MFmin of each genomic target per dose group
  # Order the genomic targets by their genic context.

  # Example data consists of 24 mouse bone marrow DNA samples imported
  # using import_mut_data() and filtered with filter_mut as in Example 4.
  # Sequenced on TS Mouse Mutagenesis Panel. Example data is
  # retrieved from MutSeqRData, an ExperimentHub data package.
  library(ExperimentHub)
  eh <- ExperimentHub()
  example_data <- eh[["EH9861"]]

  mf <- calculate_mf(
    mutation_data = example_data,
    cols_to_group = c("sample", "label"),
    retain_metadata_cols = c("dose_group", "genic_context")
  )
  # Define the order of the genomic targets
  label_order <- mf %>%
    dplyr::arrange(genic_context) %>%
    dplyr::pull(label) %>%
    unique()
  # Calculate the mean MF per dose_group for each target.
  mean <- mf %>%
    dplyr::group_by(dose_group, label) %>%
    dplyr::summarise(mean = mean(mf_min))
  # Set the order of each column
  mean$dose_group <- factor(mean$dose_group,
    levels = c(
      "Control",
      "Low",
      "Medium",
      "High"
    )
  )
  mean$label <- factor(mean$label,
    levels = label_order
  )
  # Plot
  plot <- plot_radar(
    mf_data = mean,
    response_col = "mean",
    label_col = "label",
    facet_col = "dose_group",
    indiv_y = FALSE
  )
}
```
