# Join Sample Metadata

This function imports the sample metadata and joins it with the mutation
data.

## Usage

``` r
import_sample_data(mutation_data, sample_data, sd_sep = "\t")
```

## Arguments

- mutation_data:

  A data frame containing mutation data.

- sample_data:

  The path to the file containing the sample metadata. Alternatively, a
  data frame can be provided directly.

- sd_sep:

  The separator used in the sample metadata file. Default is tab (`\t`).

## Value

A data frame that combines the mutation data with the sample metadata.
