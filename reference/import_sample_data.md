# Import Sample Metadata

This function imports sample metadata from a file or accepts a data
frame directly, and performs basic validation checks.

## Usage

``` r
import_sample_data(sample_data, sd_sep = "\t")
```

## Arguments

- sample_data:

  The path to the file containing the sample metadata, or a data frame
  provided directly.

- sd_sep:

  The separator used in the sample metadata file. Default is tab (`\t`).

## Value

A validated data frame containing sample metadata, including a required
column named `sample`.
