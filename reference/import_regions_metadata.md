# Join Regions Metadata

This function imports the regions metadata and joins it with the
mutation data.

## Usage

``` r
import_regions_metadata(
  mutation_granges,
  regions,
  rg_sep,
  is_0_based_rg,
  padding
)
```

## Arguments

- mutation_granges:

  A data frame containing mutation data.

- regions:

  The path to the file containing the regions metadata. Alternatively, a
  data frame can be provided directly.

- rg_sep:

  The separator used in the regions metadata file. Default is tab
  (`\t`).

- is_0_based_rg:

  A logical value indicating whether the regions file is 0-based (TRUE)
  or 1-based (FALSE). Default is FALSE.

- padding:

  An integer value indicating the number of base pairs to pad the
  regions on either side. Default is 0.

## Value

A GRanges object that combines the mutation data with the regions
metadata.
