# Derive the alternate allele depth from AD

Used when the mutation data has no `alt_depth` (VD) column. The
alternate allele depth is the sum of the AD values after the first, and
0 for records without an alternate allele.

## Usage

``` r
alt_depth_from_ad(data, ad_columns)
```

## Arguments

- data:

  The mutation data, with `variation_type` and `AD_` columns.

- ad_columns:

  The `AD_` column names, in allele order.

## Value

A numeric vector of alternate allele depths, or `NULL` when AD has no
alternate allele column.
