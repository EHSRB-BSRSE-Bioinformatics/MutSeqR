# A list of reference contexts at different resolutions

A list of reference contexts at different resolutions

## Usage

``` r
context_list
```

## Format

A list with corresponding values

## Value

A named list whose values are either `NA` for the `none` and `type`
resolutions or character vectors of reference bases or trinucleotide
contexts.

## Examples

``` r
context_list[["base_6"]]
#> [1] "C" "T"
context_list[["base_12"]]
#> [1] "A" "C" "G" "T"
```
