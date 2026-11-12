# Subtype Resolution Dictionary

Associates the subtype resolution names with their corresponding column
names in the mutation data.

## Usage

``` r
subtype_dict
```

## Format

A vector with corresponding values

## Value

A named character vector mapping subtype resolutions to mutation-data
column names, with `NA` for the `none` resolution.

## Examples

``` r
subtype_dict["base_96"]
#>                            base_96 
#> "normalized_context_with_mutation" 
subtype_dict["type"]
#>             type 
#> "variation_type" 
```
