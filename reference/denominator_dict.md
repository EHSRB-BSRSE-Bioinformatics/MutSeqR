# Values used for denominators in frequency calculations

These values are used to cross reference base substitution types to
their appropriate denominators for calculations. That is", "for example,
the 6 base substitution frequency should be subsetted based on the
normalized_ref column which would contain only T or C (i.e., the
pyrimidine context for base substitutions).

## Usage

``` r
denominator_dict
```

## Format

A vector with corresponding values

## Examples

``` r
denominator_dict["base_96"]
#>              base_96 
#> "normalized_context" 
denominator_dict["type"]
#> type 
#>   NA 
```
