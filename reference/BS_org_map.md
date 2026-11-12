# BS genome organism dictionary

A dictionary keyed by BSgenome organism identifiers, with corresponding
scientific organism names and common-name aliases as values.

## Usage

``` r
BS_org_map
```

## Format

A named list with BSgenome organism identifiers as keys and character
vectors of scientific names and common-name aliases as values.

## Value

A named list mapping BSgenome organism identifiers to character vectors
of scientific names and common-name aliases.

## Examples

``` r
BS_org_map["Hsapiens"]
#> $Hsapiens
#> [1] "homo sapiens"         "homo sapiens sapiens" "human"               
#> 
```
