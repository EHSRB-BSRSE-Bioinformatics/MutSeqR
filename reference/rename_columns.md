# Map column names of mutation data to default column names. A utility function that renames columns of mutation data to default columns names.

Map column names of mutation data to default column names. A utility
function that renames columns of mutation data to default columns names.

## Usage

``` r
rename_columns(data, column_map = op$column)
```

## Arguments

- data:

  mutation data

- column_map:

  a list that maps synonymous column names to their default.

## Value

the mutation data with column names changed to match default.

## Examples

``` r
df <- data.frame(
  chromosome = c("chr1", "chr2", "chr3"),
  pos = c(100, 200, 300),
  end = c(100, 200, 300),
  sample_id = c("S1", "S2", "S3"),
  reference = c("G", "C", "T"),
  alternate = c("A", "T", "G")
)
renamed_data <- rename_columns(df, column_map = op$column)
#> Expected 'alt' but found 'alternate', renaming it.
#> Expected 'contig' but found 'chromosome', renaming it.
#> Expected 'ref' but found 'reference', renaming it.
#> Expected 'sample' but found 'sample_id', renaming it.
#> Expected 'start' but found 'pos', renaming it.
```
