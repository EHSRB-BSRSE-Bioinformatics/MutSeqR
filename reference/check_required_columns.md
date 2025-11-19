# Check that all required columns are present before proceeding with the function

A utility function that will check that all required columns are
present.

## Usage

``` r
check_required_columns(data, required_columns)
```

## Arguments

- data:

  mutation data

- required_columns:

  a list of required column names.

## Value

an error

## Examples

``` r
df <- data.frame(
  contig = c("chr1", "chr2", "chr3"),
  start = c(100, 200, 300),
  end = c(100, 200, 300),
  sample = c("S1", "S2", "S3"),
  ref = c("G", "C", "T"),
  alt = c("A", "T", "G")
)
check_required_columns(df, required_columns = op$base_required_mut_cols)
#>   contig start end sample ref alt
#> 1   chr1   100 100     S1   G   A
#> 2   chr2   200 200     S2   C   T
#> 3   chr3   300 300     S3   T   G
```
