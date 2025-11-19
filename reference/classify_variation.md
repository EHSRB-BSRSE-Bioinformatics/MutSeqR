# classify_variation

Classify the variation type of a mutation based on its ref and alt
values.

## Usage

``` r
classify_variation(ref, alt)
```

## Arguments

- ref:

  The reference allele.

- alt:

  The alternate allele.

## Value

A character indicating the type of variation.

## Examples

``` r
df <- data.frame(
  ref = c("A", "CAGT", "GCC", "T", "ACG", "C", "G", "T", "A"),
  alt = c("R", "TGA", "G", "TC", "TAC", "C", "<DEL>", "G", "???")
)
df$variation_type <- mapply(classify_variation, df$ref, df$alt)
df
#>    ref   alt variation_type
#> 1    A     R      ambiguous
#> 2 CAGT   TGA        complex
#> 3  GCC     G       deletion
#> 4    T    TC      insertion
#> 5  ACG   TAC            mnv
#> 6    C     C     no_variant
#> 7    G <DEL>             sv
#> 8    T     G            snv
#> 9    A   ???  uncategorized
```
