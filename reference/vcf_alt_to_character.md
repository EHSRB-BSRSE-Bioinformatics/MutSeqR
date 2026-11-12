# Convert VCF alternate alleles to character values

Converts the list-like ALT field from a VCF object without dispatching
through Bioconductor's
[`as.data.frame()`](https://rdrr.io/pkg/BiocGenerics/man/as.data.frame.html)
methods. The conversion preserves one output value per VCF record.
Multiple alternate alleles in a record are joined with commas, and
records without an alternate allele are represented by `NA_character_`.

## Usage

``` r
vcf_alt_to_character(alt)
```

## Arguments

- alt:

  A list-like collection of alternate alleles, such as the value
  returned by
  [`VariantAnnotation::alt()`](https://rdrr.io/pkg/VariantAnnotation/man/VCF-class.html).

## Value

A character vector with one element per VCF record.
