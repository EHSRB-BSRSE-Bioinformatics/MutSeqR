# Derive end coordinates from a VCF object

`VariantAnnotation::rowRanges()` defines the range width from `REF`,
which is appropriate for small variants but does not honor
structural-variant span annotations such as `INFO/END`. This helper
prefers explicit `END` values from the INFO field and falls back to
`SVLEN` for symbolic non-insertion alleles when `END` is absent.

## Usage

``` r
get_vcf_end_positions(vcf)
```

## Arguments

- vcf:

  A `VCF` object from
  [`VariantAnnotation::readVcf()`](https://rdrr.io/pkg/VariantAnnotation/man/readVcf-methods.html).

## Value

An integer vector of end coordinates aligned to the rows of `vcf`.
