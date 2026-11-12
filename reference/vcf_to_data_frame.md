# Extract the records of a single-sample VCF into a data frame

Builds one row per VCF record with the fixed fields, the FORMAT fields
of the VCF's one sample and the INFO fields. FORMAT fields with several
values per record, such as AD, are spread into one column per value
(`AD_1`, `AD_2`, ...). INFO fields that share a name with a fixed or
FORMAT column, ignoring case, are dropped in their favor.

## Usage

``` r
vcf_to_data_frame(vcf, add_chr = FALSE)
```

## Arguments

- vcf:

  A single-sample VCF read with
  [`VariantAnnotation::readVcf()`](https://rdrr.io/pkg/VariantAnnotation/man/readVcf-methods.html).

- add_chr:

  Whether to prepend "chr" to contig names missing it and rename "MT" to
  "chrM".

## Value

A data frame with one row per VCF record.
