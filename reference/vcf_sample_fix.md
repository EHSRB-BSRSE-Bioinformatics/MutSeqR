# Retrieve the sample column from VCF files

Checks to find the sample name of the vcf in the INFO field or in the
FORMAT header. Can also handle sample name synonyms.

## Usage

``` r
vcf_sample_fix(vcf)
```

## Arguments

- vcf:

  The imported VCF

## Value

The vcf with sample column name corrected
