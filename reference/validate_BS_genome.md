# Validate BSgenome Input

Internal utility function to validate the `BS_genome` argument prior to
sequence context extraction. Ensures that the provided genome is a valid
BSgenome package name and that it is installed locally.

## Usage

``` r
validate_BS_genome(BS_genome)
```

## Arguments

- BS_genome:

  A character string specifying the package name of a BSgenome object
  (e.g., `"BSgenome.Hsapiens.UCSC.hg38"`), or `NULL`.

## Value

Invisibly returns `TRUE` if validation passes; otherwise, an error is
raised.

## Details

This function performs three checks:

1.  If `BS_genome` is `NULL`, an error is thrown indicating that a
    genome must be provided when sequence context is required.

2.  If `BS_genome` is not among the available BSgenome packages, an
    error is thrown.

3.  If `BS_genome` is valid but not installed locally, an error is
    thrown with instructions to install it via
    [`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html).

This function is intended to be called only when sequence context needs
to be populated (i.e., when a `context` column is absent or incomplete).
