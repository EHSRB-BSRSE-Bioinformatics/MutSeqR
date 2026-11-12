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

If `BS_genome` is `NULL`, an error is thrown indicating that a genome
must be provided when sequence context is required. If the genome is
installed, validation passes without network access. Otherwise,
[`BSgenome::available.genomes()`](https://rdrr.io/pkg/BSgenome/man/available.genomes.html)
is queried (requires internet;) and an error is thrown with installation
instructions, or stating that the genome is not recognized.

This function is intended to be called only when sequence context needs
to be populated (i.e., when a `context` column is absent or incomplete).
