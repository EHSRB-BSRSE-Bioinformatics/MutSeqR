# Populate Sequence context

This function populates the trinucleotide context for each mutation in
the mutation data.

## Usage

``` r
populate_sequence_context(mutation_granges, BS_genome, n = 1)
```

## Arguments

- mutation_granges:

  A GRanges object containing mutation data.

- BS_genome:

  The name of the Bioconductor BSgenome package to use for retrieving
  the reference genome sequence.

- n:

  An integer value indicating the number of base pairs to include on
  either side of the mutation for context. Default is 1 (trinucleotide
  context including the mutation).

## Value

A GRanges object with an additional column for the trinucleotide
context.
