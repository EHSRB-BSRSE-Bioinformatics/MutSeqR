# Annotate CpG sites

A simple method to test whether your trinucleotide context contains a
CpG site. Vectorized version of Biostrings::vcountPattern is used.

## Usage

``` r
annotate_cpg_sites(mutation_data, motif = "CG", column_query = "context", ...)
```

## Arguments

- mutation_data:

  A dataframe or GRanges object containing the genomic regions of
  interest in which to look for CpG sites.

- motif:

  Default "CG", which returns CpG sites. You could in theory use an
  arbitrary string to look at different motifs. Use with caution. In
  this case the pattern being searched must be a column in the mutation
  data.

- column_query:

  Default "context" but can be any column in the mutation data that you
  wish to look for a motif in.

- ...:

  Additional arguments to vcountPattern()

## Value

A data frame with the same number of rows as there were ranges in the
input, but with an additional metadata column indicating CpG sites in
the target sequence of the mutation.
