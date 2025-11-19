# Get the coordinates of the CpG sites within your genomic regions

Filters the ranges of your genomic regions to find all positions with a
specific motif. The default is CpG sites, but can be customizable.

## Usage

``` r
get_cpg_regions(regions, motif = "CG")
```

## Arguments

- regions:

  A GRanges object containing the genomic regions of interest in which
  to look for CpG sites. Must have the metadata column "sequence"
  populated with the raw nucleotide sequence to search for CpGs. This
  object can be obtained using the get_seq() function.

- motif:

  Default "CG", which returns CpG sites. You could in theory use an
  arbitrary string to look at different motifs. Use with caution.

## Value

A GRanges object where each range is a CpG site (a subset of ranges from
the larger object provided to the function).
