# Get the reverse complement of a DNA or RNA sequence.

Get the reverse complement of a DNA or RNA sequence.

## Usage

``` r
reverseComplement(
  x,
  content = c("dna", "rna"),
  case = c("lower", "upper", "as is")
)
```

## Arguments

- x:

  A character vector of DNA or RNA sequences.

- content:

  c("dna", "rna") The type of sequence to be reversed.

- case:

  c("lower", "upper", "as is") The case of the output sequence.

## Value

A character vector of the reverse complement sequences.

## Details

This file is part of the source code for SPGS: an R package for
identifying statistical patterns in genomic sequences. Copyright (C)
2015 Universidad de Chile and INRIA-Chile A copy of Version 2 of the GNU
Public License is available in the share/licenses/gpl-2 file in the R
installation directory or from http://www.R-project.org/Licenses/GPL-2.
reverseComplement.R
