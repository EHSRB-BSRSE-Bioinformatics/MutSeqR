# Write FASTA file of reference sequences.

Write FASTA file of reference sequences.

## Usage

``` r
write_reference_fasta(regions_gr, output_path = NULL)
```

## Arguments

- regions_gr:

  A GRanges object including the sequences of the reference regions
  included for the data. This can be generated from the `get_seq`
  function.

- output_path:

  The directory where the FASTA file should be written. Default is NULL,
  which will write the file to the current working directory.

## Value

Writes a FASTA reference file "reference_output.fasta". If multiple
ranges are included in the GRanges object, the sequences will be written
to a single FASTA file. Sequences names will be the seqnames (contig) of
the range.

## Details

Generate an arbitrary multi-sequence FASTA file from GRanges including
the reference sequences.

## Examples

``` r
# Write FASTA files for the 20 genomic target sequences
# of TwinStrand's Mouse Mutagenesis Panel.
output_path <- tempdir()
rg <- get_seq("TSpanel_mouse")
#> Loading reference genome: BSgenome.Mmusculus.UCSC.mm10.
write_reference_fasta(rg, output_path = output_path)
```
