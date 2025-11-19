# Write mutation_data to a VCF file

Export your mutation_data to a VCF file for downstream applications.

## Usage

``` r
write_vcf_from_mut(mutation_data, output_path = NULL)
```

## Arguments

- mutation_data:

  A data frame of a GRanges object containing your mutation data. This
  can be the output of `import_mut_data`, `import_vcf_data`, or
  `filter_mut.` Coordinates must be 1-based. Required columns are
  "contig", "start", "end", "ref", "alt", "sample", "alt_depth",
  "total_depth", and "ref_depth". Additional columns are allowed.

- output_path:

  The directory where the VCF file should be written. Default is NULL,
  which will write the file to the current working directory.

## Value

Writes a VCF file of mutations "mutation_output.vcf".

## Examples

``` r
if (requireNamespace("MutSeqRData", quietly = TRUE)) {
  # Example data consists of 24 mouse bone marrow DNA samples imported
  # using import_mut_data() and filtered with filter_mut as in Example 4.
  # Sequenced on TS Mouse Mutagenesis Panel. Example data is
  # retrieved from MutSeqRData, an ExperimentHub data package.
  library(ExperimentHub)
  eh <- ExperimentHub()
  example_data <- eh[["EH9861"]]

  write_vcf_from_mut(mutation_data = example_data, output_path = tempdir())
}
```
