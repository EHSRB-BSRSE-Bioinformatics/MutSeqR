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
# Example data  consists of 24 mouse bone marrow
# samples exposed to three doses of BaP alongside vehicle controls.
# Libraries were sequenced with Duplex Sequencing using
# the TwinStrand Mouse Mutagenesis Panel which consists of 20 2.4kb
# targets = 48kb of sequence. Example data can be retrieved from
# MutSeqRData, an ExperimentHub data package:
## library(ExperimentHub)
## eh <- ExperimentHub()
## query(eh, "MutSeqRData")
# The data is a subset of variants from the target chr1
# from samples of the high dose group (50mg).
example_data <- readRDS(system.file("extdata", "Example_files",
                                    "variants_subset_d50_chr1.rds",
                                     package = "MutSeqR")
)
# Export mutation data of the four samples to a multi-sample VCF file.
write_vcf_from_mut(mutation_data = example_data, output_path = tempdir())
```
