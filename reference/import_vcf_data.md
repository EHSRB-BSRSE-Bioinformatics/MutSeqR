# Import a VCF file

The function reads VCF file(s) and extracts the data into a dataframe.

## Usage

``` r
import_vcf_data(
  vcf_file,
  sample_data = NULL,
  sd_sep = "\t",
  regions = NULL,
  rg_sep = "\t",
  is_0_based_rg = FALSE,
  padding = 0,
  BS_genome = NULL,
  output_granges = FALSE
)
```

## Arguments

- vcf_file:

  The path to the .vcf (.gvcf, gzip, bgzip) to be imported. If you
  specify a directory, the function will attempt to read all files in
  the directory and combine them into a single table. VCF files should
  follow the VCF specifications, version 4.5. Multisample VCF files are
  not supported; VCF files must contain one sample each. Required fields
  are listed in details.

- sample_data:

  An optional file containing additional sample metadata (dose,
  timepoint, etc.). This can be a data frame or a file path. Metadata
  will be joined with the mutation data based on the sample column.
  Required columns are `sample` and any additional columns you wish to
  include.

- sd_sep:

  The delimiter for importing sample metadata tables. Default is
  tab-delimited.

- regions:

  An optional file containing metadata of genomic regions. Region
  metadata will be joined with mutation data and variants will be
  checked for overlap with the regions. `regions` can be either a file
  path, a data frame, or a GRanges object. File paths will be read using
  the rg_sep. Users can also choose from the built-in TwinStrand's
  Mutagenesis Panels by inputting "TSpanel_human", "TSpanel_mouse", or
  "TSpanel_rat". Required columns for the regions file are "contig",
  "start", and "end". For a GRanges object, the required columns are
  "seqnames", "start", and "end". Default is NULL.

- rg_sep:

  The delimiter for importing the custom_regions. The default is
  tab-delimited "\t".

- is_0_based_rg:

  A logical variable. Indicates whether the position coordinates in
  `regions` are 0 based (TRUE) or 1 based (FALSE). If TRUE, positions
  will be converted to 1-based (start + 1). Need not be supplied for
  TSpanels. Default is TRUE.

- padding:

  Extend the range of your regions in both directions by the given
  amount. Ex. Structural variants and indels may start outside of the
  regions. Adjust the `padding` to include these variants in your
  region's ranges.

- BS_genome:

  The pkgname of a BS genome. A BS genome must be installed prior to
  import to populate the context column (trinucleotide context for each
  position). Only required if data does not already include a context
  column. Please install the appropriate BS genome using
  BiocManager::install("pkgname") where pkgname is the name of the
  BSgenome package. The pkgname can be found using the find_BS_genome()
  function, which requires the species and assembly version. Ex.
  "BSgenome.Hsapiens.UCSC.hg38" \| "BSgenome.Hsapiens.UCSC.hg19" \|
  "BSgenome.Mmusculus.UCSC.mm10" \| "BSgenome.Mmusculus.UCSC.mm39" \|
  "BSgenome.Rnorvegicus.UCSC.rn6"

- output_granges:

  `TRUE` or `FALSE`; whether you want the mutation data to output as a
  GRanges object. Default output is as a dataframe.

## Value

A table where each row is a mutation, and columns indicate the location,
type, and other data. If `output_granges` is set to TRUE, the mutation
data will be returned as a GRanges object, otherwise mutation data is
returned as a dataframe.

Output Column Definitions:

- `short_ref`: The reference base at the start position.

- `normalized_ref`: The short_ref in C/T-base notation for this position
  (e.g. A -\> T, G -\> C).

- `context` The trinucleotide context at this position. Consists of the
  reference base and the two flanking bases (e.g. TAC).

- `normalized_context`: The trinucleotide context in C/T base notation
  for this position (e.g. TAG -\> CTA).

- `variation_type` The type of variant (snv, mnv, insertion, deletion,
  complex, sv, no_variant, ambiguous, uncategorized).

- `subtype` The substitution type for the snv variant (12-base spectrum;
  e.g. A\>C).

- `normalized_subtype` The C/T-based substitution type for the snv
  variant (6-base spectrum; e.g. A\>C -\> T\>G).

- `context_with_mutation`: The substitution type for the snv variant
  including the two flanking nucleotides (192-trinucleotide spectrum;
  e.g. `T[A>C]G`)

- `normalized_context_with_mutation`: The C/T-based substitution type
  for the snv variant including the two flanking nucleotides (96-base
  spectrum e.g. `T[A>C]G` -\> `C[T>G]A`).

- `nchar_ref`: The length (in bp) of the reference allele.

- `nchar_alt`: The length (in bp) of the alternate allele.

- `varlen`: The length (in bp) of the variant.

- `ref_depth`: The depth of the reference allele. Calculated as
  `total_depth` - `alt_depth`, if applicable.

- `vaf` : The variant allele fraction. Calculated as
  `alt_depth`/`total_depth`.

- `gc_content`: % GC of the trinucleotide context at this position.

- `is_known`: TRUE or FALSE. Flags known variants (ID != ".").

- `row_has_duplicate`: TRUE or FALSE. Flags rows whose position is the
  same as that of at least one other row for the same sample.

- `filter_mut` : A logical value, initially set to FALSE that indicates
  to calculte_mf() if the variant should be excluded from mutation
  counts. See the filter_mut function for more detail.

## Details

The required fields are:

**FIXED FIELDS**

- `CHROM`: The name of the reference sequence. Equivalent to `contig`.

- `POS`: The 1-based start position of the feature. Equivalent to
  `start`.

- `REF`: The reference allele at this position.

- `ALT`: The left-aligned, normalized, alternate allele at this
  position. Multiple alt alleles called for a single position should be
  represented as separate rows in the table.

**INFO FIELDS**

- `END`: The half-open end position of the feature.

- `sample`: An identifying field for your samples; either in the INFO
  field or as the header to the FORMAT field.

**SUGGESTED FIELDS**

The following **FORMAT** fields are not required, but are recommended
for full package functionality:

- `AD`: The allelic depths for the reference and alternate allele in the
  order listed. The sum of AD is equivalent to the `total_depth` (read
  depth at this position excluding N-calls).

- `DP`: The read depth at this position (including N-calls). Equivalent
  to `depth`. Note that in many VCF files, the DP field is defined as
  `total_depth`. However, in most cases, the DP field includes N-calls.

- `VD`: The read depth supporting the alternate allele. If not included,
  the function will add this column, assuming a value of 1. Equivalent
  to `alt_depth`.

We recommend that files include a record for every sequenced position,
regardless of whether a variant was called, along with the `AD` for each
record. This enables site-specific depth calculations required for some
downstream analyses. AD is used to calculate the `total_depth` (the read
depth excluding No-calls). If AD is not available, the `DP` field will
be used as the `total_depth`.

## Examples

``` r
# Mutation data is just for example purposes. It does not reflect real data
file <- system.file("extdata", "Example_files", 
                   "simple_vcf_data.vcf", package = "MutSeqR")
# Import the data
imported_example_data <- import_vcf_data(
 vcf_file = file,
BS_genome = find_BS_genome("mouse", "mm10"))
#> 'getOption("repos")' replaces Bioconductor standard repositories, see
#> 'help("repositories", package = "BiocManager")' for details.
#> Replacement repositories:
#>     CRAN: https://cran.rstudio.com
#> Selected reference genome: BSgenome.Mmusculus.UCSC.mm10
#> Reference genome is already installed.
#> Once installed, supply 'BSgenome.Mmusculus.UCSC.mm10' as the BS_genome parameter.
#> 'getOption("repos")' replaces Bioconductor standard repositories, see
#> 'help("repositories", package = "BiocManager")' for details.
#> Replacement repositories:
#>     CRAN: https://cran.rstudio.com
#> Warning: info fields with no header: sample
#> Expected 'alt' but found 'alt.value', renaming it.
#> Expected 'alt_depth' but found 'VD', renaming it.
#> Loading reference genome: BSgenome.Mmusculus.UCSC.mm10.
#> Retrieving context sequences from BSgenome
```
