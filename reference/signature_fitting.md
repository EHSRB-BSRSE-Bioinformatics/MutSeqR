# Run COSMIC signatures comparison using SigProfilerAssignment

Run COSMIC signatures comparison using SigProfilerAssignment

## Usage

``` r
signature_fitting(
  mutation_data,
  project_name = "Default",
  project_genome = "GRCh38",
  env_name = "MutSeqR",
  group = "sample",
  output_path = NULL,
  python_version,
  dry_run = FALSE
)
```

## Arguments

- mutation_data:

  A data frame containing mutation data.

- project_name:

  The name of the project. This is used to format the data into required
  .txt format for SigProfiler tools.

- project_genome:

  The reference genome to use. On first use, the function will install
  the genome using SigProfilerMatrixGenerator::install. e.x. GRCh37,
  GRCH38, mm10, mm9, rn6

- env_name:

  The name of the virtual environment. This will be created on first
  use.

- group:

  The column in the mutation data used to aggregate groups. Signature
  assignment will be performed on each group separately.

- output_path:

  The filepath to the directory in which the output folder will be
  created to store results. Default is NULL. This will store results in
  the current working directory.

- python_version:

  The version of python installed on the user's computer.

- dry_run:

  Logical. If `TRUE`, validate and format the mutation data, then return
  the SigProfiler input data without configuring Python, installing
  dependencies or reference genomes, or running SigProfiler. Default is
  `FALSE`.

## Value

If `dry_run = TRUE`, a data frame formatted for SigProfiler input.
Otherwise, creates a subfolder "SigProfiler" in the output directory
with SigProfiler tools results. For a complete breakdown of the results,
see the Readme file for MutSeqR. Most relevant results are stored in
SigProfiler \> group \> matrices \> output \> Assignment_Solution \>
Activities \> SampleReconstruction \> WebPNGs. These plots show a
summary of the signature assignment results for each group. In each
plot, the top left panel represents the base_96 mutation count for the
group. The bottom left panel represents the reconstructed profile. Below
the reconstruction are the solution statistics that indicate the
goodness of fit of the reconstructed profile to the observed profile.
(Recommended cosine similarity \> 0.9). The panels on the right
represent the SBS signatures that contribute to the reconstructed
profile. The signature name and its contribution % are shown in the
panel. A high contribution means a high association of the signature
with the group's mutation spectra.

## Details

Assign COSMIC SBS signatures to mutation data using
SigProfilerAssignment. Data is cleaned and formatted for input into
SigProfiler tools. This function will create a virtual environment using
reticulate to run python, as this is a requirement for the SigProfiler
suite of tools. Note that it will also install several python
dependencies using a conda virtual environment on first use. Please be
aware of the implications of this. For advanced use, it is suggested to
use the SigProfiler python tools directly in python as described in
their respective documentation. Users must have python installed on
their computer to use this function.

Mutation data will be filtered to only include SNVs. Variants flagged by
the filter_mut column will be excluded.

## Examples

``` r
example_data <- readRDS(system.file(
  "extdata", "Example_files", "filtered_simple_mutation_data.rds",
  package = "MutSeqR"
))
sigprofiler_input <- signature_fitting(
  mutation_data = example_data,
  project_name = "Example",
  project_genome = "mm10",
  group = "sample",
  dry_run = TRUE
)
#> Creating cleaned data for input into SigProfiler...
head(sigprofiler_input)
#>   Project        Samples ID Genome mut_type chrom pos_start pos_end ref alt
#> 1 Example sample_sample1  .   mm10      SNP     1       102     102   C   A
#> 2 Example sample_sample1  .   mm10      SNP     1       104     104   T   C
#> 3 Example sample_sample2  .   mm10      SNP     2       202     202   T   A
#> 4 Example sample_sample2  .   mm10      SNP     2       206     206   G   T
#> 5 Example sample_sample3  .   mm10      SNP     1       106     106   G   A
#>      Type
#> 1 SOMATIC
#> 2 SOMATIC
#> 3 SOMATIC
#> 4 SOMATIC
#> 5 SOMATIC
```
