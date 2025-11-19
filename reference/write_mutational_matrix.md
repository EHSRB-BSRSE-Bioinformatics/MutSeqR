# Write a Mutational Matrix to input into the sigprofiler web application

Creates a .txt file from mutation data that can be used for mutational
signatures analysis using the SigProfiler web application. Can handle
group analyses (ex dose, tissue, etc). Currently only supports SBS
matrices i.e. snvs.

## Usage

``` r
write_mutational_matrix(
  mutation_data,
  group = "dose",
  subtype_resolution = "base_96",
  mf_type = "min",
  output_path = NULL
)
```

## Arguments

- mutation_data:

  The object containing the mutation data. The output of
  import_mut_data() or import_vcf_data().

- group:

  The column in the mutation data used to aggregate groups (e.g.,
  sample, tissue, dose).

- subtype_resolution:

  The resolution of the mutation subtypes. Options are "base_6" or
  "base_96". Default is "base_96".

- mf_type:

  The mutation counting method to use. Options are "min" or "max".
  Default is "min".

- output_path:

  The path to save the output file. If not provided, the file will be
  saved in the current working directory. Default is NULL.

## Value

a .txt file that can be uploaded to the SigProfiler Assignment web
application (https://cancer.sanger.ac.uk/signatures/assignment/) as a
"Mutational Matrix".

## Details

Mutations will be be filtered for SNVs. Mutations flagged in
`filter_mut` will be excluded from the output. Mutations will be summed
across the groups specified in the `group` argument.

## Examples

``` r
if (requireNamespace("MutSeqRData", quietly = TRUE)) {
  library(ExperimentHub)
  eh <- ExperimentHub()
  example_data <- eh[["EH9861"]]
  temp_output <- tempdir()

  write_mutational_matrix(
    mutation_data = example_data,
    group = "dose_group",
    subtype_resolution = "base_96",
    mf_type = "min",
    output_path = temp_output
  )
  list.files(temp_output)
  # The file is saved in the temporary directory
  # To view the file, use the following code:
  ## output_file <- file.path(temp_output, "dose_group_base_96_mutational_matrix.txt")
  ## file.show(output_file)
}
```
