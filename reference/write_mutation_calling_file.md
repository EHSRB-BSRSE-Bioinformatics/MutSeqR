# Write the mutation calling file to input into the SigProfiler Assignment web application.

Creates a .txt file from mutation data that can be used for mutational
signatures analysis using the SigProfiler Assignment web application.
Currently only supports SBS analysis i.e. snvs.

## Usage

``` r
write_mutation_calling_file(
  mutation_data,
  project_name = "Example",
  project_genome = "GRCm38",
  output_path = NULL
)
```

## Arguments

- mutation_data:

  The object containing the mutation data. The output of
  import_mut_data() or import_vcf_data().

- project_name:

  The name of the project. Default is "Example".

- project_genome:

  The reference genome to use. (e.g., Human: GRCh38, Mouse mm10: GRCm38)

- output_path:

  The path to save the output file. If NULL, files will be saved in the
  current working directory. Default is NULL.

## Value

a .txt file that can be uploaded to the SigProfiler Assignment web
application (https://cancer.sanger.ac.uk/signatures/assignment/) as a
"Mutational calling file".

## Details

Mutations will be be filtered for SNVs. Mutations flagged in
`filter_mut` will be excluded from the output.

## Examples

``` r
if (requireNamespace("MutSeqRData", quietly = TRUE)) {
  library(ExperimentHub)
  eh <- ExperimentHub()
  example_data <- eh[["EH9861"]]

  temp_output <- tempdir()
  write_mutation_calling_file(
    mutation_data = example_data,
    project_name = "Example",
    project_genome = "GRCm38",
    output_path = temp_output
  )
  list.files(temp_output)
  # The file is saved in the temporary directory
  # To view the file, use the following code:
  ## output_file <- file.path(temp_output, "mutation_calling_file.txt")
  ## file.show(output_file
}
```
