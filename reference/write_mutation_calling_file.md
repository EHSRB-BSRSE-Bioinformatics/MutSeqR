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
 write_mutation_calling_file(
    mutation_data = example_data,
    project_name = "Example",
    project_genome = "GRCm38",
    output_path = tempdir()
  )
  list.files(tempdir())
#>  [1] "bslib-e392982257a0a7e1a4ae522df2f1f23b"                                                        
#>  [2] "downlit"                                                                                       
#>  [3] "file255744d7c73a"                                                                              
#>  [4] "file255749206024"                                                                              
#>  [5] "file255763f87a68"                                                                              
#>  [6] "file2557e240741"                                                                               
#>  [7] "libloc_177_2cc935f9605fe05a.rds"                                                               
#>  [8] "libloc_182_8ba9a5ce909cb448.rds"                                                               
#>  [9] "libloc_182_d3198563f59e5c50.rds"                                                               
#> [10] "libloc_191_2359b704beb3739.rds"                                                                
#> [11] "mutation_calling_file.txt"                                                                     
#> [12] "repos_https%3A%2F%2Fbioconductor.org%2Fpackages%2F3.23%2Fdata%2Fannotation%2Fsrc%2Fcontrib.rds"
#> [13] "temp_libpath2557288a0290"                                                                      
#> [14] "test_list.xlsx"                                                                                
#> [15] "test_single.xlsx"                                                                              
  # The file is saved in the temporary directory
  # To view the file, use the following code:
  ## output_file <- file.path(tempdir(), "mutation_calling_file.txt")
  ## file.show(output_file)
```
