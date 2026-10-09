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
#>  [1] "bslib-06dfa953f10b44af4b6ca2dc86bfa180"                                                        
#>  [2] "downlit"                                                                                       
#>  [3] "file21df114437c5"                                                                              
#>  [4] "file21df35e60300"                                                                              
#>  [5] "file21df4a34f7e3"                                                                              
#>  [6] "file21df4d78e56f"                                                                              
#>  [7] "libloc_187_10f1dc1f01615a9f.rds"                                                               
#>  [8] "libloc_192_4da063ac4afb603c.rds"                                                               
#>  [9] "libloc_192_5fc0888d6df7bedc.rds"                                                               
#> [10] "libloc_201_893962b1c5113ec4.rds"                                                               
#> [11] "mutation_calling_file.txt"                                                                     
#> [12] "repos_https%3A%2F%2Fbioconductor.org%2Fpackages%2F3.23%2Fdata%2Fannotation%2Fsrc%2Fcontrib.rds"
#> [13] "temp_libpath21df27f322a8"                                                                      
#> [14] "test_list.xlsx"                                                                                
#> [15] "test_single.xlsx"                                                                              
  # The file is saved in the temporary directory
  # To view the file, use the following code:
  ## output_file <- file.path(tempdir(), "mutation_calling_file.txt")
  ## file.show(output_file)
```
