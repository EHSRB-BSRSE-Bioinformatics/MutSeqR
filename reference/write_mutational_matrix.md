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

  write_mutational_matrix(
    mutation_data = example_data,
    group = "sample",
    subtype_resolution = "base_96",
    mf_type = "min",
    output_path = tempdir()
  )
#> Warning: No depth data provided. Mutation frequencies will not be calculated. Mutation subtype proportions will not be normalized to the total_depth.
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
#> [13] "sample_base_96_mutational_matrix.txt"                                                          
#> [14] "temp_libpath21df27f322a8"                                                                      
#> [15] "test_list.xlsx"                                                                                
#> [16] "test_single.xlsx"                                                                              
  # The file is saved in the temporary directory
  # To view the file, use the following code:
  ## output_file <- file.path(tempdir(), "sample_base_96_mutational_matrix.txt")
  ## file.show(output_file)
```
