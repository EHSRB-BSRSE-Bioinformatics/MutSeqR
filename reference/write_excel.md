# Write Excel tables

Writes data to an Excel file.

## Usage

``` r
write_excel(data, output_path = NULL, workbook_name, model_results = FALSE)
```

## Arguments

- data:

  A data frame or a list of data frames. If a data frame, it will be
  written to a single sheet in the Excel workbook. If a list, each data
  frame will be written to a separate sheet in the Excel workbook. Data
  may also be the output to model_mf, in which case set
  `model_results = TRUE`.

- output_path:

  The directory where the Excel file should be written. Default is NULL,
  which will write the file to the current working directory.

- workbook_name:

  The file name for the Excel file.

- model_results:

  A logical value indicating whether the data is the output of model_mf.
  Default is FALSE. If TRUE, the function will grab the model_data,
  point_estimates, and pairwise_comparisons data frames from the
  model_mf output and write them to separate sheets in the Excel
  workbook.

## Value

A saved Excel workbook.

## Examples

``` r
# Example data consists of 24 mouse bone marrow DNA samples imported
# using import_mut_data(), filtered with filter_mut, and summarized
# using calculate_mf().
outputpath <- tempdir()
mf_example <- readRDS(system.file("extdata/Example_files/mf_data_global.rds",
  package = "MutSeqR"
))
mf_example2 <- readRDS(system.file("extdata/Example_files/mf_data_6.rds",
  package = "MutSeqR"
))
mf_example3 <- readRDS(system.file("extdata/Example_files/mf_data_6_sample.rds",
  package = "MutSeqR"
))
list <- list(mf_example, mf_example2, mf_example3)
names(list) <- c("Global MF", "Base 6 Spectra", "Base 6 Sample Spectra")

# save a single data frame to an Excel file
write_excel(mf_example, output_path = outputpath, workbook_name = "test_single")

# save a list of data frames to an Excel file
write_excel(list, output_path = outputpath, workbook_name = "test_list")
```
