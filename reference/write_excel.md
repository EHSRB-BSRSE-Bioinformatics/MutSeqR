# Write results to Excel tables

Writes data to an Excel file.

## Usage

``` r
write_excel(data, output_path = NULL, workbook_name, model_results = FALSE)
```

## Arguments

- data:

  A data frame, a list of data frames, or model_mf output.

- output_path:

  Directory to write to. Defaults to current working directory.

- workbook_name:

  Filename (without extension).

- model_results:

  Logical. Set to TRUE if data is output from model_mf.

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
mf_example3 <- readRDS(
 system.file("extdata/Example_files/mf_data_6_sample.rds",
     package = "MutSeqR"
 )
)
list <- list(mf_example, mf_example2, mf_example3)
names(list) <- c("Global MF", "Base 6 Spectra", "Base 6 Sample Spectra")

# save a single data frame to an Excel file
write_excel(
 mf_example,
 output_path = outputpath,
 workbook_name = "test_single"
)
#> Saved: /tmp/RtmpTC4ccn/test_single.xlsx

# save a list of data frames to an Excel file
write_excel(list, output_path = outputpath, workbook_name = "test_list")
#> Saved: /tmp/RtmpTC4ccn/test_list.xlsx
```
