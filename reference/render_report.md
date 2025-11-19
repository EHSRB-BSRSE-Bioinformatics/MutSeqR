# Read configuration file and render R Markdown document

This function reads a configuration file in YAML format, extracts the
parameters, and renders an R Markdown document using the specified
parameters.

## Usage

``` r
render_report(
  config_filepath,
  output_file = "./MutSeqR_Summary_Report.html",
  output_format = "html_document"
)
```

## Arguments

- config_filepath:

  The path to the configuration file.

- output_file:

  The name of the output file. Will be saved to the outputdir in config
  params.

- output_format:

  The format of the output file. Options are "html_document" (default),
  "pdf_document", or "all".

## Value

A rendered R Markdown document.

## Examples

``` r
# Step 1: Copy the example configuration file to your working directory
## config <- system.file("extdata", "inputs", "summary_config.yaml", package = "MutSeqR")
## file.copy(from = config, to = "your/working/directory/summary_config.yaml")
# Step 2: Edit the configuration file with your inputs
# Step 3: Render the report
## render_report(config_filepath = "your/working/directory/summary_config.yaml",
##  output_file = "MutSeqR_Summary_Report.html",
##  output_format = "html_document")
```
