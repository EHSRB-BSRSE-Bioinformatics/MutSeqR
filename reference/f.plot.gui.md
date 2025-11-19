# Manages plotting for PROAST

Runs through the plotting functions depending on data type and plot
type.

## Usage

``` r
f.plot.gui(
  ans.all,
  HTML = FALSE,
  model.summ = TRUE,
  display_plots = TRUE,
  .proast_env = NULL,
  output_type = NULL,
  filename = NULL,
  interactive_mode = TRUE,
  knitting = FALSE,
  return_plots = FALSE
)
```

## Arguments

- ans.all:

  The proast object that gets passed to all functions.

- HTML:

  Keep FALSE

- model.summ:

  Keep TRUE

- display_plots:

  A logical variable - whether we want to display the plots or not.

- .proast_env:

  environment

- output_type:

  The format that you wish to export the plots as.

- filename:

  The name of the file to be read.

- interactive_mode:

  A TRUE/FALSE value specifying whether you want to run interactively
  (i.e., TRUE, the default) or using command-line mode (i.e., FALSE,
  non-interactive). If FALSE, you must provide all other parameters.

- knitting:

  A TRUE/FALSE value specifying whether you are knitting a document. If
  TRUE, the function will adjust the plot display accordingly.

- return_plots:

  A logical variable indicating whether you want to return the plots as
  a list (TRUE) or not (FALSE, the default). If TRUE, the function will
  return a list of recorded plots.

## Value

Either the proast object (default) or a list of recorded plots if
return_plots = TRUE.
