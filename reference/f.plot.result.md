# Plot the PROAST results

Independently generate the model plots from the raw results.

## Usage

``` r
f.plot.result(
  proast_results_list,
  output_path = NULL,
  output_type = "svg",
  prefix = NULL,
  model_averaging = FALSE
)
```

## Arguments

- proast_results_list:

  The raw results list. This is the output of
  [`f.proast`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/f.proast.md)

- output_path:

  The file path to the output directory. If the output_path is NULL, it
  will save it to the working directory. If the output_path doesn't
  exist, it will be created.

- output_type:

  The file type to export the plots. Options are 'svg', 'jpeg', 'pdf',
  'png', 'tiff', or 'none'. If "none", the plots will be displayed to
  the graphics window, recorded with recordPlot(), and returned as a
  list.

- prefix:

  A custom prefix to append to the file names. Default is "PROAST\_".

- model_averaging:

  A logical variable indicating whether you want to generate the model
  averaging figure (TRUE) or the plots of the individual models (FALSE).
  You plot one or the other, not both. Plotting the model averaging
  figure will require the function to re-run the bootstrapping so it
  might take a while. You may think this seems rather inefficient. Well,
  it is, but I'm too tired to fix it, so we all just have to deal with
  it for now.

## Value

Generates plots. Either saves them to an output path or records them and
returns them as a list.
