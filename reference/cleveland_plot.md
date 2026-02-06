# Cleveland Plot

Make a Cleveland plot for the PROAST results.

## Usage

``` r
cleveland_plot(results, covariate_col = NULL, output_path = NULL)
```

## Arguments

- results:

  PROAST results object.

- covariate_col:

  Covariate column name.

- output_path:

  Output path for the plot. If the output_path doesn't exist, it will be
  created. If NULL, the plots will not be exported.

## Value

A single ggplot object with facets for each response.
