# Prepare precalculated depth tables for the summary report

Validates context-resolved precalculated depth inputs and derives the
global, base_6, and base_96 tables used by the summary report. Explicit
inputs at a requested resolution take precedence over derived tables.
When an input is omitted, the finest available context-resolved table is
used to derive it.

## Usage

``` r
prepare_report_depth_data(
  base192 = NULL,
  base96 = NULL,
  global = NULL,
  base6 = NULL
)
```

## Arguments

- base192:

  Optional base_192 depth data. Must be a data frame accepted by
  [`write_depth_data()`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/write_depth_data.md).

- base96:

  Optional base_96 depth data. Must be a data frame accepted by
  [`write_depth_data()`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/write_depth_data.md).

- global:

  Optional global depth data with `sample` and `group_depth` columns.

- base6:

  Optional base_6 depth data. Must be a data frame accepted by
  [`write_depth_data()`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/write_depth_data.md).

## Value

A list with `global`, `base6`, `base96`, and `base192` elements. The
context-resolved tables include `group_depth` in addition to their
context and `subtype_depth` columns. Elements are `NULL` when neither an
explicit input nor a compatible higher-resolution input is available.
