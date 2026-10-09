# Prepare depth tables from a single resolution-specific input

Reads and validates a precalculated depth data frame or delimited file,
detects its subtype resolution from the denominator context column, and
derives all lower-resolution depth tables that can be recovered from
that input.

## Usage

``` r
write_depth_data(depth_data, d_sep = "\t", group_cols = "sample")
```

## Arguments

- depth_data:

  A data frame or path to a delimited file. It must contain the columns
  named in `group_cols`, `subtype_depth`, and exactly one context column
  from `denominator_dict` for `base_6`, `base_12`, `base_96`, or
  `base_192`.

- d_sep:

  A single-character delimiter used to read `depth_data` when it is a
  file path.

- group_cols:

  Character vector of columns that identify independent groups in
  `depth_data`. Defaults to `"sample"`; use other columns when depth
  data is already aggregated to those groups.

## Value

A named list of the depth tables at available subtype resolutions:
"base192", "base96", "base12", "base6", "global" (none). Each subtype
table contains the `group_cols`, its `denominator_dict` context column,
`subtype_depth`, and `group_depth`; the global table contains
`group_cols` and `group_depth`.

## Details

The input must contain the grouping columns, `subtype_depth`, and the
appropriate context column. Use
[`MutSeqR::denominator_dict`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/denominator_dict.md)
to see the context column name associated with each subtype resolution.
Each group must have exactly one row for every context in
[`MutSeqR::context_list`](https://ehsrb-bsrse-bioinformatics.github.io/MutSeqR/reference/context_list.md)
at the given resolution; contexts with zero depth must be included
explicitly. Duplicate group-context pairs, unknown contexts, missing
contexts, or non-finite per-group depth totals produce an error.

*Base192* input produces Base96, Base12, Base6, and global depth tables.
*Base96* input produces Base6 and global tables, but cannot produce
Base12 because its normalized contexts have already combined
reverse-complement strands. *Base12* input produces Base6 and global
tables. *Base6* input produces global depth.
