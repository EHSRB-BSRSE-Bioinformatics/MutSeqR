# Calculate mutation frequency

Calculates mutation frequencies for arbitrary groupings and creates a
new dataframe with the results. Mutation frequency is calculated by
dividing the sum of mutations by the sum of the total_depth for a given
group (mutations/bp). The operation is run using both the minimum and
maximum independent mutation counting methods.

## Usage

``` r
calculate_mf(
  mutation_data,
  cols_to_group = "sample",
  subtype_resolution = "none",
  variant_types = c("snv", "deletion", "insertion", "complex", "mnv", "sv", "ambiguous",
    "uncategorized"),
  calculate_depth = TRUE,
  correct_depth = TRUE,
  correct_depth_by_indel_priority = FALSE,
  precalc_depth_data = NULL,
  d_sep = "\t",
  summary = TRUE,
  retain_metadata_cols = NULL
)
```

## Arguments

- mutation_data:

  The data frame (or GRanges) to be processed containing mutation data.
  Required columns are listed in details.

- cols_to_group:

  A vector of grouping variables. This should be the groups of interest
  that you want to calculate a frequency for. For instance, getting the
  frequency by `"sample"`. Other options might include an experimental
  group Ex. `"dose"` or a locus Ex. `c("sample", "locus")`. All listed
  variables must be a column in the mutation_data. Do not include
  mutation subtype columns in this field. Please refer to
  subtype_resolution to group by subtype as the calculation will differ.

- subtype_resolution:

  The degree at which to resolve the mutation subtypes when calculating
  frequencies. Mutation frequency will be calculated across all
  `col_to_groups` for each mutation subtype given the desired
  resolution. Subtype proportions will also be calculated. Options are
  "none", "type", "base_6", "base_12", "base_96", and "base_192". See
  details for definitions.

- variant_types:

  Use this parameter to choose which variation types to include in the
  mutation counts. Provide a character vector of the variation types
  that you want to include. Alternatively, provide a character vector of
  the variation types that you want to exclude preceded by "-". Options
  are: "snv", "complex", "deletion", "insertion", "mnv", "sv",
  "ambiguous", "uncategorized". Ex. inclusion: "snv", exclusion: "-snv".
  Default includes all variants. For `calculate_depth = TRUE`:
  Regardless of whether or not a variant is included in the mutation
  counts, the total_depth for that position will be counted.

- calculate_depth:

  A logical variable, whether to calculate the per-group total_depth
  from the `mutation_data`. If set to TRUE, the mutation data must
  contain a total_depth value for every sequenced base (including
  variants AND no-variant calls). If set to FALSE, pre-calculated
  per-group `total_depth` values may be supplied at the desired
  subtype_resolution using the precalc_depth_data parameter.
  Alternatively, if no per-group `total_depth` is available, per-group
  mutation counts will be calculated, but mutation frequency will not.
  In such cases, mutation subtype proportions will not be normalized to
  the `total_depth`.

- correct_depth:

  A logical value. If TRUE, the function will correct the `total_depth`
  column in `mutation_data` in order to prevent double-counting the
  `total_depth` values for the same genomic position. For rows with the
  same sample, contig, and start values, the `total_depth` will be
  retained for only one row. All other rows in the group will have their
  `total_depth` set to 0. The default is TRUE.

- correct_depth_by_indel_priority:

  A logical value. If TRUE, during depth correction, should there be
  different `total_depth` values within a group of rows with the same
  sample, contig, and start values, the `total_depth` value for the row
  with the highest priority `variation_type` will be retained, while the
  other rows will have their `total_depth` set to 0. `variation_type`
  priority order is: deletion, complex, insertion, snv, mnv, sv,
  uncategorised, ambiguous, no_variant. If FALSE, the `total_depth`
  value for the first row in the group will be retained, while the other
  rows will have their `total_depth` set to 0. The default is FALSE.

- precalc_depth_data:

  A data frame or a file path to a text file containing pre-calculated
  depth values. See details for required columns.

- d_sep:

  The delimiter used in the precalc_depth_data, if applicable. Default
  is tab-delimited.

- summary:

  A logical variable, whether to return a summary table (i.e., where
  only relevant columns for frequencies and groupings are returned).
  Setting this to false returns all columns in the original
  mutation_data, which might make plotting more difficult, but may
  provide additional flexibility to power users. Default is TRUE.

- retain_metadata_cols:

  a character vector that contains the names of the metadata columns
  that you would like to retain in the summary table. This may be useful
  for plotting your summary data. Ex. retain the "dose" column when
  summarising by "sample".

## Value

A data frame with the mutation frequency calculated. If summary is set
to TRUE, the data frame will be a summary table with the mutation
frequency calculated for each group. If summary is set to FALSE, the
mutation frequency will be appended to each row of the original
mutation_data.

- `sum_min`: The sum of all mutations within the group, calculated using
  the "min" method for mutation counting. All identical mutations within
  a sample are assumed to be the result of clonal expansion and are only
  counted once.

- `sum_max`: The sum of all mutations within the group, calculated using
  the "max" method for mutaiton counting. All identical mutations within
  a sample are assumed to be idenpendant mutational events and are
  included in the mutation frequency calculation.

- `group_depth`: The total_depth summed across groups.

- `subtype_depth`: The total_depth summed across groups for a given
  sequence context. Used for calculating snv subtype frequencies.

- `mf_min`: The mutation frequency calculated using the "min" method for
  mutation counting. mf_min = sum_min / depth.

- `mf_max`: The mutation frequency calculated using the "max" method for
  mutation counting. mf_max = sum_max / depth.

- `proportion_min`: The proportion of each mutation subtype within the
  group, normalized to the depth. Calculated using the "min" method.
  This is only calculated if `subtype_resolution` is not "none". If no
  depth is calculated or provided, proportion is calculated without
  normalization to the depth.

- `proportion_max`: The proportion of each mutation subtype within the
  group, normalized to its read depth. Calculated using the "max"
  method. This is only calculated if `subtype_resolution` is not "none".
  If no depth is calculated or provided, proportion is calculated
  without normalization to the depth.

## Details

**Required columns:**

- `contig`: (or `seqnames`) The reference sequence name.

- `start`: 1-based start position of the feature.

- `alt_depth`: The read depth supporting the alternate allele.

- `variation_type`: The category to which this variant is assigned.

- subtype_col: The column containing the mutation subtype. This column
  depends on the `subtype_resolution` parameter.

- reference context: The column containing the referene base(s) for the
  mutation. This column depends on the `subtype_resolution` parameter.

- cols to group: all columns across which you want to calculate the
  mutation frequency. Ex. `c("tissue", "dose")`. These columns should be
  listed in cols_to_group.

It is also required to include the total_depth column if you are
calculating depth from the mutation data. If you are using precalculated
depth data, the total_depth column is not required.

MF is total mutation count divided by total depth. If grouping by an
experimental variable note that MF is NOT the arithmetic mean of sample
MFs.

**Subtype Resolutions:**

- "none" calculates mutation frequencies across all selected grouping
  columns.

- "type" calculates mutation frequencies across all selected grouping
  columns for each `variation_type` seperately; snv, mnv, deletion,
  insertion, complex, sv, ambiguous, uncategorized.

- "base_6" calculates mutation frequencies across all selected grouping
  columns for each variation_type with snv mutations separated by
  `normalized_subtype`; C\>A, C\>G, C\>T, T\>A, T\>C, T\>G. The
  reference context is `normalized_ref`.

- "base_12" calculates mutation frequencies across all selected grouping
  columns for each variation_type with snv mutations separated by
  `subtype`; A\>C, A\>G, A\>T, C\>A, C\>G, C\>T, G\>A, G\>C, G\>T, T\>A,
  T\>C, T\>G. The reference context is `short_ref`.

- "base_96" calculates mutation frequencies across all selected grouping
  columns for each variation_type with snv mutations separated by
  `normalized_context_with_mutation`, i.e. the 96-base trinucleotide
  context. Ex. A\[C\>T\]A. The reference context is
  `normalized_context`.

- "base_192" calculates mutation frequencies across all selected
  grouping columns for each variation_type with snv mutations separated
  by `context_with_mutation`, i.e. the 192-base trinucleotide context.
  Ex A\[G\>A\]A. The reference context is `context`.

**Subtype depth:** For SNV subtypes, the total_depth is summed based on
the sequence context in which the SNV subtype occurs. Ex. for base_6,
the two possible reference bases are C or T; hence, the total_depth is
summed seperately for C:G positions and T:A positions. The MF for C\>T
mutations is calculated as the total number of C\>T mutations /
total_depth for C and G positions (sum / subtype_depth). Non-SNV
mutation subtypes will be caluclated as their sum / group_depth, since
they can occur in the context of any nucleotide.

**retain_metadata_cols at subtype_resolution:** The summary table uses a
pre-defined list of possible subtypes for each resolution. If a
particular subtype within a given group is not recorded in the mutation
data, the summary table will have no frame of reference for populating
the metadata_cols. Thus, for subtypes that do not occur in the mutation
data for a given group, the corresponding metadata_col will be NA.

**Precalculated depth:** The precalculated depth data should include
depth calculated per your grouping variable, resolved to the nucleotide
context that supports the `subtype_resolution`. Required columns are:

- `cols_to_group`: (or `seqnames`) The grouping column(s).

- a context column: The nucleotide context, depending on the subtype
  resolution. Must be one column from `denominator_dict`. This is not
  required for subtype resolutions "none" or "type".

- `subtype_depth` or `group_depth`: The total_depth summed across groups
  (group_depth) and groups & contexts (subtype_depth). group_depth is
  used for subtype_resolution = "none" and "type", while subtype_depth
  is used for all SNV resolutions.

Generally, depth should be precalculated to the grouping variables in
`cols_to_group`. However, the function is capable of aggregating depth
supplied per sample across a lower resolution grouping variable. The
grouping variable is mapped to sample from the mutation data. This
mapping requires each sample to identify exactly one group value.
Mapping to higher- resolution grouping variables is not possible. For
example, you may calculate mutation frequencies per "dose", and supply
precalc_depth_data as group_depth per sample, assuming that each sample
belongs to exactly one dose. However, you cannot calculate mutation
frequencies per contig using this precalc_depth_data, since each sample
maps to more than one contig.

precalc_depth_data may be supplied at a higher context resolution than
is required by subtype_resolution. The function will automatically
aggregate depth across possible contexts to support lower-level
resolutions. As such, precalc_depth_data with subtype_depth resolved to
the 64 trinucleotide context (base_192; "context") can be used to
calculate mutation frequencies at any given subtype_resolution. You can
access a list of context values for each subtype resolution using
`MutSeqR::context_list$your_subtype_resolution`.

**Variant filtering:** Variants flagged as TRUE in the `filter_mut`
column will be excluded from the mutation counts. However, the
total_depth of these variants will be included in the group/subtype
depths if calculating depth.

**Depth correction** is important for preventing double-counting of
reads in mutation data when summing the total_depth across samples or
other groups. Generally, when several mutations have been detected at
the same genomic position, within a sample, the total_depth value will
be the same for all of them. However, in some datasets, whenever a
deletion is detected, the data may contain an additional row with the
same genomic position calling a "no_variant". The total_depth will
differ between the deletion and the no_variant. In these cases,
correct_depth_by_indel_priority == TRUE will ensure that the total_depth
value for the deletion is retained, while the total_depth value for the
no_variant is removed.

## Examples

``` r
# Mutation data is just for example purposes. It does not reflect real data
mutation_data <- readRDS(system.file("extdata", "Example_files",
                                     "filtered_simple_mutation_data.rds",
                                     package = "MutSeqR"))
# Calculate mutation frequency by sample.
# Calculate depth from the mutation data (default)
# Correct the Depth (default) with indel priority (set)
mf_example <- calculate_mf(
  mutation_data = mutation_data,
  cols_to_group = "sample",
  correct_depth_by_indel_priority = TRUE
)
#> Performing internal depth correction to prevent double-counting
#> Internal depth correction complete.
#> Joining with `by = join_by(sample)`
```
