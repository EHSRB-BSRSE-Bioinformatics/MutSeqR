# Hierarchical Clustering

perform hierarchical clustering of samples based on the mutation
spectra.

## Usage

``` r
cluster_spectra(
  mf_data = mf_data,
  group_col = "sample",
  response_col = "proportion_min",
  subtype_col = "normalized_subtype",
  dist = "cosine",
  cluster_method = "ward.D"
)
```

## Arguments

- mf_data:

  A data frame containing the mutation data. This data must include a
  column containing the mutation subtypes, a column containing the
  sample/cohort names, and a column containing the response variable.

- group_col:

  The name of the column in data that contains the sample/cohort names.

- response_col:

  The name of the column in data that contains the response variable.
  Typical response variables can be the subtype mf, proportion, or
  count.

- subtype_col:

  The name of the column in data that contains the mutation subtypes.

- dist:

  the distance measure to be used. This must be one of "cosine",
  "euclidean", "maximum", "manhattan","canberra", "binary" or
  "minkowski". See [dist](https://rdrr.io/r/stats/dist.html) for
  details.

- cluster_method:

  The agglomeration method to be used. See
  [hclust](https://rdrr.io/r/stats/hclust.html) for details.

## Value

A dendrogram object representing the hierarchical clustering of the
samples.

## Details

The cosine distance measure represents the inverted cosine similarity
between samples:

\\\text{Cosine Dissimilarity} = 1 - \frac{\mathbf{A} \cdot
\mathbf{B}}{\\ \mathbf{A} \\ \cdot \\ \mathbf{B} \\}\\

This equation calculates the cosine dissimilarity between two vectors A
and B.

Leaves are sorted using dendsort, if installed, otherwise leaves are
unsorted.
