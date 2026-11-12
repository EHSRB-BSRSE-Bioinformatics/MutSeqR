# Estimate the total depth of each record

Estimates `total_depth` (the read depth excluding N-calls) from the
depth columns available: `depth` - `no_calls`, then the sum of AD, then
DP. A record whose AD lacks the reference depth, or the alternate depth
of a variant, takes DP instead of the sum of what AD has, since that sum
misses reads. Callers such as VarDict leave the reference depth of MNV
and complex calls missing.

## Usage

``` r
estimate_total_depth(data)
```

## Arguments

- data:

  The mutation data, with `variation_type` and any of `depth`,
  `no_calls` and `AD_` columns.

## Value

A numeric vector of total depths, or `NULL` when no depth column is
available.
