# Pipe operator

See `magrittr::%>%` for details.

## Usage

``` r
lhs %>% rhs
```

## Arguments

- lhs:

  A value or the magrittr placeholder.

- rhs:

  A function call using the magrittr semantics.

## Value

The result of calling `rhs(lhs)`.

## Examples

``` r
df <- data.frame(x = 1:5, y = rnorm(5))
df %>% dplyr::mutate(z = x + y)
#>   x           y         z
#> 1 1 -0.05260191 0.9473981
#> 2 2  0.54299634 2.5429963
#> 3 3 -0.91407483 2.0859252
#> 4 4  0.46815442 4.4681544
#> 5 5  0.36295126 5.3629513
df %>% head(3) %>% summary()
#>        x             y          
#>  Min.   :1.0   Min.   :-0.9141  
#>  1st Qu.:1.5   1st Qu.:-0.4833  
#>  Median :2.0   Median :-0.0526  
#>  Mean   :2.0   Mean   :-0.1412  
#>  3rd Qu.:2.5   3rd Qu.: 0.2452  
#>  Max.   :3.0   Max.   : 0.5430  
```
