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
#>   x          y        z
#> 1 1  0.3833950 1.383395
#> 2 2 -0.2161434 1.783857
#> 3 3  0.2922274 3.292227
#> 4 4 -1.0101391 2.989861
#> 5 5 -0.5412957 4.458704
df %>% head(3) %>% summary()
#>        x             y           
#>  Min.   :1.0   Min.   :-0.21614  
#>  1st Qu.:1.5   1st Qu.: 0.03804  
#>  Median :2.0   Median : 0.29223  
#>  Mean   :2.0   Mean   : 0.15316  
#>  3rd Qu.:2.5   3rd Qu.: 0.33781  
#>  Max.   :3.0   Max.   : 0.38339  
```
