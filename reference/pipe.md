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
#>   x            y          z
#> 1 1 -1.400043517 -0.4000435
#> 2 2  0.255317055  2.2553171
#> 3 3 -2.437263611  0.5627364
#> 4 4 -0.005571287  3.9944287
#> 5 5  0.621552721  5.6215527
df %>% head(3) %>% summary()
#>        x             y          
#>  Min.   :1.0   Min.   :-2.4373  
#>  1st Qu.:1.5   1st Qu.:-1.9187  
#>  Median :2.0   Median :-1.4000  
#>  Mean   :2.0   Mean   :-1.1940  
#>  3rd Qu.:2.5   3rd Qu.:-0.5724  
#>  Max.   :3.0   Max.   : 0.2553  
```
