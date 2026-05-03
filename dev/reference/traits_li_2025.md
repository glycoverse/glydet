# Get Traits in Li et al. 2025

These traits are the ones used by Li et al. 2025
(https://doi.org/10.1038/s41467-025-57633-9). These traits were created
based only on glycan compositions, not structures. Here the traits are
remapped to use glycan structures, so they are not exactly the same as
the ones in the paper. Generally we don't recommend using these traits
because they capture too little information. We include this function
because this paper is very influential in the field.

## Usage

``` r
traits_li_2025()
```

## Value

A named list of derived traits.

## Examples

``` r
traits_li_2025()[1:5]
#> $S0
#> prop(nS == 0, na_action = "keep")
#> 
#> $S
#> prop(nS > 0, na_action = "keep")
#> 
#> $S1
#> prop(nS == 1, na_action = "keep")
#> 
#> $S2
#> prop(nS == 2, na_action = "keep")
#> 
#> $S3
#> prop(nS == 3, na_action = "keep")
#> 
```
