# Get Traits in Fu et al. 2026

These traits are the ones used by Fu et al. 2026
(https://doi.org/10.1038/s41467-026-68579-x). It is much like
[`traits_detailed()`](https://glycoverse.github.io/glydet/dev/reference/traits_detailed.md),
but doesn't differentiate core- and arm-fucosylation, and doesn't
include the sialic acid linkage traits. Also, many traits are too
specific to be useful and interpretable.

## Usage

``` r
traits_fu_2026()
```

## Value

A named list of derived traits.

## Examples

``` r
traits_fu_2026()[1:5]
#> $TM
#> prop(Tp == "highmannose", na_action = "keep")
#> 
#> $THy
#> prop(Tp == "hybrid", na_action = "keep")
#> 
#> $TC
#> prop(Tp == "complex", na_action = "keep")
#> 
#> $MM
#> wmean(nM, within = (Tp == "highmannose"), na_action = "keep")
#> 
#> $CA2
#> prop(nA == 2, na_action = "keep")
#> 
```
