# Get Traits in Clerc et al. 2018

These traits are the ones used by Clerc et al. 2018
(https://doi.org/10.1053/j.gastro.2018.05.030). We generally don't
recommend using these traits because they are either redundant or
missing some important traits. We include this function because this
paper is very influential in the field.

## Usage

``` r
traits_clerc_2018(sia_link = FALSE)
```

## Arguments

- sia_link:

  A boolean indicating whether to include sialic acid linkage traits.
  Default is `FALSE`.

## Value

A named list of derived traits.

## Usage of sialic acid linkage traits

To use these sialic acid linkage traits,
[`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
of the input `GlycomicSE` or `GlycoproteomicSE` must have the following
columns:

- `nE`: Number of a2,6-linked sialic acids

- `nL`: Number of a2,3-linked sialic acids

Note that you have to add these two columns even if the
`glycan_structure` column has intact linkages. This is because by
convention all traits work with glycan structures with "basic" structure
levels (i.e., with generic monosaccharides like "Hex" and "HexNAc" and
no linkages specified).

## Examples

``` r
traits_clerc_2018()[1:5]
#> $CA1
#> prop(nA == 1, na_action = "keep")
#> 
#> $CA2
#> prop(nA == 2, na_action = "keep")
#> 
#> $CA3
#> prop(nA == 3, na_action = "keep")
#> 
#> $CA4
#> prop(nA == 4, na_action = "keep")
#> 
#> $TC
#> prop(Tp == "complex", na_action = "keep")
#> 
```
