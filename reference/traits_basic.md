# Get Basic Derived Traits

These derived traits are the most basic and commonly used derived
traits. They describe global properties of a glycome including the type
of glycans, fucosylation level, sialylation level, galactosylation
level, and branching level.

## Usage

``` r
traits_basic(sia_link = FALSE)

basic_traits(sia_link = FALSE)
```

## Arguments

- sia_link:

  A boolean indicating whether to include sialic acid linkage traits.
  Default is `FALSE`.

## Value

A named list of derived traits.

## Descriptions of traits

The explanations of the derived traits are as follows:

- `TM`: Proportion of highmannose glycans

- `TH`: Proportion of hybrid glycans

- `TC`: Proportion of complex glycans

- `MM`: Average number of mannoses within highmannose glycans

- `CA2`: Proportion of bi-antennary glycans within complex glycans

- `CA3`: Proportion of tri-antennary glycans within complex glycans

- `CA4`: Proportion of tetra-antennary glycans within complex glycans

- `TF`: Proportion of fucosylated glycans

- `TFc`: Proportion of core-fucosylated glycans

- `TFa`: Proportion of arm-fucosylated glycans

- `TB`: Proportion of glycans with bisecting GlcNAc

- `GS`: Average degree of sialylation per galactose

- `AG`: Average degree of galactosylation per antenna

- `TS`: Proportion of sialylated glycans

Four additional sialic acid linkage traits are included if
`sia_link = TRUE`.

- `GE`: Average degree of a2,6-linked sialylation per galactose

- `GL`: Average degree of a2,3-linked sialylation per galactose

- `TE`: Proportion of a2,6-linked sialylated glycans

- `TL`: Proportion of a2,3-linked sialylated glycans

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
traits_basic()
#> $TM
#> prop(Tp == "highmannose", na_action = "keep")
#> 
#> $TH
#> prop(Tp == "hybrid", na_action = "keep")
#> 
#> $TC
#> prop(Tp == "complex", na_action = "keep")
#> 
#> $MM
#> wmean(nM, within = (Tp == "highmannose"), na_action = "keep")
#> 
#> $CA2
#> prop(nA == 2, within = (Tp == "complex"), na_action = "keep")
#> 
#> $CA3
#> prop(nA == 3, within = (Tp == "complex"), na_action = "keep")
#> 
#> $CA4
#> prop(nA == 4, within = (Tp == "complex"), na_action = "keep")
#> 
#> $TF
#> prop(nF > 0, na_action = "keep")
#> 
#> $TFc
#> prop(nFc > 0, na_action = "keep")
#> 
#> $TFa
#> prop(nFa > 0, na_action = "keep")
#> 
#> $TB
#> prop(B, na_action = "keep")
#> 
#> $GS
#> wmean(nS/nG, na_action = "keep")
#> 
#> $AG
#> wmean(nG/nA, na_action = "keep")
#> 
#> $TS
#> prop(nS > 0, na_action = "keep")
#> 
```
