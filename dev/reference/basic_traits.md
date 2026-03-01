# Get Basic Derived Traits

These derived traits are the most basic and commonly used derived
traits. They describe global properties of a glycome including the type
of glycans, fucosylation level, sialylation level, galactosylation
level, and branching level.

## Usage

``` r
basic_traits()
```

## Value

A named list of derived traits.

## Details

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

## Examples

``` r
basic_traits()
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
