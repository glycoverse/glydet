# Get All Derived Traits

This function returns a named list of all derived traits. Compared to
[`basic_traits()`](https://glycoverse.github.io/glydet/reference/basic_traits.md),
this function includes derived traits with more detailed `within`
conditions.

## Usage

``` r
all_traits()
```

## Value

A named list of derived traits.

## Details

The explanations of the derived traits are as follows:

- `A1F`: Proportion of fucosylated glycans within mono-antennary glycans

- `A2F`: Proportion of fucosylated glycans within bi-antennary glycans

- `A3F`: Proportion of fucosylated glycans within tri-antennary glycans

- `A4F`: Proportion of fucosylated glycans within tetra-antennary
  glycans

- `A1Fc`: Proportion of core-fucosylated glycans within mono-antennary
  glycans

- `A2Fc`: Proportion of core-fucosylated glycans within bi-antennary
  glycans

- `A3Fc`: Proportion of core-fucosylated glycans within tri-antennary
  glycans

- `A4Fc`: Proportion of core-fucosylated glycans within tetra-antennary
  glycans

- `A1Fa`: Proportion of arm-fucosylated glycans within mono-antennary
  glycans

- `A2Fa`: Proportion of arm-fucosylated glycans within bi-antennary
  glycans

- `A3Fa`: Proportion of arm-fucosylated glycans within tri-antennary
  glycans

- `A4Fa`: Proportion of arm-fucosylated glycans within tetra-antennary
  glycans

- `A1SFa`: Proportion of arm-fucosylated glycans within sialylated
  mono-antennary glycans

- `A2SFa`: Proportion of arm-fucosylated glycans within sialylated
  bi-antennary glycans

- `A3SFa`: Proportion of arm-fucosylated glycans within sialylated
  tri-antennary glycans

- `A4SFa`: Proportion of arm-fucosylated glycans within sialylated
  tetra-antennary glycans

- `A1S0Fa`: Proportion of arm-fucosylated glycans within asialylated
  mono-antennary glycans

- `A2S0Fa`: Proportion of arm-fucosylated glycans within asialylated
  bi-antennary glycans

- `A3S0Fa`: Proportion of arm-fucosylated glycans within asialylated
  tri-antennary glycans

- `A4S0Fa`: Proportion of arm-fucosylated glycans within asialylated
  tetra-antennary glycans

- `A1B`: Proportion of bisecting glycans within mono-antennary glycans

- `A2B`: Proportion of bisecting glycans within bi-antennary glycans

- `A3B`: Proportion of bisecting glycans within tri-antennary glycans

- `A4B`: Proportion of bisecting glycans within tetra-antennary glycans

- `A1FcB`: Proportion of bisecting glycans within core-fucosylated
  mono-antennary glycans

- `A2FcB`: Proportion of bisecting glycans within core-fucosylated
  bi-antennary glycans

- `A3FcB`: Proportion of bisecting glycans within core-fucosylated
  tri-antennary glycans

- `A4FcB`: Proportion of bisecting glycans within core-fucosylated
  tetra-antennary glycans

- `A1Fc0B`: Proportion of bisecting glycans within a-core-fucosylated
  mono-antennary glycans

- `A2Fc0B`: Proportion of bisecting glycans within a-core-fucosylated
  bi-antennary glycans

- `A3Fc0B`: Proportion of bisecting glycans within a-core-fucosylated
  tri-antennary glycans

- `A4Fc0B`: Proportion of bisecting glycans within a-core-fucosylated
  tetra-antennary glycans

- `A1G`: Average degree of galactosylation per antenna within
  mono-antennary glycans

- `A2G`: Average degree of galactosylation per antenna within
  bi-antennary glycans

- `A3G`: Average degree of galactosylation per antenna within
  tri-antennary glycans

- `A4G`: Average degree of galactosylation per antenna within
  tetra-antennary glycans

- `A1Gt`: Average degree of terminal galactose per antenna within
  mono-antennary glycans

- `A2Gt`: Average degree of terminal galactose per antenna within
  bi-antennary glycans

- `A3Gt`: Average degree of terminal galactose per antenna within
  tri-antennary glycans

- `A4Gt`: Average degree of terminal galactose per antenna within
  tetra-antennary glycans

- `A1S`: Average degree of sialylation per antenna within mono-antennary
  glycans

- `A2S`: Average degree of sialylation per antenna within bi-antennary
  glycans

- `A3S`: Average degree of sialylation per antenna within tri-antennary
  glycans

- `A4S`: Average degree of sialylation per antenna within
  tetra-antennary glycans

- `A1GS`: Average degree of sialylation per galactose within
  mono-antennary glycans

- `A2GS`: Average degree of sialylation per galactose within
  bi-antennary glycans

- `A3GS`: Average degree of sialylation per galactose within
  tri-antennary glycans

- `A4GS`: Average degree of sialylation per galactose within
  tetra-antennary glycans

## Examples

``` r
all_traits()
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
#> $A1F
#> prop(nF > 0, within = (nA == 1), na_action = "keep")
#> 
#> $A2F
#> prop(nF > 0, within = (nA == 2), na_action = "keep")
#> 
#> $A3F
#> prop(nF > 0, within = (nA == 3), na_action = "keep")
#> 
#> $A4F
#> prop(nF > 0, within = (nA == 4), na_action = "keep")
#> 
#> $A1Fc
#> prop(nFc > 0, within = (nA == 1), na_action = "keep")
#> 
#> $A2Fc
#> prop(nFc > 0, within = (nA == 2), na_action = "keep")
#> 
#> $A3Fc
#> prop(nFc > 0, within = (nA == 3), na_action = "keep")
#> 
#> $A4Fc
#> prop(nFc > 0, within = (nA == 4), na_action = "keep")
#> 
#> $A1Fa
#> prop(nFa > 0, within = (nA == 1), na_action = "keep")
#> 
#> $A2Fa
#> prop(nFa > 0, within = (nA == 2), na_action = "keep")
#> 
#> $A3Fa
#> prop(nFa > 0, within = (nA == 3), na_action = "keep")
#> 
#> $A4Fa
#> prop(nFa > 0, within = (nA == 4), na_action = "keep")
#> 
#> $A1SFa
#> prop(nFa > 0, within = (nA == 1 & nS > 0), na_action = "keep")
#> 
#> $A2SFa
#> prop(nFa > 0, within = (nA == 2 & nS > 0), na_action = "keep")
#> 
#> $A3SFa
#> prop(nFa > 0, within = (nA == 3 & nS > 0), na_action = "keep")
#> 
#> $A4SFa
#> prop(nFa > 0, within = (nA == 4 & nS > 0), na_action = "keep")
#> 
#> $A1S0Fa
#> prop(nFa > 0, within = (nA == 1 & nS == 0), na_action = "keep")
#> 
#> $A2S0Fa
#> prop(nFa > 0, within = (nA == 2 & nS == 0), na_action = "keep")
#> 
#> $A3S0Fa
#> prop(nFa > 0, within = (nA == 3 & nS == 0), na_action = "keep")
#> 
#> $A4S0Fa
#> prop(nFa > 0, within = (nA == 4 & nS == 0), na_action = "keep")
#> 
#> $A1B
#> prop(B, within = (nA == 1), na_action = "keep")
#> 
#> $A2B
#> prop(B, within = (nA == 2), na_action = "keep")
#> 
#> $A3B
#> prop(B, within = (nA == 3), na_action = "keep")
#> 
#> $A4B
#> prop(B, within = (nA == 4), na_action = "keep")
#> 
#> $A1FcB
#> prop(B, within = (nA == 1 & nFc > 0), na_action = "keep")
#> 
#> $A2FcB
#> prop(B, within = (nA == 2 & nFc > 0), na_action = "keep")
#> 
#> $A3FcB
#> prop(B, within = (nA == 3 & nFc > 0), na_action = "keep")
#> 
#> $A4FcB
#> prop(B, within = (nA == 4 & nFc > 0), na_action = "keep")
#> 
#> $A1Fc0B
#> prop(B, within = (nA == 1 & nFc == 0), na_action = "keep")
#> 
#> $A2Fc0B
#> prop(B, within = (nA == 2 & nFc == 0), na_action = "keep")
#> 
#> $A3Fc0B
#> prop(B, within = (nA == 3 & nFc == 0), na_action = "keep")
#> 
#> $A4Fc0B
#> prop(B, within = (nA == 4 & nFc == 0), na_action = "keep")
#> 
#> $A1G
#> wmean(nG/nA, within = (nA == 1), na_action = "keep")
#> 
#> $A2G
#> wmean(nG/nA, within = (nA == 2), na_action = "keep")
#> 
#> $A3G
#> wmean(nG/nA, within = (nA == 3), na_action = "keep")
#> 
#> $A4G
#> wmean(nG/nA, within = (nA == 4), na_action = "keep")
#> 
#> $A1Gt
#> wmean(nGt/nA, within = (nA == 1), na_action = "keep")
#> 
#> $A2Gt
#> wmean(nGt/nA, within = (nA == 2), na_action = "keep")
#> 
#> $A3Gt
#> wmean(nGt/nA, within = (nA == 3), na_action = "keep")
#> 
#> $A4Gt
#> wmean(nGt/nA, within = (nA == 4), na_action = "keep")
#> 
#> $A1S
#> wmean(nS/nA, within = (nA == 1), na_action = "keep")
#> 
#> $A2S
#> wmean(nS/nA, within = (nA == 2), na_action = "keep")
#> 
#> $A3S
#> wmean(nS/nA, within = (nA == 3), na_action = "keep")
#> 
#> $A4S
#> wmean(nS/nA, within = (nA == 4), na_action = "keep")
#> 
#> $A1GS
#> wmean(nS/nG, within = (nA == 1), na_action = "keep")
#> 
#> $A2GS
#> wmean(nS/nG, within = (nA == 2), na_action = "keep")
#> 
#> $A3GS
#> wmean(nS/nG, within = (nA == 3), na_action = "keep")
#> 
#> $A4GS
#> wmean(nS/nG, within = (nA == 4), na_action = "keep")
#> 
```
