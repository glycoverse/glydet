#' Get Basic Derived Traits
#'
#' These derived traits are the most basic and commonly used derived traits.
#' They describe global properties of a glycome including the type of glycans,
#' fucosylation level, sialylation level, galactosylation level,
#' and branching level.
#'
#' @details
#' The explanations of the derived traits are as follows:
#'
#' - `TM`: Proportion of highmannose glycans
#' - `TH`: Proportion of hybrid glycans
#' - `TC`: Proportion of complex glycans
#' - `MM`: Average number of mannoses within highmannose glycans
#' - `CA2`: Proportion of bi-antennary glycans within complex glycans
#' - `CA3`: Proportion of tri-antennary glycans within complex glycans
#' - `CA4`: Proportion of tetra-antennary glycans within complex glycans
#' - `TF`: Proportion of fucosylated glycans
#' - `TFc`: Proportion of core-fucosylated glycans
#' - `TFa`: Proportion of arm-fucosylated glycans
#' - `TB`: Proportion of glycans with bisecting GlcNAc
#' - `GS`: Average degree of sialylation per galactose
#' - `AG`: Average degree of galactosylation per antenna
#' - `TS`: Proportion of sialylated glycans
#'
#' @returns
#' A named list of derived traits.
#'
#' @examples
#' basic_traits()
#'
#' @export
basic_traits <- function() {
  list(
    # Proportion of highmannose glycans
    TM = prop(Tp == "highmannose"),
    # Proportion of hybrid glycans
    TH = prop(Tp == "hybrid"),
    # Proportion of complex glycans
    TC = prop(Tp == "complex"),
    # Average number of mannoses within highmannose glycans
    MM = wmean(nM, within = (Tp == "highmannose")),
    # Proportion of bi-antennary glycans within complex glycans
    CA2 = prop(nA == 2, within = (Tp == "complex")),
    # Proportion of tri-antennary glycans within complex glycans
    CA3 = prop(nA == 3, within = (Tp == "complex")),
    # Proportion of tetra-antennary glycans within complex glycans
    CA4 = prop(nA == 4, within = (Tp == "complex")),
    # Proportion of fucosylated glycans
    TF = prop(nF > 0),
    # Proportion of core-fucosylated glycans
    TFc = prop(nFc > 0),
    # Proportion of arm-fucosylated glycans
    TFa = prop(nFa > 0),
    # Proportion of glycans with bisecting GlcNAc
    TB = prop(B),
    # Average degree of sialylation per galactose
    GS = wmean(nS / nG),
    # Average degree of galactosylation per antenna
    AG = wmean(nG / nA),
    # Proportion of sialylated glycans
    TS = prop(nS > 0)
  )
}

#' Get All Derived Traits
#'
#' This function returns a named list of all derived traits.
#' Compared to [basic_traits()], this function includes derived traits with more
#' detailed `within` conditions.
#'
#' @details
#' The explanations of the derived traits are as follows:
#'
#' - `A1F`: Proportion of fucosylated glycans within mono-antennary glycans
#' - `A2F`: Proportion of fucosylated glycans within bi-antennary glycans
#' - `A3F`: Proportion of fucosylated glycans within tri-antennary glycans
#' - `A4F`: Proportion of fucosylated glycans within tetra-antennary glycans
#' - `A1Fc`: Proportion of core-fucosylated glycans within mono-antennary glycans
#' - `A2Fc`: Proportion of core-fucosylated glycans within bi-antennary glycans
#' - `A3Fc`: Proportion of core-fucosylated glycans within tri-antennary glycans
#' - `A4Fc`: Proportion of core-fucosylated glycans within tetra-antennary glycans
#' - `A1Fa`: Proportion of arm-fucosylated glycans within mono-antennary glycans
#' - `A2Fa`: Proportion of arm-fucosylated glycans within bi-antennary glycans
#' - `A3Fa`: Proportion of arm-fucosylated glycans within tri-antennary glycans
#' - `A4Fa`: Proportion of arm-fucosylated glycans within tetra-antennary glycans
#' - `A1SFa`: Proportion of arm-fucosylated glycans within sialylated mono-antennary glycans
#' - `A2SFa`: Proportion of arm-fucosylated glycans within sialylated bi-antennary glycans
#' - `A3SFa`: Proportion of arm-fucosylated glycans within sialylated tri-antennary glycans
#' - `A4SFa`: Proportion of arm-fucosylated glycans within sialylated tetra-antennary glycans
#' - `A1S0Fa`: Proportion of arm-fucosylated glycans within asialylated mono-antennary glycans
#' - `A2S0Fa`: Proportion of arm-fucosylated glycans within asialylated bi-antennary glycans
#' - `A3S0Fa`: Proportion of arm-fucosylated glycans within asialylated tri-antennary glycans
#' - `A4S0Fa`: Proportion of arm-fucosylated glycans within asialylated tetra-antennary glycans
#' - `A1B`: Proportion of bisecting glycans within mono-antennary glycans
#' - `A2B`: Proportion of bisecting glycans within bi-antennary glycans
#' - `A3B`: Proportion of bisecting glycans within tri-antennary glycans
#' - `A4B`: Proportion of bisecting glycans within tetra-antennary glycans
#' - `A1FcB`: Proportion of bisecting glycans within core-fucosylated mono-antennary glycans
#' - `A2FcB`: Proportion of bisecting glycans within core-fucosylated bi-antennary glycans
#' - `A3FcB`: Proportion of bisecting glycans within core-fucosylated tri-antennary glycans
#' - `A4FcB`: Proportion of bisecting glycans within core-fucosylated tetra-antennary glycans
#' - `A1Fc0B`: Proportion of bisecting glycans within a-core-fucosylated mono-antennary glycans
#' - `A2Fc0B`: Proportion of bisecting glycans within a-core-fucosylated bi-antennary glycans
#' - `A3Fc0B`: Proportion of bisecting glycans within a-core-fucosylated tri-antennary glycans
#' - `A4Fc0B`: Proportion of bisecting glycans within a-core-fucosylated tetra-antennary glycans
#' - `A1G`: Average degree of galactosylation per antenna within mono-antennary glycans
#' - `A2G`: Average degree of galactosylation per antenna within bi-antennary glycans
#' - `A3G`: Average degree of galactosylation per antenna within tri-antennary glycans
#' - `A4G`: Average degree of galactosylation per antenna within tetra-antennary glycans
#' - `A1Gt`: Average degree of terminal galactose per antenna within mono-antennary glycans
#' - `A2Gt`: Average degree of terminal galactose per antenna within bi-antennary glycans
#' - `A3Gt`: Average degree of terminal galactose per antenna within tri-antennary glycans
#' - `A4Gt`: Average degree of terminal galactose per antenna within tetra-antennary glycans
#' - `A1S`: Average degree of sialylation per antenna within mono-antennary glycans
#' - `A2S`: Average degree of sialylation per antenna within bi-antennary glycans
#' - `A3S`: Average degree of sialylation per antenna within tri-antennary glycans
#' - `A4S`: Average degree of sialylation per antenna within tetra-antennary glycans
#' - `A1GS`: Average degree of sialylation per galactose within mono-antennary glycans
#' - `A2GS`: Average degree of sialylation per galactose within bi-antennary glycans
#' - `A3GS`: Average degree of sialylation per galactose within tri-antennary glycans
#' - `A4GS`: Average degree of sialylation per galactose within tetra-antennary glycans
#'
#' @returns
#' A named list of derived traits.
#'
#' @examples
#' all_traits()
#'
#' @export
all_traits <- function() {
  additional_traits <- list(

    # ===== Advanced Fucosylation Traits =====

    # Proportion of fucosylated glycans within mono-antennary glycans
    A1F = prop(nF > 0, within = (nA == 1)),
    # Proportion of fucosylated glycans within bi-antennary glycans
    A2F = prop(nF > 0, within = (nA == 2)),
    # Proportion of fucosylated glycans within tri-antennary glycans
    A3F = prop(nF > 0, within = (nA == 3)),
    # Proportion of fucosylated glycans within tetra-antennary glycans
    A4F = prop(nF > 0, within = (nA == 4)),

    # ===== Advanced Core-Fucosylation Traits =====

    # Proportion of core-fucosylated glycans within mono-antennary glycans
    A1Fc = prop(nFc > 0, within = (nA == 1)),
    # Proportion of core-fucosylated glycans within bi-antennary glycans
    A2Fc = prop(nFc > 0, within = (nA == 2)),
    # Proportion of core-fucosylated glycans within tri-antennary glycans
    A3Fc = prop(nFc > 0, within = (nA == 3)),
    # Proportion of core-fucosylated glycans within tetra-antennary glycans
    A4Fc = prop(nFc > 0, within = (nA == 4)),

    # ===== Advanced Arm-Fucosylation Traits =====

    # Proportion of arm-fucosylated glycans within mono-antennary glycans
    A1Fa = prop(nFa > 0, within = (nA == 1)),
    # Proportion of arm-fucosylated glycans within bi-antennary glycans
    A2Fa = prop(nFa > 0, within = (nA == 2)),
    # Proportion of arm-fucosylated glycans within tri-antennary glycans
    A3Fa = prop(nFa > 0, within = (nA == 3)),
    # Proportion of arm-fucosylated glycans within tetra-antennary glycans
    A4Fa = prop(nFa > 0, within = (nA == 4)),

    # Proportion of arm-fucosylated glycans within sialylated mono-antennary glycans
    A1SFa = prop(nFa > 0, within = (nA == 1 & nS > 0)),
    # Proportion of arm-fucosylated glycans within sialylated bi-antennary glycans
    A2SFa = prop(nFa > 0, within = (nA == 2 & nS > 0)),
    # Proportion of arm-fucosylated glycans within sialylated tri-antennary glycans
    A3SFa = prop(nFa > 0, within = (nA == 3 & nS > 0)),
    # Proportion of arm-fucosylated glycans within sialylated tetra-antennary glycans
    A4SFa = prop(nFa > 0, within = (nA == 4 & nS > 0)),

    # Proportion of arm-fucosylated glycans within asialylated mono-antennary glycans
    A1S0Fa = prop(nFa > 0, within = (nA == 1 & nS == 0)),
    # Proportion of arm-fucosylated glycans within asialylated bi-antennary glycans
    A2S0Fa = prop(nFa > 0, within = (nA == 2 & nS == 0)),
    # Proportion of arm-fucosylated glycans within asialylated tri-antennary glycans
    A3S0Fa = prop(nFa > 0, within = (nA == 3 & nS == 0)),
    # Proportion of arm-fucosylated glycans within asialylated tetra-antennary glycans
    A4S0Fa = prop(nFa > 0, within = (nA == 4 & nS == 0)),

    # ===== Advanced Bisecting Traits =====

    # Proportion of bisecting glycans within mono-antennary glycans
    A1B = prop(B, within = (nA == 1)),
    # Proportion of bisecting glycans within bi-antennary glycans
    A2B = prop(B, within = (nA == 2)),
    # Proportion of bisecting glycans within tri-antennary glycans
    A3B = prop(B, within = (nA == 3)),
    # Proportion of bisecting glycans within tetra-antennary glycans
    A4B = prop(B, within = (nA == 4)),

    # Proportion of bisecting glycans within core-fucosylated mono-antennary glycans
    A1FcB = prop(B, within = (nA == 1 & nFc > 0)),
    # Proportion of bisecting glycans within core-fucosylated bi-antennary glycans
    A2FcB = prop(B, within = (nA == 2 & nFc > 0)),
    # Proportion of bisecting glycans within core-fucosylated tri-antennary glycans
    A3FcB = prop(B, within = (nA == 3 & nFc > 0)),
    # Proportion of bisecting glycans within core-fucosylated tetra-antennary glycans
    A4FcB = prop(B, within = (nA == 4 & nFc > 0)),

    # Propotion of bisecting glycans within a-core-fucosylated mono-antennary glycans
    A1Fc0B = prop(B, within = (nA == 1 & nFc == 0)),
    # Propotion of bisecting glycans within a-core-fucosylated bi-antennary glycans
    A2Fc0B = prop(B, within = (nA == 2 & nFc == 0)),
    # Propotion of bisecting glycans within a-core-fucosylated tri-antennary glycans
    A3Fc0B = prop(B, within = (nA == 3 & nFc == 0)),
    # Propotion of bisecting glycans within a-core-fucosylated tetra-antennary glycans
    A4Fc0B = prop(B, within = (nA == 4 & nFc == 0)),

    # ===== Advanced Galactosylation Traits =====

    # Average degree of galactosylation per antenna within mono-antennary glycans
    A1G = wmean(nG / nA, within = (nA == 1)),
    # Average degree of galactosylation per antenna within bi-antennary glycans
    A2G = wmean(nG / nA, within = (nA == 2)),
    # Average degree of galactosylation per antenna within tri-antennary glycans
    A3G = wmean(nG / nA, within = (nA == 3)),
    # Average degree of galactosylation per antenna within tetra-antennary glycans
    A4G = wmean(nG / nA, within = (nA == 4)),

    # Average degree of terminal galactose per antenna within mono-antennary glycans
    A1Gt = wmean(nGt / nA, within = (nA == 1)),
    # Average degree of terminal galactose per antenna within bi-antennary glycans
    A2Gt = wmean(nGt / nA, within = (nA == 2)),
    # Average degree of terminal galactose per antenna within tri-antennary glycans
    A3Gt = wmean(nGt / nA, within = (nA == 3)),
    # Average degree of terminal galactose per antenna within tetra-antennary glycans
    A4Gt = wmean(nGt / nA, within = (nA == 4)),

    # ===== Advanced Sialylation Traits =====

    # Average degree of sialylation per antenna within mono-antennary glycans
    A1S = wmean(nS / nA, within = (nA == 1)),
    # Average degree of sialylation per antenna within bi-antennary glycans
    A2S = wmean(nS / nA, within = (nA == 2)),
    # Average degree of sialylation per antenna within tri-antennary glycans
    A3S = wmean(nS / nA, within = (nA == 3)),
    # Average degree of sialylation per antenna within tetra-antennary glycans
    A4S = wmean(nS / nA, within = (nA == 4)),

    # Average degree of sialylation per galactose within mono-antennary glycans
    A1GS = wmean(nS / nG, within = (nA == 1)),
    # Average degree of sialylation per galactose within bi-antennary glycans
    A2GS = wmean(nS / nG, within = (nA == 2)),
    # Average degree of sialylation per galactose within tri-antennary glycans
    A3GS = wmean(nS / nG, within = (nA == 3)),
    # Average degree of sialylation per galactose within tetra-antennary glycans
    A4GS = wmean(nS / nG, within = (nA == 4))
  )

  c(basic_traits(), additional_traits)
}

# To avoid note about global variables in R CMD check
Tp <- nM <- nA <- nF <- nFc <- nFa <- nG <- nS <- nGt <- nT <- B <- NULL