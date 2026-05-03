#' Get Basic Derived Traits
#'
#' These derived traits are the most basic and commonly used derived traits.
#' They describe global properties of a glycome including the type of glycans,
#' fucosylation level, sialylation level, galactosylation level,
#' and branching level.
#'
#' @details
#' # Descriptions of traits
#'
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
#' Four additional sialic acid linkage traits are included if `sia_link = TRUE`.
#'
#' - `GE`: Average degree of a2,6-linked sialylation per galactose
#' - `GL`: Average degree of a2,3-linked sialylation per galactose
#' - `TE`: Proportion of a2,6-linked sialylated glycans
#' - `TL`: Proportion of a2,3-linked sialylated glycans
#'
#' # Usage of sialic acid linkage traits
#'
#' To use these sialic acid linkage traits,
#' `var_info` of the input `glyexp::experiment()` must have the following columns:
#'
#' - `nE`: Number of a2,6-linked sialic acids
#' - `nL`: Number of a2,3-linked sialic acids
#'
#' Note that you have to add these two columns even if the `glycan_structure` column has intact linkages.
#' This is because by convention all traits work with glycan structures with "basic" structure levels
#' (i.e., with generic monosaccharides like "Hex" and "HexNAc" and no linkages specified).
#'
#' @param sia_link A boolean indicating whether to include sialic acid linkage traits. Default is `FALSE`.
#'
#' @returns
#' A named list of derived traits.
#'
#' @examples
#' basic_traits()
#'
#' @export
basic_traits <- function(sia_link = FALSE) {
  checkmate::assert_flag(sia_link)
  traits <- list(
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
  if (sia_link) {
    cli::cli_alert_info("Please ensure that {.field nE} and {.field nL} are in {.field var_info}. See {.code ?basic_traits} for details.")
    sia_traits <- list(
      # Average degree of a2,6-linked sialylation per galactose
      GE = wmean(nE / nG),
      # Average degree of a2,3-linked sialylation per galactose
      GL = wmean(nL / nG),
      # Proportion of a2,6-linked sialylated glycans
      TE = prop(nE > 0),
      # Proportion of a2,3-linked sialylated glycans
      TL = prop(nL > 0)
    )
    traits <- c(traits, sia_traits)
  }
  traits
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
#' These additional sialic acid linkage traits are included if `sia_link = TRUE`.
#'
#' - `A1E`: Average degree of a2,6-linked sialylation per antenna within mono-antennary glycans
#' - `A2E`: Average degree of a2,6-linked sialylation per antenna within bi-antennary glycans
#' - `A3E`: Average degree of a2,6-linked sialylation per antenna within tri-antennary glycans
#' - `A4E`: Average degree of a2,6-linked sialylation per antenna within tetra-antennary glycans
#' - `A1L`: Average degree of a2,3-linked sialylation per antenna within mono-antennary glycans
#' - `A2L`: Average degree of a2,3-linked sialylation per antenna within bi-antennary glycans
#' - `A3L`: Average degree of a2,3-linked sialylation per antenna within tri-antennary glycans
#' - `A4L`: Average degree of a2,3-linked sialylation per antenna within tetra-antennary glycans
#' - `A1GE`: Average degree of a2,6-linked sialylation per galactose within mono-antennary glycans
#' - `A2GE`: Average degree of a2,6-linked sialylation per galactose within bi-antennary glycans
#' - `A3GE`: Average degree of a2,6-linked sialylation per galactose within tri-antennary glycans
#' - `A4GE`: Average degree of a2,6-linked sialylation per galactose within tetra-antennary glycans
#' - `A1GL`: Average degree of a2,3-linked sialylation per galactose within mono-antennary glycans
#' - `A2GL`: Average degree of a2,3-linked sialylation per galactose within bi-antennary glycans
#' - `A3GL`: Average degree of a2,3-linked sialylation per galactose within tri-antennary glycans
#' - `A4GL`: Average degree of a2,3-linked sialylation per galactose within tetra-antennary glycans
#'
#' @inheritSection basic_traits Usage of sialic acid linkage traits
#' @inheritParams basic_traits
#' @returns
#' A named list of derived traits.
#'
#' @examples
#' all_traits()
#'
#' @export
all_traits <- function(sia_link = FALSE) {
  checkmate::assert_flag(sia_link)
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

  traits <- c(suppressMessages(basic_traits(sia_link)), additional_traits)

  if (sia_link) {
    cli::cli_alert_info("Please ensure that {.field nE} and {.field nL} are in {.field var_info}. See {.code ?all_traits} for details.")
    sia_traits <- list(
      # Average degree of a2,6-linked sialylation per antenna within mono-antennary glycans
      A1E = wmean(nE / nA, within = (nA == 1)),
      # Average degree of a2,6-linked sialylation per antenna within bi-antennary glycans
      A2E = wmean(nE / nA, within = (nA == 2)),
      # Average degree of a2,6-linked sialylation per antenna within tri-antennary glycans
      A3E = wmean(nE / nA, within = (nA == 3)),
      # Average degree of a2,6-linked sialylation per antenna within tetra-antennary glycans
      A4E = wmean(nE / nA, within = (nA == 4)),
      # Average degree of a2,3-linked sialylation per antenna within mono-antennary glycans
      A1L = wmean(nL / nA, within = (nA == 1)),
      # Average degree of a2,3-linked sialylation per antenna within bi-antennary glycans
      A2L = wmean(nL / nA, within = (nA == 2)),
      # Average degree of a2,3-linked sialylation per antenna within tri-antennary glycans
      A3L = wmean(nL / nA, within = (nA == 3)),
      # Average degree of a2,3-linked sialylation per antenna within tetra-antennary glycans
      A4L = wmean(nL / nA, within = (nA == 4)),
      # Average degree of a2,6-linked sialylation per galactose within mono-antennary glycans
      A1GE = wmean(nE / nG, within = (nA == 1)),
      # Average degree of a2,6-linked sialylation per galactose within bi-antennary glycans
      A2GE = wmean(nE / nG, within = (nA == 2)),
      # Average degree of a2,6-linked sialylation per galactose within tri-antennary glycans
      A3GE = wmean(nE / nG, within = (nA == 3)),
      # Average degree of a2,6-linked sialylation per galactose within tetra-antennary glycans
      A4GE = wmean(nE / nG, within = (nA == 4)),
      # Average degree of a2,3-linked sialylation per galactose within mono-antennary glycans
      A1GL = wmean(nL / nG, within = (nA == 1)),
      # Average degree of a2,3-linked sialylation per galactose within bi-antennary glycans
      A2GL = wmean(nL / nG, within = (nA == 2)),
      # Average degree of a2,3-linked sialylation per galactose within tri-antennary glycans
      A3GL = wmean(nL / nG, within = (nA == 3)),
      # Average degree of a2,3-linked sialylation per galactose within tetra-antennary glycans
      A4GL = wmean(nL / nG, within = (nA == 4))
    )
    traits <- c(traits, sia_traits)
  }
  traits
}

#' Get Traits in Clerc et al. 2018
#'
#' These traits are the ones used by Clerc et al. 2018 (https://doi.org/10.1053/j.gastro.2018.05.030).
#' We generally don't recommend using these traits because they are either redundant
#' or missing some important traits.
#' We include this function because this paper is very influential in the field.
#'
#' @inheritSection basic_traits Usage of sialic acid linkage traits
#'
#' @inheritParams basic_traits
#' @returns A named list of derived traits.
#' @examples
#' traits_clerc_2018()[1:5]
#' @export
traits_clerc_2018 <- function(sia_link = FALSE) {
  traits <- list(
    CA1 = prop(nA == 1),
    CA2 = prop(nA == 2),
    CA3 = prop(nA == 3),
    CA4 = prop(nA == 4),
    TC = prop(Tp == "complex"),
    TM = prop(Tp == "highmannose"),
    THy = prop(Tp == "hybrid"),
    MM = wmean(nM, within = (Tp == "highmannose")),
    MHy = ratio(Tp == "highmannose", Tp == "hybrid"),
    TA2FS0 = prop(nA == 2 & nFc == 0 & nS == 0),
    CF = prop(nF > 0, within = (Tp == "complex")),
    CFa = prop(nFa > 0, within = (Tp == "complex")),
    A2Fa = prop(nFa > 0, within = nA == 2),
    A3Fa = prop(nFa > 0, within = nA == 3),
    A1F0 = prop(nA == 1 & nF == 0),
    A2F = prop(nF > 0, within = (nA == 2)),
    A3F = prop(nF > 0, within = (nA == 3)),
    A2SF = prop(nF > 0, within = (nA == 2 & nS > 0)),
    A2S0F = prop(nF > 0, within = (nA == 2 & nS == 0)),
    CB = prop(B, within = (Tp == "complex")),
    A2B = prop(B, within = (nA == 2)),
    A2FB = prop(B, within = (nA == 2 & nF > 0)),
    A2F0B = prop(B, within = (nA == 2 & nF == 0)),
    A2SB = prop(B, within = (nA == 2 & nS > 0)),
    A2S0B = prop(B, within = (nA == 2 & nS == 0)),
    A2FSB = prop(B, within = (nA == 2 & nF > 0 & nS > 0)),
    A2F0SB = prop(B, within = (nA == 2 & nF == 0 & nS > 0)),
    A2FS0B = prop(B, within = (nA == 2 & nF > 0 & nS == 0)),
    A2F0S0B = prop(B, within = (nA == 2 & nF == 0 & nS == 0)),
    CG = prop(nG > 0, within = (Tp == "complex")),
    A2G = wmean(nG, within = (nA == 2)),
    A2FG = wmean(nG, within = (nA == 2 & nF > 0)),
    A2F0G = wmean(nG, within = (nA == 2 & nF == 0)),
    A2SG = wmean(nG, within = (nA == 2 & nS > 0)),
    A2S0G = wmean(nG, within = (nA == 2 & nS == 0)),
    A2FSG = wmean(nG, within = (nA == 2 & nF > 0 & nS > 0)),
    A2FS0G = wmean(nG, within = (nA == 2 & nF > 0 & nS == 0)),
    A2F0SG = wmean(nG, within = (nA == 2 & nF == 0 & nS > 0)),
    A2F0S0G = wmean(nG, within = (nA == 2 & nF == 0 & nS == 0)),
    CS = prop(nS > 0, within = (Tp == "complex")),
    A2S = wmean(nS, within = (nA == 2)),
    A2FS = wmean(nS, within = (nA == 2 & nF > 0)),
    A2F0S = wmean(nS, within = (nA == 2 & nF == 0)),
    A2GS = wmean(nS / nG, within = (nA == 2 & nG > 0)),
    A3GS = wmean(nS / nG, within = (nA == 3 & nG > 0)),
    A2FGS = wmean(nS, within = (nA == 2 & nF > 0 & nG > 0)),
    A3FGS = wmean(nS, within = (nA == 3 & nF > 0 & nG > 0)),
    A2F0GS = wmean(nS, within = (nA == 2 & nF == 0 & nG > 0)),
    A3F0GS = wmean(nS, within = (nA == 3 & nF == 0 & nG > 0)),
    A4F0GS = wmean(nS, within = (nA == 4 & nF == 0 & nG > 0))
  )
  if (sia_link) {
    cli::cli_alert_info("Please ensure that {.field nE} and {.field nL} are in {.field var_info}. See {.code ?traits_florent_2018} for details.")
    sia_traits <- list(
      A2L0F = prop(nF > 0, within = (nA == 2 & nL == 0)),
      A3L0F = prop(nF > 0, within = (nA == 3 & nL == 0)),
      A2E0F = prop(nF > 0, within = (nA == 2 & nE == 0)),
      A2LF = prop(nF > 0, within = (nA == 2 & nL > 0)),
      A3LF = prop(nF > 0, within = (nA == 3 & nL > 0)),
      A1EF = prop(nF > 0, within = (nA == 1 & nE > 0)),
      A2EF = prop(nF > 0, within = (nA == 2 & nE > 0)),
      A3EF = prop(nF > 0, within = (nA == 3 & nE > 0)),
      A2L = wmean(nL, within = (nA == 2)),
      A2FL = wmean(nL, within = (nA == 2 & nF > 0)),
      A2F0L = wmean(nL, within = (nA == 2 & nF == 0)),
      A2GL = wmean(nL / nG, within = (nA == 2 & nG > 0)),
      A3GL = wmean(nL / nG, within = (nA == 3 & nG > 0)),
      A2FGL = wmean(nL / nG, within = (nA == 2 & nF > 0 & nG > 0)),
      A3FGL = wmean(nL / nG, within = (nA == 3 & nF > 0 & nG > 0)),
      A2F0GL = wmean(nL / nG, within = (nA == 2 & nF == 0 & nG > 0)),
      A3F0GL = wmean(nL / nG, within = (nA == 3 & nF == 0 & nG > 0)),
      A4F0GL = wmean(nL / nG, within = (nA == 4 & nF == 0 & nG > 0)),
      A2E = wmean(nE, within = (nA == 2)),
      A2FE = wmean(nE, within = (nA == 2 & nF > 0)),
      A2F0E = wmean(nE, within = (nA == 2 & nF == 0)),
      A2GE = wmean(nE / nG, within = (nA == 2 & nG > 0)),
      A3GE = wmean(nE / nG, within = (nA == 3 & nG > 0)),
      A2FGE = wmean(nE / nG, within = (nA == 2 & nF > 0 & nG > 0)),
      A3FGE = wmean(nE / nG, within = (nA == 3 & nF > 0 & nG > 0)),
      A2F0GE = wmean(nE / nG, within = (nA == 2 & nF == 0 & nG > 0)),
      A3F0GE = wmean(nE / nG, within = (nA == 3 & nF == 0 & nG > 0)),
      A4F0GE = wmean(nE / nG, within = (nA == 4 & nF == 0 & nG > 0))
    )
    traits <- c(traits, sia_traits)
  }
  traits
}

# To avoid note about global variables in R CMD check
Tp <- nM <- nA <- nF <- nFc <- nFa <- nG <- nS <- nE <- nL <- nGt <- nT <- B <- NULL