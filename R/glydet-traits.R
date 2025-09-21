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
#' - `SG`: Average degree of sialylation per galactose
#' - `GA`: Average degree of galactosylation per antenna
#' - `TS`: Proportion of sialylated glycans
#'
#' @returns
#' A named list of derived traits.
#'
#' @export
basic_traits <- function() {
  list(
    # Proportion of highmannose glycans
    TM = prop(T == "highmannose"),
    # Proportion of hybrid glycans
    TH = prop(T == "hybrid"),
    # Proportion of complex glycans
    TC = prop(T == "complex"),
    # Average number of mannoses within highmannose glycans
    MM = wmean(nM, within = (T == "highmannose")),
    # Proportion of bi-antennary glycans within complex glycans
    CA2 = prop(nA == 2, within = (T == "complex")),
    # Proportion of tri-antennary glycans within complex glycans
    CA3 = prop(nA == 3, within = (T == "complex")),
    # Proportion of tetra-antennary glycans within complex glycans
    CA4 = prop(nA == 4, within = (T == "complex")),
    # Proportion of fucosylated glycans
    TF = prop(nF > 0),
    # Proportion of core-fucosylated glycans
    TFc = prop(nFc > 0),
    # Proportion of arm-fucosylated glycans
    TFa = prop(nFa > 0),
    # Proportion of glycans with bisecting GlcNAc
    TB = prop(B),
    # Average degree of sialylation per galactose
    SG = wmean(nS / nG),
    # Average degree of galactosylation per antenna
    GA = wmean(nG / nA),
    # Proportion of sialylated glycans
    TS = prop(nS > 0)
  )
}

# To avoid note about global variables in R CMD check
T <- nM <- nA <- nF <- nFc <- nFa <- nG <- nS <- nGt <- nT <- B <- NULL