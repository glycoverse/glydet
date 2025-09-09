#' Get All Derived Traits
#'
#' This function returns a named list of all derived traits:
#' - "TM": `prop(T == "highmannose")`
#' - "TH": `prop(T == "hybrid")`
#' - "TC": `prop(T == "complex")`
#' - "MM": `wmean(nM, within = (T == "highmannose"))`
#' - "CA2": `prop(nA == 2, within = (T == "complex"))`
#' - "CA3": `prop(nA == 3, within = (T == "complex"))`
#' - "CA4": `prop(nA == 4, within = (T == "complex"))`
#' - "TF": `prop((nFc + nFa) > 0)`
#' - "TFc": `prop(nFc > 0)`
#' - "TFa": `prop(nFa > 0)`
#' - "TB": `prop(B)`
#' - "SG": `wmean(nS / nG)`
#' - "GA": `wmean(nG / nA)`
#' - "TS": `prop(nS > 0)`
#'
#' @returns
#' A named list of derived traits.
#'
#' @export
all_traits <- function() {
  list(
    # Proportion of highmannose glycans
    TM = prop(.data$T == "highmannose"),
    # Proportion of hybrid glycans
    TH = prop(.data$T == "hybrid"),
    # Proportion of complex glycans
    TC = prop(.data$T == "complex"),
    # Average number of mannoses within highmannose glycans
    MM = wmean(.data$nM, within = (.data$T == "highmannose")),
    # Proportion of bi-antennary glycans within complex glycans
    CA2 = prop(.data$nA == 2, within = (.data$T == "complex")),
    # Proportion of tri-antennary glycans within complex glycans
    CA3 = prop(.data$nA == 3, within = (.data$T == "complex")),
    # Proportion of tetra-antennary glycans within complex glycans
    CA4 = prop(.data$nA == 4, within = (.data$T == "complex")),
    # Proportion of fucosylated glycans
    TF = prop((.data$nFc + .data$nFa) > 0),
    # Proportion of core-fucosylated glycans
    TFc = prop(.data$nFc > 0),
    # Proportion of arm-fucosylated glycans
    TFa = prop(.data$nFa > 0),
    # Proportion of glycans with bisecting GlcNAc
    TB = prop(.data$B),
    # Average degree of sialylation per galactose
    SG = wmean(.data$nS / .data$nG),
    # Average degree of galactosylation per antenna
    GA = wmean(.data$nG / .data$nA),
    # Proportion of sialylated glycans
    TS = prop(.data$nS > 0)
  )
}