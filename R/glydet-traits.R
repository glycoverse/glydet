all_traits <- function() {
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
    TF = prop((nFc + nFa) > 0),
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