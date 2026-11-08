scale_test_exp <- function(proteomics = FALSE, legacy = FALSE) {
  structures <- glyparse::parse_iupac_condensed(c("Man(b1-4)GlcNAc", "GlcNAc"))
  info <- tibble::tibble(
    variable = c("V1", "V2"),
    glycan_structure = structures,
    glycan_composition = glyrepr::as_glycan_composition(structures),
    selected = c(TRUE, FALSE),
    count = c(2, 4)
  )
  if (proteomics) {
    info$protein <- "P1"
    info$protein_site <- 10L
  }
  exp <- glyexp::experiment(
    matrix(c(1, 3, 2, 6), 2, dimnames = list(c("V1", "V2"), c("S1", "S2"))),
    tibble::tibble(sample = c("S1", "S2"), group = factor(c("A", "B"))),
    info,
    exp_type = if (proteomics) "glycoproteomics" else "glycomics",
    glycan_type = "N"
  )
  exp$meta_data$scale_test_note <- "preserved"
  if (legacy) exp else as_test_se(exp)
}

set_test_scale <- function(exp, scale) {
  if (inherits(exp, "glyexp_experiment")) {
    exp$meta_data["scale"] <- list(scale)
  } else {
    S4Vectors::metadata(exp)["scale"] <- list(scale)
  }
  exp
}
