test_that("explain_trait works for proportion traits", {
  # Basic proportion traits from all_traits()
  expect_equal(
    explain_trait(all_traits()$TM),
    "Proportion of high-mannose glycans among all glycans."
  )
  
  expect_equal(
    explain_trait(all_traits()$TH),
    "Proportion of hybrid glycans among all glycans."
  )
  
  expect_equal(
    explain_trait(all_traits()$TC),
    "Proportion of complex glycans among all glycans."
  )
  
  expect_equal(
    explain_trait(all_traits()$TFc),
    "Proportion of core-fucosylated glycans among all glycans."
  )
  
  expect_equal(
    explain_trait(all_traits()$TFa),
    "Proportion of arm-fucosylated glycans among all glycans."
  )
  
  expect_equal(
    explain_trait(all_traits()$TS),
    "Proportion of sialylated glycans among all glycans."
  )
  
  expect_equal(
    explain_trait(all_traits()$TB),
    "Proportion of glycans with bisecting GlcNAc among all glycans."
  )
  
  # Proportion traits with within condition
  expect_equal(
    explain_trait(all_traits()$CA2),
    "Proportion of bi-antennary glycans within complex glycans."
  )
  
  expect_equal(
    explain_trait(all_traits()$CA3),
    "Proportion of tri-antennary glycans within complex glycans."
  )
  
  expect_equal(
    explain_trait(all_traits()$CA4),
    "Proportion of tetra-antennary glycans within complex glycans."
  )
})

test_that("explain_trait works for custom proportion traits", {
  # Simple conditions
  expect_equal(
    explain_trait(prop(nFc > 0)),
    "Proportion of core-fucosylated glycans among all glycans."
  )
  
  # Proportion with within condition
  expect_equal(
    explain_trait(prop(nFc > 0, within = (T == "complex"))),
    "Proportion of core-fucosylated glycans within complex glycans."
  )
  
  # Combined conditions with logical operators
  expect_equal(
    explain_trait(prop(nS > 0 & nFc > 0)),
    "Proportion of sialylated glycans with core-fucosylation among all glycans."
  )
  
  expect_equal(
    explain_trait(prop(nFc > 0 | nFa > 0)),
    "Proportion of core-fucosylated or arm-fucosylated glycans among all glycans."
  )
  
  # Negation
  expect_equal(
    explain_trait(prop(!B)),
    "Proportion of glycans without bisecting GlcNAc among all glycans."
  )
})

test_that("explain_trait works for ratio traits", {
  # Simple ratio
  expect_equal(
    explain_trait(ratio(T == "complex", T == "hybrid")),
    "Ratio of complex glycans to hybrid glycans among all glycans."
  )
  
  # Ratio with within condition
  expect_equal(
    explain_trait(ratio(nFc > 0, nFc == 0, within = (T == "complex"))),
    "Ratio of core-fucosylated glycans to non-core-fucosylated glycans within complex glycans."
  )
  
  # Ratio of bisecting vs non-bisecting
  expect_equal(
    explain_trait(ratio(B, !B)),
    "Ratio of glycans with bisecting GlcNAc to glycans without bisecting GlcNAc among all glycans."
  )
  
  # Test the fixed ambiguity issue
  expect_equal(
    explain_trait(prop(nFc > 0, within = nA == 2 & B)),
    "Proportion of core-fucosylated glycans within bi-antennary glycans with bisecting GlcNAc."
  )
})

test_that("explain_trait works for weighted-mean traits", {
  # Basic weighted means from all_traits()
  expect_equal(
    explain_trait(all_traits()$MM),
    "Abundance-weighted mean of mannose count within high-mannose glycans."
  )
  
  expect_equal(
    explain_trait(all_traits()$SG),
    "Abundance-weighted mean of degree of sialylation per galactose among all glycans."
  )
  
  expect_equal(
    explain_trait(all_traits()$GA),
    "Abundance-weighted mean of degree of galactosylation per antenna among all glycans."
  )
  
  # Simple weighted mean
  expect_equal(
    explain_trait(wmean(nS)),
    "Abundance-weighted mean of sialic acid count among all glycans."
  )
  
  # Weighted mean with within condition
  expect_equal(
    explain_trait(wmean(nA, within = (T == "complex"))),
    "Abundance-weighted mean of antenna count within complex glycans."
  )
  
  # Arithmetic expressions
  expect_equal(
    explain_trait(wmean(nS / nG)),
    "Abundance-weighted mean of degree of sialylation per galactose among all glycans."
  )
  
  # Test the fixed grammar issue
  expect_equal(
    explain_trait(wmean(nS / nG, within = nA == 4 & nFc > 0)),
    "Abundance-weighted mean of degree of sialylation per galactose within tetra-antennary glycans with core-fucosylation."
  )
  
  expect_equal(
    explain_trait(wmean(nG / nA)),
    "Abundance-weighted mean of degree of galactosylation per antenna among all glycans."
  )
})

test_that("explain_trait handles complex expressions with fallback", {
  # Complex expression should fallback to generic description
  complex_trait <- prop(nS > 2 & (nFc + nFa) >= 1 & nA %in% c(3, 4))
  result <- explain_trait(complex_trait)
  expect_true(stringr::str_detect(result, "satisfying"))
  expect_true(stringr::str_detect(result, "among all glycans"))
})

test_that("explain_trait errors for non-trait objects", {
  expect_error(
    explain_trait("not a trait"),
    "glydet trait"
  )
  
  expect_error(
    explain_trait(function(x) x),
    "glydet trait"
  )
  
  expect_error(
    explain_trait(list(a = 1)),
    "glydet trait"
  )
})

test_that("explain_trait handles different na_action values", {
  # The na_action shouldn't affect the explanation
  expect_equal(
    explain_trait(prop(nFc > 0, na_action = "zero")),
    "Proportion of core-fucosylated glycans among all glycans."
  )
  
  expect_equal(
    explain_trait(prop(nFc > 0, na_action = "keep")),
    "Proportion of core-fucosylated glycans among all glycans."
  )
})
