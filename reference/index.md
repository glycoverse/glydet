# Package index

## Main functions

Entry points for the package.

- [`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md)
  : Calculate Derived Traits
- [`derive_traits_()`](https://glycoverse.github.io/glydet/reference/derive_traits_.md)
  : Calculate Derived Traits from Tidy Data
- [`quantify_motifs()`](https://glycoverse.github.io/glydet/reference/quantify_motifs.md)
  : Quantify Motifs in an Experiment

## Out-of-box trait sets

Curated sets of derived traits for out-of-box use.

- [`traits_basic()`](https://glycoverse.github.io/glydet/reference/traits_basic.md)
  [`basic_traits()`](https://glycoverse.github.io/glydet/reference/traits_basic.md)
  : Get Basic Derived Traits
- [`traits_clerc_2018()`](https://glycoverse.github.io/glydet/reference/traits_clerc_2018.md)
  : Get Traits in Clerc et al. 2018
- [`traits_detailed()`](https://glycoverse.github.io/glydet/reference/traits_detailed.md)
  [`all_traits()`](https://glycoverse.github.io/glydet/reference/traits_detailed.md)
  : Get Detailed Derived Traits
- [`traits_fu_2026()`](https://glycoverse.github.io/glydet/reference/traits_fu_2026.md)
  : Get Traits in Fu et al. 2026
- [`traits_li_2025()`](https://glycoverse.github.io/glydet/reference/traits_li_2025.md)
  : Get Traits in Li et al. 2025

## Utility functions

Helper functions for trait interpretation and construction.

- [`explain_trait()`](https://glycoverse.github.io/glydet/reference/explain_trait.md)
  : Explain a Derived Trait
- [`make_trait()`](https://glycoverse.github.io/glydet/reference/make_trait.md)
  **\[experimental\]** : Use a Large Language Model (LLM) to create a
  derived trait function

## Trait factories

Functions for constructing custom derived traits.

- [`prop()`](https://glycoverse.github.io/glydet/reference/prop.md) :
  Create a Proportion Trait
- [`ratio()`](https://glycoverse.github.io/glydet/reference/ratio.md) :
  Create a Ratio Trait
- [`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md) :
  Create a Weighted-Mean Trait
- [`total()`](https://glycoverse.github.io/glydet/reference/total.md) :
  Create a Total Abundance Trait
- [`wsum()`](https://glycoverse.github.io/glydet/reference/wsum.md) :
  Create a Weighted Sum Trait

## Meta-property utilities

Functions for working with glycan meta-properties.

- [`add_meta_properties()`](https://glycoverse.github.io/glydet/reference/add_meta_properties.md)
  : Add Meta-Properties to Experiment
- [`all_mp_fns()`](https://glycoverse.github.io/glydet/reference/all_mp_fns.md)
  : Get All Meta-Property Functions
- [`get_meta_properties()`](https://glycoverse.github.io/glydet/reference/get_meta_properties.md)
  : Get Meta-Properties of Glycans

## Meta-property functions

Functions for calculating specific glycan meta-properties.

- [`n_glycan_type()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)
  [`has_bisecting()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)
  [`n_antennae()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)
  [`n_fuc()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)
  [`n_core_fuc()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)
  [`n_arm_fuc()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)
  [`n_gal()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)
  [`n_terminal_gal()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)
  [`n_sia()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)
  [`n_man()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)
  : Determine N-Glycan Key Properties
