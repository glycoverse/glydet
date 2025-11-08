# Defining Custom Traits

This vignette provides a comprehensive guide for defining custom derived
traits using `glydet`. Before proceeding, please ensure you have read
the “Get Started with glydet” vignette and are familiar with the
fundamental concepts of derived traits and meta-properties.

``` r
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(glydet)
library(glyexp)
#> 
#> Attaching package: 'glyexp'
#> The following object is masked from 'package:dplyr':
#> 
#>     select_var
library(glyclean)
#> 
#> Attaching package: 'glyclean'
#> The following object is masked from 'package:stats':
#> 
#>     aggregate
library(glyrepr)

exp <- auto_clean(real_experiment)
#> ℹ Normalizing data (Median)
#> ✔ Normalizing data (Median) [146ms]
#> 
#> ℹ Removing variables with >50% missing values
#> ✔ Removing variables with >50% missing values [27ms]
#> 
#> ℹ Imputing missing values
#> ℹ Sample size <= 30, using sample minimum imputation
#> ℹ Imputing missing values✔ Imputing missing values [25ms]
#> 
#> ℹ Aggregating data
#> ✔ Aggregating data [1.1s]
#> 
#> ℹ Normalizing data again
#> ✔ Normalizing data again [18ms]
```

## Custom Traits

`glydet` provides three trait factory functions for creating custom
derived traits:

- [`prop()`](https://glycoverse.github.io/glydet/reference/prop.md) for
  calculating the abundance proportion of a specific glycan subset
- [`ratio()`](https://glycoverse.github.io/glydet/reference/ratio.md)
  for computing the abundance ratio between two glycan subsets
- [`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md)
  for calculating the weighted mean of a quantitative property, weighted
  by glycan abundances

All derived traits in `glydet` are constructed using these three factory
functions. The definitions of all built-in traits demonstrate their
usage:

- `TM`: `prop(Tp == "highmannose")`
- `TH`: `prop(Tp == "hybrid")`
- `TC`: `prop(Tp == "complex")`
- `MM`: `wmean(nM, within = (Tp == "highmannose"))`
- `CA2`: `prop(nA == 2, within = (Tp == "complex"))`
- `CA3`: `prop(nA == 3, within = (Tp == "complex"))`
- `CA4`: `prop(nA == 4, within = (Tp == "complex"))`
- `TF`: `prop(nF > 0)`
- `TFc`: `prop(nFc > 0)`
- `TFa`: `prop(nFa > 0)`
- `TB`: `prop(B)`
- `GS`: `wmean(nS / nG)`
- `AG`: `wmean(nG / nA)`
- `TS`: `prop(nS > 0)`

These definitions utilize meta-properties as building blocks. For
instance, `T` represents glycan type, `nM` denotes the number of mannose
residues, and so forth. To retrieve all available built-in
meta-properties, use `names(all_mp_fns())`.

The complete list of built-in meta-properties is provided below:

| Name  | Function                                                                             | Type    | Description                                                                      |
|-------|--------------------------------------------------------------------------------------|---------|----------------------------------------------------------------------------------|
| `T`   | [`n_glycan_type()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)  | factor  | Type of the glycan, either “complex”, “hybrid”, “highmannose”, or “pausimannose” |
| `B`   | [`has_bisecting()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)  | logical | Whether the glycan has a bisecting GlcNAc                                        |
| `nA`  | [`n_antennae()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)     | integer | Number of antennae                                                               |
| `nF`  | [`n_fuc()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)          | integer | Number of fucoses                                                                |
| `nFc` | [`n_core_fuc()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)     | integer | Number of core fucoses                                                           |
| `nFa` | [`n_arm_fuc()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)      | integer | Number of arm fucoses                                                            |
| `nG`  | [`n_gal()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)          | integer | Number of galactoses                                                             |
| `nGt` | [`n_terminal_gal()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md) | integer | Number of terminal galactoses                                                    |
| `nS`  | [`n_sia()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)          | integer | Number of sialic acids                                                           |
| `nM`  | [`n_man()`](https://glycoverse.github.io/glydet/reference/n_glycan_type.md)          | integer | Number of mannoses                                                               |

The following sections provide detailed explanations of each factory
function.

### Trait Factories

#### `prop()`

[`prop()`](https://glycoverse.github.io/glydet/reference/prop.md)
creates proportion traits, which represent the most common type of
derived traits. These traits calculate the relative abundance of a
specific glycan subset within a defined population. Examples include the
proportion of core-fucosylated glycans within all glycans, or the
proportion of tetra-antennary glycans within the complex glycan subset.

[`prop()`](https://glycoverse.github.io/glydet/reference/prop.md)
accepts an expression that evaluates to a logical vector. Both built-in
and custom meta-properties (covered later) can be referenced within the
expression.

For example, the proportion of core-fucosylated glycans within all
glycans is defined as:

``` r
prop(nFc > 0)
```

Since `nFc` is an integer meta-property representing the number of core
fucoses, the expression `nFc > 0` evaluates to a logical vector,
creating a valid trait definition.

Consider this simpler example:

``` r
prop(B)
```

This demonstrates a straightforward trait definition. Since `B`
(bisecting GlcNAc presence) is already a logical meta-property,
`prop(B)` constitutes a valid trait definition.

A more complex example demonstrates compound logical operations:

``` r
prop(nS > 0 & nFa > 0)
```

This trait calculates the proportion of glycans containing both sialic
acid and arm fucose. Both `nS > 0` and `nFa > 0` represent logical
expressions, combined using the logical AND operator (`&`). Any R
logical operator (including `|`, `!`, etc.) can be utilized within these
expressions.

Arithmetic calculations can also be incorporated within expressions. For
example, the `CF` trait is defined as:

``` r
prop(nF > 0)
```

We know that `nFc + nFa` is equivalent to `nF`, so the above definition
is equivalent to:

``` r
prop((nFc + nFa) > 0)
```

As a general principle, any R expression that evaluates to a logical
vector is valid for use with
[`prop()`](https://glycoverse.github.io/glydet/reference/prop.md).

The traits described above calculate proportions relative to the entire
glycan population. To customize the denominator (reference population),
the `within` parameter can be employed. For example, calculating the
proportion of bi-antennary glycans within the **complex** glycan subset
(`CA2`) requires this approach.

The `within` parameter accepts an expression that evaluates to a logical
vector, following the same syntax as the primary parameter.

The `CA2` trait is defined as follows:

``` r
prop(nA == 2, within = (Tp == "complex"))
```

Note that parentheses around `Tp == "complex"` are optional.
Consequently, this alternative definition is equally valid:

``` r
prop(nA == 2, within = Tp == "complex")
```

Combining these two parameters enables the creation of sophisticated
trait definitions. For instance, the proportion of sialylated glycans
within the subset of core-fucosylated tetra-antennary glycans is defined
as:

``` r
prop(nS > 0, within = (nFc > 0 & nA == 4))
```

Having mastered the usage of
[`prop()`](https://glycoverse.github.io/glydet/reference/prop.md),
examine the definitions of built-in proportion traits to reinforce your
understanding.

Practice exercises:

1.  Define a trait calculating the proportion of glycans containing
    terminal galactose.
2.  Define a trait calculating the proportion of glycans with terminal
    galactose but no sialic acid.
3.  Define a trait calculating the proportion of glycans with exactly
    two antennae and no bisecting GlcNAc.
4.  Define a trait calculating the proportion of glycans with bisecting
    GlcNAc within the bi-antennary glycan subset.

Solutions are provided at the end of this vignette.

#### `ratio()`

[`ratio()`](https://glycoverse.github.io/glydet/reference/ratio.md)
creates ratio traits that represent the quotient of total abundances
between two glycan groups. Examples include the ratio of complex to
hybrid glycans, or the ratio of bisecting to non-bisecting glycans.

[`ratio()`](https://glycoverse.github.io/glydet/reference/ratio.md)
accepts two expressions that evaluate to logical vectors, similar to
[`prop()`](https://glycoverse.github.io/glydet/reference/prop.md).

For example, the ratio of complex to hybrid glycans is defined as:

``` r
ratio(Tp == "complex", Tp == "hybrid")
```

The first expression defines the numerator, while the second expression
defines the denominator.

The `within` parameter can be employed to apply restrictions to both
numerator and denominator. For example, the ratio of bisecting to
non-bisecting glycans within the bi-antennary subset is defined as:

``` r
ratio(B, !B, within = (nA == 2))
```

This represents syntactic sugar for the more verbose:

``` r
ratio(B & (nA == 2), (!B) & (nA == 2))
```

The `within` parameter provides clearer semantics and reduces
redundancy.

Note that
[`prop()`](https://glycoverse.github.io/glydet/reference/prop.md)
represents a special case of
[`ratio()`](https://glycoverse.github.io/glydet/reference/ratio.md),
where `prop(cond, within)` is mathematically equivalent to
`ratio(cond & within, within)`. The following two definitions are
functionally identical:

``` r
# using prop()
prop(nFc > 0, within = (Tp == "complex"))

# using ratio()
ratio(nFc > 0 & Tp == "complex", Tp == "complex")
```

However,
[`prop()`](https://glycoverse.github.io/glydet/reference/prop.md) is
recommended when appropriate for enhanced readability, as it provides
more intuitive semantics that align with natural language. The
[`prop()`](https://glycoverse.github.io/glydet/reference/prop.md)
example above can be interpreted as “the proportion of core-fucosylated
glycans within all complex glycans,” whereas the
[`ratio()`](https://glycoverse.github.io/glydet/reference/ratio.md)
version requires additional cognitive processing.

Practice exercises:

1.  Define a trait calculating the ratio of core-fucosylated to
    non-core-fucosylated glycans.
2.  Define a trait calculating the ratio of tetra-antennary to
    tri-antennary glycans within the complex glycan subset.
3.  Redefine the `CA2` trait using
    [`ratio()`](https://glycoverse.github.io/glydet/reference/ratio.md).

#### `wmean()`

[`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md)
creates weighted-mean traits, which represent the most sophisticated yet
powerful type of derived traits. These traits calculate the
abundance-weighted average of quantitative properties across glycan
populations. Examples include the average number of antennae across all
glycans, or the average degree of sialylation per galactose across the
entire glycan repertoire.

Before exploring
[`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md)
usage, it is essential to understand the weighted-mean concept.

Consider three glycans: `G1`, `G2`, and `G3`. Their respective antennae
counts are 1, 2, and 3, with relative abundances of 50%, 20%, and 30%.
The weighted-mean number of antennae is calculated as:

``` r
(1 * 0.5 + 2 * 0.2 + 3 * 0.3) / (0.5 + 0.2 + 0.3)
#> [1] 1.8
```

This value represents the abundance-weighted average degree of branching
across the glycan population.

This concept provides substantial analytical power, enabling the
calculation of numerous biologically meaningful traits.

Consider another example: the average degree of sialylation per antenna
across all glycans. Using the same three glycans (`G1`, `G2`, and `G3`)
with antennae counts of 1, 2, and 3, sialic acid counts of 1, 1, and 3,
respectively, and relative abundances of 50%, 20%, and 30%.

First, calculate the sialylation degree per antenna for each glycan
(sialic acid count divided by antenna count):

``` r
c(1/1, 1/2, 3/3)
#> [1] 1.0 0.5 1.0
```

Subsequently, incorporate the abundance weighting:

``` r
(1/1 * 0.5 + 1/2 * 0.2 + 3/3 * 0.3) / (0.5 + 0.2 + 0.3)
#> [1] 0.9
```

This yields the abundance-weighted average degree of sialylation per
antenna across the glycan population.

Having established the weighted-mean concept and mastered
[`prop()`](https://glycoverse.github.io/glydet/reference/prop.md) and
[`ratio()`](https://glycoverse.github.io/glydet/reference/ratio.md)
fundamentals, we can now examine
[`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md)
implementation.

[`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md)
accepts expressions that evaluate to numeric vectors, contrasting with
the logical vectors required by
[`prop()`](https://glycoverse.github.io/glydet/reference/prop.md) and
[`ratio()`](https://glycoverse.github.io/glydet/reference/ratio.md).

For example, the average number of antennae across all glycans is
defined as:

``` r
wmean(nA)
```

The average degree of sialylation per antenna across all glycans is
defined as:

``` r
wmean(nS / nA)
```

The `within` parameter can be utilized to restrict the weighted-mean
calculation to specific glycan subsets. For example, the average number
of antennae within complex glycans is defined as:

``` r
wmean(nA, within = (Tp == "complex"))
```

Practice exercises:

1.  Define a trait calculating the average degree of sialylation per
    antenna.
2.  Define a trait calculating the average number of arm fucoses.

#### `total()`

[`total()`](https://glycoverse.github.io/glydet/reference/total.md)
creates total abundance traits, which is the sum of the abundances of a
group of glycans. This is the simplest type of derived traits in
`glydet`.

For example, the total abundance of all complex glycans is defined as:

``` r
total(Tp == "complex")
```

There is no `within` parameter for
[`total()`](https://glycoverse.github.io/glydet/reference/total.md), as
you can always add a restriction to the expression by using `&`.

``` r
total(Tp == "complex" & nA == 4)
```

#### `wsum()`

[`wsum()`](https://glycoverse.github.io/glydet/reference/wsum.md) is the
cousin of
[`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md), but
instead of calculating the average value, it calculates the sum of the
values.

Let’s see an example to understand the difference between
[`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md) and
[`wsum()`](https://glycoverse.github.io/glydet/reference/wsum.md). Say
we have three glycans: `G1`, `G2`, and `G3` in two samples: `S1` and
`S2`.

- Sample S1: `G1` (10), `G2` (20), `G3` (30) (intensity values)
- Sample S2: `G1` (20), `G2` (40), `G3` (60) (intensity values)

``` r
expr_mat <- matrix(c(10, 20, 30, 20, 40, 60), nrow = 3)
rownames(expr_mat) <- c("G1", "G2", "G3")
colnames(expr_mat) <- c("S1", "S2")
expr_mat
#>    S1 S2
#> G1 10 20
#> G2 20 40
#> G3 30 60
```

For simplicity, we don’t use any meta-properties here, just the number
1.

``` r
trait1 <- wmean(1)
trait1(expr_mat, tibble::tibble())
#> [1] 1 1
```

``` r
trait2 <- wsum(1)
trait2(expr_mat, tibble::tibble())
#> [1]  60 120
```

See the difference?
[`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md)
calculates the average value of some numeric property, so as long as the
ratio of each glycan’s abundance is the same, the result is always the
same. But
[`wsum()`](https://glycoverse.github.io/glydet/reference/wsum.md) is
sensitive to the abundance of each glycan, so the result is different.

### Using Custom Traits

Having learned how to define custom traits, we can now implement them.

[`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md)
includes a `trait_fns` parameter that accepts a named list of derived
trait functions. Note that the three factory functions described above
return functions as their output.

``` r
class(prop(nFc > 0))
#> [1] "glydet_prop"  "glydet_trait"
```

R supports functional programming paradigms, treating functions as
first-class objects. Functions can be passed as arguments to other
functions and returned from function calls.

While understanding this concept enhances R proficiency, it is not
prerequisite for usage—simply define traits as a named list and provide
it to the `trait_fns` parameter.

The following example examines sialylation degree within glycan subsets
of varying antenna counts:

``` r
my_traits <- list(
  A2S = wmean(nS / nA, within = (nA == 2)),
  A3S = wmean(nS / nA, within = (nA == 3)),
  A4S = wmean(nS / nA, within = (nA == 4))
)
derive_traits(exp, trait_fns = my_traits)
```

The identifiers “A2S”, “A3S”, and “A4S” represent the derived trait
names.

When `trait_fns` is omitted,
[`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md)
internally invokes
[`basic_traits()`](https://glycoverse.github.io/glydet/reference/basic_traits.md)
to utilize built-in trait functions.

### Validating Trait Definitions

How can you ensure your trait definitions are accurate and meaningful?
The
[`explain_trait()`](https://glycoverse.github.io/glydet/reference/explain_trait.md)
function provides an intuitive way to validate your trait definitions by
generating human-readable explanations.

This function helps you:

- Verify that your trait logic matches your intended analysis
- Understand complex trait expressions at a glance  
- Debug and refine trait definitions before deployment

Let’s examine the `A2S` trait defined earlier as an example:

``` r
explain_trait(wmean(nS / nA, within = (nA == 2)))
#> [1] "Abundance-weighted mean of degree of sialylation per antenna within bi-antennary glycans."
```

The function interprets your trait expression and returns a clear,
natural language description of what the trait calculates. This
validation step is particularly valuable when working with complex trait
definitions or when collaborating with team members who need to
understand your analytical approach.

## Custom Meta-Properties

The preceding examples utilized exclusively built-in meta-properties. As
the building blocks of derived traits, understanding how to define
custom meta-properties will bring you more flexibility.

`glydet` provides two means of defining custom meta-properties:

- By providing meta-property functions to the `mp_fns` parameter of
  [`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md).
- By using columns in the variable information tibble as
  meta-properties.

### Defining Custom Meta-Properties with Functions

Meta-property functions must accept a
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
vector and return a vector of corresponding meta-property values.

For example, the built-in meta-property `B` is defined as follows:

``` r
# Simplified version for illustration purposes
function(glycans) {
  motif <- "HexNAc(??-?)Hex(??-?)HexNAc(??-?)HexNAc(??-"
  glymotif::have_motif(glycans, motif, alignment = "core")
}
```

This implementation utilizes
[`glymotif::have_motif()`](https://glycoverse.github.io/glymotif/reference/have_motif.html)
to determine motif presence within glycan structures. Meta-property
function definitions offer considerable flexibility. Any tools from the
`glycoverse` ecosystem or broader R environment can be employed,
provided the function accepts a
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
vector and returns an atomic vector. All built-in meta-properties
utilize the `glymotif` package for motif matching and counting
operations. For detailed information, consult the [glymotif
documentation](https://glycoverse.github.io/glymotif/).

The following example demonstrates custom meta-property creation:

``` r
my_mp_fns <- list(
  # Number of Lewis antigens
  nLe = ~ glymotif::count_motif(.x, "Hex(??-?)[dHex(??-?)]HexNAc(??-"),
  # Number of poly-LacNAc units
  nPl = ~ glymotif::count_motif(.x, "Hex(??-?)HexNAc(??-?)Hex(??-?)HexNAc(??-")
)
```

Note that purrr-style lambda functions provide a concise syntax for
meta-property function definitions.

Traits can be defined using these new meta-properties:

``` r
my_traits <- list(
  # Average number of Lewis antigens per antenna
  LeA = wmean(nLe / nA),
  # Average number of poly-LacNAc units per glycan
  Pl = wmean(nPl)
)
```

These expressions incorporate both custom meta-properties (`nLe` and
`nPl`) and the built-in meta-property `nA`.

Implementation of these custom meta-properties and traits proceeds as
follows:

``` r
derive_traits(exp, trait_fns = my_traits, mp_fns = c(my_mp_fns, all_mp_fns()))
#> 
#> ── Traitproteomics Experiment ──────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 548 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: protein <chr>, protein_site <int>, trait <chr>, gene <chr>
```

Ensure that custom meta-properties are combined with built-in
meta-properties when the latter are required within trait expressions.

Let’s see another example about sialic acid linkage types. Neu5Ac can
have two different linkage types: a2-3 and a2-6. We can define two
custom meta-properties to count the number of sialic acids with each
linkage type:

``` r
sia_mp_fns <- list(
  # Number of a2-3 sialic acids
  nL = ~ glymotif::count_motif(.x, "NeuAc(a2-3)Hex(??-"),
  # Number of a2-6 sialic acids
  nE = ~ glymotif::count_motif(.x, "NeuAc(a2-6)Hex(??-")
)
```

And define two traits to calculate the degree of a2-3 and a2-6
sialylation per galactose:

``` r
sia_traits <- list(
  # Average degree of a2-3 sialylation per galactose
  GL = wmean(nL / nG),
  # Average degree of a2-6 sialylation per galactose
  GE = wmean(nE / nG)
)
```

And calculate the traits:

``` r
derive_traits(exp, trait_fns = sia_traits, mp_fns = c(sia_mp_fns, all_mp_fns()))
#> 
#> ── Traitproteomics Experiment ──────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 548 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: protein <chr>, protein_site <int>, trait <chr>, gene <chr>
```

This is an example of how you can violate the [ambiguity
assumption](https://glycoverse.github.io/glydet/articles/glydet.html#working-with-structural-ambiguity)
of built-in meta-properties and derived traits, by introducing more
sophisticated and linkage-awared meta-properties.

### Defining Custom Meta-Properties with Columns

The first method is handy if all information you need is already encoded
in glycan structures. However, there are situations that you need to
rely on other meta-data. In this case, you can use columns in the
variable information tibble as meta-properties directly, by specifying
the `mp_cols` parameter.

For example, let’s use the sialic acid linkage type example above. You
might have used special derivatization methods to differentiate the two
linkage types. This information might end up in a column in the variable
information tibble, not directly in the glycan structures.

``` r
# Here we assume all sialic acids are a2-6
exp2 <- exp |>
  mutate_var(
    n_a26_sia = count_mono(glycan_structure, "NeuAc"),
    n_a23_sia = 0L
  )
```

``` r
exp2 |>
  get_var_info() |>
  filter(n_a26_sia > 0) |>
  pull(glycan_structure)
#> <glycan_structure[2526]>
#> [1] NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)[NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-
#> [2] NeuAc(??-?)Hex(??-?)HexNAc(??-?)[HexNAc(??-?)]Hex(??-?)[Hex(??-?)Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-
#> [3] NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)[Hex(??-?)HexNAc(??-?)Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-
#> [4] NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)[NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-
#> [5] NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)[HexNAc(??-?)Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-
#> [6] NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)[NeuAc(??-?)Hex(??-?)[dHex(??-?)]HexNAc(??-?)Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-
#> [7] NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)[NeuAc(??-?)Hex(??-?)[dHex(??-?)]HexNAc(??-?)Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-
#> [8] NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)[Hex(??-?)]Hex(??-?)HexNAc(??-?)[dHex(??-?)]HexNAc(??-
#> [9] NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)[Hex(??-?)HexNAc(??-?)Hex(??-?)]Hex(??-?)HexNAc(??-?)[dHex(??-?)]HexNAc(??-
#> [10] NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)[NeuAc(??-?)Hex(??-?)HexNAc(??-?)Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-
#> ... (2516 more not shown)
#> # Unique structures: 526
```

See? No linkage information about sialic acids is in the glycan
structures. But we have two columns: `n_a26_sia` and `n_a23_sia` storing
this information.

Now we can define two traits to calculate the degree of a2-3 and a2-6
sialylation per galactose:

``` r
sia_traits <- list(
  # Average degree of a2-3 sialylation per galactose
  GL = wmean(nL / nG),
  # Average degree of a2-6 sialylation per galactose
  GE = wmean(nE / nG)
)
```

And use `mp_cols` to tell
[`derive_traits()`](https://glycoverse.github.io/glydet/reference/derive_traits.md)
to use these columns as meta-properties:

``` r
derive_traits(exp2, trait_fns = sia_traits, mp_cols = c(nL = "n_a23_sia", nE = "n_a26_sia"))
#> 
#> ── Traitproteomics Experiment ──────────────────────────────────────────────────
#> ℹ Expression matrix: 12 samples, 548 variables
#> ℹ Sample information fields: group <fct>
#> ℹ Variable information fields: protein <chr>, protein_site <int>, trait <chr>, gene <chr>, n_a23_sia <int>
```

## Exercise Solutions

**[`prop()`](https://glycoverse.github.io/glydet/reference/prop.md)**

1.  `prop(nGt > 0)`
2.  `prop(nGt > 0 & nS == 0)`
3.  `prop(nA == 2 & !B)`
4.  `prop(B, within = (nA == 2))`

**[`ratio()`](https://glycoverse.github.io/glydet/reference/ratio.md)**

1.  `ratio(nFc > 0, nFc == 0)`
2.  `ratio(nA == 4, nA == 3, within = (Tp == "complex"))`
3.  `ratio(nA == 2 & Tp == "complex", Tp == "complex")`

**[`wmean()`](https://glycoverse.github.io/glydet/reference/wmean.md)**

1.  `wmean(nS / nA)`
2.  `wmean(nFa)`
