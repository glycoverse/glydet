# Determine N-Glycan Key Properties

These functions check key properties of an N-glycan:

- `n_glycan_type()`: Determine the N-glycan type.

- `has_bisecting()`: Check if the glycan has a bisecting GlcNAc.

- `n_antennae()`: Count the number of antennae.

- `n_fuc()`: Count the number of fucoses.

- `n_core_fuc()`: Count the number of core fucoses.

- `n_arm_fuc()`: Count the number of arm fucoses.

- `n_gal()`: Count the number of galactoses.

- `n_terminal_gal()`: Count the number of terminal galactoses.

- `n_sia()`: Count the number of sialic acids.

- `n_man()`: Count the number of mannoses.

All functions assume the glycans are N-glycans without validation, thus
may return meaningless values for non-N-glycans. Therefore, please make
sure to pass in N-glycans only.

All functions put minimum requirement on the glycans, i.e. they work
with glycans with generic monosaccharides (e.g. "Hex", "HexNAc") and no
linkage information. This type of structures are common in
glycoproteomics and glycomics studies.

## Usage

``` r
n_glycan_type(glycans)

has_bisecting(glycans)

n_antennae(glycans)

n_fuc(glycans)

n_core_fuc(glycans)

n_arm_fuc(glycans)

n_gal(glycans)

n_terminal_gal(glycans)

n_sia(glycans)

n_man(glycans)
```

## Arguments

- glycans:

  A
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
  vector.

## Value

- `n_glycan_type()`: A factor vector indicating the N-glycan type,
  either "highmannose", "hybrid", "complex", or "paucimannose".

- `has_bisecting()`: A logical vector indicating if the glycan has a
  bisecting GlcNAc.

- `n_antennae()`: An integer vector indicating the number of antennae.

- `n_fuc()`: An integer vector indicating the number of fucoses.

- `n_core_fuc()`: An integer vector indicating the number of core
  fucoses.

- `n_arm_fuc()`: An integer vector indicating the number of arm fucoses.

- `n_gal()`: An integer vector indicating the number of galactoses.

- `n_terminal_gal()`: An integer vector indicating the number of
  terminal galactoses.

- `n_sia()`: An integer vector indicating the number of sialic acids.

- `n_man()`: An integer vector indicating the number of mannoses.

## `n_glycan_type()`: N-Glycan Types

Four types of N-glycans are recognized: high mannose, hybrid, complex,
and paucimannose. For more information about N-glycan types, see
[Essentials of
Glycobiology](https://www.ncbi.nlm.nih.gov/books/NBK579964/#_s9_2_).

## `has_bisecting()`: Bisecting GlcNAc

Bisecting GlcNAc is a GlcNAc residue attached to the core mannose of
N-glycans.

         Man
            \
    GlcNAc - Man - GlcNAc - GlcNAc -
    ~~~~~~  /
         Man

## `n_antennae()`: Number of Antennae

The number of antennae is the number of branching GlcNAc to the core
mannoses.

## `n_fuc()`: Number of Fucoses

Number of fucoses. This function assumes the fucose is a dHex.

## `n_core_fuc()`: Number of Core Fucoses

Core fucoses are those fucose residues attached to the core GlcNAc of an
N-glycan.

    Man             Fuc  <- core fucose
       \             |
        Man - GlcNAc - GlcNAc -
       /
    Man

## `n_arm_fuc()`: Number of Arm Fucoses

Arm focuses are those focuse residues attached to the branching GlcNAc
of an N-glycan.

     Fuc  <- arm fucose
      |
    GlcNAc - Man
                \
                 Man - GlcNAc - GlcNAc -
                /
    GlcNAc - Man

## `n_gal()`: Number of Galactoses

This function seems useless and silly. It is, if you have a
well-structured glycan with concrete monosaccharides. However, if you
only have "Hex" or "H" at hand, it is tricky to know how many of them
are "Gal" and how many are "Man". This function makes a simply
assumption that all the rightmost "H" in a "H-H-N-H" unit is a
galactose. The two "H" on the left are mannoses of the N-glycan core.
The "N" is a GlcNAc attached to one core mannose.

## `n_terminal_gal()`: Number of Terminal Galactoses

Terminal galactoses are those galactose residues on the non-reducing end
without sialic acid capping.

             Gal - GlcNAc - Man
             ~~~               \
         terminal Gal           Man - GlcNAc - GlcNAc -
                               /
    Neu5Ac - Gal - GlcNAc - Man
             ~~~
       not terminal Gal

## `n_sia()`: Number of Sialic Acids

Number of sialic acids (Neu5Ac). Neu5Gc is not counted.

## `n_man()`: Number of Mannoses

Number of mannoses. This function assumes the Hex of the N-glycan core
is mannoses. Also, for high-mannose and paucimannose glycans, all Hex
are mannoses. Finally, for hybrid glycans, all the rightmost (the side
without branching HexNAc) are mannoses.
