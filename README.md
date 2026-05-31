# GEIO

Research-alpha R package for local equilibrium-index diagnostics in
small parameterized finite-game systems.

GEIO attempts to recover equilibria of F_theta(x) = 0 from multiple
starts. For each recovered regular equilibrium, it computes
sign(det(D_x F_theta(x_star))) and a singular-value near-singularity
diagnostic.

## Scope

GEIO tracks how equilibrium structure changes as game parameters move —
which equilibria persist, which collapse, and where the system approaches
singularity. Built-in examples cover two-action coordination games
(Stag Hunt, Battle of the Sexes, general bimatrix). Custom equilibrium
systems are also supported. Exhaustive equilibrium enumeration is not
guaranteed for arbitrary games; the index-sum consistency check applies
to selected generic built-in examples.

## Development disclosure

Development used AI coding assistance (Anthropic Claude, OpenAI GPT)
under continuous human direction. All algorithmic design, methodology
selection, and validation decisions were made by the human author.

## Installation

```r
install.packages("GEIO")
```

Requires GALAHAD (>= 2.0.0).
