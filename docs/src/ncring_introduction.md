```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```
# Introduction

Noncommutative rings work mostly like the commutative rings described in
[Ring functionality](@ref), with a few differences described here.

Their parents belong to the abstract type `NCRing` and their elements to
`NCRingElem`. As `Ring <: NCRing` and `RingElem <: NCRingElem`, these types
cover all rings, commutative or not; "noncommutative" means "not necessarily
commutative". The union type `NCRingElement` adds Julia's number types to
`NCRingElem`, as `RingElement` does for `RingElem`.

AbstractAlgebra provides univariate polynomials over a noncommutative ring,
free associative algebras and matrix algebras, each described on its own page
in this section.

## Exact division

As `g*a` and `a*g` may differ, exact division comes in a left and a right
version; `divexact` is only available for commutative rings.

```@docs
divexact_left
divexact_right
```
