```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# Quotient modules

AbstractAlgebra allows the construction of quotient modules/spaces of
AbstractAlgebra modules over euclidean domains. These are given as the
quotient of a module by a submodule of that module.

We define two quotient modules to be equal if they are quotients of the
same module $M$ by two equal submodules.

## Generic quotient module type

AbstractAlgebra implements the generic quotient module type
`Generic.QuotientModule{T}` where `T` is the element type of the base ring,
in `src/generic/QuotientModule.jl`.

Elements of generic quotient modules have type `Generic.QuotientModuleElem{T}`.

## Abstract types

Quotient module types belong to the `FPModule{T}` abstract type and their
elements to `FPModuleElem{T}`.

## Constructors

```@docs
quo(M::FPModule{T}, v::Generic.Submodule{T}) where T <: RingElement
```

Note that a preimage of the canonical projection can be obtained using the
preimage function described in the section on module homomorphisms. Note that
a preimage element of the canonical projection is not unique and has no special
properties.

**Examples**

```jldoctest
julia> M = free_module(ZZ, 2)
Free module of rank 2 over integers

julia> m = M([ZZ(1), ZZ(2)])
(1, 2)

julia> N, f = sub(M, [m])
(Submodule over integers with 1 generator and no relations, Hom: N -> M)

julia> Q, g = quo(M, N)
(Quotient module over integers with 1 generator and no relations, Hom: M -> Q)

julia> p = M([ZZ(3), ZZ(1)])
(3, 1)

julia> v2 = g(p)
(-5)

julia> V = vector_space(QQ, 2)
Vector space of dimension 2 over rationals

julia> m = V([QQ(1), QQ(2)])
(1//1, 2//1)

julia> N, f = sub(V, [m])
(Subspace over rationals with 1 generator and no relations, Hom: N -> V)

julia> Q, g = quo(V, N)
(Quotient space over rationals with 1 generator and no relations, Hom: V -> Q)

```

## Functionality for submodules

In addition to the Module interface, AbstractAlgebra submodules implement the
following functionality.

### Basic manipulation

```@docs
supermodule(M::Generic.QuotientModule{T}) where T <: RingElement

dim(N::Generic.QuotientModule{T}) where T <: FieldElement
```

**Examples**

```jldoctest
julia> M = free_module(ZZ, 2)
Free module of rank 2 over integers

julia> m = M([ZZ(2), ZZ(3)])
(2, 3)

julia> N, g = sub(M, [m])
(Submodule over integers with 1 generator and no relations, Hom: N -> M)

julia> Q, h = quo(M, N)
(Quotient module over integers with 2 generators and relations:
[2 3], Hom: M -> Q)

julia> supermodule(Q) == M
true

julia> V = vector_space(QQ, 2)
Vector space of dimension 2 over rationals

julia> m = V([QQ(1), QQ(2)])
(1//1, 2//1)

julia> N, f = sub(V, [m])
(Subspace over rationals with 1 generator and no relations, Hom: N -> V)

julia> Q, g = quo(V, N)
(Quotient space over rationals with 1 generator and no relations, Hom: V -> Q)

julia> dim(V)
2

julia> dim(Q)
1

```

