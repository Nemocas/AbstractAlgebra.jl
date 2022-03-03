```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Submodules

AbstractAlgebra allows the construction of submodules/subvector spaces of
AbstractAlgebra modules over euclidean domains. These are given as the
submodule generated by a finite list of elements in the original module.

We define two submodules to be equal if they are (transitively) submodules
of the same module $M$ and their generators generate the same set of elements.

## Generic submodule type

AbstractAlgebra implements a generic submodule type `Generic.Submodule{T}`
where `T` is the element type of the base ring in `src/generic/Submodule.jl`.
See `src/generic/GenericTypes.jl` for more details of the type definition.

Elements of a generic submodule have type `Generic.SubmoduleElem{T}`.

## Abstract types

Submodule types belong to the abstract type `FPModule{T}` and their elements
to `FPModuleElem{T}`.

## Constructors

```@docs
sub(::FPModule{T}, ::Vector{FPModuleElem{T}}) where T <: RingElement
```

```@docs
sub(::FPModule{T}, ::Vector{Generic.Submodule{T}}) where T <: RingElement
```

Note that the preimage of the canonical injection can be obtained using the
preimage function described in the section on module homomorphisms. As the
canonical injection is injective, this is unique.

**Examples**

```jldoctest
julia> M = FreeModule(ZZ, 2)
Free module of rank 2 over Integers

julia> m = M([ZZ(1), ZZ(2)])
(1, 2)

julia> n = M([ZZ(2), ZZ(-1)])
(2, -1)

julia> N, f = sub(M, [m, n])
(Submodule over Integers with 2 generators and no relations
, Module homomorphism with
Domain: Submodule over Integers with 2 generators and no relations

Codomain: Free module of rank 2 over Integers)

julia> v = N([ZZ(3), ZZ(4)])
(3, 4)

julia> v2 = f(v)
(3, 26)

julia> V = VectorSpace(QQ, 2)
Vector space of dimension 2 over Rationals

julia> m = V([QQ(1), QQ(2)])
(1//1, 2//1)

julia> n = V([QQ(2), QQ(-1)])
(2//1, -1//1)

julia> N, f = sub(V, [m, n])
(Subspace over Rationals with 2 generators and no relations
, Module homomorphism with
Domain: Subspace over Rationals with 2 generators and no relations

Codomain: Vector space of dimension 2 over Rationals)

```

## Functionality for submodules

In addition to the Module interface, AbstractAlgebra submodules implement the
following functionality.

### Basic manipulation

```@docs
supermodule(::Generic.Submodule{T}) where T <: RingElement
```

```@docs
issubmodule(::FPModule{T}, ::FPModule{T}) where T <: RingElement
```

```@docs
iscompatible(::FPModule{T}, ::FPModule{T}) where T <: RingElement
```

```@docs
dim(N::Generic.Submodule{T}) where T <: FieldElement
```

**Examples**

```jldoctest
julia> M = FreeModule(ZZ, 2)
Free module of rank 2 over Integers

julia> m = M([ZZ(2), ZZ(3)])
(2, 3)

julia> n = M([ZZ(1), ZZ(4)])
(1, 4)

julia> N1, = sub(M, [m, n])
(Submodule over Integers with 2 generators and no relations
, Module homomorphism with
Domain: Submodule over Integers with 2 generators and no relations

Codomain: Free module of rank 2 over Integers)

julia> N2, = sub(M, [m])
(Submodule over Integers with 1 generator and no relations
, Module homomorphism with
Domain: Submodule over Integers with 1 generator and no relations

Codomain: Free module of rank 2 over Integers)

julia> supermodule(N1) == M
true

julia> iscompatible(N1, N2)
(true, Free module of rank 2 over Integers)

julia> issubmodule(N1, M)
false


julia> V = VectorSpace(QQ, 2)
Vector space of dimension 2 over Rationals

julia> m = V([QQ(2), QQ(3)])
(2//1, 3//1)

julia> N, = sub(V, [m])
(Subspace over Rationals with 1 generator and no relations
, Module homomorphism with
Domain: Subspace over Rationals with 1 generator and no relations

Codomain: Vector space of dimension 2 over Rationals)

julia> dim(V)
2

julia> dim(N)
1

```

### Intersection

```@docs
Base.intersect(M::FPModule{T}, N::FPModule{T}) where
T <: RingElement
```

**Examples**

```jldoctest
julia> M = FreeModule(ZZ, 2)
Free module of rank 2 over Integers

julia> m = M([ZZ(2), ZZ(3)])
(2, 3)

julia> n = M([ZZ(1), ZZ(4)])
(1, 4)

julia> N1 = sub(M, [m, n])
(Submodule over Integers with 2 generators and no relations
, Module homomorphism with
Domain: Submodule over Integers with 2 generators and no relations

Codomain: Free module of rank 2 over Integers)

julia> N2 = sub(M, [m])
(Submodule over Integers with 1 generator and no relations
, Module homomorphism with
Domain: Submodule over Integers with 1 generator and no relations

Codomain: Free module of rank 2 over Integers)

julia> I = intersect(N1, N2)
Union{AbstractAlgebra.Generic.ModuleHomomorphism{BigInt}, AbstractAlgebra.Generic.Submodule{BigInt}}[]
```