```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Quotient modules

AbstractAlgebra allows the construction of quotient modules/spaces of
AbstractAlgebra modules over euclidean domains. These are given as the
quotient of a module by a submodule of that module.

We define two quotient modules to be equal if they are quotients of the
same module $M$ by two equal submodules.

As well as implementing the entire Module interface, AbstractAlgebra
submodules also provide the following interface.

## Constructors

```@docs
QuotientModule(M::AbstractAlgebra.FPModule{T}, v::Generic.Submodule{T}) where T <: RingElement
```

```@docs
QuotientSpace(M::AbstractAlgebra.FPModule{T}, v::Generic.Submodule{T}) where T <: FieldElement

```

**Examples**

```jldoctest
julia> M = FreeModule(ZZ, 2)

julia> m = M([ZZ(1), ZZ(2)])

julia> N, f = Submodule(M, [m])

julia> Q, g = QuotientModule(M, N)

julia> p = Q([ZZ(3)])

julia> v2 = g(p)

julia> V = VectorSpace(QQ, 2)

julia> m = V([QQ(1), QQ(2)])

julia> N, f = Subspace(V, [m])

julia> Q, g = QuotientSpace(V, N)

```

## Functionality for submodules

In addition to the Module interface, AbstractAlgebra submodules implement the
following functionality.

### Basic manipulation

```@docs
supermodule(M::AbstractAlgebra.Generic.QuotientModule{T}) where T <: RingElement
```

**Examples**

```jldoctest
julia> M = FreeModule(ZZ, 2)

julia> m = M([ZZ(2), ZZ(3)])

julia> N, g = Submodule(M, [m])

julia> Q, h = QuotientModule(M, N)

julia> supermodule(Q) == M

```

