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

```julia
M = FreeModule(ZZ, 2)

m = M([ZZ(1), ZZ(2)])

N, f = Submodule(M, [m])
Q, g = QuotientModule(M, N)

p = Q([ZZ(3)])
v2 = g(p)

V = VectorSpace(QQ, 2)

m = V([QQ(1), QQ(2)])

N, f = Subspace(V, [m])
Q, g = QuotientSpace(V, N)
```

## Functionality for submodules

In addition to the Module interface, AbstractAlgebra submodules implement the
following functionality.

### Basic manipulation

```@docs
supermodule(M::AbstractAlgebra.Generic.QuotientModule{T}) where T <: RingElement
```

**Examples**

```julia
M = FreeModule(ZZ, 2)
m = M([ZZ(2), ZZ(3)])
N, g = Submodule(M, [m])
Q, h = QuotientModule(M, N)

supermodule(Q) == M
```

