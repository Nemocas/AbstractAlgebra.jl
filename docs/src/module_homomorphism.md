```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# Module Homomorphisms

Abstract Algebra provides homomorphisms of finitely presented modules.

## Generic module homomorphism types

AbstractAlgebra defines two module homomorphism types, namely
`Generic.ModuleHomomorphism` and `Generic.ModuleIsomorphism`. Functionality
for these is implemented in `src/generic/ModuleHomomorphism.jl`.

## Abstract types

The `Generic.ModuleHomomorphism` and `Generic.ModuleIsomorphism` types inherit
from `Map(FPModuleHomomorphism)`.

## Generic functionality

The following generic functionality is provided for module homomorphisms.

### Constructors

Homomorphisms of AbstractAlgebra modules, $f : R^s \to R^t$, can be represented by
$s\times t$ matrices over $R$.

```@docs
ModuleHomomorphism(M1::FPModule{T}, M2::FPModule{T}, m::MatElem{T}) where T <: RingElement
```

```@docs
ModuleIsomorphism(M1::FPModule{T}, M2::FPModule{T}, m::MatElem{T}) where T <: RingElement
```

**Examples**

```jldoctest
julia> M = free_module(ZZ, 2)
Free module of rank 2 over integers

julia> f = ModuleHomomorphism(M, M, matrix(ZZ, 2, 2, [1, 2, 3, 4]))
Module homomorphism
  from free module of rank 2 over integers
  to free module of rank 2 over integers

julia> m = M([ZZ(1), ZZ(2)])
(1, 2)

julia> f(m)
(7, 10)

```

They can also be created by giving images (in the codomain) of the generators of the domain:
```@doc
ModuleHomomorphism(M1::FPModule{T}, M2::FPModule{T}, v::Vector{<:FPModuleElem{T}}) where T <: RingElement
```

### Kernels

```@docs
kernel(f::Map(FPModuleHomomorphism))
```

**Examples**

```jldoctest
julia> M = free_module(ZZ, 3)
Free module of rank 3 over integers

julia> m = M([ZZ(1), ZZ(2), ZZ(3)])
(1, 2, 3)

julia> S, f = sub(M, [m])
(Submodule over integers with 1 generator and no relations, Hom: S -> M)

julia> Q, g = quo(M, S)
(Quotient module over integers with 2 generators and no relations, Hom: M -> Q)

julia> kernel(g)
(Submodule over integers with 1 generator and no relations, Hom: submodule over integers with 1 generator and no relations -> M)

```

### Images

```@docs
image(::Map(FPModuleHomomorphism))
```

```julia
M = free_module(ZZ, 3)

m = M([ZZ(1), ZZ(2), ZZ(3)])

S, f = sub(M, [m])
Q, g = quo(M, S)
K, k = kernel(g)

image(compose(k, g))
```

### Preimages

```@docs
preimage(::Map(FPModuleHomomorphism), ::FPModuleElem{T}) where T <: RingElement
has_preimage_with_preimage(::Map(FPModuleHomomorphism), ::FPModuleElem{T}) where T <: RingElement
```

```julia
M = free_module(ZZ, 3)

m = M([ZZ(1), ZZ(2), ZZ(3)])

S, f = sub(M, [m])
Q, g = quo(M, S)

m = rand(M, -10:10)
n = g(m)

p = preimage(g, n)
```

### Inverses

Module isomorphisms can be cheaply inverted.

```@docs
Base.inv(::Map(Generic.ModuleIsomorphism))
```

```julia
M = free_module(ZZ, 2)
N = matrix(ZZ, 2, 2, BigInt[1, 0, 0, 1])
f = ModuleIsomorphism(M, M, N)

g = inv(f)
```
