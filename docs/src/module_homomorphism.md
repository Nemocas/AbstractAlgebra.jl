```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Module Homomorphisms

Abstract Algebra provides homomorphisms of finitely presented modules.

## Generic functionality

The following generic functionality is provided for module homomorphisms.

### Constructors

Homomorphisms of AbstractAlgebra modules, $f : R^s \to R^t$, can be represented by
$s\times t$ matrices over $R$.

```@docs
Generic.ModuleHomomorphism(M1::AbstractAlgebra.FPModule{T}, M2::AbstractAlgebra.FPModule{T}, m::AbstractAlgebra.MatElem{T}) where T <: RingElement
```

```@docs
Generic.ModuleIsomorphism(M1::AbstractAlgebra.FPModule{T}, M2::AbstractAlgebra.FPModule{
T}, m::AbstractAlgebra.MatElem{T}) where T <: RingElement
```

**Examples**

```jldoctest
julia> M = FreeModule(ZZ, 2)
Free module of rank 2 over Integers

julia> f = ModuleHomomorphism(M, M, matrix(ZZ, 2, 2, [1, 2, 3, 4]))
Module homomorphism with
Domain: Free module of rank 2 over Integers
Codomain: Free module of rank 2 over Integers

julia> m = M([ZZ(1), ZZ(2)])
(1, 2)

julia> f(m)
(7, 10)

```

### Kernels

```@docs
kernel(f::Map(AbstractAlgebra.FPModuleHomomorphism))
```

**Examples**

```jldoctest
julia> M = FreeModule(ZZ, 3)
Free module of rank 3 over Integers

julia> m = M([ZZ(1), ZZ(2), ZZ(3)])
(1, 2, 3)

julia> S, f = sub(M, [m])
(Submodule over Integers with 1 generator and no relations
, Module homomorphism with
Domain: Submodule over Integers with 1 generator and no relations

Codomain: Free module of rank 3 over Integers)

julia> Q, g = quo(M, S)
(Quotient module over Integers with 2 generators and no relations
, Module homomorphism with
Domain: Free module of rank 3 over Integers
Codomain: Quotient module over Integers with 2 generators and no relations
)

julia> kernel(g)
(Submodule over Integers with 1 generator and no relations
, Module homomorphism with
Domain: Submodule over Integers with 1 generator and no relations

Codomain: Free module of rank 3 over Integers)

```

### Images

```@docs
image(::Map(AbstractAlgebra.FPModuleHomomorphism))
```

```julia
M = FreeModule(ZZ, 3)

m = M([ZZ(1), ZZ(2), ZZ(3)])

S, f = sub(M, [m])
Q, g = quo(M, S)
K, k = kernel(g)

image(compose(k, g))
```

### Preimages

```@docs
preimage(::Map(AbstractAlgebra.FPModuleHomomorphism), ::AbstractAlgebra.FPModuleElem{T}) where T <: RingElement
```

```julia
M = FreeModule(ZZ, 3)

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
inv(::Map(Generic.ModuleIsomorphism))
```

```julia
M = FreeModule(ZZ, 2)
N = matrix(ZZ, 2, 2, BigInt[1, 0, 0, 1])
f = ModuleIsomorphism(M, M, N)

g = inv(f)
```
