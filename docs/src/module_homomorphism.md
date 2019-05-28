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

**Examples**

```julia
M = FreeModule(ZZ, 2)
f = ModuleHomomorphism(M, M, matrix(ZZ, 2, 2, [1, 2, 3, 4]))

m = M([ZZ(1), ZZ(2)])

f(m)
```

### Kernels

```@docs
kernel(f::Generic.ModuleHomomorphism{T}) where T <: RingElement
```

**Examples**

```julia
M = FreeModule(ZZ, 3)

m = M([ZZ(1), ZZ(2), ZZ(3)])

S, f = Submodule(M, [m])
Q, g = QuotientModule(M, S)

kernel(g)
```

### Images

```@docs
image(f::Generic.ModuleHomomorphism{T}) where T <: RingElement
```

```julia
M = FreeModule(ZZ, 3)

m = M([ZZ(1), ZZ(2), ZZ(3)])

S, f = Submodule(M, [m])
Q, g = QuotientModule(M, S)
K, k = kernel(g)

image(compose(k, g))
```
