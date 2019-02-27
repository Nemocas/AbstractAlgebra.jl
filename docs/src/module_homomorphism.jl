## Module Homomorphisms

Homomorphisms of AbstractAlgebra modules, $f : R^s \to R^t$, can be represented by
$s\times t$ matrices over $R$.

```@docs
ModuleHomomorphism(M1::AbstractAlgebra.Module{T}, M2::AbstractAlgebra.Module{T}, m::AbstractA
lgebra.MatElem{T}) where T <: Union{RingElement, NCRingElem}
```

**Examples**

```julia
M = FreeModule(ZZ, 2)
f = ModuleHomomorphism(M, M, matrix(ZZ, 2, 2, [1, 2, 3, 4]))

m = M([ZZ(1), ZZ(2)])

f(m)
```

