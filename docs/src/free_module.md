# Free Modules

AbstractAlgebra allows the construction of the free module of any rank over any
ring. By default the system considers the free module of a given rank over a
given ring to be unique.

## Types and parents

AbstractAlgebra provides the type `FreeModule{T}` for free modules, where `T`
is the type of the elements of the ring $R$ over which the module is built.
The type `FreeModule{T}` belongs to `AbstractAlgebra.Module{T}`.

The free module of a given rank over a given ring is made unique on the
system by caching them (unless an optional `cache` parameter is set to
`false`).

See `src/generic/GenericTypes.jl` for an example of how to implement such a
cache (which usually makes use of a dictionary).

## Functionality for free modules

### Basic manipulation

```@docs
rank{T <: RingElem}(M::FreeModule{T})
```

**Examples**

```julia
M = FreeModule(ZZ, 3)

rank(M)
```

### Morphisms

Morphisms of free modules, $f : R^s \to R^t$, can be represented by $s\times t$
matrices over $R$. Note that elements of a free module are interpreted as row
vectors in this context and if $f_1$ and $f_2$ are composable morphisms then
$f_1(f_2(v))$ corresponds to multiplying the row vector $v$ on the right by 
$M_1\times M_2$ where $M_i$ is the matrix corresponding to the morphism $f_i$.

```@docs
FreeModuleMorphism(M1::FreeModule{T}, M2::FreeModule{T}, m::AbstractAlgebra.MatElem{T}) where T <: Union{RingElement, NCRingElem}
```

**Examples**

```julia
M = FreeModule(ZZ, 2)
f = FreeModuleMorphism(M, M, matrix(ZZ, 2, 2, [1, 2, 3, 4]))

m = M([ZZ(1), ZZ(2)])

f(m)
```
 


