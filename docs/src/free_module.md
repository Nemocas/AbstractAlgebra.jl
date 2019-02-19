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

