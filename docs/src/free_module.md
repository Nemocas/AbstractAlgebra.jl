# Free Modules and Vector Spaces

AbstractAlgebra allows the construction of the free module of any rank over any
Euclidean ring and the vector space of any dimension over a field. By default
the system considers the free module of a given rank over a given ring or
vector space of given dimension over a field to be unique.

## Types and parents

AbstractAlgebra provides the type `FreeModule{T}` for free modules, where `T`
is the type of the elements of the ring $R$ over which the module is built.
The type `FreeModule{T}` belongs to `AbstractAlgebra.FPModule{T}`.

Vector spaces are simply free modules over a field.

The free module of a given rank over a given ring is made unique on the
system by caching them (unless an optional `cache` parameter is set to
`false`).

See `src/generic/GenericTypes.jl` for an example of how to implement such a
cache (which usually makes use of a dictionary).

## Functionality for free modules

As well as implementing the entire module interface, free modules provide the
following functionality.

### Constructors

```@docs
FreeModule(R::AbstractAlgebra.Ring, rank::Int)
VectorSpace(F::AbstractAlgebra.Field, dim::Int)
```

Construct the free module/vector space of given rank/dimension.

**Examples**

```julia
M = FreeModule(ZZ, 3)
V = VectorSpace(QQ, 2)
```

### Basic manipulation

```julia
rank(M::Generic.FreeModule{T}) where T <: AbstractAlgebra.RingElem
dim(V::Generic.FreeModule{T}) where T <: AbstractAlgebra.FieldElem
```

**Examples**

```julia
M = FreeModule(ZZ, 3)
V = VectorSpace(QQ, 2)

rank(M)
dim(V)
```


