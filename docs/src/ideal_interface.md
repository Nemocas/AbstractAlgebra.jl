```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# Ideal Interface

AbstractAlgebra.jl generic code makes use of a standardised set of functions which it
expects to be implemented by anyone implementing ideals for AbstractAlgebra rings. 
Here we document this interface. All libraries which want to make use of the generic
capabilities of AbstractAlgebra.jl must supply all of the required functionality for their ideals.
There are already many helper methods in AbstractAlgebra.jl for the methods mentioned below.

In addition to the required functions, there are also optional functions which can be
provided for certain types of ideals e.g., for ideals of polynomial rings. If implemented,
these allow the generic code to provide additional functionality for those ideals, or in
some cases, to select more efficient algorithms.

## Types and parents

New ideal types should come with the following type information:

```julia
ideal_type(::Type{NewRing}) = NewIdealType 
base_ring_type(::Type{NewIdeal}) = NewRingType
parent_type(::Type{NewIdeal{T}}) = DefaultIdealSet{T}
```

However, new implementations of ideals needn't necessarily supply new types and could just extend
the existing functionality for new rings as AbstractAlgebra.jl provides a generic ideal type
based on Julia arrays which is implemented in `src/generic/Ideal.jl`. For information 
about implementing new rings, see the [Ring interface](@ref "Ring Interface").

## Required functionality for ideals

In the following, we list all the functions that are required to be provided for ideals
in AbstractAlgebra.jl or by external libraries wanting to use AbstractAlgebra.jl.

We give this interface for fictitious type `NewIdeal` and `Ring` or `NewRing` for the type of the base ring
object `R`, and `RingElem` for the type of the elements of the ring.
We assume that the function

```julia
ideal(R::Ring, xs::Vector{U})
```

with `U === elem_type(Ring)` and `xs` a list of generators,
is implemented by anyone implementing ideals for AbstractAlgebra rings. 
Additionally, the following constructors are already implemented generically:

```julia
ideal(R::Ring, x::U)
ideal(xs::Vector{U}) = ideal(parent(xs[1]), xs)
ideal(x::U) = ideal(parent(x), x)
*(x::RingElem, R::Ring)
*(R::Ring, x::RingElem)
```

An implementation of an Ideal subtype should also provide the
following methods:

```julia
base_ring(I::NewIdeal)
```
```julia
gen(I::NewIdeal, k::Int)
```
```julia
gens(I::NewIdeal)
```
```julia
ngens(I::NewIdeal)
```

## Optional functionality for ideals

Some functionality is difficult or impossible to implement for all ideals.
If it is provided, additional functionality or performance may become available. Here
is a list of all functions that are considered optional and can't be relied on by
generic functions in the AbstractAlgebra Ideal interface.

It may be that no algorithm, or no efficient algorithm is known to implement these
functions. As these functions are optional, they do not need to exist. Julia will
already inform the user that the function has not been implemented if it is called but
doesn't exist.

```julia
in(v::RingElem, I::NewIdeal)
```
```julia
issubset(I::NewIdeal, J::NewIdeal)
```
```julia
iszero(I::NewIdeal)
```
```julia
+(I::T, J::T) where {T <: NewIdeal}
```
```julia
*(I::T, J::T) where {T <: NewIdeal}
```
```julia
intersect(I::T, J::T) where {T <: NewIdeal}
```
```julia
==(I::T, J::T) where {T <: NewIdeal}
```
