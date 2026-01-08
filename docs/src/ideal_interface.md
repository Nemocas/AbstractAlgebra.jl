```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# Ideal Interface

AbstractAlgebra.jl generic code makes use of a standardised set of functions which it
expects to be implemented by anyone implementing ideals for commutative rings.
Here we document this interface. All libraries which want to make use of the generic
capabilities of AbstractAlgebra.jl must supply all of the required functionality for their ideals.
There are already many helper methods in AbstractAlgebra.jl for the methods mentioned below.

In addition to the required functions, there are also optional functions which can be
provided for certain types of ideals e.g., for ideals of polynomial rings. If implemented,
these allow the generic code to provide additional functionality for those ideals, or in
some cases, to select more efficient algorithms.


## Types and parents

Below we describe this interface for a fictitious type `NewIdeal` representing
ideals over a base ring of type `NewRing`, with element type `NewRingElem`. To
make use of the functionality described on this page, `NewIdeal` must be a
subtype of `Ideal{NewRingElem}`. To inform the system about this relationship,
it is necessary to provide the following method:
```julia
ideal_type(::Type{NewRing}) = NewIdeal
```julia
The system automatically provides the following reverse method:
```
base_ring_type(::Type{NewIdeal}) = NewRing
```

For ideals of a Euclidean domain, it may also be possibly to opt into using the existing
functionality  which is implemented in `src/generic/Ideal.jl`. In that case you would
essentially defined `const NewIdeal = Generic.Ideal{NewRingElem}`. For more information
about implementing new rings, see the [Ring interface](@ref "Ring Interface").

## Required functionality for ideals

In the following, we list all the functions that are required to be provided for ideals
in AbstractAlgebra.jl or by external libraries wanting to use AbstractAlgebra.jl.

To facilitate construction of new ideals, implementations must provide a method with signature
```julia
ideal(R::NewRing, xs::Vector{NewRingElem})
```
Here `xs` is a list of generators, and `NewRingElem === elem_type(NewRing)` holds.

With this in place, the following additional ideal constructors will automatically work via
generic implementations:
```julia
ideal(R::NewRing, x::RingElement...) = ideal(R, [x...])
ideal(x::RingElement, y::RingElement...) = ideal(parent(x), x, y...)
ideal(xs::Vector{NewRingElem}) = ideal(parent(xs[1]), xs)
*(x::NewRingElem, R::NewRing) = ideal(R, x)
*(R::NewRing, x::NewRingElem) = ideal(R, x)
```
In addition sums and products of ideals can be formed:
```julia
+(I::T, J::T) where {T <: NewIdeal}
*(I::T, J::T) where {T <: NewIdeal}
```

An implementation of an `Ideal` subtype must also provide the following methods:
```julia
base_ring(I::NewIdeal)
gen(I::NewIdeal, k::Int)
gens(I::NewIdeal)
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

The following method have no generic implementation and only work when explicitly
implemented.
```julia
in(v::NewRingElem, I::NewIdeal)
intersect(I::T, J::T) where {T <: NewIdeal}
```

If a method for `in` as above is provided, then the following automatically works:
```julia
issubset(I::NewIdeal, J::NewIdeal)
```

If a method for `in` as above is provided (e.g. indirectly by providing method for `in`),
then the following automatically works:
```julia
==(I::T, J::T) where {T <: NewIdeal}
```
Note that implementing `==` for a Julia type means that we have to provide
a matching `hash` method which preserves the invariant that `I == J` implies `hash(I) == hash(J)`.
We provide such a method but by necessity it is very conservative and hence does
not provide good hashing. You may wish to implement a better `hash` methods.

The following method is implemented generically via the ideal generators.
```julia
iszero(I::Ideal) = all(iszero, gens(I))
```
