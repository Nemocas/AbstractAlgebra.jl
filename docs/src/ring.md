# Ring functionality

AbstractAlgebra has both commutative and noncommutative rings. Together we
refer to them below as rings.

## Abstract types for rings

All commutative ring types in AbstractAlgebra belong to the `Ring` abstract
type and commutative ring elements belong to the `RingElem` abstract type.

Noncommutative ring types belong to the `NCRing` abstract type and their
elements to `NCRingElem`.

As Julia types cannot belong to our `RingElem` type hierarchy, we also
provide the union type `RingElement` which includes `RingElem` in union with
the Julia types `Integer`, `Rational` and `AbstractFloat`.

Similarly `NCRingElement` includes the Julia types just mentioned in union
with `NCRingElem`.

Note that

```julia
Ring <: NCRing
RingElem <: NCRingElem
RingElement <: NCRingElement
```

## Functions for types and parents of rings

```julia
parent_type(::Type{T}) where T <: NCRingElement
elem_type(::Type{T}) where T <: NCRing
```

Return the type of the parent (resp. element) type corresponding to the given
ring element (resp. parent) type.

```julia
base_ring(R::NCRing)
base_ring(a::NCRingElement)
```

For generic ring constructions over a base ring (e.g. polynomials over a
coefficient ring), return the parent object of that base ring.


```julia
parent(a::NCRingElement)
```

Return the parent of the given ring element.

```julia
isdomain_type(::Type{T}) where T <: NCRingElement
isexact_type(::Type{T}) where T <: NCRingElement
```

Return true if the given ring element type can only belong to elements of an
integral domain or exact ring respectively. (An exact ring is one whose
elements are represented exactly in the system without approximation.)

The following function is implemented where mathematically and algorithmically
possible.

```julia
characteristic(R::NCRing)
```

## Constructors

If `R` is a parent object of a ring in AbstractAlgebra, it can always be used
to construct certain objects in that ring.

```julia
(R::NCRing)() # constructs zero
(R::NCRing)(c::Integer)
(R::NCRing)(c::elem_type(R))
(R::NCRing{T})(a::T) where T <: RingElement
```

## Basic functions

All rings in AbstractAlgebra are expected to implement basic ring operations,
unary minus, binary addition, subtraction and multiplication, equality testing,
powering.

In addition, the following are implemented for parents/elements just as they
would be in Julia for types/objects.

```julia
zero(R::NCRing)
one(R::NCRing)
iszero(a::NCRingElement)
isone(a::NCRingElement)
```

In addition, the following is implemented where it is
mathematically/algorithmically viable to do so.

```julia
isunit(a::NCRingElement)
```

The following standard Julia functions are also implemented for all ring
elements.

```julia
hash(f::RingElement, h::UInt)
deepcopy_internal(a::RingElement, dict::ObjectIdDict)
show(io::IO, R::NCRing)
show(io::IO, a::NCRingElement)
```

## Basic functionality for inexact rings only

By default, inexact ring elements in AbstractAlgebra compare equal if they are
the same to the minimum precision of the two elements. However, we also provide
the following more strict notion of equality, which also requires the
precisions to be the same.

```julia
isequal(a::T, b::T) where T <: NCRingElement
```

For floating point and ball arithmetic it is sometimes useful to be able to
check if two elements are approximately equal, e.g. to suppress numerical noise
in comparisons. For this, the following are provided.

```julia
isapprox(a::T, b::T; atol::Real=sqrt(eps())) where T <: RingElement
```

Similarly, for a parameterised ring with type `MyElem{T}` over such an inexact
ring we have the following.

```julia
isapprox(a::MyElem{T}, b::T; atol::Real=sqrt(eps())) where T <: RingElement
isapprox(a::T, b::MyElem{T}; atol::Real=sqrt(eps())) where T <: RingElement
```

These notionally perform a coercion into the parameterised ring before doing
the approximate equality test.

## Basic functionality for commutative rings only

```julia
divexact(a::T, b::T) where T <: RingElement
inv(a::T)
```

Return `a/b` or `1/a` respectively, where the slash here refers to the
mathematical notion of division in the ring, not Julia's floating point
division operator.

## Basic functionality for noncommutative rings only

```julia
divexact_left(a::T, b::T) where T <: NCRingElement
divexact_right(a::T, b::T) where T <: NCRingElement
```

As per `divexact` above, except that division by `b` happens on the left or
right, respectively, of `a`.

## Unsafe ring operators

To speed up polynomial arithmetic, various unsafe operators are provided, which
mutate the output rather than create a new object.

```julia
zero!(a::NCRingElement)
mul!(a::T, b::T, c::T) where T <: NCRingElement
add!(a::T, b::T, c::T) where T <: NCRingElement
addeq!(a::T, b::T) where T <: NCRingElement
addmul!(a::T, b::T, c::T, t::T) where T <: NCRingElement
```

In each case the mutated object is the leftmost parameter.

The `addeq!(a, b)` operation does the same thing as `add!(a, a, b)`. The
optional `addmul!(a, b, c, t)` operation does the same thing as
`mul!(t, b, c); addeq!(a, t)` where `t` is a temporary which can be mutated so
that an addition allocation is not needed.

## Random generation

The Julia random interface is implemented for all ring parents (instead of
for types). The exact interface differs depending on the ring, but the
parameters supplied are usually ranges, e.g. `-1:10` for the range of allowed
degrees for a univariate polynomial.

```julia
rand(R::NCRing, v...)
```

## Factorization

For commutative rings supporting factorization and irreducibility testing, the
following optional functions may be implemented.

```julia
isirreducible(a::T) where T <: RingElement
issquarefree(a::T) where T <: RingElement
```

Decide whether `a` is irreducible or squarefree, respectively.

```julia
factor(a::T) where T <: RingElement
factor_squarefree(a::T) where T <: RingElement
```

Return a factorization into irreducible or squarefree elements, respectively.
The return is an object of type `Fac{T}`.

```@docs
Fac
unit(a::Fac)
evaluate(a::Fac)
getindex(a::Fac, b)
setindex!(a::Fac{Int}, c::Int, b::Int)
```

