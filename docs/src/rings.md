# Ring Interface

AbstractAlgebra.jl generic code makes use of a standardised set of functions which it
expects to be implemented for all rings. Here we document this interface. All libraries
which want to make use of the generic capabilities of AbstractAlgebra.jl must supply
all of the required functionality for their rings.

In addition to the required functions, there are also optional functions which can be
provided for certain types of rings, e.g. GCD domains or fields, etc. If implemented,
these allow the generic code to provide additional functionality for those rings, or in
some cases, to select more efficient algorithms.

## Types

Most rings must supply two types:
  - a type for the parent object (representing the ring itself)
  - a type for elements of that ring

For example, the generic univariate polynomial type in AbstractAlgebra.jl provides two
types in generic/GenericTypes.jl:

  - `Generic.PolyRing{T}` for the parent objects
  - `Generic.Poly{T}` for the actual polynomials

The parent type must belong to `AbstractAlgebra.Ring` and the element type must belong
to `AbstractAlgebra.RingElem`. Of course, the types may belong to these abstract types
transitively, e.g. `Poly{T}` actually belongs to `AbstractAlgebra.PolyElem{T}` which in
turn belongs to `AbstractAlgebra.RingElem`.

For parameterised rings, we advise that the types of both the parent objects and
element objects to be parameterised by the types of the elements of the base ring
(see the function `base_ring` below for a definition).

There can be variations on this theme: e.g. in some areas of mathematics there is a
notion of a coefficient domain, in which case it may make sense to parameterise all
types by the type of elements of this coefficient domain. But note that this may have
implications for the ad hoc operators one might like to explicitly implement.

## RingElement type union

Because of its lack of multiple inheritance, Julia does not allow Julia Base
types to belong to `AbstractAlgebra.RingElem`. To allow us to work equally with
AbstractAlgebra and Julia types that represent elements of rings we define a
union type `AbstractAlgebra.RingElement` in `src/julia/JuliaTypes`.

So far, in addition to `AbstractAlgebra.RingElem` the  union type
`AbstractAlgebra.RingElement` includes the Julia types `Integer`, `Rational`
and `AbstractFloat`.

Most of the generic code in AbstractAlgebra makes use of the union type
`AbstractAlgebra.RingElement` instead of `AbstractAlgebra.RingElem` so that the
generic functions also accept the Julia Base ring types.

One must be careful when defining ad hoc binary operations for ring element
types. It is often necessary to define separate versions of the functions for
`AbstractAlgebra.RingElem` then for each of the Julia types separately in
order to avoid ambiguity warnings.

Note that even though `AbstractAlgebra.RingElement` is a union type we still
have the following inclusion

```julia
AbstractAlgebra.RingElement <: AbstractAlgebra.NCRingElement
```

## Parent object caches

In many cases, it is desirable to have only one object in the system to represent each
ring. This means that if the same ring is constructed twice, elements of the two rings
will be compatible as far as arithmetic is concerned.

In order to facilitate this, global caches of rings are stored in AbstractAlgebra.jl,
usually implemented using dictionaries. For example, the `Generic.PolyRing` parent
objects are looked up in a dictionary `PolyID` to see if they have been previously
defined.

Whether these global caches are provided or not, depends on both mathematical and
algorithmic considerations. E.g. in the case of number fields, it isn't desirable to
identify all number fields with the same defining polynomial, as they may be considered
with distinct embeddings into one another. In other cases, identifying whether two rings
are the same may be prohibitively expensive. Generally, it may only make sense
algorithmically to identify two rings if they were constructed from identical data.

If a global cache is provided, it must be optionally possible to construct the parent
objects without caching. This is done by passing a boolean value `cached` to the inner
constructor of the parent object. See `generic/GenericTypes.jl` for examples of how to
construct and handle such caches.

## Required functions for all rings

In the following, we list all the functions that are required to be provided for rings
in AbstractAlgebra.jl or by external libraries wanting to use AbstractAlgebra.jl.

We give this interface for fictitious types `MyParent` for the type of the ring parent
object `R` and `MyElem` for the type of the elements of the ring.

Note that generic functions in AbstractAlgebra.jl may not rely on the existence of
functions that are not documented here. If they do, those functions will only be
available for rings that implement that additional functionality, and should be
documented as such.

### Data type and parent object methods

```julia
parent_type(::Type{MyElem})
```

Return the type of the corresponding parent object for the given element type. For
example, `parent_type(Generic.Poly{T})` will return `Generic.PolyRing{T}`.

```julia
elem_type(::Type{MyParent})
```

Return the type of the elements of the ring whose parent object has the given type.
This is the inverse of the `parent_type` function, i.e. `elem_type(Generic.PolyRing{T})`
will return `Generic.Poly{T}`.

```julia
base_ring(R::MyParent)
```

Given a parent object `R`, representing a ring, this function returns the parent object
of any base ring that parameterises this ring. For example, the base ring of the ring
of polynomials over the integers would be the integer ring.

If the ring is not parameterised by another ring, this function must return `Union{}`.

Note that there is a distinction between a base ring and other kinds of parameters. For
example, in the ring $\mathbb{Z}/n\mathbb{Z}$, the modulus $n$ is a parameter, but the
only base ring is $\mathbb{Z}$. We consider the ring $\mathbb{Z}/n\mathbb{Z}$ to have
been constructed from the base ring $\mathbb{Z}$ by taking its quotient by a (principal)
ideal.

```julia
parent(f::MyElem)
```

Return the parent object of the given element, i.e. return the ring to which the given
element belongs.

This is usually stored in a field `parent` in each ring element. (If the parent objects
have `mutable struct` types, the internal overhead here is just an additional machine
pointer stored in each element of the ring.)

For some element types it isn't necessary to append the parent object as a field of
every element. This is the case when the parent object can be reconstructed just given
the type of the elements. For example, this is the case for the ring of integers and
in fact for any ring element type that isn't parameterised or generic in any way.

```julia
isdomain_type(::Type{MyElem})
```

Return `true` if every element of the given element type (which may be parameterised
or an abstract type) necessarily has a parent that is an integral domain, otherwise
if this cannot be guaranteed, the function returns `false`.

For example, if `MyElem` was the type of elements of generic residue rings of a
polynomial ring, the answer to the question would depend on the modulus of the residue
ring. Therefore `isdomain_type` would have to return `false`, since we cannot guarantee
that we are dealing with elements of an integral domain in general. But if the given
element type was for rational integers, the answer would be `true`, since every rational
integer has as parent the ring of rational integers, which is an integral domain.

Note that this function depends only on the type of an element and cannot access
information about the object itself, or its parent.

```julia
isexact_type(::Type{MyElem})
```

Return `true` if every element of the given type is represented exactly. For example,
$p$-adic numbers, real and complex floating point numbers and power series are not
exact, as we can only represent them in general with finite truncations. Similarly
polynomials and matrices over inexact element types are themselves inexact.

Integers, rationals, finite fields and polynomials and matrices over them are always
exact.

Note that `MyElem` may be parameterised or an abstract type, in which case every
element of every type represented by `MyElem` must be exact, otherwise the function
must return `false`.

```julia
Base.hash(f::MyElem, h::UInt)
```

Return a hash for the object $f$ of type `UInt`. This is used as a hopefully cheap way
to distinguish objects that differ arithmetically.

If the object has components, e.g. the coefficients of a polynomial or elements of a
matrix, these should be hashed recursively, passing the same parameter `h` to all
levels. Each component should then be xor'd with `h` before combining the individual
component hashes to give the final hash.

The hash functions in AbstractAlgebra.jl usually start from some fixed 64 bit
hexadecimal  value that has been picked at random by the library author for that type.
That is then truncated to fit a `UInt` (in case the latter is not 64 bits). This ensures
that objects that are the same arithmetically (or that have the same components), but
have different types (or structures), are unlikely to hash to the same value.

```julia
deepcopy_internal(f::MyElem, dict::ObjectIdDict)
```

Return a copy of the given element, recursively copying all components of the object.

Obviously the parent, if it is stored in the element, should not be copied. The new
element should have precisely the same parent as the old object.

For types that cannot self-reference themselves anywhere internally, the `dict` argument
may be ignored.

In the case that internal self-references are possible, please consult the Julia
documentation on how to implement `deepcopy_internal`.

### Constructors

Outer constructors for most AbstractAlgebra types are provided by overloading the call
syntax for parent objects.

If `R` is a parent object for a given ring we provide the following constructors.

```julia
(R::MyParent)()
```

Return the zero object of the given ring.

```julia
(R::MyParent)(a::Integer)
```

Coerce the given integer into the given ring.

```julia
(R::MyParent)(a::MyElem)
```

If $a$ belongs to the given ring, the function returns it (without making a copy).
Otherwise an error is thrown.

For parameterised rings we also require a function to coerce from the base ring into
the parent ring.

```julia
(R::MyParent{T})(a::T) where T <: AbstractAlgebra.RingElem
```

Coerce $a$ into the ring $R$ if $a$ belongs to the base ring of $R$.

### Basic manipulation of rings and elements

```julia
zero(R::MyParent)
```

Return the zero element of the given ring.

```julia
one(R::MyParent)
```

Return the multiplicative identity of the given ring.

```julia
iszero(f::MyElem)
```

Return `true` if the given element is the zero element of the ring it belongs to.

```julia
isone(f::MyElem)
```

Return `true` if the given element is the multiplicative identity of the ring it belongs
to.

### Canonicalisation

```julia
canonical_unit(f::MyElem)
```

When fractions are created with two elements of the given type, it is nice to be able
to represent them in some kind of canonical form. This is of course not always possible.
But for example, fractions of integers can be canonicalised by first removing any common
factors of the numerator and denominator, then making the denominator positive.

In AbstractAlgebra.jl, the denominator would be made positive by dividing both the
numerator and denominator by the canonical unit of the denominator. For a negative
denominator, this would be $-1$.

For elements of a field, `canonical_unit` simply returns the element itself. In general,
`canonical_unit` of an invertible element should be that element. Finally, if $a = ub$
we should have the identity `canonical_unit(a) = canonical_unit(u)*canonical_unit(b)`.

For some rings, it is completely impractical to implement this function, in which case
it may return $1$ in the given ring. The function must however always exist, and always
return an element of the ring.

### String I/O

```julia
show(io::IO, R::MyParent)
```

This should print an English description of the parent ring (to the given IO object).
If the ring is parameterised, it can call the corresponding `show` function for any
rings it depends on.

```julia
show(io::IO, f::MyElem)
```

This should print a human readable, textual representation of the object (to the given
IO object). It can recursively call the corresponding `show` functions for any of its
components.

!!! note

    The functionality of the function `needs_parentheses` has been replaced by
    `expressify` and `needs_parentheses` will be removed in future versions.

It may be necessary in some cases to print parentheses around components of $f$ or to
print signs of components. For these, the following functions will exist for each
component or component type.

```julia
needs_parentheses(f::MyElem)
```

Should return `true` if parentheses are needed around this object when printed, e.g. as
a coefficient of a polynomial. As an example, non-constant polynomials would need such
parentheses if used as coefficients of another polynomial.

Note that since this approach quickly leads to unnecessary parentheses, the
expression method below is preferred.

### Expressions

To obtain best results when printing composed types derived from other types, e.g., polynomials,
the following method should be implemented.

```julia
expressify(f::MyElem; context = nothing)
```

which must return either `Expr`, `Symbol`, `Integer` or `String`. In case one
implements `expressify`, one can define the following show methods for `MyElem`:

```julia
function Base.show(io::IO, a::MyElem)
  show_via_expressify(io, a)
end

function Base.show(io::IO, mi::MIME"text/plain", a::MyElem)
  show_via_expressify(io, mi, a)
end

function Base.show(io::IO, mi::MIME"text/latex", a::MyElem)
  show_via_expressify(io, mi, a)
end

function Base.show(io::IO, mi::MIME"text/html", a::MyElem)
  show_via_expressify(io, mi, a)
end
```

As an example, assume that an object `f` of type `MyElem` has two components
`f.a` and `f.b` of integer type, which should be printed as `a^b`, this can be
implemented as

```julia
expressify(f::MyElem; context = nothing) = Expr(:call, :^, f.a, f.b)
```

If `f.a` and `f.b` themselves are objects that can be expressified, this can
be implemented as

```julia
function expressify(f::MyElem; context = nothing)
  return Expr(:call, :^, expressify(f.a, context = context),
                         expressify(f.b, context = context))
end
```

### Unary operations

```julia
-(f::MyElem)
```

Return $-f$.

### Binary operations

```julia
+(f::MyElem, g::MyElem)
-(f::MyElem, g::MyElem)
*(f::MyElem, g::MyElem)
```

Return $f + g$, $f - g$ or $fg$, respectively.

### Comparison

```
==(f::MyElem, g::MyElem)
```

Return `true` if $f$ and $g$ are arithmetically equal. In the case where the two
elements are inexact, the function returns `true` if they agree to the minimum precision
of the two.

```
isequal(f::MyElem, g::MyElem)
```

For exact rings, this should return the same thing as `==` above. For inexact rings,
this returns `true` only if the two elements are arithmetically equal and have the same
precision.

### Powering

```julia
^(f::MyElem, e::Int)
```

Return $f^e$. The function should throw a `DomainError()` if negative exponents don't
make sense but are passed to the function.

### Exact division

```julia
divexact(f::MyElem, g::MyElem)
```

Return $f/g$, though note that Julia uses `/` for floating point division. Here we
mean exact division in the ring, i.e. return $q$ such that $f = gq$. A `DivideError()`
should be thrown if $g$ is zero. If no exact quotient exists or an impossible inverse
is unavoidably encountered, an error should be thrown.

### Inverse

```julia
inv(f::MyElem)
```

Return the inverse of $f$, i.e. $1/f$, though note that Julia uses `/` for floating
point division. Here we mean exact division in the ring.

A fallback for this function is provided in terms of `divexact` so an implementation
can be omitted if preferred.

### Unsafe operators

To speed up polynomial and matrix arithmetic, it sometimes makes sense to mutate values
in place rather than replace them with a newly created object every time they are
modified.

For this purpose, certain mutating operators are required. In order to support immutable
types (struct in Julia) and systems that don't have in-place operators, all unsafe
operators must return the (ostensibly) mutated value. Only the returned value is used
in computations, so this lifts the requirement that the unsafe operators actually
mutate the value.

Note the exclamation point is a convention, which indicates that the object may be
mutated in-place.

To make use of these functions, one must be certain that no other references are held
to the object being mutated, otherwise those values will also be changed!

The results of `deepcopy` and all arithmetic operations, including powering and division
can be assumed to be new objects without other references being held, as can objects
returned from constructors.

Note that `R(a)` where `R` is the ring `a` belongs to, does not create a new value. For
this case, use `deepcopy(a)`.

```julia
zero!(f::MyElem)
```

Set the value $f$ to zero in place. Return the mutated value.

```julia
mul!(c::MyElem, a::MyElem, b::MyElem)
```

Set $c$ to the value $ab$ in place. Return the mutated value. Aliasing is permitted.

```julia
add!(c::MyElem, a::MyElem, b::MyElem)
```

Set $c$ to the value $a + b$ in place. Return the mutated value. Aliasing is permitted.

```julia
addeq!(a::MyElem, b::MyElem)
```

Set $a$ to $a + b$ in place. Return the mutated value. Aliasing is permitted.

### Random generation

The random functions are only used for test code to generate test data. They therefore
don't need to provide any guarantees on uniformity, and in fact, test values that are
known to be a good source of corner cases can be supplied.

```julia
rand(R::MyParent, v...)
```

Return a random element in the given ring of the specified size.

There can be as many arguments as is necessary to specify the size of the test example
which is being produced.

### Promotion rules

In order for AbstractAlgebra to be able to automatically coerce up towers of rings,
certain promotion rules must be defined. For every ring, one wants to be able to coerce
integers into the ring. And for any ring constructed over a base ring, one would like to
be able to coerce from the base ring into the ring.

The promotion rules look a bit different depending on whether the element type is
parameterised or not and whether it is built on a base ring.

For ring element types `MyElem` that are neither parameterised nor built over a base
ring, the promotion rules can be defined as follows:

```julia
promote_rule(::Type{MyElem}, ::Type{T}) where {T <: Integer} = MyElem
```

For ring element types `MyElem` that aren't parameterised, but which have a base ring
with concrete element type `T` the promotion rules can be defined as follows:

```julia
promote_rule(::Type{MyElem}, ::Type{U}) where U <: Integer = MyElem
```

```julia
promote_rule(::Type{MyElem}, ::Type{T}) = MyElem
```

For ring element types `MyElem{T}` that are parameterised by the type of elements of
the base ring, the promotion rules can be defined as follows:

```julia
promote_rule(::Type{MyElem{T}}, ::Type{MyElem{T}}) where T <: RingElement = MyElem{T}
```

```julia
function promote_rule(::Type{MyElem{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? MyElem{T} : Union{}
end
```

## Required functionality for inexact rings

### Approximation (floating point and ball arithmetic only)

```julia
isapprox(f::MyElem, g::MyElem; atol::Real=sqrt(eps()))
```

This is used by test code that uses rings involving floating point or ball arithmetic.
The function should return `true` if all components of $f$ and $g$ are equal to
within the square root of the Julia epsilon, since numerical noise may make an exact
comparison impossible.

For parameterised rings over an inexact ring, we also require the following ad hoc
approximation functionality.

```julia
isapprox(f::MyElem{T}, g::T; atol::Real=sqrt(eps())) where T <: AbstractAlgebra.RingElem
```

```julia
isapprox(f::T, g::MyElem{T}; atol::Real=sqrt(eps())) where T <: AbstractAlgebra.RingElem
```

These notionally coerce the element of the base ring into the parameterised ring and do
a full comparison.

## Optional functionality

Some functionality is difficult or impossible to implement for all rings in the system.
If it is provided, additional functionality or performance may become available. Here
is a list of all functions that are considered optional and can't be relied on by
generic functions in the AbstractAlgebra Ring interface.

It may be that no algorithm, or no efficient algorithm is known to implement these
functions. As these functions are optional, they do not need to exist. Julia will
already inform the user that the function has not been implemented if it is called but
doesn't exist.

### Optional basic manipulation functionality

```julia
isunit(f::MyElem)
```

Return `true` if the given element is a unit in the ring it belongs to.

```julia
characteristic(R::MyParent)
```

Return the characteristic of the ring. The function should not be defined if
it is not possible to unconditionally give the characteristic as the function
is used in some generic code for correctness, but will always take the safe
path if the function is not defined.

### Optional binary ad hoc operators

By default, ad hoc operations are handled by AbstractAlgebra.jl if they are not defined
explicitly, by coercing both operands into the same ring and then performing the
required operation.

In some cases, e.g. for matrices, this leads to very inefficient behaviour. In such
cases, it is advised to implement some of these operators explicitly.

It can occasionally be worth adding a separate set of ad hoc binary operators for the
type `Int`, if this can be done more efficiently than for arbitrary Julia Integer types.

```julia
+(f::MyElem, c::Integer)
-(f::MyElem, c::Integer)
*(f::MyElem, c::Integer)
```

```julia
+(c::Integer, f::MyElem)
-(c::Integer, f::MyElem)
*(c::Integer, f::MyElem)
```

For parameterised types, it is also sometimes more performant to provide explicit ad
hoc operators with elements of the base ring.

```julia
+(f::MyElem{T}, c::T) where T <: AbstractAlgebra.RingElem
-(f::MyElem{T}, c::T) where T <: AbstractAlgebra.RingElem
*(f::MyElem{T}, c::T) where T <: AbstractAlgebra.RingElem
```

```julia
+(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem
-(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem
*(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem
```

### Optional ad hoc comparisons

```julia
==(f::MyElem, c::Integer)
```

```julia
==(c::Integer, f::MyElem)
```

```julia
==(f::MyElem{T}, c:T) where T <: AbstractAlgebra.RingElem
```

```julia
==(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem
```

### Optional ad hoc exact division functions

```julia
divexact(a::MyElem{T}, b::T) where T <: AbstractAlgebra.RingElem
```

```julia
divexact(a::MyElem, b::Integer)
```

### Optional powering functions

```julia
^(f::MyElem, e::BigInt)
```

In case $f$ cannot explode in size when powered by a very large integer, and it is
practical to do so, one may provide this function to support powering with `BigInt`
exponents (or for external modules, any other big integer type).

### Optional unsafe operators

```julia
addmul!(c::MyElem, a::MyElem, b::MyElem, t::MyElem)
```

Set $c = c + ab$ in-place. Return the mutated value. The value $t$ should be a temporary
of the same type as $a$, $b$ and $c$, which can be used arbitrarily by the
implementation to speed up the computation. Aliasing between $a$, $b$ and $c$ is
permitted.
