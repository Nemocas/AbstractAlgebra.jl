```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```
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

The parent type must belong to `Ring` and the element type must belong
to `RingElem`. Of course, the types may belong to these abstract types
transitively, e.g. `Poly{T}` actually belongs to `PolyRingElem{T}` which in
turn belongs to `RingElem`.

For parameterised rings, we advise that the types of both the parent objects and
element objects to be parameterised by the types of the elements of the base ring
(see [`base_ring`](@ref)).

There can be variations on this theme: e.g. in some areas of mathematics there is a
notion of a coefficient domain, in which case it may make sense to parameterise all
types by the type of elements of this coefficient domain. But note that this may have
implications for the ad hoc operators one might like to explicitly implement.

## RingElement type union

Julia's own number types cannot be subtypes of `RingElem`, so AbstractAlgebra
defines the union type `RingElement`; see [Abstract types for rings](@ref ring-abstract-types).
Generic code should accept `RingElement` rather than `RingElem`, so that it
also works for Julia's number types.

!!! note

    One must be careful when defining ad hoc binary operations for ring element
    types. It is often necessary to define separate versions of the functions for
    `RingElem` then for each of the Julia types separately in
    order to avoid ambiguity warnings.

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
constructor of the parent object. See `src/generic/GenericTypes.jl` for examples of how to
construct and handle such caches.

## Required functions for all rings

In the following, we list all the functions that are required to be provided for rings
in AbstractAlgebra.jl or by external libraries wanting to use AbstractAlgebra.jl.

We give this interface for fictitious types `MyParent` for the type of the ring parent
object `R` and `MyElem` for the type of the elements of the ring.

!!! note

    Generic functions in AbstractAlgebra.jl may not rely on the existence of
    functions that are not documented here. If they do, those functions will only be
    available for rings that implement that additional functionality, and should be
    documented as such.

### Data type and parent object methods

- [`parent_type(::Type{MyElem})`](@ref parent_type)
- [`elem_type(::Type{MyParent})`](@ref elem_type)
- [`base_ring_type(::Type{MyParent})`](@ref base_ring_type): return `Union{}`
  if the ring is not parameterised by another ring.
- [`base_ring(R::MyParent)`](@ref base_ring): only if the ring is
  parameterised by another ring.
- [`coefficient_ring(R::MyParent)`](@ref coefficient_ring) and
  [`coefficient_ring_type(::Type{MyParent})`](@ref coefficient_ring_type): only
  if there is a well-defined notion of a coefficient ring, e.g. for polynomial
  rings.
- [`parent(f::MyElem)`](@ref parent): usually returns a field `parent` stored
  in each element. Element types whose parent can be reconstructed from the
  type alone, e.g. types that are not parameterised, need not store it.
- [`is_domain_type(::Type{MyElem})`](@ref is_domain_type) and
  [`is_exact_type(::Type{MyElem})`](@ref is_exact_type): these see only the
  type, never an element or its parent.
- `Base.hash(f::MyElem, h::UInt)`: hash the components of `f`, e.g. the
  coefficients of a polynomial, recursively, passing the same `h` to all
  levels, and xor each component hash with `h` before combining them. The hash
  functions in AbstractAlgebra start from a 64 bit constant picked at random
  for each type and truncated to `UInt`, so that arithmetically equal objects
  of different types are unlikely to hash to the same value.
- `deepcopy_internal(f::MyElem, dict::IdDict)`: copy all components of `f`
  recursively, but not its parent: the copy must have the identical parent.
  For types that cannot reference themselves internally, `dict` may be
  ignored; otherwise see the Julia documentation of `deepcopy_internal`.

### Constructors

Elements are constructed by calling the parent object; see
[Constructors](@ref ring-constructors) for what each call returns.

- `(R::MyParent)()`
- `(R::MyParent)(a::Integer)`
- `(R::MyParent)(a::MyElem)`: return `a` itself, without copying it, if its
  parent is `R`; otherwise throw an exception.
- `(R::MyParent{T})(a::T) where T <: RingElem`: for parameterised rings,
  coerce `a` from the base ring.

### Basic manipulation of rings and elements

- [`zero(R::MyParent)`](@ref zero)
- [`one(R::MyParent)`](@ref one)
- [`iszero(f::MyElem)`](@ref iszero)
- [`isone(f::MyElem)`](@ref isone)

### Canonicalisation

- [`canonical_unit(f::MyElem)`](@ref canonical_unit). This must exist and
  return an element of the ring. Where no useful normalisation is practical,
  return `one(parent(f))`; otherwise the result must have the properties
  listed in the docstring.

### String I/O

- `show(io::IO, R::MyParent)`: print an English description of the ring. For a
  parameterised ring, it can call `show` on the rings it depends on.
- `show(io::IO, f::MyElem)`: print a human readable, textual representation of
  `f`. It can call `show` on the components of `f`.

### Expressions

To obtain best results when printing composed types derived from other types, e.g., polynomials,
the following method should be implemented.

```@docs
expressify
@enable_all_show_via_expressify
```

`@enable_all_show_via_expressify MyElem` defines the following show methods for
the specified type `MyElem`:

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

As noted above, expressify should return an `Expr`, `Symbol`, `Integer` or
`String`. The rendering of such expressions with a particular MIME type to an
output context is controlled by the following rules which are subject to change
slightly in future versions of AbstracAlgebra.

`Integer`: The printing of integers is straightforward and automatically
includes transformations such as `1 + (-2)*x => 1 - 2*x` as this is cumbersome
to implement per-type.

`Symbol`: Since variable names are stored as mere symbols in AbstractAlgebra,
some transformations related to subscripts are applied to symbols automatically
in latex output. The `\operatorname{` in the following table is actually
replaced with the more portable `\mathop{\mathrm{`.

expressify             | latex output
:----------------------|:-----------------------------
`Symbol("a")`          | `a`
`Symbol("α")`          | `{\alpha}`
`Symbol("x1")`         | `\operatorname{x1}`
`Symbol("xy_1")`       | `\operatorname{xy}_{1}`
`Symbol("sin")`        | `\operatorname{sin}`
`Symbol("sin_cos")`    | `\operatorname{sin\_cos}`
`Symbol("sin_1")`      | `\operatorname{sin}_{1}`
`Symbol("sin_cos_1")`  | `\operatorname{sin\_cos}_{1}`
`Symbol("αaβb_1_2")`   | `\operatorname{{\alpha}a{\beta}b}_{1,2}`

`Expr`: These are the most versatile as the `Expr` objects themselves contain
a symbolic head and any number of arguments. What looks like `f(a,b)` in textual
output is `Expr(:call, :f, :a, :b)` under the hood. AbstractAlgebra currently
contains the following printing rules for such expressions.

expressify                 | output    | latex notes
:--------------------------|:----------|:------------
`Expr(:call, :+, a, b)`    | `a + b`   |
`Expr(:call, :*, a, b)`    | `a*b`     | one space for implied multiplication
`Expr(:call, :cdot, a, b)` | `a * b`   | a real `\cdot` is used
`Expr(:call, :^, a, b)`    | `a^b`     | may include some courtesy parentheses
`Expr(:call, ://, a, b)`   | `a//b`    | will create a fraction box
`Expr(:call, :/, a, b)`    | `a/b`     | will not create a fraction box
`Expr(:call, a, b, c)`     | `a(b, c)` |
`Expr(:ref, a, b, c)`      | `a[b, c]` |
`Expr(:vcat, a, b)`        | `[a; b]`  | actually vertical
`Expr(:vect, a, b)`        | `[a, b]`  |
`Expr(:tuple, a, b)`       | `(a, b)`  |
`Expr(:list, a, b)`        | `{a, b}`  |
`Expr(:series, a, b)`      | `a, b`    |
`Expr(:sequence, a, b)`    | `ab`      |
`Expr(:row, a, b)`         | `a b`     | combine with `:vcat` to make matrices
`Expr(:hcat, a, b)`        | `a b`     |

`String`: Strings are printed verbatim and should only be used as a last resort
as they provide absolutely no precedence information on their contents.


### Unary operations

- `-(f::MyElem)`

### Binary operations

- `+(f::MyElem, g::MyElem)`
- `-(f::MyElem, g::MyElem)`
- `*(f::MyElem, g::MyElem)`

### Comparison

- `==(f::MyElem, g::MyElem)`: for inexact elements, return `true` if `f` and
  `g` agree to the minimum precision of the two.
- `isequal(f::MyElem, g::MyElem)`: the same as `==` for exact rings; for
  inexact rings, additionally require equal precision.

### Powering

- `^(f::MyElem, e::Int)`: throw a `DomainError` if `e` is negative and that
  makes no sense for `f`.

### Exact division

- [`divexact(f::MyElem, g::MyElem; check::Bool=true)`](@ref divexact). Throw a
  `DivideError` if `g` is zero. With `check=false` the exactness check may be
  skipped; what happens on an inexact division is then undefined, including
  meaningless results, hanging or crashing.

### Inverse

- [`inv(f::MyElem)`](@ref Base.inv(::RingElem)). Optional: the default calls
  `divexact(one(parent(f)), f)`.

### Random generation

- `rand(R::MyParent, v...)`: return a random element of `R`, of a size given by
  the arguments `v`. It is only used to generate test data, so it need not be
  uniform, and values known to be good sources of corner cases are welcome.
  See the [Random interface](@ref).

### Promotion rules

AbstractAlgebra currently has a very simple coercion model. With few exceptions
only simple coercions are supported. For example if $x \in \mathbb{Z}$ and
$y \in \mathbb{Z}[x]$ then $x + y$ can be computed by coercing $x$ into
the same ring as $y$ and then adding in that ring.

Complex coercions such as adding elements of $\mathbb{Q}$ and $\mathbb{Z}[x]$
are not supported, as this would require finding and creating a common
overring in which the elements could be added.

AbstractAlgebra supports simple coercions by overloading parent object call
syntax `R(x)` to coerce the object `x` into the ring `R`. However, to coerce
elements up a tower of rings, one needs to also have a promotion system
similar to Julia's type promotion system.

As for Julia, AbstractAlgebra's promotion system only specifies what happens
to types. It is the coercions themselves that must deal with the mathematical
situation at the level of rings, including checking that the object can even
be coerced into the given ring.

Note that `AbstractAlgebra.promote_rule` is a separate function from
`Base.promote_rule` and is the one that must be extended.

We now describe the required AbstractAlgebra type promotion rules.

For every ring, one wants to be able to coerce integers into the ring. And for
any ring constructed over a base ring, one would like to be able to coerce from
the base ring into the ring.

The required promotion rules to support this look a bit different depending on
whether the element type is parameterised or not and whether it is built on a
base ring.

For ring element types `MyElem` that are neither parameterised nor built over a
base ring, the promotion rules can be defined as follows:

```julia
AbstractAlgebra.promote_rule(::Type{MyElem}, ::Type{T}) where {T <: Integer} = MyElem
```

For ring element types `MyElem` that aren't parameterised, but which have a
base ring with concrete element type `T` the promotion rules can be defined as
follows:

```julia
AbstractAlgebra.promote_rule(::Type{MyElem}, ::Type{U}) where U <: Integer = MyElem
```

```julia
AbstractAlgebra.promote_rule(::Type{MyElem}, ::Type{T}) = MyElem
```

For ring element types `MyElem{T}` that are parameterised by the type of
elements of the base ring, the promotion rules can be defined as follows:

```julia
AbstractAlgebra.promote_rule(::Type{MyElem{T}}, ::Type{MyElem{T}}) where T <: RingElement = MyElem{T}
```

```julia
function AbstractAlgebra.promote_rule(::Type{MyElem{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   AbstractAlgebra.promote_rule(T, U) == T ? MyElem{T} : Union{}
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
isapprox(f::MyElem{T}, g::T; atol::Real=sqrt(eps())) where T <: RingElem
```

```julia
isapprox(f::T, g::MyElem{T}; atol::Real=sqrt(eps())) where T <: RingElem
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

### Optional unsafe operators

The various operators described in [Unsafe ring operators](@ref) such as
`add!` and `mul!` have default implementations which are not faster than their
regular safe counterparts. Implementors may wish to implement some or all of
them for their rings. Note that in general only the variants with the most
arguments needs to be implemented. E.g. for `add!` only `add!(z,a,b)` has to be
implemented for any new ring type, as `add!(a,b)` delegates to `add!(a,a,b)`.

### Optional basic manipulation functionality

- [`is_unit(f::MyElem)`](@ref is_unit)
- [`is_zero_divisor(f::MyElem)`](@ref is_zero_divisor): without it, the total
  ring of fractions may not be usable over the ring. All fields have a
  fallback.
- [`characteristic(R::MyParent)`](@ref characteristic): define it only if the
  characteristic can be given unconditionally; otherwise the default raises an
  exception.

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
+(f::MyElem{T}, c::T) where T <: RingElem
-(f::MyElem{T}, c::T) where T <: RingElem
*(f::MyElem{T}, c::T) where T <: RingElem
```

```julia
+(c::T, f::MyElem{T}) where T <: RingElem
-(c::T, f::MyElem{T}) where T <: RingElem
*(c::T, f::MyElem{T}) where T <: RingElem
```

### Optional ad hoc comparisons

```julia
==(f::MyElem, c::Integer)
```

```julia
==(c::Integer, f::MyElem)
```

```julia
==(f::MyElem{T}, c:T) where T <: RingElem
```

```julia
==(c::T, f::MyElem{T}) where T <: RingElem
```

### Optional ad hoc exact division functions

```julia
divexact(a::MyElem{T}, b::T) where T <: RingElem
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

## Minimal example of ring implementation

Here is a minimal example of implementing the Ring Interface for a constant
polynomial type (i.e. polynomials of degree less than one).

```jldoctest ConstPoly
# ConstPoly.jl : Implements constant polynomials

using AbstractAlgebra

using Random: Random, SamplerTrivial
using AbstractAlgebra.RandomExtensions: RandomExtensions, Make2, AbstractRNG

import AbstractAlgebra: parent_type, elem_type, base_ring, base_ring_type, parent, is_domain_type,
       is_exact_type, canonical_unit, isequal, divexact, zero!, mul!, add!,
       get_cached!, is_unit, characteristic, Ring, RingElem, expressify,
       @show_name, @show_special, is_terse, pretty, terse, Lowercase,
       promote_rule

import Base: show, +, -, *, ^, ==, inv, isone, iszero, one, zero, rand,
             deepcopy_internal, hash

@attributes mutable struct ConstPolyRing{T <: RingElement} <: Ring
   base_ring::Ring

   function ConstPolyRing{T}(R::Ring, cached::Bool) where T <: RingElement
      return get_cached!(ConstPolyID, R, cached) do
         new{T}(R)
      end::ConstPolyRing{T}
   end
end

const ConstPolyID = AbstractAlgebra.CacheDictType{Ring, ConstPolyRing}()

mutable struct ConstPoly{T <: RingElement} <: RingElem
   c::T
   parent::ConstPolyRing{T}

   function ConstPoly{T}(c::T) where T <: RingElement
      return new(c)
   end
end

# Data type and parent object methods

parent_type(::Type{ConstPoly{T}}) where T <: RingElement = ConstPolyRing{T}

elem_type(::Type{ConstPolyRing{T}}) where T <: RingElement = ConstPoly{T}

base_ring_type(::Type{ConstPolyRing{T}}) where T <: RingElement = parent_type(T)

base_ring(R::ConstPolyRing) = R.base_ring::base_ring_type(R)

parent(f::ConstPoly) = f.parent

is_domain_type(::Type{ConstPoly{T}}) where T <: RingElement = is_domain_type(T)

is_exact_type(::Type{ConstPoly{T}}) where T <: RingElement = is_exact_type(T)

function hash(f::ConstPoly, h::UInt)
   r = 0x65125ab8e0cd44ca
   return xor(r, hash(f.c, h))
end

function deepcopy_internal(f::ConstPoly{T}, dict::IdDict) where T <: RingElement
   r = ConstPoly{T}(deepcopy_internal(f.c, dict))
   r.parent = f.parent # parent should not be deepcopied
   return r
end

# Basic manipulation

zero(R::ConstPolyRing) = R()

one(R::ConstPolyRing) = R(1)

iszero(f::ConstPoly) = iszero(f.c)

isone(f::ConstPoly) = isone(f.c)

is_unit(f::ConstPoly) = is_unit(f.c)

characteristic(R::ConstPolyRing) = characteristic(base_ring(R))

# Canonical unit

canonical_unit(f::ConstPoly) = canonical_unit(f.c)

# String I/O

function show(io::IO, R::ConstPolyRing)
   @show_name(io, R)
   @show_special(io, R)
   print(io, "Constant polynomials")
   if !is_terse(io)
     io = pretty(io)
     print(terse(io), " over ", Lowercase(), base_ring(R))
   end
end

function show(io::IO, f::ConstPoly)
   print(io, f.c)
end

# Expressification (optional)

function expressify(R::ConstPolyRing; context = nothing)
   return Expr(:sequence, Expr(:text, "Constant polynomials over "),
                          expressify(base_ring(R), context = context))
end

function expressify(f::ConstPoly; context = nothing)
   return expressify(f.c, context = context)
end

# Unary operations

function -(f::ConstPoly)
   R = parent(f)
   return R(-f.c)
end

# Binary operations

function +(f::ConstPoly{T}, g::ConstPoly{T}) where T <: RingElement
   check_parent(f, g)
   R = parent(f)
   return R(f.c + g.c)
end

function -(f::ConstPoly{T}, g::ConstPoly{T}) where T <: RingElement
   check_parent(f, g)
   R = parent(f)
   return R(f.c - g.c)
end

function *(f::ConstPoly{T}, g::ConstPoly{T}) where T <: RingElement
   check_parent(f, g)
   R = parent(f)
   return R(f.c*g.c)
end

# Comparison

function ==(f::ConstPoly{T}, g::ConstPoly{T}) where T <: RingElement
   check_parent(f, g)
   return f.c == g.c
end

function isequal(f::ConstPoly{T}, g::ConstPoly{T}) where T <: RingElement
   check_parent(f, g)
   return isequal(f.c, g.c)
end

# Powering need not be implemented if * is

# Exact division

function divexact(f::ConstPoly{T}, g::ConstPoly{T}; check::Bool = true) where T <: RingElement
   check_parent(f, g)
   R = parent(f)
   return R(divexact(f.c, g.c, check = check))
end

# Inverse

function inv(f::ConstPoly)
   R = parent(f)
   return R(AbstractAlgebra.inv(f.c))
end

# Unsafe operators

function zero!(f::ConstPoly)
   f.c = zero(base_ring(parent(f)))
   return f
end

function one!(f::ConstPoly)
   f.c = one(base_ring(parent(f)))
   return f
end

function mul!(f::ConstPoly{T}, g::ConstPoly{T}, h::ConstPoly{T}) where T <: RingElement
   f.c = g.c*h.c
   return f
end

function add!(f::ConstPoly{T}, g::ConstPoly{T}, h::ConstPoly{T}) where T <: RingElement
   f.c = g.c + h.c
   return f
end

# Random generation

RandomExtensions.maketype(R::ConstPolyRing, _) = elem_type(R)

rand(rng::AbstractRNG, sp::SamplerTrivial{<:Make2{ConstPoly,ConstPolyRing}}) =
        sp[][1](rand(rng, sp[][2]))

rand(rng::AbstractRNG, R::ConstPolyRing, n::AbstractUnitRange{Int}) = R(rand(rng, n))

rand(R::ConstPolyRing, n::AbstractUnitRange{Int}) = rand(Random.default_rng(), R, n)

# Promotion rules

AbstractAlgebra.promote_rule(::Type{ConstPoly{T}}, ::Type{ConstPoly{T}}) where T <: RingElement = ConstPoly{T}

function AbstractAlgebra.promote_rule(::Type{ConstPoly{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   AbstractAlgebra.promote_rule(T, U) == T ? ConstPoly{T} : Union{}
end

# Constructors

function (R::ConstPolyRing{T})() where T <: RingElement
   r = ConstPoly{T}(base_ring(R)(0))
   r.parent = R
   return r
end

function (R::ConstPolyRing{T})(c::Integer) where T <: RingElement
   r = ConstPoly{T}(base_ring(R)(c))
   r.parent = R
   return r
end

# Needed to prevent ambiguity
function (R::ConstPolyRing{T})(c::T) where T <: Integer
   r = ConstPoly{T}(base_ring(R)(c))
   r.parent = R
   return r
end

function (R::ConstPolyRing{T})(c::T) where T <: RingElement
   base_ring(R) != parent(c) && error("Unable to coerce element")
   r = ConstPoly{T}(c)
   r.parent = R
   return r
end

function (R::ConstPolyRing{T})(f::ConstPoly{T}) where T <: RingElement
   R != parent(f) && error("Unable to coerce element")
   return f
end

# Parent constructor

function constant_polynomial_ring(R::Ring, cached::Bool=true)
   T = elem_type(R)
   return ConstPolyRing{T}(R, cached)
end

# output

constant_polynomial_ring (generic function with 2 methods)
```

The above implementation of `constant_polynomial_ring` may be tested as follows.

```jldoctest ConstPoly; filter = r".*"s
using Test

function ConformanceTests.generate_element(R::ConstPolyRing{elem_type(ZZ)})
   n = rand(1:999)
   return R(rand(-n:n))
end

test_Ring_interface(constant_polynomial_ring(ZZ))

# output
Test Summary:                                                                       |  Pass  Total  Time
Ring interface for Constant polynomials over integers of type ConstPolyRing{BigInt} | 13844  13844  0.9s
```

Note that we only showed a minimal implementation of the ring interface.
Additional interfaces exists, e.g. for Euclidean rings. Additional interface
usually require implementing additional methods, and in some cases we also
provide additional conformance tests. In this case, just one necessary
method is missing.

```jldoctest ConstPoly
function Base.divrem(a::ConstPoly{elem_type(ZZ)}, b::ConstPoly{elem_type(ZZ)})
   check_parent(a, b)
   q, r = AbstractAlgebra.divrem(a.c, b.c)
   return parent(a)(q), parent(a)(r)
end

# output

```

We can test it like this.

```jldoctest ConstPoly; filter = r".*"s
test_EuclideanRing_interface(constant_polynomial_ring(ZZ))

# output
Test Summary:                                                                                 | Pass  Total  Time
Euclidean Ring interface for Constant polynomials over integers of type ConstPolyRing{BigInt} | 2212   2212  0.1s
```
