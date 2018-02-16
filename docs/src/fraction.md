```@meta
CurrentModule = AbstractAlgebra
```

# Generic fraction fields

AbstractAlgebra.jl provides a module, implemented in `src/generic/Fraction.jl` for
generic fraction fields over any gcd domain belonging to the AbstractAlgebra.jl
abstract type hierarchy.

As well as implementing the Fraction Field interface a number of generic algorithms are
implemented for fraction fields. We describe this generic functionality below.

All of the generic functionality is part of a submodule of AbstractAlgebra called
`Generic`. This is exported by default so that it is not necessary to qualify the
function names with the submodule name.

## Types and parent objects

Fractions implemented using the AbstractAlgebra generics have type `Generic.Frac{T}`
where `T` is the type of elements of the base ring. See the file
`src/generic/GenericTypes.jl` for details.

Parent objects of such fraction elements have type `Generic.FracField{T}`.

The fraction element types belong to the abstract type `AbstractAlgebra.FracElem{T}`
and the fraction field types belong to the abstract type `AbstractAlgebra.FracRing{T}`.
This enables one to write generic functions that can accept any AbstractAlgebra
fraction type.

Note that both the generic fraction field type `Generic.FracField{T}` and the abstract
type it belongs to, `AbstractAlgebra.FracField{T}` are both called `FracField`. The 
former is a (parameterised) concrete type for a fraction field over a given base ring
whose elements have type `T`. The latter is an abstract type representing all
fraction field types in AbstractAlgebra.jl, whether generic or very specialised (e.g.
supplied by a C library).

## Fraction field constructors

In order to construct fractions in AbstractAlgebra.jl, one can first construct the
fraction field itself. This is accomplished with the following constructor.

```julia
FractionField(R::AbstractAlgebra.Ring; cached::Bool = true)
```

Given a base ring `R` return the parent object of the fraction field of $R$. By default
the parent object `S` will depend only on `R` and will be cached. Setting the optional
argument `cached` to `false` will prevent the parent object `S` from being cached.

Here are some examples of creating fraction fields and making use of the
resulting parent objects to coerce various elements into the fraction field.

**Examples**

```julia
R, x = PolynomialRing(ZZ, "x")
S = FractionField(R)

f = S()
g = S(123)
h = S(BigInt(1234))
k = S(x + 1)
```

All of the examples here are generic fraction fields, but specialised implementations
of fraction fields provided by external modules will also usually provide a
`FractionField` constructor to allow creation of the fraction fields they provide.

## Basic field functionality

Fraction fields in AbstractAlgebra.jl implement the full Field interface. Of course
the entire Fraction Field interface is also implemented.

We give some examples of such functionality.

**Examples**

```julia
R, x = PolynomialRing(QQ, "x")
S = FractionField(R)

f = S(x + 1)
g = (x^2 + x + 1)//(x^3 + 3x + 1)

h = zero(S)
k = one(S)
isone(k) == true
iszero(f) == false
m = characteristic(S)
U = base_ring(S)
V = base_ring(f)
T = parent(f)
r = deepcopy(f)
n = numerator(g)
d = denominator(g)
```

## Fraction field functionality provided by AbstractAlgebra.jl

The functionality listed below is automatically provided by AbstractAlgebra.jl for
any fraction field module that implements the full Fraction Field interface.
This includes AbstractAlgebra.jl's own generic fraction fields.

But if a C library provides all the functionality documented in the Fraction Field
interface, then all the functions described here will also be automatically supplied by
AbstractAlgebra.jl for that fraction field type.

Of course, modules are free to provide specific implementations of the functions
described here, that override the generic implementation.

### Greatest common divisor

```@docs
gcd{T <: RingElem}(::FracElem{T}, ::FracElem{T})
```

**Examples**

```julia
R, x = PolynomialRing(QQ, "x")

f = (x + 1)//(x^3 + 3x + 1)
g = (x^2 + 2x + 1)//(x^2 + x + 1)

h = gcd(f, g)
```

### Remove and valuation

When working over a Euclidean domain, it is convenient to extend valuations to the
fraction field. To facilitate this, we define the following functions.

```@docs
remove{T <: RingElem}(::FracElem{T}, ::T)
```

```@docs
valuation{T <: RingElem}(::FracElem{T}, ::T)
```

**Examples**

```julia
R, x = PolynomialRing(ZZ, "x")

f = (x + 1)//(x^3 + 3x + 1)
g = (x^2 + 1)//(x^2 + x + 1)

v, q = remove(f^3*g, x + 1)
v = valuation(f^3*g, x + 1)
```

