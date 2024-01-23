```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Total ring of fractions

AbstractAlgebra.jl provides a module, implemented in
`src/generic/TotalFraction.jl`, for the total ring of fractions of a ring.

The total ring of fractions of a ring `R` is the localisation of `R` at the
non-zero divisors of `R`, the latter being a multiplicative subset of `R`.

There are no restrictions on the ring except the function `is_zero_divisor`
must be defined and effective for `R`.

In particular, we do not assume that all elements of `R` which are not zero
divisors are units in `R`. This has the effect of making exact division
impossible generically in the total ring of fractions of `R`.

This in turn limits the usefulness of the total ring of fractions as a ring
in AbstractAlgebra as a great deal of generic code relies on `divexact`.
Should this be a limitation, the user can define their own `divexact`
function for the total ring of fractions in question.

Note that in most cases `a*inv(b)` is not a sufficient definition of
`divexact(a, b)` due to the possibility that `b` is not a unit in the total
ring of fractions.

It is also possible to construct a total ring of fractions of `R` without the
`is_zero_divisor` function existing for `R`, but some functions such as
`is_unit`, `inv`, `rand` and ad hoc arithmetic operations involving rational
numbers are not available for the total ring of fractions. One must also
construct fractions using the option `check=false` and it is one's own
responsibility to check that the denominator is not a zero divisor.

Note that although the total ring of fractions of an integral domain `R` is
mathematically the same thing as the fraction field of `R`, these will be
different objects in AbstractAlgebra and have different types.

## Generic total ring of fraction types

AbstractAlgebra.jl implements a generic type for elements of a total ring of
fractions, namely`Generic.TotFrac{T}` where `T` is the type of elements of the
base ring. See the file `src/generic/GenericTypes.jl` for details.

Parent objects of such elements have type `Generic.TotFracRing{T}`.

## Abstract types

The types for elements of a total ring of fractions belong directly to the
abstract type `RingElem` and the type for the total ring of fractions parent
object belongs directly to the abstract type `Ring`.

## Total ring of fractions constructors

In order to construct fractions in a total ring of fractions in
AbstractAlgebra.jl, one must first construct the parent object for the total
ring of fractions itself. This is accomplished with the following constructor.

```julia
total_ring_of_fractions(R::Ring; cached::Bool = true)
```

Given a base ring `R` return the parent object of the total ring of fractions
of $R$. By default the parent object `S` will depend only on `R` and will be
cached. Setting the optional argument `cached` to `false` will prevent the
parent object `S` from being cached.

Here are some examples of creating a total ring of fractions and making use of
the resulting parent objects to coerce various elements into the ring.

**Examples**

```jldoctest
julia> R, x = polynomial_ring(ZZ, "x")
(Univariate polynomial ring in x over integers, x)

julia> S = total_ring_of_fractions(R)
Total ring of fractions of Univariate polynomial ring in x over integers

julia> f = S()
0

julia> g = S(123)
123

julia> h = S(BigInt(1234))
1234

julia> k = S(x + 1)
x + 1
```

## Fraction constructors

One can construct fractions using the total ring of fractions parent object,
as for any ring or field.

```julia
(R::TotFracRing)() # constructs zero
(R::TotFracRing)(c::Integer)
(R::TotFracRing)(c::elem_type(R))
(R::TotFracRing{T})(a::T) where T <: RingElement
```

Although one cannot use the double slash operator `//` to construct elements
of a total ring of fractions, as no parent has been specified, one can use the
double slash operator to construct elements of a total ring of fractions so
long as one of the arguments to the double slash operator is already in the
total ring of fractions in question.

**Examples**

```jldoctest
julia> R, x = polynomial_ring(QQ, "x")
(Univariate polynomial ring in x over rationals, x)

julia> S = total_ring_of_fractions(R)
Total ring of fractions of Univariate polynomial ring in x over rationals

julia> f = S(x + 1)
x + 1

julia> f//3
(x + 1)//3

julia> 3//f
3//(x + 1)

julia> f//x
(x + 1)//x
```

## Functions for types and parents of total rings of fractions

Total rings of fractions in AbstractAlgebra.jl implement the Ring interface
except for the `divexact` function which is not generically possible to
implement.

```julia
base_ring(R::TotFracRing)
base_ring(a::TotFrac)
```

Return the base ring of which the total ring of fractions was constructed.

```julia
parent(a::TotFrac)
```

Return the total ring of fractions that the given fraction belongs to.

```julia
characteristic(R::TotFracRing)
```

Return the characteristic of the base ring of the total ring of fractions. If
the characteristic is not known an exception is raised.


**Examples**

```jldoctest
julia> R, x = polynomial_ring(QQ, "x")
(Univariate polynomial ring in x over rationals, x)

julia> S = total_ring_of_fractions(R)
Total ring of fractions of Univariate polynomial ring in x over rationals

julia> f = S(x + 1)
x + 1

julia> U = base_ring(S)
Univariate polynomial ring in x over rationals

julia> V = base_ring(f)
Univariate polynomial ring in x over rationals

julia> T = parent(f)
Total ring of fractions of Univariate polynomial ring in x over rationals

julia> m = characteristic(S)
0
```

## Total ring of fractions functions

### Basic functions

Total rings of fractions implement the Ring interface.

```julia
zero(R::TotFracRing)
one(R::TotFracRing)
iszero(a::TotFrac)
isone(a::TotFrac)
```

```julia
inv(a::T) where T <: TotFrac
```

They also implement some of the following functions which would usually be
associated with the field and fraction field interfaces.

```julia
is_unit(f::TotFrac)
```

```julia
numerator(a::TotFrac)
denominator(a::TotFrac)
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(QQ, "x")
(Univariate polynomial ring in x over rationals, x)

julia> S = total_ring_of_fractions(R)
Total ring of fractions of Univariate polynomial ring in x over rationals

julia> f = S(x + 1)
x + 1

julia> g = f//(x^3 + 3x + 1)
(x + 1)//(x^3 + 3*x + 1)

julia> h = zero(S)
0

julia> k = one(S)
1

julia> isone(k)
true

julia> iszero(f)
false

julia> r = deepcopy(f)
x + 1

julia> n = numerator(g)
x + 1

julia> d = denominator(g)
x^3 + 3*x + 1
```

### Random generation

Random fractions can be generated using `rand`. The parameters passed after the
total ring of fractions tell `rand` how to generate random elements of the base
ring.

```julia
rand(R::TotFracRing, v...)
```

**Examples**

```jldoctest; setup = :(import Random; Random.seed!(42))
julia> R, = residue_ring(ZZ, 12);

julia> K = total_ring_of_fractions(R)
Total ring of fractions of Residue ring of integers modulo 12

julia> f = rand(K, 0:11)
7//5

julia> R, x = polynomial_ring(ZZ, "x")
(Univariate polynomial ring in x over integers, x)

julia> S = total_ring_of_fractions(R)
Total ring of fractions of Univariate polynomial ring in x over integers

julia> g = rand(S, -1:3, -10:10)
(4*x + 4)//(-4*x^2 - x + 4)
```
