```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Generic univariate polynomials over a noncommutative ring

AbstractAlgebra.jl provides a module, implemented in `src/generic/NCPoly.jl` for generic
polynomials over any noncommutative ring belonging to the AbstractAlgebra abstract type
hierarchy.

As well as implementing the Univariate Polynomial interface, there are many additional
generic algorithms implemented for such polynomial rings. We describe this generic
functionality below.

All of the generic functionality is part of a submodule of AbstractAlgebra called
`Generic`. This is exported by default so that it is not necessary to qualify the
function names with the submodule name.

## Types and parent objects

Polynomials implemented using the AbstractAlgebra generics have type `Generic.NCPoly{T}`
where `T` is the type of elements of the coefficient ring. Internally they consist of
a Julia array of coefficients and some additional fields for length and a parent object,
etc. See the file `src/generic/GenericTypes.jl` for details.

Parent objects of such polynomials have type `Generic.NCPolyRing{T}`.

The string representation of the variable of the polynomial ring and the
base/coefficient ring $R$ is stored in the parent object. 

The polynomial element types belong to the abstract type `AbstractAlgebra.NCPolyElem{T}`
and the polynomial ring types belong to the abstract type
`AbstractAlgebra.NCPolyRing{T}`. This enables one to write generic functions that can
accept any AbstractAlgebra polynomial type.

Note that both the generic polynomial ring type `Generic.NCPolyRing{T}` and the abstract
type it belongs to, `AbstractAlgebra.NCPolyRing{T}` are both called `NCPolyRing`. The 
former is a (parameterised) concrete type for a polynomial ring over a given base ring
whose elements have type `T`. The latter is an abstract type representing all
polynomial ring types in AbstractAlgebra.jl, whether generic or very specialised (e.g.
supplied by a C library).

## Polynomial ring constructors

In order to construct polynomials in AbstractAlgebra.jl, one must first construct the
polynomial ring itself. This is accomplished with the following constructor.

```julia
PolynomialRing(R::AbstractAlgebra.NCRing, s::AbstractString; cached::Bool = true)
```

Given a base ring `R` and string `s` specifying how the generator (variable) should be
printed, return a tuple `S, x` representing the new polynomial ring $S = R[x]$ and the
generator $x$ of the ring. By default the parent object `S` will depend only on `R` and 
`x` and will be cached. Setting the optional argument `cached` to `false` will prevent
the parent object `S` from being cached.

A shorthand version of this function is provided: given a base ring `R`, we abbreviate
the constructor as follows.

```julia
R["x"]
```

Here are some examples of creating polynomial rings and making use of the
resulting parent objects to coerce various elements into the polynomial ring.

**Examples**

```jldoctest
julia> R = MatrixAlgebra(ZZ, 2)
Matrix Algebra of degree 2 over Integers

julia> S, x = PolynomialRing(R, "x")
(Univariate Polynomial Ring in x over Matrix Algebra of degree 2 over Integers, x)

julia> T, y = PolynomialRing(S, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Matrix Algebra of degree 2 over Integers, y)

julia> U, z = R["z"]
(Univariate Polynomial Ring in z over Matrix Algebra of degree 2 over Integers, z)

julia> f = S()
[0 0]
[0 0]

julia> g = S(123)
([123 0]
[0 123])

julia> h = T(BigInt(1234))
([1234 0]
[0 1234])

julia> k = T(x + 1)
(x+([1 0]
[0 1]))

julia> m = U(z + 1)
z+([1 0]
[0 1])

```

All of the examples here are generic polynomial rings, but specialised implementations
of polynomial rings provided by external modules will also usually provide a
`PolynomialRing` constructor to allow creation of their polynomial rings.

## Basic ring functionality

Once a polynomial ring is constructed, there are various ways to construct
polynomials in that ring.

The easiest way is simply using the generator returned by the `PolynomialRing`
constructor and build up the polynomial using basic arithmetic, as described in
the Ring interface. 

The Julia language also has special syntax for the construction of polynomials in terms
of a generator, e.g. we can write `2x` instead of `2*x`.

The polynomial rings in AbstractAlgebra.jl implement the full Ring interface. Of course
the entire Univariate Polynomial Ring interface is also implemented.

We give some examples of such functionality.

**Examples**

```jldoctest
julia> R = MatrixAlgebra(ZZ, 2)
Matrix Algebra of degree 2 over Integers

julia> S, x = PolynomialRing(R, "x")
(Univariate Polynomial Ring in x over Matrix Algebra of degree 2 over Integers, x)

julia> T, y = PolynomialRing(S, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Matrix Algebra of degree 2 over Integers, y)

julia> f = x^3 + 3x + 21
x^3+([3 0]
[0 3])*x+([21 0]
[0 21])

julia> g = (x + 1)*y^2 + 2x + 1
(x+([1 0]
[0 1]))*y^2+(([2 0]
[0 2])*x+([1 0]
[0 1]))

julia> h = zero(T)
[0 0]
[0 0]

julia> k = one(S)
([1 0]
[0 1])

julia> isone(k)
true

julia> iszero(f)
false

julia> n = length(g)
3

julia> U = base_ring(T)
Univariate Polynomial Ring in x over Matrix Algebra of degree 2 over Integers

julia> V = base_ring(y + 1)
Univariate Polynomial Ring in x over Matrix Algebra of degree 2 over Integers

julia> v = var(T)
:y

julia> U = parent(y + 1)
Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Matrix Algebra of degree 2 over Integers

julia> g == deepcopy(g)
true
```

## Polynomial functionality provided by AbstractAlgebra.jl

The functionality listed below is automatically provided by AbstractAlgebra.jl for
any polynomial module that implements the full Univariate Polynomial Ring interface
over a noncommutative ring. This includes AbstractAlgebra.jl's own generic polynomial
rings.

But if a C library provides all the functionality documented in the Univariate
Polynomial Ring interface over a noncommutative ring, then all the functions described
here will also be automatically supplied by AbstractAlgebra.jl for that polynomial type.

Of course, modules are free to provide specific implementations of the functions
described here, that override the generic implementation.

### Basic functionality

```@docs
lead(::NCPolyElem)
trail(::NCPolyElem)
```

```@docs
gen(::NCPolyRing)
```

```@docs
isgen(::NCPolyElem)
```

```@docs
isunit(::NCPolyElem)
```

```@docs
ismonomial(::NCPolyElem)
```

```@docs
isterm(::NCPolyElem)
```

**Examples**

```jldoctest
julia> R = MatrixAlgebra(ZZ, 2)
Matrix Algebra of degree 2 over Integers

julia> S, x = PolynomialRing(R, "x")
(Univariate Polynomial Ring in x over Matrix Algebra of degree 2 over Integers, x)

julia> T, y = PolynomialRing(S, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Matrix Algebra of degree 2 over Integers, y)

julia> a = zero(T)
[0 0]
[0 0]

julia> b = one(T)
([1 0]
[0 1])

julia> c = BigInt(1)*y^2 + BigInt(1)
y^2+([1 0]
[0 1])

julia> d = x*y^2 + (x + 1)*y + 3
(x)*y^2+(x+([1 0]
[0 1]))*y+([3 0]
[0 3])

julia> f = lead(d)
x

julia> y = gen(T)
y

julia> g = isgen(y)
true

julia> m = isunit(b)
true

julia> n = degree(d)
2

julia> isterm(2y^2)
true

julia> ismonomial(y^2)
true

```

### Truncation

```@docs
truncate(::NCPolyElem, ::Int)
```

```@docs
mullow(::NCPolyElem{T}, ::NCPolyElem{T}, ::Int) where T <: NCRingElem
```

**Examples**

```jldoctest
julia> R = MatrixAlgebra(ZZ, 2)
Matrix Algebra of degree 2 over Integers

julia> S, x = PolynomialRing(R, "x")
(Univariate Polynomial Ring in x over Matrix Algebra of degree 2 over Integers, x)

julia> T, y = PolynomialRing(S, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Matrix Algebra of degree 2 over Integers, y)

julia> f = x*y^2 + (x + 1)*y + 3
(x)*y^2+(x+([1 0]
[0 1]))*y+([3 0]
[0 3])

julia> g = (x + 1)*y + (x^3 + 2x + 2)
(x+([1 0]
[0 1]))*y+(x^3+([2 0]
[0 2])*x+([2 0]
[0 2]))

julia> h = truncate(f, 1)
([3 0]
[0 3])

julia> k = mullow(f, g, 4)
(x^2+x)*y^3+(x^4+([3 0]
[0 3])*x^2+([4 0]
[0 4])*x+([1 0]
[0 1]))*y^2+(x^4+x^3+([2 0]
[0 2])*x^2+([7 0]
[0 7])*x+([5 0]
[0 5]))*y+(([3 0]
[0 3])*x^3+([6 0]
[0 6])*x+([6 0]
[0 6]))

```

### Reversal

```@docs
reverse(::NCPolyElem, ::Int)
reverse(::NCPolyElem)
```

**Examples**

```jldoctest
julia> R = MatrixAlgebra(ZZ, 2)
Matrix Algebra of degree 2 over Integers

julia> S, x = PolynomialRing(R, "x")
(Univariate Polynomial Ring in x over Matrix Algebra of degree 2 over Integers, x)

julia> T, y = PolynomialRing(S, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Matrix Algebra of degree 2 over Integers, y)

julia> f = x*y^2 + (x + 1)*y + 3
(x)*y^2+(x+([1 0]
[0 1]))*y+([3 0]
[0 3])

julia> g = reverse(f, 7)
([3 0]
[0 3])*y^6+(x+([1 0]
[0 1]))*y^5+(x)*y^4

julia> h = reverse(f)
([3 0]
[0 3])*y^2+(x+([1 0]
[0 1]))*y+(x)

```

### Shifting

```@docs
shift_left(::NCPolyElem, ::Int)
```

```@docs
shift_right(::NCPolyElem, ::Int)
```

**Examples**

```jldoctest
julia> R = MatrixAlgebra(ZZ, 2)
Matrix Algebra of degree 2 over Integers

julia> S, x = PolynomialRing(R, "x")
(Univariate Polynomial Ring in x over Matrix Algebra of degree 2 over Integers, x)

julia> T, y = PolynomialRing(S, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Matrix Algebra of degree 2 over Integers, y)

julia> f = x*y^2 + (x + 1)*y + 3
(x)*y^2+(x+([1 0]
[0 1]))*y+([3 0]
[0 3])

julia> g = shift_left(f, 7)
(x)*y^9+(x+([1 0]
[0 1]))*y^8+([3 0]
[0 3])*y^7

julia> h = shift_right(f, 2)
(x)

```

### Evaluation

```@docs
evaluate{T <: NCRingElem}(::NCPolyElem{T}, ::T)
evaluate(::NCPolyElem, ::Integer)
```

We also overload the functional notation so that the polynomial $f$ can be
evaluated at $a$ by writing $f(a)$. 

**Examples**

```jldoctest
julia> R = MatrixAlgebra(ZZ, 2)
Matrix Algebra of degree 2 over Integers

julia> S, x = PolynomialRing(R, "x")
(Univariate Polynomial Ring in x over Matrix Algebra of degree 2 over Integers, x)

julia> T, y = PolynomialRing(S, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Matrix Algebra of degree 2 over Integers, y)


julia> f = x*y^2 + (x + 1)*y + 3
(x)*y^2+(x+([1 0]
[0 1]))*y+([3 0]
[0 3])

julia> k = evaluate(f, 3)
([12 0]
[0 12])*x+([6 0]
[0 6])

julia> m = evaluate(f, x^2 + 2x + 1)
x^5+([4 0]
[0 4])*x^4+([7 0]
[0 7])*x^3+([7 0]
[0 7])*x^2+([4 0]
[0 4])*x+([4 0]
[0 4])

julia> r = f(23)
([552 0]
[0 552])*x+([26 0]
[0 26])

```

### Derivative

```@docs
derivative(::NCPolyElem)
```

**Examples**

```jldoctest
julia> R = MatrixAlgebra(ZZ, 2)
Matrix Algebra of degree 2 over Integers

julia> S, x = PolynomialRing(R, "x")
(Univariate Polynomial Ring in x over Matrix Algebra of degree 2 over Integers, x)

julia> T, y = PolynomialRing(S, "y")
(Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Matrix Algebra of degree 2 over Integers, y)

julia> f = x*y^2 + (x + 1)*y + 3
(x)*y^2+(x+([1 0]
[0 1]))*y+([3 0]
[0 3])

julia> h = derivative(f)
(([2 0]
[0 2])*x)*y+(x+([1 0]
[0 1]))

```

