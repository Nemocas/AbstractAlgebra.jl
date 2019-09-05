```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Generic matrices

AbstractAlgebra.jl allows the creation of dense matrices over any computable commutative
ring $R$. Generic matrices over a commutative ring are implemented in
`src/generic/Matrix.jl`. Much of the functionality there covers both matrix spaces and
matrix algebras.

Functions specific to generic matrix algebras of $m\times m$ matrices are implemented in
`src/generic/MatrixAlgebra.jl`.

As well as implementing the entire Matrix interface, including the optional
functionality, there are many additional generic algorithms implemented for matrix
spaces. We describe this functionality below.

All of this generic functionality is part of the Generic submodule of
AbstractAlgebra.jl. This is exported by default, so it is not necessary to qualify names
of functions.

## Types and parent objects

Generic matrices in AbstractAlgebra.jl have type `Generic.MatSpaceElem{T}` for matrices
in a matrix space, or `Generic.MatAlgElem{T}` for matrices in a matrix algebra, where
`T` is the type of elements of the matrix. Internally, generic matrices are implemented
using an object wrapping a Julia two dimensional array, though they are not themselves
Julia arrays. See the file `src/generic/GenericTypes.jl` for details.

For the most part, one doesn't want to work directly with the `MatSpaceElem` type though,
but with an abstract type called `Generic.Mat` which includes `MatSpaceElem` and views
thereof.

Parents of generic matrices (matrix spaces) have type `Generic.MatSpace{T}`. Parents of
matrices in a matrix algebra have type `Generic.MatAlgebra{T}`.

The generic matrix types (matrix spaces) belong to the abstract type
`AbstractAlgebra.MatElem{T}` and the matrix space parent types belong to
`AbstractAlgebra.MatSpace{T}`. Similarly the generic matrix algebra matrix types belong
to the abstract type `AbstractAlgebra.MatAlgElem{T}` and the parent types belong to
 `AbstractAlgebra.MatAlgebra{T}` Note that both
the concrete type of a matrix space parent object and the abstract class it belongs to
have the name `MatElem`, therefore disambiguation is required to specify which is
intended. The same is true for the abstract types for matrix spaces and their elements.

The dimensions and base ring $R$ of a generic matrix are stored in its parent object,
however to allow creation of matrices without first creating the matrix space parent,
generic matrices in Julia do not contain a reference to their parent. They contain the
row and column numbers (or degree, in the case of matrix algebras) and the base ring
on a per matrix basis. The parent object can then be reconstructed from this data on
demand.

## Matrix space constructors

A matrix space in AbstractAlgebra.jl represents a collection of all matrices with
given dimensions and base ring.

In order to construct matrices in AbstractAlgebra.jl, one can first construct the
matrix space itself. This is accomplished with the following constructor. We discuss
creation of matrix algebras separately in a dedicated section elsewhere in the
documentation.

```julia
MatrixSpace(R::Ring, rows::Int, cols::Int; cache::Bool=true)
```

Construct the space of matrices with the given number of rows and columns over the
given base ring. By default such matrix spaces are cached based on the base ring and
numbers of rows and columns. If the optional named parameter `cached` is set to false,
no caching occurs.

Here are some examples of creating matrix spaces and making use of the
resulting parent objects to coerce various elements into the matrix space.

**Examples**

```jldoctest
julia> R, t = PolynomialRing(QQ, "t")
(Univariate Polynomial Ring in t over Rationals, t)

julia> S = MatrixSpace(R, 3, 3)
Matrix Space of 3 rows and 3 columns over Univariate Polynomial Ring in t over Rationals

julia> A = S()
[0//1 0//1 0//1]
[0//1 0//1 0//1]
[0//1 0//1 0//1]

julia> B = S(12)
[12//1 0//1 0//1]
[0//1 12//1 0//1]
[0//1 0//1 12//1]

julia> C = S(R(11))
[11//1 0//1 0//1]
[0//1 11//1 0//1]
[0//1 0//1 11//1]

```

We also allow matrices over a given base ring to be constructed directly (see the
Matrix interface).

## Matrix element constructors

In addition to coercing elements into a matrix space as above, we provide the
following functions for constructing explicit matrices.

Also see the Matrix interface for a list of other ways to create matrices.

```julia
R[a b c...;...]
```

Create the matrix over the base ring $R$ consisting of the given rows (separated by
semicolons). Each entry is coerced into $R$  automatically. Note that parentheses may
be placed around individual entries if the lists would otherwise be ambiguous, e.g.
`R[1 2; 2 (-3)]`.

Beware that this syntax does not support the creation of column vectors. See
the notation below for creating those.

```julia
R[a b c...]
```

Create the row vector with entries in $R$ consisting of the given entries
(separated by spaces). Each entry is coerced into $R$ automatically. Note that
parentheses may be placed around individual entries if the list would otherwise
be ambiguous, e.g. `R[1 2 (-3)]`.

```julia
R[a b c...]'
```

Create the column vector with entries in $R$ consisting of the given entries
(separated by spaces). Each entry is coerced into $R$ automatically. Observe
the dash that is used to transpose the row vector notation (for free) to turn
it into a column vector. Note that parentheses may be placed around individual
entries if the list would otherwise be ambiguous, e.g. `R[1 2 (-3)]'`.

**Examples**

```jldoctest
julia> R, t = PolynomialRing(QQ, "t")
(Univariate Polynomial Ring in t over Rationals, t)

julia> S = MatrixSpace(R, 3, 3)
Matrix Space of 3 rows and 3 columns over Univariate Polynomial Ring in t over Rationals

julia> M = R[t + 1 1; t^2 0]
[t+1//1 1//1]
[t^2 0//1]

julia> N = R[t + 1 2 t]
[t+1//1 2//1 t]

julia> P = R[1 2 t]'
[1//1]
[2//1]
[t]

```

## Submatrices

In addition to the functionality described in the Matrix interface for taking
submatrices of a matrix, the following function variant is also available.

```@docs
sub(::MatElem, ::Int, ::Int, ::Int, ::Int)
```

**Examples**

```jldoctest
julia> M = ZZ[1 2 3; 2 3 4]
[1 2 3]
[2 3 4]

julia> N = sub(M, 1, 1, 2, 2)
[1 2]
[2 3]

```

## Matrix functionality provided by AbstractAlgebra.jl

Most of the following generic functionality is available for both matrix spaces and
matrix algebras. Exceptions include functions that do not return or accept square
matrices or which cannot specify a parent. Such functions include `solve`, `kernel`,
and `nullspace` which can't be provided for matrix algebras.

For details on functionality that is provided for matrix algebras only, see the dedicated
section of the documentation.

### Basic matrix functionality

As well as the Ring and Matrix interfaces, the following functions are provided to
manipulate matrices and to set and retrieve entries and other basic data associated
with the matrices.

```@docs
dense_matrix_type(::Ring)
```

```@docs
nrows(::MatElem)
```

```@docs
ncols(::MatElem)
```

```@docs
length(::MatElem)
```

```@docs
isempty(::MatElem)
```

```@docs
change_base_ring(::AbstractAlgebra.MatElem, ::AbstractAlgebra.Ring)
```

**Examples**

```jldoctest
julia> R, t = PolynomialRing(QQ, "t")
(Univariate Polynomial Ring in t over Rationals, t)

julia> S = MatrixSpace(R, 3, 3)
Matrix Space of 3 rows and 3 columns over Univariate Polynomial Ring in t over Rationals

julia> A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
[t+1//1 t 1//1]
[t^2 t t]
[-2//1 t+2//1 t^2+t+1//1]

julia> B = S([R(2) R(3) R(1); t t + 1 t + 2; R(-1) t^2 t^3])
[2//1 3//1 1//1]
[t t+1//1 t+2//1]
[-1//1 t^2 t^3]

julia> T = dense_matrix_type(R)
AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Poly{Rational{BigInt}}}

julia> r = nrows(B)
3

julia> c = ncols(B)
3

julia> length(B)
9

julia> isempty(B)
false

julia> M = A + B
[t+3//1 t+3//1 2//1]
[t^2+t 2//1*t+1//1 2//1*t+2//1]
[-3//1 t^2+t+2//1 t^3+t^2+t+1//1]

julia> N = 2 + A
[t+3//1 t 1//1]
[t^2 t+2//1 t]
[-2//1 t+2//1 t^2+t+3//1]

julia> M1 = deepcopy(A)
[t+1//1 t 1//1]
[t^2 t t]
[-2//1 t+2//1 t^2+t+1//1]

julia> A != B
true

julia> isone(one(S))
true

julia> V = A[1:2, :]
[t+1//1 t 1//1]
[t^2 t t]

julia> W = A^3
[3//1*t^4+4//1*t^3+t^2-3//1*t-5//1 t^4+5//1*t^3+10//1*t^2+7//1*t+4//1 2//1*t^4+7//1*t^3+9//1*t^2+8//1*t+1//1]
[t^5+4//1*t^4+3//1*t^3-7//1*t^2-4//1*t 4//1*t^4+8//1*t^3+7//1*t^2+2//1*t t^5+5//1*t^4+9//1*t^3+7//1*t^2-t]
[t^5+3//1*t^4-10//1*t^2-16//1*t-2//1 t^5+6//1*t^4+12//1*t^3+11//1*t^2+5//1*t-2//1 t^6+3//1*t^5+8//1*t^4+15//1*t^3+10//1*t^2+t-5//1]

julia> Z = divexact(2*A, 2)
[t+1//1 t 1//1]
[t^2 t t]
[-2//1 t+2//1 t^2+t+1//1]

```

### Powering

```@docs
powers(::MatElem, ::Int)
```

**Examples**

```jldoctest
julia> M = ZZ[1 2 3; 2 3 4; 4 5 5]
[1 2 3]
[2 3 4]
[4 5 5]

julia> A = powers(M, 4)
5-element Array{AbstractAlgebra.Generic.MatSpaceElem{BigInt},1}:
 [1 0 0]
[0 1 0]
[0 0 1]
 [1 2 3]
[2 3 4]
[4 5 5]
 [17 23 26]
[24 33 38]
[34 48 57]
 [167 233 273]
[242 337 394]
[358 497 579]
 [1725 2398 2798]
[2492 3465 4044]
[3668 5102 5957]

```

### Gram matrix

```@docs
gram(::MatElem)
```

**Examples**

```jldoctest
julia> R, t = PolynomialRing(QQ, "t")
(Univariate Polynomial Ring in t over Rationals, t)

julia> S = MatrixSpace(R, 3, 3)
Matrix Space of 3 rows and 3 columns over Univariate Polynomial Ring in t over Rationals

julia> A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
[t+1//1 t 1//1]
[t^2 t t]
[-2//1 t+2//1 t^2+t+1//1]

julia> B = gram(A)
[2//1*t^2+2//1*t+2//1 t^3+2//1*t^2+t 2//1*t^2+t-1//1]
[t^3+2//1*t^2+t t^4+2//1*t^2 t^3+3//1*t]
[2//1*t^2+t-1//1 t^3+3//1*t t^4+2//1*t^3+4//1*t^2+6//1*t+9//1]

```

### Trace

```@docs
tr(::MatElem)
```

**Examples**

```jldoctest
julia> R, t = PolynomialRing(QQ, "t")
(Univariate Polynomial Ring in t over Rationals, t)

julia> S = MatrixSpace(R, 3, 3)
Matrix Space of 3 rows and 3 columns over Univariate Polynomial Ring in t over Rationals

julia> A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
[t+1//1 t 1//1]
[t^2 t t]
[-2//1 t+2//1 t^2+t+1//1]

julia> b = tr(A)
t^2+3//1*t+2//1

```

### Content

```@docs
content(::MatElem)
```

**Examples**

```jldoctest
julia> R, t = PolynomialRing(QQ, "t")
(Univariate Polynomial Ring in t over Rationals, t)

julia> S = MatrixSpace(R, 3, 3)
Matrix Space of 3 rows and 3 columns over Univariate Polynomial Ring in t over Rationals

julia> A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
[t+1//1 t 1//1]
[t^2 t t]
[-2//1 t+2//1 t^2+t+1//1]

julia> b = content(A)
1//1

```

### Permutation

```@docs
*(::perm, ::MatElem)
```

**Examples**

```jldoctest
julia> R, t = PolynomialRing(QQ, "t")
(Univariate Polynomial Ring in t over Rationals, t)

julia> S = MatrixSpace(R, 3, 3)
Matrix Space of 3 rows and 3 columns over Univariate Polynomial Ring in t over Rationals

julia> G = PermGroup(3)
Permutation group over 3 elements

julia> A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
[t+1//1 t 1//1]
[t^2 t t]
[-2//1 t+2//1 t^2+t+1//1]

julia> P = G([1, 3, 2])
(2,3)

julia> B = P*A
[t+1//1 t 1//1]
[-2//1 t+2//1 t^2+t+1//1]
[t^2 t t]

```

### LU factorisation

```@docs
lu{T <: FieldElem}(::MatElem{T}, ::PermGroup)
```

```@docs
fflu{T <: RingElem}(::MatElem{T}, ::PermGroup)
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(QQ, "x")
(Univariate Polynomial Ring in x over Rationals, x)

julia> K, a = NumberField(x^3 + 3x + 1, "a")
(Residue field of Univariate Polynomial Ring in x over Rationals modulo x^3+3//1*x+1//1, x)

julia> S = MatrixSpace(K, 3, 3)
Matrix Space of 3 rows and 3 columns over Residue field of Univariate Polynomial Ring in x over Rationals modulo x^3+3//1*x+1//1

julia> A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 - 2 a - 1 2a])
[0//1 2//1*x+3//1 x^2+1//1]
[x^2-2//1 x-1//1 2//1*x]
[x^2-2//1 x-1//1 2//1*x]

julia> r, P, L, U = lu(A)
(2, (1,2), [1//1 0//1 0//1]
[0//1 1//1 0//1]
[1//1 0//1 1//1], [x^2-2//1 x-1//1 2//1*x]
[0//1 2//1*x+3//1 x^2+1//1]
[0//1 0//1 0//1])

julia> r, d, P, L, U = fflu(A)
(2, 3//1*x^2-10//1*x-8//1, (1,2), [x^2-2//1 0//1 0//1]
[0//1 3//1*x^2-10//1*x-8//1 0//1]
[x^2-2//1 0//1 1//1], [x^2-2//1 x-1//1 2//1*x]
[0//1 3//1*x^2-10//1*x-8//1 -4//1*x^2-x-2//1]
[0//1 0//1 0//1])

```

### Reduced row-echelon form

```@docs
rref{T <: RingElem}(::MatElem{T})
rref{T <: FieldElem}(::MatElem{T})
```

```@docs
isrref{T <: RingElem}(::MatElem{T})
isrref{T <: FieldElem}(::MatElem{T})
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(QQ, "x")
(Univariate Polynomial Ring in x over Rationals, x)

julia> K, a = NumberField(x^3 + 3x + 1, "a")
(Residue field of Univariate Polynomial Ring in x over Rationals modulo x^3+3//1*x+1//1, x)

julia> S = MatrixSpace(K, 3, 3)
Matrix Space of 3 rows and 3 columns over Residue field of Univariate Polynomial Ring in x over Rationals modulo x^3+3//1*x+1//1

julia> M = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])
[0//1 2//1*x+3//1 x^2+1//1]
[x^2-2//1 x-1//1 2//1*x]
[x^2+3//1*x+1//1 2//1*x 1//1]

julia> r, A = rref(M)
(3, [1//1 0//1 0//1]
[0//1 1//1 0//1]
[0//1 0//1 1//1])

julia> isrref(A)
true

julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S = MatrixSpace(R, 3, 3)
Matrix Space of 3 rows and 3 columns over Univariate Polynomial Ring in x over Integers

julia> M = S([R(0) 2x + 3 x^2 + 1; x^2 - 2 x - 1 2x; x^2 + 3x + 1 2x R(1)])
[0 2*x+3 x^2+1]
[x^2-2 x-1 2*x]
[x^2+3*x+1 2*x 1]

julia> r, d, A = rref(M)
(3, -x^5-2*x^4-15*x^3-18*x^2-8*x-7, [-x^5-2*x^4-15*x^3-18*x^2-8*x-7 0 0]
[0 -x^5-2*x^4-15*x^3-18*x^2-8*x-7 0]
[0 0 -x^5-2*x^4-15*x^3-18*x^2-8*x-7])

julia> isrref(A)
true
```

### Hermite normal form

```@docs
hnf{T <: RingElem}(::MatElem{T})
hnf_with_transform{T <: RingElem}(::MatElem{T})
```

```@docs
ishnf{T <: RingElem}(::MatElem{T})
```

### Determinant

```@docs
det{T <: RingElem}(::MatElem{T})
det{T <: FieldElem}(::MatElem{T})
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(QQ, "x")
(Univariate Polynomial Ring in x over Rationals, x)

julia> K, a = NumberField(x^3 + 3x + 1, "a")
(Residue field of Univariate Polynomial Ring in x over Rationals modulo x^3+3//1*x+1//1, x)

julia> S = MatrixSpace(K, 3, 3)
Matrix Space of 3 rows and 3 columns over Residue field of Univariate Polynomial Ring in x over Rationals modulo x^3+3//1*x+1//1

julia> A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])
[0//1 2//1*x+3//1 x^2+1//1]
[x^2-2//1 x-1//1 2//1*x]
[x^2+3//1*x+1//1 2//1*x 1//1]

julia> d = det(A)
11//1*x^2-30//1*x-5//1

```

### Rank

```@docs
rank{T <: RingElem}(::MatElem{T})
rank{T <: FieldElem}(::MatElem{T})
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(QQ, "x")
(Univariate Polynomial Ring in x over Rationals, x)

julia> K, a = NumberField(x^3 + 3x + 1, "a")
(Residue field of Univariate Polynomial Ring in x over Rationals modulo x^3+3//1*x+1//1, x)

julia> S = MatrixSpace(K, 3, 3)
Matrix Space of 3 rows and 3 columns over Residue field of Univariate Polynomial Ring in x over Rationals modulo x^3+3//1*x+1//1

julia> A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])
[0//1 2//1*x+3//1 x^2+1//1]
[x^2-2//1 x-1//1 2//1*x]
[x^2+3//1*x+1//1 2//1*x 1//1]

julia> d = rank(A)
3

```

### Linear solving

```@docs
solve{T <: FieldElem}(::MatElem{T}, ::MatElem{T})
```

```@docs
solve_rational{T <: RingElem}(::MatElem{T}, ::MatElem{T})
```

```@docs
solve_left{T <: RingElem}(::MatElem{T}, ::MatElem{T})
```

```@docs
solve_triu{T <: FieldElem}(::MatElem{T}, ::MatElem{T}, ::Bool)
```

```@docs
can_solve_left_reduced_triu{T <: RingElement}(::MatElem{T}, ::MatElem{T})
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(QQ, "x")
(Univariate Polynomial Ring in x over Rationals, x)

julia> K, a = NumberField(x^3 + 3x + 1, "a")
(Residue field of Univariate Polynomial Ring in x over Rationals modulo x^3+3//1*x+1//1, x)

julia> S = MatrixSpace(K, 3, 3)
Matrix Space of 3 rows and 3 columns over Residue field of Univariate Polynomial Ring in x over Rationals modulo x^3+3//1*x+1//1

julia> U = MatrixSpace(K, 3, 1)
Matrix Space of 3 rows and 1 columns over Residue field of Univariate Polynomial Ring in x over Rationals modulo x^3+3//1*x+1//1

julia> A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])
[0//1 2//1*x+3//1 x^2+1//1]
[x^2-2//1 x-1//1 2//1*x]
[x^2+3//1*x+1//1 2//1*x 1//1]

julia> b = U([2a a + 1 (-a - 1)]')
[2//1*x]
[x+1//1]
[-x-1//1]

julia> x = solve(A, b)
[1984//7817*x^2+1573//7817*x-937//7817]
[-2085//7817*x^2+1692//7817*x+965//7817]
[-3198//7817*x^2+3540//7817*x-3525//7817]

julia> A = S([a + 1 2a + 3 a^2 + 1; K(0) a^2 - 1 2a; K(0) K(0) a])
[x+1//1 2//1*x+3//1 x^2+1//1]
[0//1 x^2-1//1 2//1*x]
[0//1 0//1 x]

julia> bb = U([2a a + 1 (-a - 1)]')
[2//1*x]
[x+1//1]
[-x-1//1]

julia> x = solve_triu(A, bb, false)
[1//3*x^2+8//3*x+13//3]
[-3//5*x^2-3//5*x-12//5]
[x^2+2//1]

julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S = MatrixSpace(R, 3, 3)
Matrix Space of 3 rows and 3 columns over Univariate Polynomial Ring in x over Integers

julia> U = MatrixSpace(R, 3, 2)
Matrix Space of 3 rows and 2 columns over Univariate Polynomial Ring in x over Integers

julia> A = S([R(0) 2x + 3 x^2 + 1; x^2 - 2 x - 1 2x; x^2 + 3x + 1 2x R(1)])
[0 2*x+3 x^2+1]
[x^2-2 x-1 2*x]
[x^2+3*x+1 2*x 1]

julia> bbb = U([2x x + 1 (-x - 1); x + 1 (-x) x^2]')
[2*x x+1]
[x+1 -x]
[-x-1 x^2]

julia> x, d = solve_rational(A, bbb)
([3*x^4-10*x^3-8*x^2-11*x-4 -x^5+3*x^4+x^3-2*x^2+3*x-1]
[-2*x^5-x^4+6*x^3+2*x+1 x^6+x^5+4*x^4+9*x^3+8*x^2+5*x+2]
[6*x^4+12*x^3+15*x^2+6*x-3 -2*x^5-4*x^4-6*x^3-9*x^2-4*x+1], x^5+2*x^4+15*x^3+18*x^2+8*x+7)

julia> S = MatrixSpace(ZZ, 3, 3)
Matrix Space of 3 rows and 3 columns over Integers

julia> T = MatrixSpace(ZZ, 3, 1)
Matrix Space of 3 rows and 1 columns over Integers

julia> A = S([BigInt(2) 3 5; 1 4 7; 9 2 2])
[2 3 5]
[1 4 7]
[9 2 2]

julia> B = T([BigInt(4), 5, 7])
[4]
[5]
[7]
```

### Inverse

```@docs
inv{T <: RingElem}(::MatElem{T})
inv{T <: FieldElem}(::MatElem{T})
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(QQ, "x")
(Univariate Polynomial Ring in x over Rationals, x)

julia> K, a = NumberField(x^3 + 3x + 1, "a")
(Residue field of Univariate Polynomial Ring in x over Rationals modulo x^3+3//1*x+1//1, x)

julia> S = MatrixSpace(K, 3, 3)
Matrix Space of 3 rows and 3 columns over Residue field of Univariate Polynomial Ring in x over Rationals modulo x^3+3//1*x+1//1

julia> A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])
[0//1 2//1*x+3//1 x^2+1//1]
[x^2-2//1 x-1//1 2//1*x]
[x^2+3//1*x+1//1 2//1*x 1//1]

julia> X = inv(A)
[-343//7817*x^2+717//7817*x-2072//7817 -4964//23451*x^2+2195//23451*x-11162//23451 -232//23451*x^2-4187//23451*x-1561//23451]
[128//7817*x^2-655//7817*x+2209//7817 599//23451*x^2-2027//23451*x-1327//23451 -1805//23451*x^2+2702//23451*x-7394//23451]
[545//7817*x^2+570//7817*x+2016//7817 -1297//23451*x^2-5516//23451*x-337//23451 8254//23451*x^2-2053//23451*x+16519//23451]

julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S = MatrixSpace(R, 3, 3)
Matrix Space of 3 rows and 3 columns over Univariate Polynomial Ring in x over Integers

julia> A = S([R(0) 2x + 3 x^2 + 1; x^2 - 2 x - 1 2x; x^2 + 3x + 1 2x R(1)])
[0 2*x+3 x^2+1]
[x^2-2 x-1 2*x]
[x^2+3*x+1 2*x 1]

julia> X, d = pseudo_inv(A)
([4*x^2-x+1 -2*x^3+3 x^3-5*x^2-5*x-1]
[-2*x^3-5*x^2-2*x-2 x^4+3*x^3+2*x^2+3*x+1 -x^4+x^2+2]
[-x^3+2*x^2+2*x-1 -2*x^3-9*x^2-11*x-3 2*x^3+3*x^2-4*x-6], -x^5-2*x^4-15*x^3-18*x^2-8*x-7)

```

### Nullspace

```@docs
nullspace{T <: RingElem}(::MatElem{T})
nullspace{T <: FieldElem}(::MatElem{T})
```

**Examples**

```jldoctest
julia> R, x = PolynomialRing(ZZ, "x")
(Univariate Polynomial Ring in x over Integers, x)

julia> S = MatrixSpace(R, 4, 4)
Matrix Space of 4 rows and 4 columns over Univariate Polynomial Ring in x over Integers

julia> M = S([-6*x^2+6*x+12 -12*x^2-21*x-15 -15*x^2+21*x+33 -21*x^2-9*x-9;
              -8*x^2+8*x+16 -16*x^2+38*x-20 90*x^2-82*x-44 60*x^2+54*x-34;
              -4*x^2+4*x+8 -8*x^2+13*x-10 35*x^2-31*x-14 22*x^2+21*x-15;
              -10*x^2+10*x+20 -20*x^2+70*x-25 150*x^2-140*x-85 105*x^2+90*x-50])
[-6*x^2+6*x+12 -12*x^2-21*x-15 -15*x^2+21*x+33 -21*x^2-9*x-9]
[-8*x^2+8*x+16 -16*x^2+38*x-20 90*x^2-82*x-44 60*x^2+54*x-34]
[-4*x^2+4*x+8 -8*x^2+13*x-10 35*x^2-31*x-14 22*x^2+21*x-15]
[-10*x^2+10*x+20 -20*x^2+70*x-25 150*x^2-140*x-85 105*x^2+90*x-50]

julia> n, N = nullspace(M)
(2, [1320*x^4-330*x^2-1320*x-1320 1056*x^4+1254*x^3+1848*x^2-66*x-330]
[-660*x^4+1320*x^3+1188*x^2-1848*x-1056 -528*x^4+132*x^3+1584*x^2+660*x-264]
[396*x^3-396*x^2-792*x 0]
[0 396*x^3-396*x^2-792*x])
```

### Kernel

```@docs
kernel{T <: RingElem}(::MatElem{T})
left_kernel{T <: RingElem}(::MatElem{T})
right_kernel{T <: RingElem}(::MatElem{T})
```

**Examples**

```jldoctest
julia> S = MatrixSpace(ZZ, 4, 4)
Matrix Space of 4 rows and 4 columns over Integers

julia> M = S([1 2 0 4;
              2 0 1 1;
              0 1 1 -1;
              2 -1 0 2])
[1 2 0 4]
[2 0 1 1]
[0 1 1 -1]
[2 -1 0 2]

julia> nr, Nr = kernel(M)
(1, [-8]
[-6]
[11]
[5])

julia> nl, Nl = left_kernel(M)
(1, [0 -1 1 1])

```

### Hessenberg form

```@docs
hessenberg{T <: RingElem}(::MatElem{T})
```

```@docs
ishessenberg{T <: RingElem}(::MatElem{T})
```

**Examples**

```jldoctest
julia> R = ResidueRing(ZZ, 7)
Residue ring of Integers modulo 7

julia> S = MatrixSpace(R, 4, 4)
Matrix Space of 4 rows and 4 columns over Residue ring of Integers modulo 7

julia> M = S([R(1) R(2) R(4) R(3); R(2) R(5) R(1) R(0);
              R(6) R(1) R(3) R(2); R(1) R(1) R(3) R(5)])
[1 2 4 3]
[2 5 1 0]
[6 1 3 2]
[1 1 3 5]

julia> A = hessenberg(M)
[1 5 5 3]
[2 1 1 0]
[0 1 3 2]
[0 0 2 2]

julia> ishessenberg(A)
true

```

### Characteristic polynomial

```@docs
charpoly{T <: RingElem}(::Ring, ::MatElem{T})
```

**Examples**

```jldoctest
julia> R = ResidueRing(ZZ, 7)
Residue ring of Integers modulo 7

julia> S = MatrixSpace(R, 4, 4)
Matrix Space of 4 rows and 4 columns over Residue ring of Integers modulo 7

julia> T, x = PolynomialRing(R, "x")
(Univariate Polynomial Ring in x over Residue ring of Integers modulo 7, x)

julia> M = S([R(1) R(2) R(4) R(3); R(2) R(5) R(1) R(0);
              R(6) R(1) R(3) R(2); R(1) R(1) R(3) R(5)])
[1 2 4 3]
[2 5 1 0]
[6 1 3 2]
[1 1 3 5]

julia> A = charpoly(T, M)
x^4+2*x^2+6*x+2

```

### Minimal polynomial

```@docs
minpoly{T <: RingElem}(::Ring, ::MatElem{T}, ::Bool)
minpoly{T <: FieldElem}(::Ring, ::MatElem{T}, ::Bool)
```

**Examples**

```jldoctest
julia> R = GF(13)
Finite field F_13

julia> T, y = PolynomialRing(R, "y")
(Univariate Polynomial Ring in y over Finite field F_13, y)

julia> M = R[7 6 1;
             7 7 5;
             8 12 5]
[7 6 1]
[7 7 5]
[8 12 5]

julia> A = minpoly(T, M)
y^2+10*y

```

### Transforms

```@docs
similarity!{T <: RingElem}(::MatElem{T}, ::Int, ::T)
```

**Examples**

```jldoctest
julia> R = ResidueRing(ZZ, 7)
Residue ring of Integers modulo 7

julia> S = MatrixSpace(R, 4, 4)
Matrix Space of 4 rows and 4 columns over Residue ring of Integers modulo 7

julia> M = S([R(1) R(2) R(4) R(3); R(2) R(5) R(1) R(0);
              R(6) R(1) R(3) R(2); R(1) R(1) R(3) R(5)])
[1 2 4 3]
[2 5 1 0]
[6 1 3 2]
[1 1 3 5]

julia> similarity!(M, 1, R(3))

```

### (Weak) Popov form

AbstractAlgebra.jl provides algorithms for computing the (weak) Popov of a matrix with
entries in a univariate polynomial ring over a field.

```@docs
weak_popov{T <: PolyElem}(::Generic.Mat{T})
weak_popov_with_transform{T <: PolyElem}(::Generic.Mat{T})
popov{T <: PolyElem}(::Generic.Mat{T})
popov_with_transform{T <: PolyElem}(::Generic.Mat{T})
```
