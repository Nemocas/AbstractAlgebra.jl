```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# Matrix Spaces

## Matrix Space Constructors

A matrix space in AbstractAlgebra.jl represents a collection of all matrices with
given dimensions and base ring.

In order to construct matrices in AbstractAlgebra.jl, one can first construct the
matrix space itself. This is accomplished with the following constructor. We discuss
creation of matrix algebras separately in a dedicated section elsewhere in the
documentation.

```julia
matrix_space(R::Ring, rows::Int, cols::Int)
```

Construct the space of matrices with the given number of rows and columns over the
given base ring.

Here are some examples of creating matrix spaces and making use of the
resulting parent objects to coerce various elements into the matrix space.

**Examples**

```jldoctest
julia> R, t = polynomial_ring(QQ, :t)
(Univariate polynomial ring in t over rationals, t)

julia> S = matrix_space(R, 3, 3)
Matrix space of 3 rows and 3 columns
  over univariate polynomial ring in t over rationals

julia> A = S()
[0   0   0]
[0   0   0]
[0   0   0]

julia> B = S(12)
[12    0    0]
[ 0   12    0]
[ 0    0   12]

julia> C = S(R(11))
[11    0    0]
[ 0   11    0]
[ 0    0   11]

```

## Basic matrix space functionality

```@docs
number_of_rows(::MatSpace)
number_of_columns(::MatSpace)
zero(::MatSpace)
one(::MatSpace)
```

```jldoctest
julia> R, t = polynomial_ring(QQ, :t)
(Univariate polynomial ring in t over rationals, t)

julia> S = matrix_space(R, 3, 3)
Matrix space of 3 rows and 3 columns
  over univariate polynomial ring in t over rationals

julia> A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
[t + 1       t             1]
[  t^2       t             t]
[   -2   t + 2   t^2 + t + 1]

julia> B = S([R(2) R(3) R(1); t t + 1 t + 2; R(-1) t^2 t^3])
[ 2       3       1]
[ t   t + 1   t + 2]
[-1     t^2     t^3]

julia> r = number_of_rows(B)
3

julia> c = number_of_columns(B)
3

julia> length(B)
9

julia> isempty(B)
false

julia> M = A + B
[  t + 3         t + 3                   2]
[t^2 + t       2*t + 1             2*t + 2]
[     -3   t^2 + t + 2   t^3 + t^2 + t + 1]

julia> N = 2 + A
[t + 3       t             1]
[  t^2   t + 2             t]
[   -2   t + 2   t^2 + t + 3]

julia> M1 = deepcopy(A)
[t + 1       t             1]
[  t^2       t             t]
[   -2   t + 2   t^2 + t + 1]

julia> A != B
true

julia> isone(one(S))
true

julia> V = A[1:2, :]
[t + 1   t   1]
[  t^2   t   t]

julia> W = A^3
[    3*t^4 + 4*t^3 + t^2 - 3*t - 5            t^4 + 5*t^3 + 10*t^2 + 7*t + 4                 2*t^4 + 7*t^3 + 9*t^2 + 8*t + 1]
[t^5 + 4*t^4 + 3*t^3 - 7*t^2 - 4*t               4*t^4 + 8*t^3 + 7*t^2 + 2*t                 t^5 + 5*t^4 + 9*t^3 + 7*t^2 - t]
[  t^5 + 3*t^4 - 10*t^2 - 16*t - 2   t^5 + 6*t^4 + 12*t^3 + 11*t^2 + 5*t - 2   t^6 + 3*t^5 + 8*t^4 + 15*t^3 + 10*t^2 + t - 5]

julia> Z = divexact(2*A, 2)
[t + 1       t             1]
[  t^2       t             t]
[   -2   t + 2   t^2 + t + 1]

julia> M = matrix(ZZ, BigInt[2 3 0; 1 1 1])
[2   3   0]
[1   1   1]

julia> M[1, 2] = BigInt(4)
4

julia> c = M[1, 1]
2
```

## Matrix element constructors

There are a few ways to construct matrices other than by coercing elements
as shown above. The first method is from an array of elements.

This can be done with either two or one dimensional arrays.

```julia
(S::MatSpace{T})(A::Matrix{S}) where {S <: RingElement, T <: RingElement}
(S::MatRing{T})(A::Matrix{S}) where {S <: RingElement, T <: RingElement}
```

Create the matrix in the given space/algebra whose $(i, j)$ entry is given by `A[i, j]`,
where `S` is the type of elements that can be coerced into the base ring of the matrix.

```julia
(S::MyMatSpace{T})(A::Vector{S}) where {S <: RingElem, T <: RingElem}
(S::MyMatAlgebra{T})(A::Vector{S}) where {S <: RingElem, T <: RingElem}
```

Create the matrix in the given space/algebra of matrices (with dimensions $m\times n$
say), whose $(i, j)$ entry is given by `A[i*(n - 1) + j]` and where `S` is the type of
elements that can be coerced into the base ring of the matrix.

We also provide the
following syntax for constructing literal matrices (similar to how Julia
arrays can be be constructed).

```julia
R[a b c...;...]
```

Create the matrix over the base ring $R$ consisting of the given rows (separated by
semicolons). Each entry is coerced into $R$  automatically. Note that parentheses may
be placed around individual entries if the lists would otherwise be ambiguous, e.g.
`R[1 2; 2 (- 3)]`.

Also see the Matrix interface for a list of other ways to create matrices.

**Examples**

```jldoctest
julia> S = matrix_space(QQ, 2, 3)
Matrix space of 2 rows and 3 columns
  over rationals

julia> T = matrix_ring(QQ, 2)
Matrix ring of degree 2
  over rationals

julia> M1 = S(Rational{BigInt}[2 3 1; 1 0 4])
[2//1   3//1   1//1]
[1//1   0//1   4//1]

julia> M2 = S(BigInt[2 3 1; 1 0 4])
[2//1   3//1   1//1]
[1//1   0//1   4//1]

julia> M3 = S(BigInt[2, 3, 1, 1, 0, 4])
[2//1   3//1   1//1]
[1//1   0//1   4//1]

julia> N1 = T(Rational{BigInt}[2 3; 1 0])
[2//1   3//1]
[1//1   0//1]

julia> N2 = T(BigInt[2 3; 1 0])
[2//1   3//1]
[1//1   0//1]

julia> N3 = T(BigInt[2, 3, 1, 1])
[2//1   3//1]
[1//1   1//1]

julia> R, t = polynomial_ring(QQ, :t)
(Univariate polynomial ring in t over rationals, t)

julia> S = matrix_space(R, 3, 3)
Matrix space of 3 rows and 3 columns
  over univariate polynomial ring in t over rationals

julia> M = R[t + 1 1; t^2 0]
[t + 1   1]
[  t^2   0]

julia> N = R[t + 1 2 t] # create a row vector
[t + 1   2   t]

julia> P = R[1; 2; t] # create a column vector
[1]
[2]
[t]
```
