```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Matrix functionality

AbstractAlgebra.jl provides a module, implemented in `src/Matrix.jl` for
matrices over any ring belonging to the AbstractAlgebra abstract type
hierarchy. This functionality will work for any matrix type which
follows the Matrix interface.

Similarly, AbstractAlgebra.jl provides a module in `src/MatrixAlgebra.jl` for
matrix algebras over a ring.

## Generic matrix types

AbstractAlgebra.jl allows the creation of dense matrices over any computable
ring $R$. Generic matrices over a ring are implemented in
`src/generic/Matrix.jl`.

Generic matrix algebras of $m\times m$ matrices are implemented in
`src/generic/MatrixAlgebra.jl`.

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

The dimensions and base ring $R$ of a generic matrix are stored in its parent object,
however to allow creation of matrices without first creating the matrix space parent,
generic matrices in Julia do not contain a reference to their parent. They contain the
row and column numbers (or degree, in the case of matrix algebras) and the base ring
on a per matrix basis. The parent object can then be reconstructed from this data on
demand.

## Abstract types

The generic matrix types (matrix spaces) belong to the abstract type
`MatElem{T}` and the matrix space parent types belong to
`MatSpace{T}`. Similarly the generic matrix algebra matrix types belong
to the abstract type `MatAlgElem{T}` and the parent types belong to
 `MatAlgebra{T}` Note that both
the concrete type of a matrix space parent object and the abstract class it belongs to
have the name `MatElem`, therefore disambiguation is required to specify which is
intended. The same is true for the abstract types for matrix spaces and their elements.

## Matrix space constructors

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
julia> R, t = polynomial_ring(QQ, "t")
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

## Matrix element constructors

There are a few ways to construct matrices other than by coercing elements
as shown above. The first method is from an array of elements.

This can be done with either two or one dimensional arrays.

```julia
(S::MatSpace{T})(A::Matrix{S}) where {S <: RingElement, T <: RingElement}
(S::MatAlgebra{T})(A::Matrix{S}) where {S <: RingElement, T <: RingElement}
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

julia> T = MatrixAlgebra(QQ, 2)
Matrix algebra of degree 2
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

julia> R, t = polynomial_ring(QQ, "t")
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

It is also possible to create matrices (in a matrix space only) directly, without first
creating the corresponding matrix space (the inner constructor being called directly).

```julia
matrix(R::Ring, arr::Matrix{T}) where T <: RingElement
```

Given an $m\times n$ Julia matrix of entries, construct the corresponding
AbstractAlgebra.jl matrix over the given ring `R`, assuming all the entries can be
coerced into `R`.


```julia
matrix(R::Ring, r::Int, c::Int, A::Vector{T}) where T <: RingElement
```

Construct the given $r\times c$ AbstractAlgebra.jl matrix over the ring `R` whose
$(i, j)$ entry is given by `A[c*(i - 1) + j]`, assuming that all the entries can be
coerced into `R`.

```julia
zero_matrix(R::Ring, r::Int, c::Int)
```

Construct the $r\times c$ AbstractAlgebra.jl zero matrix over the ring `R`.

**Examples**

```jldoctest
julia> M = matrix(ZZ, BigInt[3 1 2; 2 0 1])
[3   1   2]
[2   0   1]

julia> N = matrix(ZZ, 3, 2, BigInt[3, 1, 2, 2, 0, 1])
[3   1]
[2   2]
[0   1]

julia> P = zero_matrix(ZZ, 3, 2)
[0   0]
[0   0]
[0   0]

julia> R = MatrixAlgebra(ZZ, 2)
Matrix algebra of degree 2
  over integers

julia> M = R()
[0   0]
[0   0]
```

## Block diagonal matrix constructors

It is also possible to create block diagonal matrices from a vector of
existing matrices. It is also possible to construct them from Julia
matrices if one supplies the base ring.

Note that if the input matrices are not square, the output matrix may
not be square.

```@docs
block_diagonal_matrix(::Vector{<:MatElem{T}}) where T <: RingElement
block_diagonal_matrix(::Ring, ::Vector{<:Matrix{T}}) where T <: RingElement
```

**Examples**

```jldoctest
julia> block_diagonal_matrix(ZZ, [[1 2; 3 4], [4 5 6; 7 8 9]])
[1   2   0   0   0]
[3   4   0   0   0]
[0   0   4   5   6]
[0   0   7   8   9]

julia> M = matrix(ZZ, [1 2; 3 4])
[1   2]
[3   4]

julia> N = matrix(ZZ, [4 5 6; 7 8 9])
[4   5   6]
[7   8   9]

julia> block_diagonal_matrix([M, N])
[1   2   0   0   0]
[3   4   0   0   0]
[0   0   4   5   6]
[0   0   7   8   9]
```

## Conversion to Julia matrices, iteration and broacasting

While `AbstractAlgebra` matrices are not instances of `AbstractArray`,
they are closely related to Julia matrices. For convenience, a `Matrix`
and an `Array` constructors taking an `AbstractAlgebra` matrix as input
are provided:

```@docs
Matrix(::MatrixElem{T}) where T <: RingElement
Array(::MatrixElem{T}) where T <: RingElement
```

Matrices also support iteration, and therefore functions accepting an iterator
can be called on them, e.g.:

```jldoctest
julia> M = matrix_space(ZZ, 2, 3); x = M(1:6)
[1   2   3]
[4   5   6]

julia> collect(x)
2×3 Matrix{BigInt}:
 1  2  3
 4  5  6

julia> Set(x)
Set{BigInt} with 6 elements:
  5
  4
  6
  2
  3
  1
```

Matrices also support broadcasting, which amounts to elementwise application
of functions to matrices:

```jldoctest
julia> k = GF(5);

julia> A = ZZ[1 2; 3 4];

julia> k.(A)
[1   2]
[3   4]

julia> 3 .* A .+ 2
[ 5    8]
[11   14]

julia> B = ZZ[3 4; 5 6];

julia> ((x, y) -> x^2 + y^2).(A, B)
[10   20]
[34   52]
```

### Views

As per Julia, AbstractAlgebra supports the construction of matrix views. 
These allow one to work with a submatrix of a given
matrix. Modifying the submatrix also modifies the original matrix.

The syntax for views is as for Julia's own views.

**Examples**

```jldoctest
julia> M = matrix(ZZ, 3, 3, BigInt[1, 2, 3, 2, 3, 4, 3, 4, 5])
[1   2   3]
[2   3   4]
[3   4   5]

julia> N1 = @view M[1:2, :]
[1   2   3]
[2   3   4]

julia> N2 = @view M[:, 1:2]
[1   2]
[2   3]
[3   4]

julia> R = N1*N2
[14   20]
[20   29]
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
number_of_rows(::MatSpace)
number_of_columns(::MatSpace)
```

```@docs
number_of_rows(::MatrixElem{T}) where T <: RingElement
number_of_columns(::MatrixElem{T}) where T <: RingElement
```

```@docs
length(::MatrixElem{T}) where T <: RingElement
```

```@docs
isempty(::MatrixElem{T}) where T <: RingElement
```

```@docs
identity_matrix(::Ring, ::Int)
```

```@docs
identity_matrix(::MatElem{T}) where T <: RingElement
```

```@docs
scalar_matrix(R::Ring, n::Int, a::RingElement)
```

```@docs
diagonal_matrix(::RingElement, ::Int, ::Int)
```

```@docs
zero(::MatSpace)
zero(::MatrixElem{T}, ::Ring) where T <: RingElement
```

```@docs
one(::MatSpace)
one(::MatElem{T}) where T <: RingElement
```

```@docs
lower_triangular_matrix(L::AbstractVector{T}) where {T <: RingElement}
```

```@docs
upper_triangular_matrix(L::AbstractVector{T}) where {T <: RingElement}
```

```@docs
strictly_lower_triangular_matrix(L::AbstractVector{T}) where {T <: RingElement}
```

```@docs
strictly_upper_triangular_matrix(L::AbstractVector{T}) where {T <: RingElement}
```

```@docs
is_lower_triangular(::MatrixElem)
```

```@docs
is_upper_triangular(::MatrixElem)
```

```@docs
is_diagonal(::MatrixElem)
```

```@docs
change_base_ring(::Ring, ::MatElem{T}) where T <: RingElement
```

```@docs
Base.map(f, ::MatrixElem{T}) where T <: RingElement
```

```@docs
Base.map!(f, ::MatrixElem{S}, ::MatrixElem{T}) where {S <: RingElement, T <: RingElement}
```

**Examples**

```jldoctest
julia> R, t = polynomial_ring(QQ, "t")
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

julia> T = dense_matrix_type(R)
AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Poly{Rational{BigInt}}}

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

### Transpose

```@docs
transpose(::MatrixElem{T}) where T <: RingElement
```

### Submatrices

Submatrices are only available for matrix spaces, not for matrix
algebras and generally only available for generic matrices built
on Julia arrays.

Submatrices return a new matrix with the same entries as the
submatrix with the given range of rows and columns. They are best
illustrated with examples.

**Examples**

```jldoctest
julia> M = matrix(ZZ, BigInt[1 2 3; 2 3 4; 3 4 5])
[1   2   3]
[2   3   4]
[3   4   5]

julia> N1 = M[1:2, :]
[1   2   3]
[2   3   4]

julia> N2 = M[:, :]
[1   2   3]
[2   3   4]
[3   4   5]

julia> N3 = M[2:3, 2:3]
[3   4]
[4   5]

```

### Elementary row and column operations

```@docs
add_column(::MatElem{T}, ::Int, ::Int, ::Int) where T <: RingElement
add_column!(::MatElem{T}, ::Int, ::Int, ::Int) where T <: RingElement
add_row(::MatElem{T}, ::Int, ::Int, ::Int) where T <: RingElement
add_row!(::MatElem{T}, ::Int, ::Int, ::Int) where T <: RingElement
multiply_column(::MatElem{T}, ::Int, ::Int) where T <: RingElement
multiply_column!(::MatElem{T}, ::Int, ::Int) where T <: RingElement
multiply_row(::MatElem{T}, ::Int, ::Int) where T <: RingElement
multiply_row!(::MatElem{T}, ::Int, ::Int) where T <: RingElement

```

**Examples**
```jldoctest
julia> M = ZZ[1 2 3; 2 3 4; 4 5 5]
[1   2   3]
[2   3   4]
[4   5   5]

julia> add_column(M, 2, 3, 1)
[ 7   2   3]
[10   3   4]
[14   5   5]

julia> add_row(M, 1, 2, 3)
[1   2   3]
[2   3   4]
[6   8   9]

julia> multiply_column(M, 2, 3)
[1   2    6]
[2   3    8]
[4   5   10]

julia> multiply_row(M, 2, 3)
[1    2    3]
[2    3    4]
[8   10   10]
```

### Swapping rows and columns

```@docs
swap_rows(a::MatrixElem{T}, i::Int, j::Int) where T <: RingElement
swap_rows!(a::MatrixElem{T}, i::Int, j::Int) where T <: RingElement
swap_cols(a::MatrixElem{T}, i::Int, j::Int) where T <: RingElement
swap_cols!(a::MatrixElem{T}, i::Int, j::Int) where T <: RingElement
```

Swap the rows of `M` in place. The function returns the mutated matrix (since
matrices are assumed to be mutable in AbstractAlgebra.jl).

### Concatenation

The following are only available for matrix spaces, not for matrix algebras.

```julia
hcat(M::T, N::T) where T <: MatElem
```

Return the horizontal concatenation of $M$ and $N$. It is assumed that the number of
rows of $M$ and $N$ are the same.

```julia
vcat(M::T, N::T) where T <: MatElem
```

Return the vertical concatenation of $M$ and $N$. It is assumed that the number of
columns of $M$ and $N$ are the same.

**Examples**

```jldoctest
julia> M = matrix(ZZ, BigInt[1 2 3; 2 3 4; 3 4 5])
[1   2   3]
[2   3   4]
[3   4   5]

julia> N = matrix(ZZ, BigInt[1 0 1; 0 1 0; 1 0 1])
[1   0   1]
[0   1   0]
[1   0   1]

julia> P = hcat(M, N)
[1   2   3   1   0   1]
[2   3   4   0   1   0]
[3   4   5   1   0   1]

julia> Q = vcat(M, N)
[1   2   3]
[2   3   4]
[3   4   5]
[1   0   1]
[0   1   0]
[1   0   1]

```

### Similar and zero

Both `similar` and `zero` construct new matrices, but
the entries are either undefined with `similar` or zero-initialized with `zero`.

```julia
similar(x::MatElem, R::Ring=base_ring(x))
zero(x::MatElem, R::Ring=base_ring(x))
```

Construct the matrix with the same dimensions as the given matrix, and the
same base ring unless explicitly specified.

```julia
similar(x::MatElem, R::Ring, r::Int, c::Int)
similar(x::MatElem, r::Int, c::Int)
zero(x::MatElem, R::Ring, r::Int, c::Int)
zero(x::MatElem, r::Int, c::Int)
```

Construct the $r\times c$ matrix with `R` as base ring (which defaults to the
base ring of the the given matrix).
If $x$ belongs to a matrix algebra and $r \neq c$, an exception is raised, and it's
also possible to specify only one `Int` as the order (e.g. `similar(x, n)`).

```julia
Base.isassigned(M::MatElem, i, j)
```

Test whether the given matrix has a value associated with indices `i` and `j`.

**Examples**

```jldoctest
julia> M = matrix(ZZ, BigInt[3 1 2; 2 0 1])
[3   1   2]
[2   0   1]

julia> isassigned(M, 1, 2)
true

julia> isassigned(M, 4, 4)
false

julia> A = similar(M)
[#undef   #undef   #undef]
[#undef   #undef   #undef]

julia> isassigned(A, 1, 2)
false

julia> B = zero(M)
[0   0   0]
[0   0   0]

julia> C = similar(M, 4, 5)
[#undef   #undef   #undef   #undef   #undef]
[#undef   #undef   #undef   #undef   #undef]
[#undef   #undef   #undef   #undef   #undef]
[#undef   #undef   #undef   #undef   #undef]

julia> base_ring(B)
Integers

julia> D = zero(M, QQ, 2, 2)
[0//1   0//1]
[0//1   0//1]

julia> base_ring(D)
Rationals
```

### Symmetry testing

```@docs
is_symmetric(::MatrixElem)
```

```@docs
is_skew_symmetric(::MatrixElem)
```

### Powering

```@docs
powers(::MatElem, ::Int)
```

### Gram matrix

```@docs
gram(::MatElem)
```

### Trace

```@docs
tr(::MatElem{T}) where T <: RingElement
```

### Content

```@docs
content(::MatElem{T}) where T <: RingElement
```

### Permutation

```@docs
*(::Perm, ::MatElem{T}) where T <: RingElement
```

### LU factorisation

```@docs
lu{T <: FieldElem}(::MatElem{T}, ::SymmetricGroup)
```

```@docs
fflu{T <: RingElem}(::MatElem{T}, ::SymmetricGroup)
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(QQ, "x")
(Univariate polynomial ring in x over rationals, x)

julia> K, = residue_field(R, x^3 + 3x + 1); a = K(x);

julia> S = matrix_space(K, 3, 3)
Matrix space of 3 rows and 3 columns
  over residue field of univariate polynomial ring modulo x^3 + 3*x + 1

julia> A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 - 2 a - 1 2a])
[      0   2*x + 3   x^2 + 1]
[x^2 - 2     x - 1       2*x]
[x^2 - 2     x - 1       2*x]

julia> r, P, L, U = lu(A)
(2, (1,2), [1 0 0; 0 1 0; 1 0 1], [x^2-2 x-1 2*x; 0 2*x+3 x^2+1; 0 0 0])

julia> r, d, P, L, U = fflu(A)
(2, 3*x^2 - 10*x - 8, (1,2), [x^2-2 0 0; 0 3*x^2-10*x-8 0; x^2-2 0 1], [x^2-2 x-1 2*x; 0 3*x^2-10*x-8 -4*x^2-x-2; 0 0 0])

```

### Reduced row-echelon form

```@docs
rref_rational{T <: RingElem}(::MatElem{T})
rref{T <: FieldElem}(::MatElem{T})
```

```@docs
is_rref{T <: RingElem}(::MatElem{T})
is_rref{T <: FieldElem}(::MatElem{T})
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(QQ, "x")
(Univariate polynomial ring in x over rationals, x)

julia> K, = residue_field(R, x^3 + 3x + 1); a = K(x);

julia> S = matrix_space(K, 3, 3)
Matrix space of 3 rows and 3 columns
  over residue field of univariate polynomial ring modulo x^3 + 3*x + 1

julia> M = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])
[            0   2*x + 3   x^2 + 1]
[      x^2 - 2     x - 1       2*x]
[x^2 + 3*x + 1       2*x         1]

julia> r, A = rref(M)
(3, [1 0 0; 0 1 0; 0 0 1])

julia> is_rref(A)
true

julia> R, x = polynomial_ring(ZZ, "x")
(Univariate polynomial ring in x over integers, x)

julia> S = matrix_space(R, 3, 3)
Matrix space of 3 rows and 3 columns
  over univariate polynomial ring in x over integers

julia> M = S([R(0) 2x + 3 x^2 + 1; x^2 - 2 x - 1 2x; x^2 + 3x + 1 2x R(1)])
[            0   2*x + 3   x^2 + 1]
[      x^2 - 2     x - 1       2*x]
[x^2 + 3*x + 1       2*x         1]

julia> r, A, d = rref_rational(M)
(3, [-x^5-2*x^4-15*x^3-18*x^2-8*x-7 0 0; 0 -x^5-2*x^4-15*x^3-18*x^2-8*x-7 0; 0 0 -x^5-2*x^4-15*x^3-18*x^2-8*x-7], -x^5 - 2*x^4 - 15*x^3 - 18*x^2 - 8*x - 7)

julia> is_rref(A)
true
```

### Determinant

```@docs
det{T <: RingElem}(::MatElem{T})
```

### Rank

```@docs
rank{T <: RingElem}(::MatElem{T})
```

### Minors

```@docs
minors(::MatElem, ::Int)
```

### Exterior power

```@docs
exterior_power(::MatElem, ::Int)
```

### Pfaffian

```@docs
pfaffian(::MatElem)
pfaffians(::MatElem, ::Int)
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(QQ, ["x$i" for i in 1:6])
(Multivariate polynomial ring in 6 variables over rationals, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}[x1, x2, x3, x4, x5, x6])

julia> M = R[0 x[1] x[2] x[3]; -x[1] 0 x[4] x[5]; -x[2] -x[4] 0 x[6]; -x[3] -x[5] -x[6] 0]
[  0    x1    x2   x3]
[-x1     0    x4   x5]
[-x2   -x4     0   x6]
[-x3   -x5   -x6    0]

julia> pfaffian(M)
x1*x6 - x2*x5 + x3*x4

julia> pfaffians(M, 2)
6-element Vector{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}:
 x1
 x2
 x4
 x3
 x5
 x6
 
```

### Linear solving

```@docs
solve{T <: FieldElem}(::MatElem{T}, ::MatElem{T})
```

```@docs
solve_rational{T <: RingElem}(::MatElem{T}, ::MatElem{T})
```

```@docs
can_solve_with_solution{T <: RingElement}(::MatElem{T}, ::MatElem{T})
```

```@docs
can_solve{T <: RingElement}(::MatElem{T}, ::MatElem{T})
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
julia> R, x = polynomial_ring(QQ, "x")
(Univariate polynomial ring in x over rationals, x)

julia> K, = residue_field(R, x^3 + 3x + 1); a = K(x);

julia> S = matrix_space(K, 3, 3)
Matrix space of 3 rows and 3 columns
  over residue field of univariate polynomial ring modulo x^3 + 3*x + 1

julia> U = matrix_space(K, 3, 1)
Matrix space of 3 rows and 1 column
  over residue field of univariate polynomial ring modulo x^3 + 3*x + 1

julia> A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])
[            0   2*x + 3   x^2 + 1]
[      x^2 - 2     x - 1       2*x]
[x^2 + 3*x + 1       2*x         1]

julia> b = U([2a, a + 1, (-a - 1)])
[   2*x]
[ x + 1]
[-x - 1]

julia> x = solve(A, b)
[  1984//7817*x^2 + 1573//7817*x - 937//7817]
[ -2085//7817*x^2 + 1692//7817*x + 965//7817]
[-3198//7817*x^2 + 3540//7817*x - 3525//7817]

julia> A = matrix(ZZ, 2, 2, [1, 2, 0, 2])
[1   2]
[0   2]

julia> b = matrix(ZZ, 2, 1, [2, 1])
[2]
[1]

julia> can_solve(A, b, side = :right)
false

julia> A = matrix(QQ, 2, 2, [3, 4, 5, 6])
[3//1   4//1]
[5//1   6//1]

julia> b = matrix(QQ, 1, 2, [2, 1])
[2//1   1//1]

julia> can_solve_with_solution(A, b; side = :left)
(true, [-7//2 5//2])

julia> A = S([a + 1 2a + 3 a^2 + 1; K(0) a^2 - 1 2a; K(0) K(0) a])
[x + 1   2*x + 3   x^2 + 1]
[    0   x^2 - 1       2*x]
[    0         0         x]

julia> bb = U([2a, a + 1, (-a - 1)])
[   2*x]
[ x + 1]
[-x - 1]

julia> x = solve_triu(A, bb, false)
[ 1//3*x^2 + 8//3*x + 13//3]
[-3//5*x^2 - 3//5*x - 12//5]
[                   x^2 + 2]

julia> R, x = polynomial_ring(ZZ, "x")
(Univariate polynomial ring in x over integers, x)

julia> S = matrix_space(R, 3, 3)
Matrix space of 3 rows and 3 columns
  over univariate polynomial ring in x over integers

julia> U = matrix_space(R, 3, 2)
Matrix space of 3 rows and 2 columns
  over univariate polynomial ring in x over integers

julia> A = S([R(0) 2x + 3 x^2 + 1; x^2 - 2 x - 1 2x; x^2 + 3x + 1 2x R(1)])
[            0   2*x + 3   x^2 + 1]
[      x^2 - 2     x - 1       2*x]
[x^2 + 3*x + 1       2*x         1]

julia> bbb = U(transpose([2x x + 1 (-x - 1); x + 1 (-x) x^2]))
[   2*x   x + 1]
[ x + 1      -x]
[-x - 1     x^2]

julia> x, d = solve_rational(A, bbb)
([3*x^4-10*x^3-8*x^2-11*x-4 -x^5+3*x^4+x^3-2*x^2+3*x-1; -2*x^5-x^4+6*x^3+2*x+1 x^6+x^5+4*x^4+9*x^3+8*x^2+5*x+2; 6*x^4+12*x^3+15*x^2+6*x-3 -2*x^5-4*x^4-6*x^3-9*x^2-4*x+1], x^5 + 2*x^4 + 15*x^3 + 18*x^2 + 8*x + 7)

julia> S = matrix_space(ZZ, 3, 3)
Matrix space of 3 rows and 3 columns
  over integers

julia> T = matrix_space(ZZ, 3, 1)
Matrix space of 3 rows and 1 column
  over integers

julia> A = S([BigInt(2) 3 5; 1 4 7; 9 2 2])
[2   3   5]
[1   4   7]
[9   2   2]

julia> B = T([BigInt(4), 5, 7])
[4]
[5]
[7]
```

### Inverse

```@docs
Base.inv{T <: RingElement}(::MatrixElem{T})
```

```@docs
is_invertible_with_inverse{T <: RingElement}(::MatrixElem{T})
```

```@docs
is_invertible{T <: RingElement}(::MatrixElem{T})
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(QQ, "x")
(Univariate polynomial ring in x over rationals, x)

julia> K, = residue_field(R, x^3 + 3x + 1); a = K(x);

julia> S = matrix_space(K, 3, 3)
Matrix space of 3 rows and 3 columns
  over residue field of univariate polynomial ring modulo x^3 + 3*x + 1

julia> A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])
[            0   2*x + 3   x^2 + 1]
[      x^2 - 2     x - 1       2*x]
[x^2 + 3*x + 1       2*x         1]

julia> X = inv(A)
[-343//7817*x^2 + 717//7817*x - 2072//7817   -4964//23451*x^2 + 2195//23451*x - 11162//23451    -232//23451*x^2 - 4187//23451*x - 1561//23451]
[ 128//7817*x^2 - 655//7817*x + 2209//7817      599//23451*x^2 - 2027//23451*x - 1327//23451   -1805//23451*x^2 + 2702//23451*x - 7394//23451]
[ 545//7817*x^2 + 570//7817*x + 2016//7817     -1297//23451*x^2 - 5516//23451*x - 337//23451   8254//23451*x^2 - 2053//23451*x + 16519//23451]

julia> is_invertible(A)
true

julia> is_invertible_with_inverse(A)
(true, [-343//7817*x^2+717//7817*x-2072//7817 -4964//23451*x^2+2195//23451*x-11162//23451 -232//23451*x^2-4187//23451*x-1561//23451; 128//7817*x^2-655//7817*x+2209//7817 599//23451*x^2-2027//23451*x-1327//23451 -1805//23451*x^2+2702//23451*x-7394//23451; 545//7817*x^2+570//7817*x+2016//7817 -1297//23451*x^2-5516//23451*x-337//23451 8254//23451*x^2-2053//23451*x+16519//23451])

julia> R, x = polynomial_ring(ZZ, "x")
(Univariate polynomial ring in x over integers, x)

julia> S = matrix_space(R, 3, 3)
Matrix space of 3 rows and 3 columns
  over univariate polynomial ring in x over integers

julia> A = S([R(0) 2x + 3 x^2 + 1; x^2 - 2 x - 1 2x; x^2 + 3x + 1 2x R(1)])
[            0   2*x + 3   x^2 + 1]
[      x^2 - 2     x - 1       2*x]
[x^2 + 3*x + 1       2*x         1]

julia> X, d = pseudo_inv(A)
([4*x^2-x+1 -2*x^3+3 x^3-5*x^2-5*x-1; -2*x^3-5*x^2-2*x-2 x^4+3*x^3+2*x^2+3*x+1 -x^4+x^2+2; -x^3+2*x^2+2*x-1 -2*x^3-9*x^2-11*x-3 2*x^3+3*x^2-4*x-6], -x^5 - 2*x^4 - 15*x^3 - 18*x^2 - 8*x - 7)

```

### Nullspace

```@docs
nullspace{T <: FieldElem}(::MatElem{T})
```

### Kernel

```@docs
kernel{T <: RingElem}(::MatElem{T})
left_kernel{T <: RingElem}(::MatElem{T})
right_kernel{T <: RingElem}(::MatElem{T})
```

**Examples**

```jldoctest
julia> S = matrix_space(ZZ, 4, 4)
Matrix space of 4 rows and 4 columns
  over integers

julia> M = S([1 2 0 4;
              2 0 1 1;
              0 1 1 -1;
              2 -1 0 2])
[1    2   0    4]
[2    0   1    1]
[0    1   1   -1]
[2   -1   0    2]

julia> nr, Nr = kernel(M)
(1, [-8; -6; 11; 5])

julia> nl, Nl = left_kernel(M)
(1, [0 -1 1 1])

```

### Hessenberg form

```@docs
hessenberg{T <: RingElem}(::MatElem{T})
```

```@docs
is_hessenberg{T <: RingElem}(::MatElem{T})
```

**Examples**

```jldoctest
julia> R, = residue_ring(ZZ, 7);

julia> S = matrix_space(R, 4, 4)
Matrix space of 4 rows and 4 columns
  over residue ring of integers modulo 7

julia> M = S([R(1) R(2) R(4) R(3); R(2) R(5) R(1) R(0);
              R(6) R(1) R(3) R(2); R(1) R(1) R(3) R(5)])
[1   2   4   3]
[2   5   1   0]
[6   1   3   2]
[1   1   3   5]

julia> A = hessenberg(M)
[1   5   5   3]
[2   1   1   0]
[0   1   3   2]
[0   0   2   2]

julia> is_hessenberg(A)
true

```

### Characteristic polynomial

```@docs
charpoly{T <: RingElem}(::PolyRing{T}, ::MatrixElem{T})
```

### Minimal polynomial

```@docs
minpoly{T <: RingElem}(::PolyRing{T}, ::MatElem{T}, ::Bool)
```

### Transforms

```@docs
similarity!{T <: RingElem}(::MatElem{T}, ::Int, ::T)
```

### Hermite normal form

```@docs
hnf{T <: RingElem}(::MatElem{T})
hnf_with_transform{T <: RingElem}(::MatElem{T})
```

```@docs
is_hnf{T <: RingElem}(::MatElem{T})
```

**Examples**

```jldoctest
julia> A = matrix(ZZ, [2 3 -1; 3 5 7; 11 1 12])
[ 2   3   -1]
[ 3   5    7]
[11   1   12]

julia> H = hnf(A)
[1   0   255]
[0   1    17]
[0   0   281]

julia> is_hnf(H)
true

julia> H, U = hnf_with_transform(A)
([1 0 255; 0 1 17; 0 0 281], [-47 28 1; -3 2 0; -52 31 1])

julia> U*A
[1   0   255]
[0   1    17]
[0   0   281]
```

### Smith normal form

```@docs
is_snf(::MatrixElem{T}) where T <: RingElement
```

```@docs
snf{T <: RingElem}(::MatElem{T})
snf_with_transform{T <: RingElem}(::MatElem{T})
```

**Examples**

```jldoctest
julia> A = matrix(ZZ, [2 3 -1; 3 5 7; 11 1 12])
[ 2   3   -1]
[ 3   5    7]
[11   1   12]

julia> S = snf(A)
[1   0     0]
[0   1     0]
[0   0   281]

julia> S, T, U = snf_with_transform(A)
([1 0 0; 0 1 0; 0 0 281], [1 0 0; 7 1 0; 229 31 1], [0 -3 26; 0 2 -17; -1 0 1])

julia> T*A*U
[1   0     0]
[0   1     0]
[0   0   281]
```

### (Weak) Popov form

AbstractAlgebra.jl provides algorithms for computing the (weak) Popov of a matrix with
entries in a univariate polynomial ring over a field.

```@docs
is_weak_popov(P::MatrixElem{T}, rank::Int) where T <: Generic.Poly
```

```@docs
weak_popov{T <: PolyRingElem}(::MatElem{T})
weak_popov_with_transform{T <: PolyRingElem}(::MatElem{T})
popov{T <: PolyRingElem}(::MatElem{T})
popov_with_transform{T <: PolyRingElem}(::MatElem{T})
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(QQ, "x");

julia> A = matrix(R, map(R, Any[1 2 3 x; x 2*x 3*x x^2; x x^2+1 x^3+x^2 x^4+x^2+1]))
[1         2           3               x]
[x       2*x         3*x             x^2]
[x   x^2 + 1   x^3 + x^2   x^4 + x^2 + 1]

julia> P = weak_popov(A)
[   1                        2                    3   x]
[   0                        0                    0   0]
[-x^3   -2*x^3 + x^2 - 2*x + 1   -2*x^3 + x^2 - 3*x   1]

julia> P, U = weak_popov_with_transform(A)
([1 2 3 x; 0 0 0 0; -x^3 -2*x^3+x^2-2*x+1 -2*x^3+x^2-3*x 1], [1 0 0; -x 1 0; -x^3-x 0 1])

julia> U*A
[   1                        2                    3   x]
[   0                        0                    0   0]
[-x^3   -2*x^3 + x^2 - 2*x + 1   -2*x^3 + x^2 - 3*x   1]
```
