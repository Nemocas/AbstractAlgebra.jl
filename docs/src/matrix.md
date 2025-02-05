```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# Matrix functionality

## Basic matrix functionality

As well as the Ring and Matrix interfaces, the following functions are provided to
manipulate matrices and to set and retrieve entries and other basic data associated
with the matrices.


It is possible to create matrices directly, without first
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

julia> R = matrix_ring(ZZ, 2)
Matrix ring of degree 2
  over integers

julia> M = R()
[0   0]
[0   0]
```

```@docs
number_of_rows(::MatrixElem{T}) where T <: RingElement
number_of_columns(::MatrixElem{T}) where T <: RingElement
length(::MatrixElem{T}) where T <: RingElement
isempty(::MatrixElem{T}) where T <: RingElement
identity_matrix(::Ring, ::Int)
identity_matrix(::MatElem{T}) where T <: RingElement
ones_matrix(::Ring, ::Int, ::Int)
scalar_matrix(R::Ring, n::Int, a::RingElement)
diagonal_matrix(::RingElement, ::Int, ::Int)
zero(::MatElem{T}, ::Ring) where T <: RingElement
one(::MatElem{T}) where T <: RingElement
transpose(::MatrixElem{T}) where T <: RingElement
tr(::MatElem{T}) where T <: RingElement
det{T <: RingElem}(::MatElem{T})
rank{T <: RingElem}(::MatElem{T})
lower_triangular_matrix(L::AbstractVector{T}) where {T <: RingElement}
upper_triangular_matrix(L::AbstractVector{T}) where {T <: RingElement}
strictly_lower_triangular_matrix(L::AbstractVector{T}) where {T <: RingElement}
strictly_upper_triangular_matrix(L::AbstractVector{T}) where {T <: RingElement}
is_lower_triangular(::MatrixElem)
is_upper_triangular(::MatrixElem)
is_diagonal(::MatrixElem)
change_base_ring(::Ring, ::MatElem{T}) where T <: RingElement
Base.map(f, ::MatrixElem{T}) where T <: RingElement
Base.map!(f, ::MatrixElem{S}, ::MatrixElem{T}) where {S <: RingElement, T <: RingElement}
```

## Inverse

```@docs
Base.inv{T <: RingElement}(::MatrixElem{T})

is_invertible_with_inverse{T <: RingElement}(::MatrixElem{T})

is_invertible{T <: RingElement}(::MatrixElem{T})
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(QQ, :x)
(Univariate polynomial ring in x over rationals, x)

julia> K, = residue_field(R, x^3 + 3x + 1); a = K(x);

julia> S = matrix_space(K, 3, 3)
Matrix space of 3 rows and 3 columns
  over residue field of R modulo x^3 + 3*x + 1

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

julia> R, x = polynomial_ring(ZZ, :x)
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

## Submatrices

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

## Elementary row and column operations

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

## Swapping rows and columns

```@docs
swap_rows(a::MatrixElem{T}, i::Int, j::Int) where T <: RingElement
swap_rows!(a::MatrixElem{T}, i::Int, j::Int) where T <: RingElement
swap_cols(a::MatrixElem{T}, i::Int, j::Int) where T <: RingElement
swap_cols!(a::MatrixElem{T}, i::Int, j::Int) where T <: RingElement
```

Swap the rows of `M` in place. The function returns the mutated matrix (since
matrices are assumed to be mutable in AbstractAlgebra.jl).

## Concatenation

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

## Linear solving

See [Linear Solving & Kernel](@ref solving_chapter)

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

## Similar and zero

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

## LU factorisation

```@docs
lu{T <: FieldElem}(::MatElem{T}, ::SymmetricGroup)

fflu{T <: RingElem}(::MatElem{T}, ::SymmetricGroup)
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(QQ, :x)
(Univariate polynomial ring in x over rationals, x)

julia> K, = residue_field(R, x^3 + 3x + 1); a = K(x);

julia> S = matrix_space(K, 3, 3)
Matrix space of 3 rows and 3 columns
  over residue field of R modulo x^3 + 3*x + 1

julia> A = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 - 2 a - 1 2a])
[      0   2*x + 3   x^2 + 1]
[x^2 - 2     x - 1       2*x]
[x^2 - 2     x - 1       2*x]

julia> r, P, L, U = lu(A)
(2, (1,2), [1 0 0; 0 1 0; 1 0 1], [x^2-2 x-1 2*x; 0 2*x+3 x^2+1; 0 0 0])

julia> r, d, P, L, U = fflu(A)
(2, 3*x^2 - 10*x - 8, (1,2), [x^2-2 0 0; 0 3*x^2-10*x-8 0; x^2-2 0 1], [x^2-2 x-1 2*x; 0 3*x^2-10*x-8 -4*x^2-x-2; 0 0 0])

```

## Reduced row-echelon form

```@docs
rref_rational{T <: RingElem}(::MatElem{T})
rref{T <: FieldElem}(::MatElem{T})

is_rref{T <: RingElem}(::MatElem{T})
is_rref{T <: FieldElem}(::MatElem{T})
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(QQ, :x)
(Univariate polynomial ring in x over rationals, x)

julia> K, = residue_field(R, x^3 + 3x + 1); a = K(x);

julia> S = matrix_space(K, 3, 3)
Matrix space of 3 rows and 3 columns
  over residue field of R modulo x^3 + 3*x + 1

julia> M = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])
[            0   2*x + 3   x^2 + 1]
[      x^2 - 2     x - 1       2*x]
[x^2 + 3*x + 1       2*x         1]

julia> r, A = rref(M)
(3, [1 0 0; 0 1 0; 0 0 1])

julia> is_rref(A)
true

julia> R, x = polynomial_ring(ZZ, :x)
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

### Symmetry testing

```@docs
is_symmetric(::MatrixElem)

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

### Content

```@docs
content(::MatElem{T}) where T <: RingElement
```

### Permutation

```@docs
*(::Perm, ::MatElem{T}) where T <: RingElement
```

### Nilpotency

```@docs
is_nilpotent(::MatrixElem{T}) where {T <: RingElement}
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

### Nullspace

```@docs
nullspace{T <: FieldElem}(::MatElem{T})
```

### Hessenberg form

```@docs
hessenberg{T <: RingElem}(::MatElem{T})

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

weak_popov{T <: PolyRingElem}(::MatElem{T})
weak_popov_with_transform{T <: PolyRingElem}(::MatElem{T})
popov{T <: PolyRingElem}(::MatElem{T})
popov_with_transform{T <: PolyRingElem}(::MatElem{T})
```

**Examples**

```jldoctest
julia> R, x = polynomial_ring(QQ, :x);

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
