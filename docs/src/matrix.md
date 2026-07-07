```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# [Matrix functionality](@id matrix_functionality_chapter)

## Basic matrix functionality

```@docs
number_of_rows(::MatElem)
number_of_columns(::MatElem)
length(::MatrixElem{T}) where T <: RingElement
zero(::MatElem{T}, ::Ring) where T <: RingElement
one(::MatElem{T}) where T <: RingElement
transpose(::MatElem)
transpose!(::MatElem)
tr(::MatElem{T}) where T <: RingElement
det(::MatElem{T}) where T <: RingElem
rank(::MatElem{T}) where T <: RingElem
change_base_ring(::Ring, ::MatElem{T}) where T <: RingElement
```

## Inverse

```@docs; canonical=false
Base.inv(::MatElem{T}) where T <: RingElement
pseudo_inv(M::MatElem{T}) where T <: RingElement
```

**Examples**

```@jldoctest
julia> M = matrix(QQ, 3, 3, [1 2 3;4 5 6;0 0 1])
[1//1   2//1   3//1]
[4//1   5//1   6//1]
[0//1   0//1   1//1]

julia> X = inv(M)
[-5//3    2//3    1//1]
[ 4//3   -1//3   -2//1]
[ 0//1    0//1    1//1]

julia> is_invertible(M)
true

julia> is_invertible_with_inverse(M)
(true, [-5//3 2//3 1; 4//3 -1//3 -2; 0 0 1])

julia> pseudo_inv(M)
([5 -2 -3; -4 1 6; 0 0 -3], -3//1)
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

## Other functionality

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
nullspace(::MatElem{T}) where T <: FieldElem
```

### Characteristic polynomial

```@docs
charpoly(::PolyRing{T}, ::MatElem{T}) where T <: RingElem
```

### Minimal polynomial

```@docs
minpoly(::PolyRing{T}, ::MatElem{T}, ::Bool) where T <: RingElem
```

### Transforms

```@docs
similarity!(::MatElem{T}, ::Int, ::T) where T <: RingElem
```
