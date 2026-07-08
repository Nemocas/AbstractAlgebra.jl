```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# [Matrix functionality](@id matrix_functionality_chapter)

## Basic matrix functionality

```@docs
zero(::MatElem{T}, ::Ring) where T <: RingElement
one(::MatElem{T}) where T <: RingElement
transpose(::MatElem)
transpose!(::MatElem)
change_base_ring(::Ring, ::MatElem{T}) where T <: RingElement
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

### Permutation

```@docs
*(::Perm, ::MatElem{T}) where T <: RingElement
```

### Transforms

```@docs
similarity!(::MatElem{T}, ::Int, ::T) where T <: RingElem
```
