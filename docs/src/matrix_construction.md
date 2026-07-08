```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# [Constructing Algebraic Matrices](@id matrix_construction)

There are several functions for creating algebraic matrices -- in each case you need to indicate the
ring to which the matrix elements belong, unless this is already unambiguously indicated
by the ring to which the arguments belong. Here we present several of the basic functions for
constructing matrices.


## Why Julia matrices are not sufficient

Unfortunately, Julia's matrices cannot be used in our context due to two independent problems:

- In empty matrices (0 rows or columns) all that is known about a Julia matrix is the type of its entries,
however for advanced algebraic types, this information is insufficient to create elements,
hence `zero(T)` or friends cannot work.
- Many Julia functions (e.g. `det`) assume that all types used embed into the real or complex numbers. For
instance, in Julia `det(ones(Int, (1,1))) == 1.0`, so the fact that the determinant is exactly the integer `1`
is lost. Furthermore, more general rings cannot be embedded into the reals at all.

For this reason, the matrix constructors below must accept (or infer) the ring in which all entries of the
matrix reside.


## Constructing matrices from entries

```@docs
matrix(R::NCRing, arr::AbstractMatrix{T}) where {T}
matrix(R::NCRing, r::Int, c::Int, arr::AbstractVecOrMat{T}) where {T}
```

Several other signatures are supported:

```julia
matrix(arr::AbstractMatrix{T}) where {T<:NCRingElement}
matrix(arr::AbstractVector{T}) where {T<:NCRingElement}
matrix(arr::AbstractVector{<:AbstractVector{T}}) where {T<:NCRingElement}
```

Nested vectors can also be used to construct matrices over a specified
base ring:

```julia
matrix(R::NCRing, arr::AbstractVector{<:AbstractVector})
```


## Constructing matrices from existing matrices


### Changing the base ring

Algebraic matrices can be converted to another base ring using `change_base_ring`.

```@docs
change_base_ring(R::NCRing, M::MatrixElem{T}) where {T <: NCRingElement}
```

The same conversion can also be performed using the following constructors:

```julia
matrix(R::NCRing, arr::MatElem)
matrix(R::NCRing, arr::MatRingElem)
```


### Copying

An independent copy of an existing algebraic matrix can be created.

```julia
matrix(mat::MatElem{T}) where {T<:NCRingElement}
```


### Zero matrices

```@docs
zero(x::MatElem{T}, R::NCRing) where {T <: NCRingElement}
```


### Identity matrices

```@docs
identity_matrix(M::MatElem{T}) where {T <: NCRingElement}
identity_matrix(M::MatElem{T}, n::Int) where {T <: NCRingElement}
one(a::MatElem{T}) where {T <: NCRingElement}
```


### Uninitialized matrices

```@docs
similar(x::MatElem, R::NCRing, r::Int, c::Int)
```


### Concatenation

Matrices can be constructed from existing matrices by concatenating them
horizontally or vertically.

```@docs
Base.hcat(A::MatElem...)
Base.vcat(A::MatElem...)
```


### Submatrices and views

Submatrices can be constructed using Julia indexing syntax. By default, this
creates a new matrix containing the selected entries.

To avoid copying entries, one can instead construct a view using Julia's
`@view` syntax. Views are particularly useful for working with small regions of
large matrices, as modifications to a view also modify the original matrix.

**Examples**

```jldoctest
julia> M = matrix(ZZ, [1 2 3; 2 3 4; 3 4 5])
[1   2   3]
[2   3   4]
[3   4   5]

julia> M[1:2, 1:2]
[1   2]
[2   3]

julia> N = @view M[1:2, 1:2]
[1   2]
[2   3]

julia> N[1, 1] = 10
10

julia> M
[10   2   3]
[ 2   3   4]
[ 3   4   5]
```


## Special constructors


### The zero matrix

```@docs
zero_matrix(R::NCRing, r::Int, c::Int)
```


### The ones matrix

```@docs
ones_matrix(R::NCRing, r::Int, c::Int)
```


### The identity matrix

```@docs
identity_matrix(R::NCRing, n::Int)
```


### Scalar matrices

```@docs
scalar_matrix(R::NCRing, n::Int, a::NCRingElement)
scalar_matrix(n::Int, a::NCRingElement)
```


### Diagonal matrices

```@docs
diagonal_matrix(x::NCRingElement, m::Int, n::Int = m)
diagonal_matrix(R::NCRing, x::AbstractVector{<:NCRingElement})
```

Diagonal matrices can also be constructed from matrix blocks. This is
an alias for `block_diagonal_matrix`.

```@docs
diagonal_matrix(V::Vector{T}) where {T <: MatElem}
```


### Block diagonal matrix constructors

```@docs
block_diagonal_matrix(V::Vector{<:MatElem{T}}) where {T <: NCRingElement}
block_diagonal_matrix(R::NCRing, V::Vector{<:Matrix{T}}) where {T <: NCRingElement}
```


### Triangular matrices

```@docs
lower_triangular_matrix(L::AbstractVector{T}) where {T <: NCRingElement}
upper_triangular_matrix(L::AbstractVector{T}) where {T <: NCRingElement}
strictly_lower_triangular_matrix(L::AbstractVector{T}) where {T <: NCRingElement}
strictly_upper_triangular_matrix(L::AbstractVector{T}) where {T <: NCRingElement}
```
