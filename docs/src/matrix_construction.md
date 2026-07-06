```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# [Constructing Algebraic Matrices](@id matrix_construction)

There are several functions for creating algebraic matrices -- in each case you need to indicate the
ring to which the matrix elements belong, unless this is already unambiguously indicated
by the ring to which the arguments belong.  Here we present several of the basic functions for
constructing matrices.


## Why Julia matrices are not sufficient

Unfortunately, Julia's matrices cannot be used in our context due to two independent problems:

- In empty matrices (0 rows or columns) all that is known about a Julia matrix is the type of its entries,
however for the complex algebraic types, this information is not sufficient to create elements,
hence `zero(T)` or friends cannot work.
- Many Julia functions (e.g. `det`) assume that all types used embed into the real or complex numbers. For
instance, in Julia `det(ones(Int, (1,1))) == 1.0`, so the fact that the determinant is exactly the integer `1`
is lost. Furthermore, more general rings cannot be embedded into the reals at all.

For this reason, the matrix constructors below must accept (or infer) the ring in which all entries of the
matrix reside.


## Generic constructors

```@docs
matrix(R::NCRing, arr::AbstractMatrix{T}) where {T}
matrix(R::NCRing, r::Int, c::Int, arr::AbstractVecOrMat{T}) where {T}
```

Several other signatures are supported. Existing matrices can be converted
to a different base ring:

```julia
matrix(R::NCRing, arr::MatElem)
matrix(R::NCRing, arr::MatRingElem)
```

We can also create a deep copy of an existing algebraic matrix:

```julia
matrix(mat::MatElem{T}) where {T<:NCRingElement}
```

The base ring can also be inferred from the entries:

```julia
matrix(arr::AbstractMatrix{T}) where {T<:NCRingElement}
matrix(arr::AbstractVector{T}) where {T<:NCRingElement}
matrix(arr::AbstractVector{<:AbstractVector{T}}) where {T<:NCRingElement}
```

Finally, nested vectors can be used to construct matrices over a specified
base ring:

```julia
matrix(R::NCRing, arr::AbstractVector{<:AbstractVector})
```


## Special constructors

### The zero matrix

```@docs
zero_matrix(R::NCRing, r::Int, c::Int)
```

The following convenience method is also supported:

```julia
zero_matrix(::Type{MatElem}, R::Ring, n::Int, m::Int)
```


### The ones matrix

```@docs
ones_matrix(R::NCRing, r::Int, c::Int)
```


### The identity matrix

```@docs
identity_matrix(R::NCRing, n::Int)
identity_matrix(M::MatElem{T}) where T <: NCRingElement
```

The following additional signatures are also supported:

```julia
identity_matrix(M::MatElem{T}, n::Int) where T <: NCRingElement
identity_matrix(::Type{MatElem}, R::Ring, n::Int)
```

### Scalar matrices

Create the $n \times n$ diagonal matrix over the ring `R`, whose diagonal entries are all
equal to the image of `a` in `R`.

```@docs
scalar_matrix(R::Ring, n::Int, a::RingElement)
```

**Examples**

```jldoctest
julia> scalar_matrix(QQ, 3, 1//2)
[1//2   0//1   0//1]
[0//1   1//2   0//1]
[0//1   0//1   1//2]
```



### Diagonal matrices

Create a matrix over the ring `R` whose entries are all zero except for
the diagonal entries. Two constructors are provided. One function creates a rectangular
diagonal matrix whose diagonal entries are all equal to a given scalar. The
other constructs a square diagonal matrix whose diagonal entries are taken from a given
vector.

```@docs
diagonal_matrix(::RingElement, ::Int, ::Int)
diagonal_matrix(::NCRing, ::AbstractVector{<:NCRingElement})
```

**Examples**

```jldoctest
julia> diagonal_matrix(ZZ(1), 3, 4)
[1   0   0   0]
[0   1   0   0]
[0   0   1   0]

julia> diagonal_matrix(ZZ, [1,2,3])
[1   0   0]
[0   2   0]
[0   0   3]
```



## Block diagonal matrix constructors

Create a block diagonal matrix from a vector of existing matrices.
If Julia matrices are provided, one also has to supply the base ring.

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


## Triangular matrices

Create (strictly) triangular matrices over the ring `R`.

Why does this not require the ring as an argument?

```@docs
lower_triangular_matrix(L::AbstractVector{T}) where {T <: RingElement}
upper_triangular_matrix(L::AbstractVector{T}) where {T <: RingElement}
strictly_lower_triangular_matrix(L::AbstractVector{T}) where {T <: RingElement}
strictly_upper_triangular_matrix(L::AbstractVector{T}) where {T <: RingElement}
```

**Examples**

```jldoctest
julia> lower_triangular_matrix([1, 2, 3])
[1   0]
[2   3]

julia> strictly_lower_triangular_matrix([5])
[0   0]
[5   0]

julia> upper_triangular_matrix([3, 2, 1])
[3   2]
[0   1]

julia> strictly_upper_triangular_matrix([-2])
[0   -2]
[0    0]
```
