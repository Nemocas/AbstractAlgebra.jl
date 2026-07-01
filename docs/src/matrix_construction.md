```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# [Constructing Matrices](@id matrix_construction)


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

Given an $r \times c$ Julia matrix, the following command constructs the corresponding
matrix over the given ring `R`, assuming all the entries can be coerced into `R`.

```julia
matrix(R::Ring, arr::Matrix{T}) where T <: RingElement
```

Likewise, we construct the $r \times c$ matrix over the ring `R` whose $(i, j)$ entry is
given by `A[c*(i - 1) + j]`, assuming that all the entries can be coerced into `R`, as
follows:

```julia
matrix(R::Ring, r::Int, c::Int, A::Vector{T}) where T <: RingElement
```

**Examples**

```jldoctest
julia> M = matrix(ZZ, BigInt[3 1 2; 2 0 1])
[3   1   2]
[2   0   1]

julia> N = matrix(ZZ, 3, 2, BigInt[3, 1, 2, 2, 0, 1])
[3   1]
[2   2]
[0   1]
```


## Special constructors

### The zero matrix

Construct the $r \times c$ matrix over the ring `R` whose entries are all zero in `R`.

```julia
zero_matrix(R::Ring, r::Int, c::Int)
```

**Examples**

```jldoctest
julia> P = zero_matrix(ZZ, 3, 2)
[0   0]
[0   0]
[0   0]
```


### The identity matrix

Create the $n \times n$ identity matrix over the ring `R`.

```@docs
identity_matrix(::Ring, ::Int)
```

**Examples**

```jldoctest
julia> identity_matrix(ZZ, 2)
[1   0]
[0   1]
```


### The ones matrix

Create the $r \times c$ matrix over the ring `R` whose entries are
all equal to the multiplicative identity of `R`.

```@docs
ones_matrix(::Ring, ::Int, ::Int)
```

**Examples**

```jldoctest
julia> I = ones_matrix(ZZ, 3, 2)
[1   1]
[1   1]
[1   1]
```


### Scalar matrices

Create the $n \times n$ matrix over the ring `R`, whose entries are all
identical to a given scalar `a` in `R`.

```@docs
scalar_matrix(R::Ring, n::Int, a::RingElement)
```

**Examples**

```jldoctest
julia> scalar_matrix(QQ, 3, 1//2)
[1//2      0      0]
[   0   1//2      0]
[   0      0   1//2]
```



### Diagonal matrices

Create a matrix over the ring `R` whose entries are all zero except for
the diagonal entries. Two constructors are provided. One creates a rectangular
diagonal matrix whose diagonal entries are all equal to a given scalar. The
other constructs a diagonal matrix whose diagonal entries are taken from a given
vector.

```@docs
diagonal_matrix(::RingElement, ::Int, ::Int)
diagonal_matrix(::NCRing, ::AbstractVector{<:NCRingElement})
```

**Examples**

```jldoctest
julia> diagonal_matrix(1, 3, 4)
[1   0   0   0]
[0   1   0   0]
[0   0   1   0]

julia> diagonal_matrix(ZZ, [1,2,3])
[1   0   0]
[0   2   0]
[0   0   3]
```



## Block diagonal matrix constructors

Create block diagonal matrices from a vector of existing matrices.
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

Create (strictly) trianngular matrices over the ring `R`.

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
[0   0]
```
