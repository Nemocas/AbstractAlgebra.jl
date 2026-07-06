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

```@docs
scalar_matrix(R::Ring, n::Int, a::RingElement)
scalar_matrix(n::Int, a::NCRingElement)
```


### Diagonal matrices

```@docs
diagonal_matrix(x::NCRingElement, m::Int, n::Int)
diagonal_matrix(x::NCRingElement, m::Int)
diagonal_matrix(R::NCRing, x::AbstractVector{<:NCRingElement})
```

Diagonal matrices can also be constructed from matrix blocks. This is
an alias for `block_diagonal_matrix`.

```@docs
diagonal_matrix(V::Vector{T}) where {T <: MatElem}
```


## Block diagonal matrix constructors

```@docs
block_diagonal_matrix(V::Vector{<:MatElem{T}}) where {T <: NCRingElement}
block_diagonal_matrix(R::NCRing, V::Vector{<:Matrix{T}}) where {T <: NCRingElement}
```


## Triangular matrices

```@docs
lower_triangular_matrix(L::AbstractVector{T}) where {T <: NCRingElement}
upper_triangular_matrix(L::AbstractVector{T}) where {T <: NCRingElement}
strictly_lower_triangular_matrix(L::AbstractVector{T}) where {T <: NCRingElement}
strictly_upper_triangular_matrix(L::AbstractVector{T}) where {T <: NCRingElement}
```
