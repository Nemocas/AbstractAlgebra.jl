```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# Matrix interface

This page describes the interface that a matrix implementation must satisfy in
order to use the generic matrix functionality provided by AbstractAlgebra.jl.

It is intended for developers implementing new matrix types. For details on
the generic dense matrix implementation provided by AbstractAlgebra itself, see
the matrix implementation documentation.

AbstractAlgebra supports two kinds of matrix objects:

- matrices in matrix spaces,
- matrices in matrix algebras.

Matrices in a matrix space form an additive group but are generally not
themselves a ring. Matrix algebras consist of square matrices and implement the
noncommutative ring interface.

All matrix objects are assumed to be mutable. This is important for performance
and for compatibility with in-place algorithms.

!!! note

    AbstractAlgebra matrices are not Julia `Matrix` objects. They store
    additional algebraic information, such as the base ring, and support
    arithmetic over arbitrary rings.


## Types and parents

Matrix space implementations use the following abstract types:

- `MatSpace{T}`: parent type for matrix spaces,
- `MatElem{T}`: abstract type for elements of matrix spaces.

Matrix algebra implementations use:

- `MatRing{T}`: abstract parent type for matrix algebras,
- `MatRingElem{T}`: abstract type for elements of matrix algebras.

The parameter `T` is the type of the entries of the matrix.

Matrix spaces and matrix algebras should be uniquely determined by their
defining data, in particular their base ring and dimensions.


## Required functionality

Suppose that:

- `R` is the base ring,
- `S` is a matrix space or matrix algebra over `R`,
- `MyMat{T}` is the matrix element type,
- `MyMatParent{T}` the parent type of the matrix space or matrix algebra.

The entry type `T` must belong to the `RingElem` hierarchy.


### Required constructors

Implementations must provide constructors from Julia arrays of entries:

```julia
(P::MyMatParent{T})(A::Matrix{U}) where {T <: NCRingElem, U <: NCRingElem}
```

Implementations must provide constructors from Julia arrays of entries:

```julia
(P::MyMatParent{T})(A::Matrix{U}) where {T <: NCRingElem, U <: NCRingElem}
```

and constructors allowing coercion into the base ring:

```julia
(P::MyMatParent{T})(A::Matrix{U}) where {T <: NCRingElem, U <: NCRingElem}
```

as well as construction from flat vectors:

```julia
(P::MyMatParent{T})(A::Vector{U}) where {T <: NCRingElem, U <: NCRingElem}
```


### Dimensions and entry access

Matrix implementations must provide:

```julia
number_of_rows(M::MyMat)
number_of_columns(M::MyMat)
```

returning the dimensions of the matrix.

Entry access and mutation are provided through:

```julia
getindex(M::MyMat, i::Int, j::Int)
setindex!(M::MyMat, x, i::Int, j::Int)
```

Matrices are mutable, and `setindex!` modifies the given matrix in place.


### Transpose

Implementations must provide:

```julia
transpose(M::MyMat)
```

returning the transpose of the matrix.


### Views

Matrix types must support views of submatrices. A view behaves like a matrix,
but modifying entries of the view modifies the original matrix.

The following behaviour is required:

- `deepcopy` of a view returns a view into a copy of the original matrix;
- `similar` applied to a view returns an ordinary matrix, not another view.


## Dense matrix type selection

Implementations can specify the preferred dense matrix representation for a
given coefficient type by defining:

```julia
dense_matrix_type(::Type{T}) where T <: NCRingElem
```

The remaining methods dispatch to this one:

```julia
dense_matrix_type(::T) where T <: NCRingElem
dense_matrix_type(::Type{S}) where S <: NCRing
dense_matrix_type(::S) where S <: NCRing
```

For example, another package may specify that matrices over a particular ring
should use a specialised representation instead of the generic dense matrix
implementation.


## Optional functionality

The generic matrix code provides fallback implementations for many operations.
Custom matrix implementations may override these methods for improved
performance.


### Submatrices

```julia
Base.getindex(M::MyMat, rows::AbstractVector{Int}, cols::AbstractVector{Int})
```

Return the corresponding submatrix.


### Row operations

```julia
swap_rows!(M::MyMat, i::Int, j::Int)
```

Swap two rows of `M` in place.


### Concatenation

```julia
hcat(M::MyMat, N::MyMat)
vcat(M::MyMat, N::MyMat)
```

Return horizontal and vertical concatenations.


### Zero tests

```julia
is_zero_entry(M::MyMat, i::Int, j::Int)
is_zero_row(M::MyMat, i::Int)
is_zero_column(M::MyMat, j::Int)
```


### Similar and zero matrices

Implementations may specialise:

```julia
similar(M::MyMat, R::Ring, r::Int, c::Int)
zero(M::MyMat, R::Ring, r::Int, c::Int)
```

The function `similar` constructs a matrix of the requested type, base ring
and dimensions, but does not initialise its entries. It is unrelated to
similarity transformations of matrices.

The function `zero` constructs the corresponding zero matrix, i.e. a matrix
whose entries are all initialised to zero.


### Assigned entries

```julia
Base.isassigned(M::MyMat, i, j)
```

Tests whether an entry exists at the given position.


### Symmetry

```julia
is_symmetric(M::MyMat)
```

Return whether the matrix is symmetric.
