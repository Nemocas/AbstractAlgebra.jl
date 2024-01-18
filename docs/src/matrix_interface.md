```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Matrix Interface

Generic matrices are supported in AbstractAlgebra.jl. Both the space of $m\times n$
matrices and the algebra (ring) of $m\times m$ matrices are supported.

As the space of $m\times n$ matrices over a commutative ring is not itself a commutative
ring, not all of the Ring interface needs to be implemented for such matrices in.

In particular, the following functions do not need to be implemented: `is_domain_type`,
and `divexact`.
The `canonical_unit` function should be implemented, but simply needs to return the
corresponding value for entry $[1, 1]$ (the function is never called on empty matrices).

For matrix algebras, all of the ring interface must be implemented.

!!! note

    AbstractAlgebra.jl matrices are not the same as Julia matrices. We store a
    base ring in our matrix and matrices are row major instead of column major in order to
    support the numerous large C libraries that use this convention.

All AbstractAlgebra.jl matrices are assumed to be mutable. This is usually critical to
performance.

## Types and parents

AbstractAlgebra provides two abstract types for matrix spaces and their elements:

  * `MatSpace{T}` is the abstract type for matrix space parent types
  * `MatElem{T}` is the abstract type for matrix types belonging to a matrix space

It also provides two abstract types for matrix algebras and their elements:

  * `MatAlgebra{T}` is the abstract type for matrix algebra parent types
  * `MatAlgElem{T}` is the abstract type for matrix types belonging to a matrix algebra

Note that these abstract types are parameterised. The type `T` should usually be the
type of elements of the matrices.

Matrix spaces and matrix algebras should be made unique on the system by caching parent
objects (unless an optional `cache` parameter is set to `false`). Matrix spaces and
algebras should at least be distinguished based on their base (coefficient) ring and the
dimensions of the matrices in the space.

See `src/generic/GenericTypes.jl` for an example of how to implement such a cache (which
usually makes use of a dictionary).

## Required functionality for matrices

In addition to the required (relevant) functionality for the Ring interface (see above),
the following functionality is required for the Matrix interface.

We suppose that `R` is a fictitious base ring (coefficient ring) and that `S` is a
space of $m\times n$ matrices over `R`, or algebra of $m\times m$ matrices with parent
object `S` of type `MyMatSpace{T}` or `MyMatAlgebra{T}`, respectively. We also assume
the matrices in the space have type `MyMat{T}`, where `T` is the type of elements of
the base (element) ring.

Of course, in practice these types may not be parameterised, but we use parameterised
types here to make the interface clearer.

Note that the type `T` must (transitively) belong to the abstract type `RingElem`.

Currently only matrices over commutative rings are supported.

### Constructors

In addition to the standard constructors, the following constructors, taking an array of
elements, must be available.

```julia
(S::MyMatSpace{T})(A::Matrix{T}) where T <: RingElem
(S::MyMatAlgebra{T})(A::Matrix{T}) where T <: RingElem
```

Create the matrix in the given space/algebra whose $(i, j)$ entry is given by `A[i, j]`.

```julia
(S::MyMatSpace{T})(A::Matrix{S}) where {S <: RingElem, T <: RingElem}
(S::MyMatAlgebra{T})(A::Matrix{S}) where {S <: RingElem, T <: RingElem}
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

It is also possible to create matrices (in a matrix space only) directly, without first
creating the corresponding matrix space (the inner constructor being called directly).
Note that to support this, matrix space parent objects don't contain a reference to
their parent. Instead, parents are constructed on-the-fly if requested. (The same
strategy is used for matrix algebras.)

```julia
matrix(R::Ring, arr::Matrix{T}) where T <: RingElem
```

Given an $m\times n$ Julia matrix of entries, construct the corresponding
AbstractAlgebra.jl matrix over the given ring `R`, assuming all the entries can be
coerced into `R`.


```julia
matrix(R::Ring, r::Int, c::Int, A::Vector{T}) where T <: RingElem
```

Construct the given $r\times c$ AbstractAlgebra.jl matrix over the ring `R` whose
$(i, j)$ entry is given by `A[c*(i - 1) + j]`, assuming that all the entries can be
coerced into `R`.

```julia
zero_matrix(R::Ring, r::Int, c::Int)
```

Construct the $r\times c$ AbstractAlgebra.jl zero matrix over the ring `R`.

### Views

Just as Julia supports views of matrices, AbstractAlgebra requires all matrix
types to support views. These allow one to work with a submatrix of a given
matrix. Modifying the submatrix also modifies the original matrix.

Note that `deepcopy` of a view type must return the same type, but it should
return a view into a `deepcopy` of the original matrix. Julia enforces this
for consistency.

To support views, generic matrices in AbstractAlgebra of type
`Generic.MatSpaceElem` have an associated `Generic.MatSpaceView` type. Both
belong to the `Generic.Mat` abstract type, so that one can work with that in
functions that can accept both views and actual matrices.

The syntax for views is as for Julia's own views.

Note that the `parent_type` function returns the same type for a view as
for the original matrix type. This could potentially cause a problem if the
`elem_type` function is applied to the return value of `parent_type` and then
used in a type assertion. For this reason, there may be some limitations on
the use of views.

The `similar` function also returns a matrix of type `MatSpaceElem` when
applied to a view, rather than another view.

### Basic manipulation of matrices

```julia
dense_matrix_type(::Type{T}) where T<:NCRingElement
dense_matrix_type(::T) where T<:NCRingElement
dense_matrix_type(::Type{S}) where S<:NCRing
dense_matrix_type(::S) where S<:NCRing
```

Return the type of dense matrices whose entries have type `T` respectively
`elem_type(S)`. It suffices to provide a method with the first signature.
For the other three signatures, the default methods dispatch to the first.
E.g. in Nemo, which depends on AbstractAlgebra, we define
`dense_matrix_type(::Type{ZZRingElem}) = ZZMatrix`.

```julia
number_of_rows(M::MyMatSpace{T}) where T <: RingElem
number_of_rows(M::MyMatAlgebra{T}) where T <: RingElem
```

Return the number of rows of matrices in the matrix space.

```julia
number_of_columns(M:MyMatSpace{T}) where T <: RingElem
number_of_columns(M:MyMatAlgebra{T}) where T <: RingElem
```

Return the number of columns of matrices in the matrix space.

```julia
number_of_rows(f::MyMat{T}) where T <: RingElem
```

Return the number of rows of the given matrix.

```julia
number_of_columns(f::MyMat{T}) where T <: RingElem
```

Return the number of columns of the given matrix.

```julia
getindex(M::MyMat{T}, r::Int, c::Int) where T <: RingElem
```

Return the $(i, j)$-th entry of the matrix $M$.

```julia
setindex!(M::MyMat{T}, d::T, r::Int, c::Int) where T <: RingElem
```

Set the $(i, j)$-th entry of the matrix $M$ to $d$, which is assumed to be in the
base ring of the matrix. The matrix must have such an entry and the matrix is mutated
in place and not returned from the function.

### Transpose

```julia
transpose(::MyMat{T}) where T <: RingElem
```

Return the transpose of the given matrix.

## Optional functionality for matrices

Especially when wrapping C libraries, some functions are best implemented directly,
rather than relying on the generic functionality. The following are all provided by
the AbstractAlgebra.jl generic code, but can optionally be implemented directly for
performance reasons.

### Optional submatrices

The following are only available for matrix spaces, not for matrix algebras.

```julia
Base.getindex(M::MyMat, rows::AbstractVector{Int}, cols::AbstractVector{Int})
```

Return a new matrix with the same entries as the submatrix with the given range of rows
and columns.

### Optional row swapping

```julia
swap_rows!(M::MyMat{T}, i::Int, j::Int) where T <: RingElem
```

Swap the rows of `M` in place. The function returns the mutated matrix (since
matrices are assumed to be mutable in AbstractAlgebra.jl).

### Optional concatenation

The following are only available for matrix spaces, not for matrix algebras.

```julia
hcat(M::MyMat{T}, N::MyMat{T}) where T <: RingElem
```

Return the horizontal concatenation of $M$ and $N$. It is assumed that the number of
rows of $M$ and $N$ are the same.

```julia
vcat(M::MyMat{T}, N::MyMat{T}) where T <: RingElem
```

Return the vertical concatenation of $M$ and $N$. It is assumed that the number of
columns of $M$ and $N$ are the same.

### Optional zero tests

The following functions are available for matrices in both matrix algebras and matrix
spaces.

```julia
is_zero_entry(M::MatrixElem{T}, i::Int, j::Int) where T <: NCRingElement
is_zero_row(M::MatrixElem{T}, i::Int) where T <: NCRingElement
is_zero_column(M::MatrixElem{T}, j::Int) where T <: NCRingElement
```

### Optional similar and zero

The following functions are available for matrices in both matrix algebras and matrix
spaces. Both `similar` and `zero` construct new matrices, with the same methods, but
the entries are either undefined with `similar` or zero-initialized with `zero`.

```julia
similar(x::MyMat{T}, R::Ring=base_ring(x)) where T <: RingElem
zero(x::MyMat{T}, R::Ring=base_ring(x)) where T <: RingElem
```

Construct the matrix with the same dimensions as the given matrix, and the
same base ring unless explicitly specified.

```julia
similar(x::MyMat{T}, R::Ring, r::Int, c::Int) where T <: RingElem
similar(x::MyMat{T}, r::Int, c::Int) where T <: RingElem
zero(x::MyMat{T}, R::Ring, r::Int, c::Int) where T <: RingElem
zero(x::MyMat{T}, r::Int, c::Int) where T <: RingElem
```

Construct the $r\times c$ matrix with `R` as base ring (which defaults to the
base ring of the the given matrix).
If $x$ belongs to a matrix algebra and $r \neq c$, an exception is raised, and it's
also possible to specify only one `Int` as the order (e.g. `similar(x, n)`).

Custom matrices and rings may choose which specific matrix type is best-suited to return for the
given ring and dimensionality. If they do not specialize these functions, the default is a
`Generic.MatSpaceElem` matrix, or `Generic.MatAlgElem` for matrix algebras.
The default implementation of `zero` calls out to `similar`, so it's generally
sufficient to specialize only `similar`.
For both `similar` and `zero`, only the most general method has to be implemented
(e.g. `similar(x::MyMat, R::Ring, r::Int, c::Int)`, as all other methods (which have
defaults) call out to this more general method.

```julia
Base.isassigned(M::MyMat, i, j)
```

Test whether the given matrix has a value associated with indices `i` and `j`.
It is recommended to overload this method for custom matrices.

### Optional symmetry test

```julia
is_symmetric(a::MatrixElem)
```

Return `true` if the given matrix is symmetric with respect to its main diagonal,
otherwise return `false`.

