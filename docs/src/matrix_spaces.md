# Matrix Interface

Generic matrices are supported in AbstractAlgebra.jl. Both the space of $m\times n$ 
matrices and the algebra (ring) of $m\times m$ matrices are supported.

As the space of $m\times n$ matrices over a commutative ring is not itself a commutative
ring, not all of the Ring interface needs to be implemented for such matrices in.

In particular, the following functions do not need to be implemented: `isdomain_type`,
`needs_parentheses`, `displayed_with_minus_in_front`, `show_minus_one` and `divexact`.
The `canonical_unit` function should be implemented, but simply needs to return the
corresponding value for entry $[1, 1]$ (the function is never called on empty matrices).

For matrix algebras, all of the ring interface must be implemented.

Note that AbstractAlgebra.jl matrices are not the same as Julia matrices. We store a
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
(S::MyMatSpace{T})(A::Array{T, 2}) where T <: AbstractAlgebra.RingElem
(S::MyMatAlgebra{T})(A::Array{T, 2}) where T <: AbstractAlgebra.RingElem
```

Create the matrix in the given space/algebra whose $(i, j)$ entry is given by `A[i, j]`.

```julia
(S::MyMatSpace{T})(A::Array{S, 2}) where {S <: AbstractAlgebra.RingElem, T <: AbstractAlgebra.RingElem}
(S::MyMatAlgebra{T})(A::Array{S, 2}) where {S <: AbstractAlgebra.RingElem, T <: AbstractAlgebra.RingElem}
```

Create the matrix in the given space/algebra whose $(i, j)$ entry is given by `A[i, j]`,
where `S` is the type of elements that can be coerced into the base ring of the matrix.

```julia
(S::MyMatSpace{T})(A::Array{S, 1}) where {S <: AbstractAlgebra.RingElem, T <: AbstractAlgebra.RingElem}
(S::MyMatAlgebra{T})(A::Array{S, 1}) where {S <: AbstractAlgebra.RingElem, T <: AbstractAlgebra.RingElem}
```

Create the matrix in the given space/algebra of matrices (with dimensions $m\times n$
say), whose $(i, j)$ entry is given by `A[i*(n - 1) + j]` and where `S` is the type of
elements that can be coerced into the base ring of the matrix.

**Examples**

```julia
S = MatrixSpace(QQ, 2, 3)
T = MatrixAlgebra(QQ, 2)

M1 = S(Rational{BigInt}[2 3 1; 1 0 4])
M2 = S(BigInt[2 3 1; 1 0 4])
M3 = S(BigInt[2, 3, 1, 1, 0, 4])

N1 = T(Rational{BigInt}[2 3; 1 0])
N2 = T(BigInt[2 3; 1 0])
N3 = T(BigInt[2, 3, 1, 1])
```

It is also possible to create matrices (in a matrix space only) directly, without first
creating the corresponding matrix space (the inner constructor being called directly).
Note that to support this, matrix space parent objects don't contain a reference to
their parent. Instead, parents are constructed on-the-fly if requested. (The same
strategy is used for matrix algebras.)

```julia
matrix(R::Ring, arr::Array{T, 2}) where T <: AbstractAlgebra.RingElem
```

Given an $m\times n$ Julia matrix of entries, construct the corresponding
AbstractAlgebra.jl matrix over the given ring `R`, assuming all the entries can be
coerced into `R`.


```julia
matrix(R::Ring, r::Int, c::Int, A::Array{T, 1}) where T <: AbstractAlgebra.RingElem
```

Construct the given $r\times c$ AbstractAlgebra.jl matrix over the ring `R` whose
$(i, j)$ entry is given by `A[c*(i - 1) + j]`, assuming that all the entries can be
coerced into `R`.

```julia
zero_matrix(R::Ring, r::Int, c::Int)
```

Construct the $r\times c$ AbstractAlgebra.jl zero matrix over the ring `R`.

```julia
identity_matrix(R::Ring, n::Int)
```

Construct the $n\times n$ AbstractAlgebra.jl identity matrix over the ring `R`.

The following functions are available for matrices in both matrix algebras and matrix
spaces.

```julia
similar(x::MyMat{T}) where T <: AbstractAlgebra.RingElem
```

Construct the zero matrix with the same dimensions and base ring as the given matrix.

```julia
similar(x::MyMat{T}, r::Int, c::Int) where T <: AbstractAlgebra.RingElem
```

Construct the $r\times c$ zero matrix with the same base ring as the given matrix. If
$x$ belongs to a matrix algebra and $r \neq c$, an exception is raised.

**Examples**

```julia
M = matrix(ZZ, BigInt[3 1 2; 2 0 1])
N = matrix(ZZ, 3, 2, BigInt[3, 1, 2, 2, 0, 1])
P = zero_matrix(ZZ, 3, 2)
Q = identity_matrix(ZZ, 4)
C = similar(P)
D = similar(Q, 4, 5)

R = MatrixAlgebra(ZZ, 2)
M = R()
F = similar(M)
```

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

Note that the parent of a view will be the same as the parent of the original
matrix. The `parent_type` function also returns the same type for a view as
for the original matrix type. This could potentially cause a problem if the
`elem_type` function is applied to the return value of `parent_type` and then
used in a type assertion. For this reason, there may be some limitations on
the use of views.

The `similar` function also returns a matrix of type `MatSpaceElem` when
applied to a view, rather than another view.

```julia
M = matrix(ZZ, 3, 3, BigInt[1, 2, 3, 2, 3, 4, 3, 4, 5])

N1 = @view M[1:2, :]
N2 = @view M[:, 1:2]

R = N1*N2
```

### Basic manipulation of matrices

```julia
dense_matrix_type(::Type{T}) where T <: AbstractAlgebra.RingElem
```

Returns the type of dense matrices whose entries have the given type. E.g.
in Nemo, which depends on AbstractAlgebra, we define
`dense_matrix_type(::Type{fmpz}) = fmpz_mat`.

```julia
nrows(f::MyMat{T}) where T <: AbstractAlgebra.RingElem
```

Return the number of rows of the given matrix.

```julia
ncols(f::MyMat{T}) where T <: AbstractAlgebra.RingElem
```

Returns the number of columns of the given matrix.

```julia
getindex(M::MyMat{T}, r::Int, c::Int) where T <: AbstractAlgebra.RingElem
```

Return the $(i, j)$-th entry of the matrix $M$.

```julia
setindex!(M::MyMat{T}, d::T, r::Int, c::Int) where T <: AbstractAlgebra.RingElem
```

Set the $(i, j)$-th entry of the matrix $M$ to $d$, which is assumed to be in the
base ring of the matrix. The matrix must have such an entry and the matrix is mutated
in place and not returned from the function.

**Examples**

```julia
M = matrix(ZZ, BigInt[2 3 0; 1 1 1])

m = nrows(M)
n = ncols(M)
M[1, 2] = BigInt(4)
c = M[1, 1]
```
### Transpose

```julia
transpose(::MyMat{T}) where T <: AbstractAlgebra.RingElem
```

Return the transpose of the given matrix.

The standard Julia tick notation can also be used for transposing a matrix.

**Examples**

```julia
R, t = PolynomialRing(QQ, "t")
S = MatrixSpace(R, 3, 3)

A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])

B = transpose(A)
C = A'
```

## Optional functionality for matrices

Especially when wrapping C libraries, some functions are best implemented directly,
rather than relying on the generic functionality. The following are all provided by
the AbstractAlgebra.jl generic code, but can optionally be implemented directly for
performance reasons.

### Optional submatrices

The following are only available for matrix spaces, not for matrix algebras.

```julia
sub(M::MyMat{T}, rows::UnitRange{Int}, cols::UnitRange{Int}) where T <: AbstractAlgebra.RingElem
```

Return a new matrix with the same entries as the submatrix with the given range of rows
and columns.

**Examples**

```julia
M = matrix(ZZ, BigInt[1 2 3; 2 3 4; 3 4 5])

N1 = M[1:2, :]
N2 = M[:, :]
N3 = M[2:3, 2:3]
```

### Optional row swapping

```julia
swap_rows!(M::MyMat{T}, i::Int, j::Int) where T <: AbstractAlgebra.RingElem
```

Swap the rows of `M` in place. The function returns the mutated matrix (since
matrices are assumed to be mutable in AbstractAlgebra.jl).

**Examples**

```julia
M = identity_matrix(ZZ, 3)

swap_rows!(M, 1, 2)
```

### Optional concatenation

The following are only available for matrix spaces, not for matrix algebras.

```julia
hcat(M::MyMat{T}, N::MyMat{T}) where T <: AbstractAlgebra.RingElem
```

Return the horizontal concatenation of $M$ and $N$. It is assumed that the number of
rows of $M$ and $N$ are the same.

```julia
vcat(M::MyMat{T}, N::MyMat{T}) where T <: AbstractAlgebra.RingElem
```

Return the vertical concatenation of $M$ and $N$. It is assumed that the number of
columns of $M$ and $N$ are the same.

**Examples**

```julia
M = matrix(ZZ, BigInt[1 2 3; 2 3 4; 3 4 5])
N = matrix(ZZ, BigInt[1 0 1; 0 1 0; 1 0 1])

P = hcat(M, N)
Q = vcat(M, N)
```
