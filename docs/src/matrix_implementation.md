```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# Matrix implementation

AbstractAlgebra.jl provides a module, implemented in `src/Matrix.jl` for
matrices over any ring belonging to the AbstractAlgebra abstract type
hierarchy. This functionality will work for any matrix type which
follows the Matrix interface.

Similarly, AbstractAlgebra.jl provides a module in `src/MatRing.jl` for
matrix algebras over a ring.

## Generic matrix types

AbstractAlgebra.jl allows the creation of dense matrices over any computable
ring $R$. Generic matrices over a ring are implemented in
`src/generic/Matrix.jl`.

Generic matrix rings of $m\times m$ matrices are implemented in
`src/generic/MatRing.jl`.

Generic matrices in AbstractAlgebra.jl have type `Generic.MatSpaceElem{T}` for matrices
in a matrix space, or `Generic.MatRingElem{T}` for matrices in a matrix algebra, where
`T` is the type of elements of the matrix. Internally, generic matrices are implemented
using an object wrapping a Julia two dimensional array, though they are not themselves
Julia arrays. See the file `src/generic/GenericTypes.jl` for details.

For the most part, one doesn't want to work directly with the `MatSpaceElem` type though,
but with an abstract type called `Generic.Mat` which includes `MatSpaceElem` and views
thereof.

Parents of generic matrices (matrix spaces) have type `MatSpace{T}`. Parents of
matrices in a matrix algebra have type `Generic.MatRing{T}`.

The dimensions and base ring $R$ of a generic matrix are stored in its parent object,
however to allow creation of matrices without first creating the matrix space parent,
generic matrices in Julia do not contain a reference to their parent. They contain the
row and column numbers (or degree, in the case of matrix algebras) and the base ring
on a per matrix basis. The parent object can then be reconstructed from this data on
demand.

## Abstract types

The generic matrix types (matrix spaces) belong to the abstract type
`MatElem{T}` and the all matrix space parents are of the concrete type
`MatSpace{T}`.
On the other hand, the generic matrix algebra matrix types belong
to the abstract type `MatRingElem{T}` and the parent types belong to the abstract
 `MatRing{T}` Note that both
the concrete type of a matrix space parent object and the abstract class it belongs to
have the name `MatElem`, therefore disambiguation is required to specify which is
intended. The same is true for the abstract types for matrix spaces and their elements.

## Conversion to Julia matrices, iteration and broadcasting

While `AbstractAlgebra` matrices are not instances of `AbstractArray`,
they are closely related to Julia matrices. For convenience, a `Matrix`
and an `Array` constructors taking an `AbstractAlgebra` matrix as input
are provided:

```@docs
Matrix(::MatrixElem{T}) where T <: RingElement
Array(::MatrixElem{T}) where T <: RingElement
```

Matrices also support iteration, and therefore functions accepting an iterator
can be called on them, e.g.:

```jldoctest
julia> M = matrix_space(ZZ, 2, 3); x = M(1:6)
[1   2   3]
[4   5   6]

julia> collect(x)
2×3 Matrix{BigInt}:
 1  2  3
 4  5  6

julia> Set(x)
Set{BigInt} with 6 elements:
  5
  4
  6
  2
  3
  1
```

Matrices also support broadcasting, which amounts to elementwise application
of functions to matrices:

```jldoctest
julia> k = GF(5);

julia> A = ZZ[1 2; 3 4];

julia> k.(A)
[1   2]
[3   4]

julia> 3 .* A .+ 2
[ 5    8]
[11   14]

julia> B = ZZ[3 4; 5 6];

julia> ((x, y) -> x^2 + y^2).(A, B)
[10   20]
[34   52]
```

```@docs
dense_matrix_type(::Ring)
```

```jldoctest
julia> R, t = polynomial_ring(QQ, :t)
(Univariate polynomial ring in t over rationals, t)

julia> T = dense_matrix_type(R)
AbstractAlgebra.Generic.MatSpaceElem{AbstractAlgebra.Generic.Poly{Rational{BigInt}}}
```
