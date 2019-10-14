```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Generic matrix algebras

AbstractAlgebra.jl allows the creation of an algebra (ring) of $m\times m$ matrices over
a computable, commutative ring.

Functions specific to generic matrix algebras of $m\times m$ matrices are implemented in
`src/generic/MatrixAlgebra.jl`. The remaining functionality is in the file
`src/generic/Matrix.jl`.

As well as implementing the entire Matrix interface, including the optional
functionality, there are many additional generic algorithms implemented for matrix
algebras.

Almost all of the functionality specified for generic matrices is available for matrix
algebras. The exceptions are functions such as `solve` and `nullspace` which may
return non-square matrices, or which don't accept square matrices.

All of the generic functionality is part of the Generic submodule of
AbstractAlgebra.jl. This is exported by default, so it is not necessary to qualify names
of functions.

## Types and parent objects

Generic matrices in AbstractAlgebra.jl have type `Generic.MatAlgElem{T}` for matrices in
a matrix algebra, where `T` is the type of elements of the matrix. Internally, generic
matrices are implemented using an object wrapping a Julia two dimensional array, though
they are not themselves Julia arrays. See the file `src/generic/GenericTypes.jl` for
details.

Parents of generic matrices in a matrix algebra have type `Generic.MatAlgebra{T}`.

Note that matrix algebras are noncommutative rings. Thus their types belong to `NCRing`
and `NCRingElem`. They cannot be used in constructions which require a commutative
ring (`Ring` and `RingElem` respectively).

The generic matrix algebra matrix types belong to the abstract type
`AbstractAlgebra.MatAlgElem{T}` and the parent types belong to
 `AbstractAlgebra.MatAlgebra{T}` Note that both of these require disambiguation from
the concrete types in `Generic` of the same name.

The degree and base ring $R$ of a generic matrix are stored in its parent object,
however to allow creation of matrices without first creating the matrix space parent,
generic matrices in Julia do not contain a reference to their parent. They contain the
row and column numbers (or degree, in the case of matrix algebras) and the base ring
on a per matrix basis. The parent object can then be reconstructed from this data on
demand.

## Matrix algebra constructors

A matrix algebra in AbstractAlgebra.jl represents a collection of all matrices with
given degree and base ring.

In order to construct matrices in AbstractAlgebra.jl, one must construct the
matrix algebra itself. This is accomplished with the following constructor.

```julia
MatrixAlgebra(R::Ring, degree::Int; cache::Bool=true)
```

Construct the algebra of matrices with the given degree over the given base ring. By
default such matrix spaces are cached based on the base ring and degree. If the
optional named parameter `cached` is set to false, no caching occurs.

Here are some examples of creating matrix algebras and making use of the
resulting parent objects to coerce various elements into the matrix algebra.

**Examples**

```jldoctest
julia> R, t = PolynomialRing(QQ, "t")
(Univariate Polynomial Ring in t over Rationals, t)

julia> S = MatrixAlgebra(R, 3)
Matrix Algebra of degree 3 over Univariate Polynomial Ring in t over Rationals

julia> A = S()
[0//1  0//1  0//1]
[0//1  0//1  0//1]
[0//1  0//1  0//1]

julia> B = S(12)
[12//1   0//1   0//1]
[ 0//1  12//1   0//1]
[ 0//1   0//1  12//1]

julia> C = S(R(11))
[11//1   0//1   0//1]
[ 0//1  11//1   0//1]
[ 0//1   0//1  11//1]

```

## Matrix algebra element constructors

The following additional constructors are provided for constructing various
kinds of matrices in a matrix algebra.

```@docs
zero(::MatAlgebra)
```

```@docs
identity_matrix(::MatAlgElem{T}) where T <: RingElement                                  ```

*Examples*

```julia
S = MatrixAlgebra(ZZ, 2)
M = zero(S)

P = identity_matrix(M)
```

## Matrix algebra functionality provided by AbstractAlgebra.jl

Most of the generic matrix functionality described in the generic matrix section of
the documentation is available for both matrix spaces and matrix algebras. Exceptions
include functions that do not return or accept square matrices or which cannot specify
a parent. Such functions include `solve` and `nullspace` which can't be provided for
matrix algebras.

In addition to the functionality described for matrix spaces, matrix algebras support
all noncommutative ring operations, and matrix algebras can be used as a base ring
for other generic constructs that accept a *noncommutative* base ring (`NCRing`).

In this section we describe functionality provided for matrix algebras only.

### Basic matrix functionality

As well as the Ring and Matrix interfaces, the following functions are provided to
manipulate matrices.

```@docs
degree(::Generic.MatAlgElem)
```

**Examples**

```jldoctest
julia> R, t = PolynomialRing(QQ, "t")
(Univariate Polynomial Ring in t over Rationals, t)

julia> S = MatrixAlgebra(R, 3)
Matrix Algebra of degree 3 over Univariate Polynomial Ring in t over Rationals

julia> A = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])
[t+1//1       t        1//1]
[   t^2       t           t]
[ -2//1  t+2//1  t^2+t+1//1]

julia> n = degree(A)
3

```
