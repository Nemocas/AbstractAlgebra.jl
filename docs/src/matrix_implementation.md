```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# Matrix implementation

This page gives an overview of how matrices are implemented in
AbstractAlgebra.jl. It is intended for contributors who want to understand the
source layout and the relation between the matrix interface, the generic dense
matrix implementation, and matrix algebras.

For the list of functions that a new matrix type has to provide, see the
matrix interface documentation.

## Source layout

The matrix code is split into interface-level functionality and concrete
generic implementations.

The file `src/Matrix.jl` contains functionality for matrices belonging to
matrix spaces. This includes parent and element methods, constructors,
indexing, iteration, views, arithmetic, basic predicates, block constructions,
concatenation, linear algebra routines, and generic fallback algorithms.

The file `src/MatRing.jl` contains functionality for square matrices viewed as
elements of a matrix algebra. Matrix algebras are rings, and this file provides
the corresponding ring operations and matrix-algebra-specific functionality.

The file `src/generic/Matrix.jl` contains the dense generic implementation for
matrices in matrix spaces. These matrices store their entries in a Julia
two-dimensional array internally, but they are not themselves Julia arrays.

The file `src/generic/MatRing.jl` contains the dense generic implementation for
matrix algebras. A generic matrix algebra element wraps a generic dense matrix
and equips it with the parent and ring structure of a matrix algebra.

Additional algorithms are implemented in separate files. For example,
`src/Matrix-Strassen.jl` contains Strassen-style and block matrix algorithms,
while `src/MatrixNormalForms.jl` contains normal form functionality such as
echelon, Hermite, and Howell forms.

## Matrix spaces vs matrix algebras

AbstractAlgebra distinguishes between matrices in matrix spaces and matrices in
matrix algebras.

A matrix space contains matrices of fixed size $m \times n$ over a base ring.
Its elements belong to the abstract type `MatElem{T}`, where `T` is the type of
the matrix entries. The parent object has type `MatSpace{T}`.

A matrix algebra contains square matrices of fixed degree $n$ over a base ring,
with multiplication as the ring multiplication. Its elements belong to the
abstract type `MatRingElem{T}`. Matrix algebra parent objects belong to
`MatRing{T}`. Since matrix multiplication is generally noncommutative, matrix
algebras belong to the noncommutative ring hierarchy.

## Implementation of generic matrix spaces

The generic dense matrix implementation is provided by the `Generic` module.

Matrices in generic matrix spaces have type `Generic.MatSpaceElem{T}`. Views of
such matrices have type `Generic.MatSpaceView{T}`. Both are subtypes of the
abstract type `Generic.Mat{T}`, so generic code can dispatch on `Generic.Mat`
when it should accept both ordinary dense matrices and views.

Generic matrix-space elements and views store their entries in Julia array-like
storage internally. The base ring is stored with the matrix object, and the
number of rows and columns is obtained from the internal storage.

## Implementation of generic matrix algebras

Generic matrix algebra elements have type `Generic.MatRingElem{T}`. Internally,
such an element stores a generic matrix and uses it for entry access, mutation,
arithmetic, and linear algebra operations.

The parent of a generic matrix algebra element has type `Generic.MatRing{T}`.
It stores the base ring and the degree of the algebra.

Generic matrix algebras reuse much of the matrix functionality from
`src/Matrix.jl` and `src/generic/Matrix.jl`, but they additionally satisfy the
ring interface. In particular, they provide ring operations such as `zero`,
`one`, multiplication, promotion, exact division where available, and parent
construction through `matrix_ring`.

## Parents

For generic matrices, the parent is not stored as a direct reference in every
matrix element. Instead, enough data is stored on the element to reconstruct the
parent on demand.

For a matrix in a matrix space, the parent is reconstructed from the base ring
and the number of rows and columns. Calling `parent(M)` returns the corresponding
matrix space.

For a matrix algebra element, the parent is reconstructed from the base ring and
the degree of the underlying square matrix. Calling `parent(A)` returns the
corresponding matrix algebra.

This design allows matrices to be constructed directly, without first explicitly
constructing the parent object, while still supporting the parent-based model
used throughout AbstractAlgebra.

## Relation to the matrix interface

The matrix interface documentation describes what a matrix implementation must
provide: required constructors, entry access, mutation, dimensions, views, and
optional methods that can be implemented for performance.

This page instead describes the implementation supplied by AbstractAlgebra
itself and how the source files fit together. In particular:

* use the matrix interface page when implementing a new matrix type;
* use this page when navigating the AbstractAlgebra matrix source code;
* use the user-facing matrix pages for construction, manipulation, properties,
  transformations, matrix spaces, and matrix algebras.
