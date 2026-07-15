```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# Matrix algebras

A matrix algebra (or matrix ring) represents the collection of all square matrices
of a fixed number of rows, referred to as the degree of the matrix algebra, over a
fixed base ring. Matrix algebras are parent objects; their elements are the
corresponding matrices.

Matrix algebras are noncommutative rings. They support the usual ring
operations, such as addition and multiplication, and can be used in
constructions accepting noncommutative base rings.

Although elements of matrix algebras mostly behave like matrices, their parent object
is a matrix algebra rather than a matrix space. Therefore, some functionality
available for matrices in matrix spaces is not available for matrix algebra
elements, in particular functions which do not preserve square matrices (e.g.
`nullspace`).


## Creating matrix algebras

```@docs
matrix_ring(R::NCRing, n::Int)
```


## Properties of matrix algebras

```@docs
degree(a::MatRing)
```

The base ring can be queried using `base_ring`.


## Creating elements of a matrix algebra


### Calling a matrix algebra

```@docs
Generic.MatRing
```


### Special elements

As for other rings, the additive and multiplicative identities of a matrix
algebra can be constructed using `zero` and `one`.

```jldoctest
julia> S = matrix_ring(ZZ, 2)
Matrix ring of degree 2
  over integers

julia> zero(S)
[0   0]
[0   0]

julia> one(S)
[1   0]
[0   1]
```


## Properties of matrix algebra elements

```@docs
parent(a::MatRingElem{T}) where T <: NCRingElement
degree(a::MatRingElem{T}) where T <: NCRingElement
```
