```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# [Matrix predicates](@id matrix_predicates)

This page collects predicates, i.e. functions returning `true` or `false`,
for testing whether matrices satisfy certain structural or normal-form conditions.


## Basic predicates

Matrices support `iszero` and `isone` for testing whether a matrix
is the zero matrix or the identity matrix, respectively.

```@docs
isempty(A::MatrixElem{T}) where {T <: NCRingElement}
Base.isassigned(A::MatrixElem{T}, i::Int, j::Int) where {T <: NCRingElement}
is_zero_row(M::Union{Matrix,MatrixElem}, i::Int)
is_zero_column(M::Union{Matrix,MatrixElem}, j::Int)
```


## Triangular and diagonal matrices

```@docs
is_lower_triangular(A::MatElem)
is_upper_triangular(A::MatElem)
is_diagonal(A::MatElem)
is_hessenberg(A::MatElem{T}) where {T <: RingElement}
```


## Invertibility

```@docs
is_invertible_with_inverse(A::MatrixElem{T}; side::Symbol = :left) where {T <: RingElement}
is_invertible(A::MatElem{T}) where {T <: RingElement}
```


## Symmetry

```@docs
is_symmetric(A::MatElem)
is_skew_symmetric(A::MatElem)
is_alternating(A::MatElem)
```


## Nilpotency

```@docs
is_nilpotent(A::MatElem{T}) where {T <: RingElement}
```


## Normal forms

```@docs
is_rref(A::MatrixElem{T}) where {T <: RingElement}
is_hnf(A::MatElem{T}) where {T <: RingElement}
is_snf(A::MatElem{T}) where {T <: RingElement}
is_weak_popov(P::MatrixElem{T}, rank::Int) where {T <: PolyRingElem}
is_popov(P::MatrixElem{T}, rank::Int) where {T <: PolyRingElem}
```
