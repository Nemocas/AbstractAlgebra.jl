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
isempty(m::MatrixElem{T}) where {T <: NCRingElement}
Base.isassigned(m::MatrixElem{T}, i::Int, j::Int) where {T <: NCRingElement}
is_zero_row(m::Union{Matrix,MatrixElem}, i::Int)
is_zero_column(m::Union{Matrix,MatrixElem}, j::Int)
```


## Triangular and diagonal matrices

```@docs
is_lower_triangular(m::MatElem)
is_upper_triangular(m::MatElem)
is_diagonal(m::MatElem)
is_hessenberg(m::MatElem{T}) where {T <: RingElement}
```


## Invertibility

```@docs
is_invertible_with_inverse(m::MatrixElem{T}; side::Symbol = :left) where {T <: RingElement}
is_invertible(m::MatElem{T}) where {T <: RingElement}
```


## Symmetry

```@docs
is_symmetric(m::MatElem)
is_skew_symmetric(m::MatElem)
is_alternating(m::MatElem)
```


## Nilpotency

```@docs
is_nilpotent(m::MatElem{T}) where {T <: RingElement}
```


## Normal forms

```@docs
is_rref(m::MatrixElem{T}) where {T <: RingElement}
is_hnf(m::MatElem{T}) where {T <: RingElement}
is_snf(m::MatElem{T}) where {T <: RingElement}
is_weak_popov(m::MatrixElem{T}, rank::Int) where {T <: PolyRingElem}
is_popov(m::MatrixElem{T}, rank::Int) where {T <: PolyRingElem}
```
