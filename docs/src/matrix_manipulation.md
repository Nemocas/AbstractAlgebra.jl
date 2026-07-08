```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# [Manipulating algebraic matrices](@id matrix_manipulation)

This page describes functions for modifying or rearranging algebraic matrices.
Many operations are available both as in-place and non-mutating variants.


## Entry-wise operations

```@docs
map_entries(f, a::MatElem{T}) where T <: NCRingElement
map_entries!(f, dst::MatElem{T}, src::MatElem{U}) where {T <: NCRingElement, U <: NCRingElement}
Base.map(f, a::MatrixElem{T}) where T <: NCRingElement
Base.map!(f, dst::MatrixElem{T}, src::MatrixElem{U}) where {T <: NCRingElement, U <: NCRingElement}
```


## Elementary row and column operations

```@docs
add_column(a::MatrixElem{T}, s::RingElement, i::Int, j::Int, rows = 1:nrows(a)) where T <: RingElement
add_column!(a::MatrixElem{T}, s::RingElement, i::Int, j::Int, rows = 1:nrows(a)) where T <: RingElement
add_row(a::MatrixElem{T}, s::RingElement, i::Int, j::Int, cols = 1:ncols(a)) where T <: RingElement
add_row!(a::MatrixElem{T}, s::RingElement, i::Int, j::Int, cols = 1:ncols(a)) where T <: RingElement
multiply_column(a::MatrixElem{T}, s::RingElement, i::Int, rows = 1:nrows(a)) where T <: RingElement
multiply_column!(a::MatrixElem{T}, s::RingElement, i::Int, rows = 1:nrows(a)) where T <: RingElement
multiply_row(a::MatrixElem{T}, s::RingElement, i::Int, cols = 1:ncols(a)) where T <: RingElement
multiply_row!(a::MatrixElem{T}, s::RingElement, i::Int, cols = 1:ncols(a)) where T <: RingElement
```


## Row and column permutations

```@docs
*(P::Perm, x::MatrixElem{T}) where T <: NCRingElement
*(x::MatrixElem{T}, P::Perm) where T <: NCRingElement
swap_rows(a::MatElem{T}, i::Int, j::Int) where T <: NCRingElement
swap_rows!(a::MatElem{T}, i::Int, j::Int) where T <: NCRingElement
swap_cols(a::MatElem{T}, i::Int, j::Int) where T <: NCRingElement
swap_cols!(a::MatElem{T}, i::Int, j::Int) where T <: NCRingElement
reverse_rows(a::MatElem{T}) where T <: NCRingElement
reverse_rows!(a::MatElem{T}) where T <: NCRingElement
reverse_cols(a::MatElem{T}) where T <: NCRingElement
reverse_cols!(a::MatElem{T}) where T <: NCRingElement
```


## Transposition

Matrices can be transposed either by creating a new matrix or by modifying an
existing square matrix in place.

```@docs
transpose(x::MatElem)
transpose!(x::MatElem)
```
