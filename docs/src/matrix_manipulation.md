```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# [Manipulating matrices](@id matrix_manipulation)

This page describes functions for modifying or rearranging matrices.
Many operations are available both as in-place and non-mutating variants.


## Entry-wise operations

```@docs
map_entries(f, M::MatElem{T}) where T <: NCRingElement
map_entries!(f, dst::MatElem{T}, src::MatElem{U}) where {T <: NCRingElement, U <: NCRingElement}
Base.map(f, M::MatrixElem{T}) where T <: NCRingElement
Base.map!(f, dst::MatrixElem{T}, src::MatrixElem{U}) where {T <: NCRingElement, U <: NCRingElement}
```


## Elementary row and column operations

```@docs
add_column(M::MatrixElem{T}, s::RingElement, i::Int, j::Int, rows = 1:nrows(M)) where T <: RingElement
add_column!(M::MatrixElem{T}, s::RingElement, i::Int, j::Int, rows = 1:nrows(M)) where T <: RingElement
add_row(M::MatrixElem{T}, s::RingElement, i::Int, j::Int, cols = 1:ncols(M)) where T <: RingElement
add_row!(M::MatrixElem{T}, s::RingElement, i::Int, j::Int, cols = 1:ncols(M)) where T <: RingElement
multiply_column(M::MatrixElem{T}, s::RingElement, i::Int, rows = 1:nrows(M)) where T <: RingElement
multiply_column!(M::MatrixElem{T}, s::RingElement, i::Int, rows = 1:nrows(M)) where T <: RingElement
multiply_row(M::MatrixElem{T}, s::RingElement, i::Int, cols = 1:ncols(M)) where T <: RingElement
multiply_row!(M::MatrixElem{T}, s::RingElement, i::Int, cols = 1:ncols(M)) where T <: RingElement
```


## Row and column permutations

```@docs
*(P::Perm, M::MatrixElem{T}) where T <: NCRingElement
*(M::MatrixElem{T}, P::Perm) where T <: NCRingElement
swap_rows(M::MatElem{T}, i::Int, j::Int) where T <: NCRingElement
swap_rows!(M::MatElem{T}, i::Int, j::Int) where T <: NCRingElement
swap_cols(M::MatElem{T}, i::Int, j::Int) where T <: NCRingElement
swap_cols!(M::MatElem{T}, i::Int, j::Int) where T <: NCRingElement
reverse_rows(M::MatElem{T}) where T <: NCRingElement
reverse_rows!(M::MatElem{T}) where T <: NCRingElement
reverse_cols(M::MatElem{T}) where T <: NCRingElement
reverse_cols!(M::MatElem{T}) where T <: NCRingElement
```


## Transposition

Matrices can be transposed either by creating a new matrix or by modifying an
existing square matrix in place.

```@docs
transpose(M::MatElem)
transpose!(M::MatElem)
```
