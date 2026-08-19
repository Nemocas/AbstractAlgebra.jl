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
map_entries(f, m::MatElem{T}) where T <: NCRingElement
map_entries!(f, dst::MatElem{T}, src::MatElem{U}) where {T <: NCRingElement, U <: NCRingElement}
Base.map(f, m::MatrixElem{T}) where T <: NCRingElement
Base.map!(f, dst::MatrixElem{T}, src::MatrixElem{U}) where {T <: NCRingElement, U <: NCRingElement}
```


## Elementary row and column operations

```@docs
add_column(m::MatrixElem{T}, s::RingElement, i::Int, j::Int, rows = 1:nrows(m)) where T <: RingElement
add_column!(m::MatrixElem{T}, s::RingElement, i::Int, j::Int, rows = 1:nrows(m)) where T <: RingElement
add_row(m::MatrixElem{T}, s::RingElement, i::Int, j::Int, cols = 1:ncols(m)) where T <: RingElement
add_row!(m::MatrixElem{T}, s::RingElement, i::Int, j::Int, cols = 1:ncols(m)) where T <: RingElement
multiply_column(m::MatrixElem{T}, s::RingElement, i::Int, rows = 1:nrows(m)) where T <: RingElement
multiply_column!(m::MatrixElem{T}, s::RingElement, i::Int, rows = 1:nrows(m)) where T <: RingElement
multiply_row(m::MatrixElem{T}, s::RingElement, i::Int, cols = 1:ncols(m)) where T <: RingElement
multiply_row!(m::MatrixElem{T}, s::RingElement, i::Int, cols = 1:ncols(m)) where T <: RingElement
```


## Row and column permutations

```@docs
*(P::Perm, m::MatrixElem{T}) where T <: NCRingElement
*(m::MatrixElem{T}, P::Perm) where T <: NCRingElement
swap_rows(m::MatElem{T}, i::Int, j::Int) where T <: NCRingElement
swap_rows!(m::MatElem{T}, i::Int, j::Int) where T <: NCRingElement
swap_cols(m::MatElem{T}, i::Int, j::Int) where T <: NCRingElement
swap_cols!(m::MatElem{T}, i::Int, j::Int) where T <: NCRingElement
reverse_rows(m::MatElem{T}) where T <: NCRingElement
reverse_rows!(m::MatElem{T}) where T <: NCRingElement
reverse_cols(m::MatElem{T}) where T <: NCRingElement
reverse_cols!(m::MatElem{T}) where T <: NCRingElement
```


## Transposition

Matrices can be transposed either by creating a new matrix or by modifying an
existing square matrix in place.

```@docs
transpose(m::MatElem)
transpose!(m::MatElem)
```
