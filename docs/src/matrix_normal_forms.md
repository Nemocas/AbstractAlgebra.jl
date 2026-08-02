```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# [Matrix normal forms](@id matrix_normal_forms)

This page collects functionality for transforming matrices into normal
forms or related decompositions. Some algorithms also return transformation
matrices which certify the result.


## LU factorisation

```@docs
lu(m::MatrixElem{T}, P = SymmetricGroup(nrows(m))) where {T <: FieldElement}
fflu(m::MatrixElem{T}, P = SymmetricGroup(nrows(m))) where {T <: RingElement}
```


## Reduced row-echelon form

```@docs
rref_rational(m::MatrixElem{T}) where {T <: RingElement}
rref(m::MatrixElem{T}) where {T <: FieldElement}
```


## Hessenberg form and similarity transformations

```@docs
hessenberg(m::MatElem{T}) where {T <: RingElement}
similarity!(m::MatrixElem{T}, r::Int, d::T) where {T <: RingElement}
```


## Hermite normal form

```@docs
hnf(m::MatElem{T}) where {T <: RingElement}
hnf_with_transform(m::MatElem{T}) where {T <: RingElement}
```


## Smith normal form

```@docs
snf(A::MatElem{T}) where {T <: RingElement}
snf_with_transform(m::MatElem{T}) where {T <: RingElement}
```


## Popov forms

```@docs
weak_popov(m::MatElem{T}) where {T <: PolyRingElem}
weak_popov_with_transform(m::MatElem{T}) where {T <: PolyRingElem}
popov(m::MatElem{T}) where {T <: PolyRingElem}
popov_with_transform(m::MatElem{T}) where {T <: PolyRingElem}
```
