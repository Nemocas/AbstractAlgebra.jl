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
lu(A::MatrixElem{T}, P = SymmetricGroup(nrows(A))) where {T <: FieldElement}
fflu(A::MatrixElem{T}, P = SymmetricGroup(nrows(A))) where {T <: RingElement}
```


## Reduced row-echelon form

```@docs
rref_rational(M::MatrixElem{T}) where {T <: RingElement}
rref(M::MatrixElem{T}) where {T <: FieldElement}
```


## Hessenberg form and similarity transformations

```@docs
hessenberg(A::MatElem{T}) where {T <: RingElement}
similarity!(A::MatrixElem{T}, r::Int, d::T) where {T <: RingElement}
```


## Hermite normal form

```@docs
hnf(A::MatElem{T}) where {T <: RingElement}
hnf_with_transform(A::MatElem{T}) where {T <: RingElement}
```


## Smith normal form

```@docs
snf(A::MatElem{T}) where {T <: RingElement}
snf_with_transform(A::MatElem{T}) where {T <: RingElement}
```


## Popov forms

```@docs
weak_popov(A::MatElem{T}) where {T <: PolyRingElem}
weak_popov_with_transform(A::MatElem{T}) where {T <: PolyRingElem}
popov(A::MatElem{T}) where {T <: PolyRingElem}
popov_with_transform(A::MatElem{T}) where {T <: PolyRingElem}
```
