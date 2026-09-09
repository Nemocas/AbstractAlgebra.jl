```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# [Matrix properties](@id matrix_properties)

This page collects functions which compute values, invariants,
or associated objects from matrices, rather than testing whether matrices
satisfy certain conditions. Examples include determinants, ranks, inverses,
nullspaces, polynomials, and related constructions.


## Basic properties

```@docs
number_of_rows(A::MatElem)
number_of_columns(A::MatElem)
length(A::MatrixElem{T}) where T <: NCRingElement
```


## Trace, determinant and rank

```@docs
tr(A::MatElem{T}) where T <: NCRingElement
det(A::MatElem{T}) where {T <: RingElement}
rank(A::MatElem{T}) where T <: RingElem
```


## Inverses

```@docs
Base.inv(A::MatElem{T}) where {T <: RingElement}
pseudo_inv(A::MatElem{T}) where {T <: RingElement}
```



## Characteristic and minimal polynomials

```@docs
charpoly(S::PolyRing{T}, A::MatElem{T}) where {T <: RingElement}
minpoly(S::PolyRing{T}, A::MatElem{T}) where {T <: RingElement}
```


## Powers

```@docs
powers(A::MatElem, d::Int)
```


## Gram matrices

```@docs
gram(A::MatElem)
```


## Content

```@docs
content(A::MatrixElem{T}) where T <: RingElement
```


## Minors and exterior powers

```@docs
minors(A::MatElem, k::Int)
minors_with_position(A::MatElem, k::Int)
minors_iterator(A::MatElem, k::Int)
minors_iterator_with_position(A::MatElem, k::Int)
exterior_power(A::MatElem, k::Int)
```


## Pfaffians

```@docs
pfaffian(A::MatElem)
pfaffians(A::MatElem, k::Int)
```
