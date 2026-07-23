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
number_of_rows(a::MatElem)
number_of_columns(a::MatElem)
length(a::MatrixElem{T}) where T <: NCRingElement
```


## Trace, determinant and rank

```@docs
tr(x::MatElem{T}) where T <: NCRingElement
det(M::MatElem{T}) where {T <: RingElement}
rank(::MatElem{T}) where T <: RingElem
```


## Inverses

```@docs
Base.inv(M::MatElem{T}) where {T <: RingElement}
pseudo_inv(M::MatElem{T}) where {T <: RingElement}
```



## Characteristic and minimal polynomials

```@docs
charpoly(S::PolyRing{T}, Y::MatElem{T}) where {T <: RingElement}
minpoly(S::PolyRing{T}, M::MatElem{T}) where {T <: RingElement}
```


## Powers

```@docs
powers(::MatElem, ::Int)
```


## Gram matrices

```@docs
gram(x::MatElem)
```


## Content

```@docs
content(::MatElem{T}) where T <: RingElement
```


## Minors and exterior powers

```@docs
minors(A::MatElem, k::Int)
minors_with_position(A::MatElem, k::Int)
minors_iterator(A::MatElem, k::Int)
minors_iterator_with_position(M::MatElem, k::Int)
exterior_power(A::MatElem, k::Int)
```


## Pfaffians

```@docs
pfaffian(M::MatElem)
pfaffians(M::MatElem, k::Int)
```
