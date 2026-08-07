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
number_of_rows(m::MatElem)
number_of_columns(m::MatElem)
length(m::MatrixElem{T}) where T <: NCRingElement
```


## Trace, determinant and rank

```@docs
tr(m::MatElem{T}) where T <: NCRingElement
det(m::MatElem{T}) where {T <: RingElement}
rank(m::MatElem{T}) where T <: RingElem
```


## Inverses

```@docs
Base.inv(m::MatElem{T}) where {T <: RingElement}
pseudo_inv(m::MatElem{T}) where {T <: RingElement}
```



## Characteristic and minimal polynomials

```@docs
charpoly(S::PolyRing{T}, m::MatElem{T}) where {T <: RingElement}
minpoly(S::PolyRing{T}, m::MatElem{T}) where {T <: RingElement}
```


## Powers

```@docs
powers(m::MatElem, ::Int)
```


## Gram matrices

```@docs
gram(m::MatElem)
```


## Content

```@docs
content(m::MatElem{T}) where T <: RingElement
```


## Minors and exterior powers

```@docs
minors(m::MatElem, k::Int)
minors_with_position(m::MatElem, k::Int)
minors_iterator(m::MatElem, k::Int)
minors_iterator_with_position(m::MatElem, k::Int)
exterior_power(m::MatElem, k::Int)
```


## Pfaffians

```@docs
pfaffian(m::MatElem)
pfaffians(m::MatElem, k::Int)
```
