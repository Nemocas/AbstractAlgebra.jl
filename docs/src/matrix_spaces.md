```@meta
CurrentModule = AbstractAlgebra
CollapsedDocStrings = true
DocTestSetup = AbstractAlgebra.doctestsetup()
```

# Matrix spaces

A matrix space represents the collection of all matrices with a fixed number
of rows and columns over a fixed base ring. Matrix spaces are parent objects;
their elements are the corresponding algebraic matrices.


## Creating matrix spaces

```@docs
matrix_space(R::NCRing, r::Int, c::Int)
```


## Properties of matrix spaces

```@docs
number_of_rows(s::MatSpace)
number_of_columns(s::MatSpace)
```

The base ring and the vector space dimension of the matrix space can be queried
with `base_space` and `vector_space_dim`, respectively.


## Creating elements of a matrix space

### Calling a matrix space

```@docs
MatSpace
```

### Special elements

```@docs
zero(s::MatSpace)
one(s::MatSpace)
```

### Recovering the parent object

```@docs
parent(M::MatElem)
```
