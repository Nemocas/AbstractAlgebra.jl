```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Finitely presented modules

AbstractAlgebra allows the construction of finitely presented modules (i.e.
with finitely many generators and relations), starting from free modules.

The generic code provided by AbstractAlgebra will only work for modules over
euclidean domains.

Free modules can be built over both commutative and noncommutative rings. Other
types of module are restricted to fields and euclidean rings.

## Abstract types

AbstractAlgebra provides two abstract types for finitely presented modules and
their elements:

  * `FPModule{T}` is the abstract type for finitely presented module parent
types
  * `FPModuleElem{T}` is the abstract type for finitely presented module
element types

Note that the abstract types are parameterised. The type `T` should usually be
the type of elements of the ring the module is over.

## Module functions

All finitely presented modules over a Euclidean domain implement the following
functions.

### Basic functions

```julia
zero(M::FPModule)
```

```julia
iszero(m::FPModuleElem{T}) where T <: RingElement
```

Return `true` if the given module element is zero.

```julia
number_of_generators(M::FPModule{T}) where T <: RingElement
```

Return the number of generators of the module $M$ in its current representation.

```julia
gen(M::FPModule{T}, i::Int) where T <: RingElement
```

Return the $i$-th generator (indexed from $1$) of the module $M$.

```julia
gens(M::FPModule{T}) where T <: RingElement
```

Return a Julia array of the generators of the module $M$.

```julia
rels(M::FPModule{T}) where T <: RingElement
```

Return a Julia vector of all the relations between the generators of `M`. Each
relation is given as an AbstractAlgebra row matrix.

**Examples**

```jldoctest
julia> M = FreeModule(QQ, 2)
Vector space of dimension 2 over rationals

julia> n = number_of_generators(M)
2

julia> G = gens(M)
2-element Vector{AbstractAlgebra.Generic.FreeModuleElem{Rational{BigInt}}}:
 (1//1, 0//1)
 (0//1, 1//1)

julia> R = rels(M)
AbstractAlgebra.Generic.MatSpaceElem{Rational{BigInt}}[]

julia> g1 = gen(M, 1)
(1//1, 0//1)

julia> !iszero(g1)
true

julia> M = FreeModule(QQ, 2)
Vector space of dimension 2 over rationals

julia> z = zero(M)
(0//1, 0//1)

julia> iszero(z)
true
```

### Element constructors

We can construct elements of a module $M$ by specifying linear combinations
of the generators of $M$. This is done by passing a vector of ring elements.

```julia
(M::FPModule{T})(v::Vector{T}) where T <: RingElement
```

Construct the element of the module $M$ corresponding to $\sum_i g[i]v[i]$
where $g[i]$ are the generators of the module $M$. The resulting element
will lie in the module $M$.

### Coercions

Given a module $M$ and an element $n$ of a module $N$, it is possible to
coerce $n$ into $M$ using the notation $M(n)$ in certain circumstances.

In particular the element $n$ will be automatically coerced along any canonical
injection of a submodule map and along any canonical projection of a quotient
map. There must be a path from $N$ to $M$ along such maps.

**Examples**

```julia
F = FreeModule(ZZ, 3)

S1, f = sub(F, [rand(F, -10:10)])

S, g = sub(F, [rand(F, -10:10)])
Q, h = quo(F, S)

m = rand(S1, -10:10)
n = Q(m)
```

### Arithmetic operators

Elements of a module can be added, subtracted or multiplied by an element of
the ring the module is defined over and compared for equality.

In the case of a noncommutative ring, both left and right scalar multiplication
are defined.

### Basic manipulation

```julia
zero(M::FPModule)
```

**Examples**

```jldoctest
julia> M = FreeModule(QQ, 2)
Vector space of dimension 2 over rationals

julia> z = zero(M)
(0//1, 0//1)
```

### Element indexing

```@docs
Base.getindex(m::FPModuleElem{T}) where T <: RingElement
```

**Examples**

```jldoctest
julia> F = FreeModule(ZZ, 3)
Free module of rank 3 over integers

julia> m = F(BigInt[2, -5, 4])
(2, -5, 4)

julia> m[1]
2
```

### Module comparison

```@docs
==(::FPModule{T}, ::FPModule{T}) where T <: RingElement
```

**Examples**

```jldoctest
julia> M = FreeModule(QQ, 2)
Vector space of dimension 2 over rationals

julia> M == M
true

```

### Isomorphism

```@docs
is_isomorphic(::FPModule{T}, ::FPModule{T}) where T <: RingElement
```

!!! note

    Note that this function relies on the Smith normal form over the base ring of
    the modules being able to be made unique. This is true for Euclidean domains
    for which `divrem` has a fixed choice of quotient and remainder, but it will
    not in general be true for Euclidean rings that are not domains.

**Examples**

```jldoctest; setup = :(import Random; Random.seed!(42))
julia> M = FreeModule(ZZ, 3)
Free module of rank 3 over integers

julia> m1 = rand(M, -10:10)
(3, -1, 0)

julia> m2 = rand(M, -10:10)
(4, 4, -7)

julia> S, f = sub(M, [m1, m2])
(Submodule over Integers with 2 generators and no relations, Hom: submodule over Integers with 2 generators and no relations -> free module of rank 3 over integers)

julia> I, g = image(f)
(Submodule over Integers with 2 generators and no relations, Hom: submodule over Integers with 2 generators and no relations -> free module of rank 3 over integers)

julia> is_isomorphic(S, I)
true

```

### Invariant Factor Decomposition

For modules over a euclidean domain one can take the invariant factor
decomposition to determine the structure of the module. The invariant factors
are unique up to multiplication by a unit, and even unique if a 
`canonical_unit` is available for the ring that canonicalises elements.

```@docs
snf(::FPModule{T}) where T <: RingElement
invariant_factors(::FPModule{T}) where T <: RingElement
```

**Examples**

```jldoctest; setup = :(import Random; Random.seed!(42))
julia> M = FreeModule(ZZ, 3)
Free module of rank 3 over integers

julia> m1 = rand(M, -10:10)
(3, -1, 0)

julia> m2 = rand(M, -10:10)
(4, 4, -7)

julia> S, f = sub(M, [m1, m2])
(Submodule over Integers with 2 generators and no relations, Hom: submodule over Integers with 2 generators and no relations -> free module of rank 3 over integers)

julia> Q, g = quo(M, S)
(Quotient module over Integers with 2 generators and relations:
[16 -21], Hom: free module of rank 3 over integers -> quotient module over Integers with 2 generators and relations:
[16 -21])

julia> I, f = snf(Q)
(Invariant factor decomposed module over Integers with invariant factors BigInt[0], Hom: invariant factor decomposed module over Integers with invariant factors BigInt[0] -> quotient module over Integers with 2 generators and relations:
[16 -21])

julia> invs = invariant_factors(Q)
1-element Vector{BigInt}:
 0

```


