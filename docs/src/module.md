```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Module Interface

Note: the module infrastructure in AbstractAlgebra should be considered
experimental at this stage. This means that the interface may change in
the future.

AbstractAlgebra allows the construction of finitely presented modules (i.e.
with finitely many generators and relations), starting from free modules. The
generic code provided by AbstractAlgebra will only work for modules over
euclidean domains, however there is nothing preventing a library from
implementing more general modules using the same interface.

All finitely presented module types in AbstractAlgebra follow the following
interface.

Free modules can be built over both commutative and noncommutative rings. Other
types of module are restricted to fields and euclidean rings.

## Types and parents

AbstractAlgebra provides two abstract types for finitely presented modules and
their elements:

  * `FPModule{T}` is the abstract type for finitely presented module parent
types
  * `FPModuleElem{T}` is the abstract type for finitely presented module
element types

Note that the abstract types are parameterised. The type `T` should usually be
the type of elements of the ring the module is over.

## Required functionality for modules

We suppose that `R` is a fictitious base ring and that `S` is a module over `R` with
parent object `S` of type `MyModule{T}`. We also assume the elements in the module have
type `MyModuleElem{T}`, where `T` is the type of elements of the ring the module is
over.

Of course, in practice these types may not be parameterised, but we use parameterised
types here to make the interface clearer.

Note that the type `T` must (transitively) belong to the abstract type `RingElement`
or `NCRingElem`.

We describe the functionality below for modules over commutative rings, i.e. with
element type belonging to `RingElement`, however similar constructors should be
available for element types belonging to `NCRingElem` instead, for free modules over
a noncommutative ring.

Although not part of the module interface, implementations of modules that wish to
follow our interface should use the same function names for submodules, quotient
modules, direct sums and module homomorphisms if they wish to remain compatible
with our module generics in the future.

### Basic manipulation

```julia
iszero(m::MyModuleElem{T}) where T <: RingElement
```

Return `true` if the given module element is zero.

```julia
ngens(M::MyModule{T}) where T <: RingElement
```

Return the number of generators of the module $M$ in its current representation.

```julia
gen(M::MyModule{T}, i::Int) where T <: RingElement
```

Return the $i$-th generator (indexed from $1$) of the module $M$.

```julia
gens(M::MyModule{T}) where T <: RingElement
```

Return a Julia array of the generators of the module $M$.

```julia
rels(M::MyModule{T}) where T <: RingElement
```

Return a Julia vector of all the relations between the generators of `M`. Each
relation is given as an AbstractAlgebra row matrix.

**Examples**

```jldoctest
julia> M = FreeModule(QQ, 2)
Vector space of dimension 2 over Rationals

julia> n = ngens(M)
2

julia> G = gens(M)
2-element Array{AbstractAlgebra.Generic.FreeModuleElem{Rational{BigInt}},1}:
 (1//1, 0//1)
 (0//1, 1//1)

julia> R = rels(M)
0-element Array{AbstractAlgebra.Generic.MatSpaceElem{Rational{BigInt}},1}

julia> g1 = gen(M, 1)
(1//1, 0//1)

julia> !iszero(g1)
true

```

### Element constructors

We can construct elements of a module $M$ by specifying linear combinations
of the generators of $M$. This is done by passing a vector of ring elements.

```julia
(M::AbstractAlgebra.Module{T})(v::Vector{T}) where T <: RingElement
```

Construct the element of the module $M$ corrsponding to $\sum_i g[i]v[i]$
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

## Generic functionality provided

AbstractAlgebra provides the following functionality for all module types
that implement the interface above. Of course, this functionality can also
be provided by special implementations if desired.

### Basic manipulation

```@docs
zero(M::AbstractAlgebra.FPModule{T}) where T <: RingElement
```

**Examples**

```jldoctest
julia> M = FreeModule(QQ, 2)
Vector space of dimension 2 over Rationals

julia> z = zero(M)
(0//1, 0//1)
```

### Element indexing

```@docs
Base.getindex(m::AbstractAlgebra.FPModuleElem{T}) where T <: RingElement
```

**Examples**

```jldoctest
julia> F = FreeModule(ZZ, 3)
Free module of rank 3 over Integers

julia> m = F(BigInt[2, -5, 4])
(2, -5, 4)

julia> m[1]
2
```

### Comparison

```@docs
==(::AbstractAlgebra.FPModule{T}, ::AbstractAlgebra.FPModule{T}) where T <: RingElement
```

**Examples**

```jldoctest
julia> M = FreeModule(QQ, 2)
Vector space of dimension 2 over Rationals

julia> M == M
true

```

### Isomorphism

```@docs
isisomorphic(::AbstractAlgebra.FPModule{T}, ::AbstractAlgebra.FPModule{T}) where T <: RingElement
```

Note that this function relies on the Smith normal form over the base ring of
the modules being able to be made unique. This is true for Euclidean domains
for which `divrem` has a fixed choice of quotient and remainder, but it will
not in general be true for Euclidean rings that are not domains.

**Examples**

```repl
julia> M = FreeModule(ZZ, 3)
Free module of rank 3 over Integers

julia> m1 = rand(M, -10:10)
(0, -8, -8)

julia> m2 = rand(M, -10:10)
(-7, -5, -10)

julia> S, f = sub(M, [m1, m2])
(Submodule over Integers with 2 generators and no relations
, Module homomorphism with
Domain: Submodule over Integers with 2 generators and no relations

Codomain: Free module of rank 3 over Integers)

julia> I, g = image(f)
(Submodule over Integers with 2 generators and no relations
, Module homomorphism with
Domain: Submodule over Integers with 2 generators and no relations

Codomain: Free module of rank 3 over Integers)

julia> isisomorphic(S, I)
true

```

### Invariant Factor Decomposition

For modules over a euclidean domain one can take the invariant factor
decomposition to determine the structure of the module. The invariant factors
are unique up to multiplication by a unit, and even unique if a 
`canonical_unit` is available for the ring that canonicalises elements.

```@docs
snf(::AbstractAlgebra.FPModule{T}) where T <: RingElement
invariant_factors(::AbstractAlgebra.FPModule{T}) where T <: RingElement
```

**Examples**

```repl
julia> M = FreeModule(ZZ, 3)
Free module of rank 3 over Integers

julia> m1 = rand(M, -10:10)
(9, 7, 7)

julia> m2 = rand(M, -10:10)
(-6, 2, -8)

julia> S, f = sub(M, [m1, m2])
(Submodule over Integers with 2 generators and no relations
, Module homomorphism with
Domain: Submodule over Integers with 2 generators and no relations

Codomain: Free module of rank 3 over Integers)

julia> Q, g = quo(M, S)
(Quotient module over Integers with 3 generators and relations:
[3 9 -1], [0 20 -10], Module homomorphism with
Domain: Free module of rank 3 over Integers
Codomain: Quotient module over Integers with 3 generators and relations:
[3 9 -1], [0 20 -10])

julia> I, f = snf(Q)
(Invariant factor decomposed module over Integers with invariant factors BigInt[10, 0], Module homomorphism with
Domain: Invariant factor decomposed module over Integers with invariant factors BigInt[10, 0]
Codomain: Quotient module over Integers with 3 generators and relations:
[3 9 -1], [0 20 -10])

julia> invs = invariant_factors(Q)
2-element Array{BigInt,1}:
 10
  0

```
