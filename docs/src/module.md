# Module Interface

AbstractAlgebra allows the construction of finitely presented modules (i.e.
with finitely many generators and relations), starting from free modules.

All module types in AbstractAlgebra follow the following interface.

Free modules can be built over both commutative and noncommutative rings. Other
types of module are restricted to fields and euclidean rings.

## Types and parents

AbstractAlgebra provides two abstract types for modules and their elements:

  * `Module{T}` is the abstract type for module parent types
  * `ModuleElem{T}` is the abstract type for module element types

Note that the abstract types are parameterised. The type `T` should usually be the type
of elements of the ring the module is over.

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

### Basic manipulation

```julia
ngens(M::MyModule{T}) where T <: RingElement
```

Return the number of generators of the module $M$.

```julia
gen(M::MyModule{T}, i::Int) where T <: RingElement
```

Return the $i$-th generator (indexed from $1$) of the module $M$.

```julia
gens(M::MyModule{T}) where T <: RingElement
```

Return an array of the generators of the module $M$.

**Examples**

```julia
M = FreeModule(QQ, 2)

n = ngens(M)
G = gens(M)
g1 = gen(M, 1)
```

### Element constructors

We can construct elements of a module $M$ by specifying linear combinations
of the generators of $M$. This is done by passing a vector of ring elements.

```julia
(M::AbstractAlgebra.Module{T})(v::Vector{T})
```

Construct the element of the module $M$ corrsponding to $\sum_i g[i]*v[i]$
where $g[i]$ are the generators of the module $M$. The resulting element
will lie in the module $M$.

### Arithmetic operators

Elements of a module can be added, subtracted or multiplied by an element of
the ring the module is defined over.

In the case of a noncommutative ring, both left and right scalar multiplication
are defined.

