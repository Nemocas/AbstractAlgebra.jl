```@meta
CurrentModule = AbstractAlgebra
```

# Real field

AbstractAlgebra.jl provides a module, implemented in `src/julia/Float.jl` for
making Julia `BigFloat`s conform to the AbstractAlgebra.jl Field interface.

In addition to providing a parent object `RealField` for Julia `BigFloat`s, we
implement any additional functionality required by AbstractAlgebra.jl.

Because `BigFloat` cannot be directly included in the AbstractAlgebra.jl
abstract type hierarchy, we achieve integration of Julia `BigFloat`s by
introducing a type union, called `FieldElement`, which is a union of
`AbstractAlgebra.FieldElem` and a number of Julia types, including `BigFloat`.
Everywhere that `FieldElem` is notionally used in AbstractAlgebra.jl, we are in fact
using `FieldElement`, with additional care being taken to avoid ambiguities.

The details of how this is done are technical, and we refer the reader to the
implementation for details. For most intents and purposes, one can think of the Julia
`BigFloat` type as belonging to `AbstractAlgebra.FieldElem`.

## Types and parent objects

Reals have type `BigFloat`, as in Julia itself. We simply supplement the
functionality for this type as required for computer algebra.

The parent objects of such integers has type `Floats{BigFloat}`.

For convenience, we also make `Float64` a part of the AbstractAlgebra.jl type
hierarchy and its parent object (accessible as `RDF`) has type `Floats{Float64}`.

## Rational constructors

In order to construct reals in AbstractAlgebra.jl, one can first construct the
real field itself. This is accomplished using the following constructor.

```julia
Floats{BigFloat}()
```

This gives the unique object of type `Floats{BigFloat}` representing the field of
reals in AbstractAlgebra.jl.

In practice, one simply uses `RealField` which is assigned to be the return value of the
above constructor. There is no need to call the constructor in practice.

Here are some examples of creating the real field and making use of the
resulting parent object to coerce various elements into the field.

**Examples**

```julia
RR = RealField

f = RR()
g = RR(123)
h = RR(BigInt(1234))
k = RR(12//7)
m = RR(2.3)
```

## Basic field functionality

The real field in AbstractAlgebra.jl implements the full Field interface.

We give some examples of such functionality.

**Examples**

```julia
RR = RealField
f = RR(12//7)

h = zero(RR)
k = one(RR)
isone(k) == true
iszero(f) == false
U = base_ring(RR)
V = base_ring(f)
T = parent(f)
f == deepcopy(f)
g = f + 12
m = inv(g)
```
