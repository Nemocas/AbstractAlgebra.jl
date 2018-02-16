```@meta
CurrentModule = AbstractAlgebra
```

# Rational field

AbstractAlgebra.jl provides a module, implemented in `src/julia/Rational.jl` for
making Julia `Rational{BigInt}`s conform to the AbstractAlgebra.jl Field interface.

In addition to providing a parent object `QQ` for Julia `Rational{BigInt}`s, we
implement any additional functionality required by AbstractAlgebra.jl.

Because `Rational{BigInt}` cannot be directly included in the AbstractAlgebra.jl
abstract type hierarchy, we achieve integration of Julia `Rational{BigInt}`s by
introducing a type union, called `FieldElement`, which is a union of
`AbstractAlgebra.FieldElem` and a number of Julia types, including `Rational{BigInt}`.
Everywhere that `FieldElem` is notionally used in AbstractAlgebra.jl, we are in fact
using `FieldElement`, with additional care being taken to avoid ambiguities.

The details of how this is done are technical, and we refer the reader to the
implementation for details. For most intents and purposes, one can think of the Julia
`Rational{BigInt}` type as belonging to `AbstractAlgebra.FieldElem`.

One other technicality is that Julia defines certain functions for `Rational{BigInt}`,
such as `sqrt` and `exp` differently to what AbstractAlgebra.jl requires. To get around
this, we redefine these functions internally to AbstractAlgebra.jl, without redefining
them for users of AbstractAlgebra.jl. This allows the internals of AbstractAlgebra.jl
to function correctly, without broadcasting pirate definitions of already defined Julia
functions to the world.

To access the internal definitions, one can use `AbstractAlgebra.sqrt` and
`AbstractAlgebra.exp`, etc.

## Types and parent objects

Rationals have type `Rational{BigInt}`, as in Julia itself. We simply supplement the
functionality for this type as required for computer algebra.

The parent objects of such integers has type `Rationals{BigInt}`.

For convenience, we also make `Rational{Int}` a part of the AbstractAlgebra.jl type
hierarchy and its parent object (accessible as `qq`) has type `Rationals{Int}`. But we
caution that this type is not particularly useful as a model of the rationals and may
not function as expected within AbstractAlgebra.jl.

## Rational constructors

In order to construct rationals in AbstractAlgebra.jl, one can first construct the
rational field itself. This is accomplished using either of the following constructors.

```julia
FractionField(R::Integers{BigInt})
```

```julia
Rationals{BigInt}()
```

This gives the unique object of type `Rationals{BigInt}` representing the field of
rationals in AbstractAlgebra.jl.

In practice, one simply uses `QQ` which is assigned to be the return value of the
above constructor. There is no need to call the constructor in practice.

Here are some examples of creating the rational field and making use of the
resulting parent object to coerce various elements into the field.

**Examples**

```julia
f = QQ()
g = QQ(123)
h = QQ(BigInt(1234))
k = QQ(BigInt(12), BigInt(7))

QQ == FractionField(ZZ)
```

## Basic field functionality

The rational field in AbstractAlgebra.jl implements the full Field and Fraction Field
interfaces.

We give some examples of such functionality.

**Examples**

```julia
f = QQ(12, 7)

h = zero(QQ)
k = one(QQ)
isone(k) == true
iszero(f) == false
U = base_ring(QQ)
V = base_ring(f)
T = parent(f)
f == deepcopy(f)
g = f + 12
r = ZZ(12)//ZZ(7)
n = numerator(r)
```

## Rational functionality provided by AbstractAlgebra.jl

The functionality below supplements that provided by Julia itself for its
`Rational{BigInt}` type.

### Square root

```@docs
AbstractAlgebra.sqrt(a::Rational{BigInt})
```

```@docs
AbstractAlgebra.exp(a::Rational{BigInt})
```

**Examples**

```julia
d = AbstractAlgebra.sqrt(ZZ(36)//ZZ(25))
m = AbstractAlgebra.exp(ZZ(0)//ZZ(1))
```


