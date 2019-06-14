```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Integer ring

AbstractAlgebra.jl provides a module, implemented in `src/julia/Integer.jl` for
making Julia `BigInt`s conform to the AbstractAlgebra.jl Ring interface.

In addition to providing a parent object `ZZ` for Julia `BigInt`s, we implement
any additional functionality required by AbstractAlgebra.jl.

Because `BigInt` cannot be directly included in the AbstractAlgebra.jl abstract type
hierarchy, we achieve integration of Julia `BigInt`s by introducing a type union, called
`RingElement`, which is a union of `AbstractAlgebra.RingElem` and a number of Julia
types, including `BigInt`. Everywhere that `RingElem` is notionally used in
AbstractAlgebra.jl, we are in fact using `RingElement`, with additional care being taken
to avoid ambiguities.

The details of how this is done are technical, and we refer the reader to the
implementation for details. For most intents and purposes, one can think of the Julia
`BigInt` type as belonging to `AbstractAlgebra.RingElem`.

One other technicality is that Julia defines certain functions for `BigInt`, such as
`sqrt` and `exp` differently to what AbstractAlgebra.jl requires. To get around this,
we redefine these functions internally to AbstractAlgebra.jl, without redefining them
for users of AbstractAlgebra.jl. This allows the internals of AbstractAlgebra.jl to
function correctly, without broadcasting pirate definitions of already defined Julia
functions to the world.

To access the internal definitions, one can use `AbstractAlgebra.sqrt` and
`AbstractAlgebra.exp`, etc.

## Types and parent objects

Integers have type `BigInt`, as in Julia itself. We simply supplement the functionality
for this type as required for computer algebra.

The parent objects of such integers has type `Integers{BigInt}`.

For convenience, we also make `Int` a part of the AbstractAlgebra.jl type hierarchy
and its parent object (accessible as `zz`) has type `Integers{Int}`. But we caution
that this type is not particularly useful as a model of the integers and may not
function as expected within AbstractAlgebra.jl.

## Integer constructors

In order to construct integers in AbstractAlgebra.jl, one can first construct the
integer ring itself. This is accomplished using the following constructor.

```julia
Integers{BigInt}()
```

This gives the unique object of type `Integers{BigInt}` representing the ring of
integers in AbstractAlgebra.jl.

In practice, one simply uses `ZZ` which is assigned to be the return value of the
above constructor. There is no need to call the constructor in practice.

Here are some examples of creating the integer ring and making use of the
resulting parent object to coerce various elements into the ring.

**Examples**

```jldoctest
julia> f = ZZ()
0

julia> g = ZZ(123)
123

julia> h = ZZ(BigInt(1234))
1234

```

## Basic ring functionality

The integer ring in AbstractAlgebra.jl implements the full Ring interface and the 
Euclidean Ring interface.

We give some examples of such functionality.

**Examples**

```jldoctest
julia> f = ZZ(12)
12

julia> h = zero(ZZ)
0

julia> k = one(ZZ)
1

julia> isone(k)
true

julia> iszero(f)
false

julia> T = parent(f)
Integers

julia> f == deepcopy(f)
true

julia> g = f + 12
24

julia> h = powmod(f, 12, ZZ(17))
4

julia> flag, q = divides(f, ZZ(3))
(true, 4)

```

## Integer functionality provided by AbstractAlgebra.jl

The functionality below supplements that provided by Julia itself for its `BigInt` type.

### Basic functionality

```@docs
isunit(::Integer)
```

**Examples**

```jldoctest
julia> r = ZZ(-1)
-1

julia> isunit(r)
true

```

### Square root

```@docs
AbstractAlgebra.sqrt(a::BigInt)
```

```@docs
AbstractAlgebra.exp(a::BigInt)
```

**Examples**

```jldoctest
julia> d = AbstractAlgebra.sqrt(ZZ(36))
6

julia> m = AbstractAlgebra.exp(ZZ(0))
1

```
### Coprime bases

```@docs
ppio(a::BigInt, b::BigInt)
```

**Examples**

```jldoctest
julia> c, n = ppio(ZZ(12), ZZ(26))
(4, 3)

```


