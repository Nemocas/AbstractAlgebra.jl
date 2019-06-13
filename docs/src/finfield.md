```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Finite fields

AbstractAlgebra.jl provides a module, implemented in `src/julia/GF.jl` for
finite fields. The module is a naive implementation that supports only fields of degree
$1$ (prime fields). They are modelled as $\mathbb{Z}/p\mathbb{Z}$ for $p$ a prime.

## Types and parent objects

Finite fields have type `GFField{T}` where `T` is either `Int` or `BigInt`.

Elements of such a finite field have type `gfelem{T}`.

## Finite field constructors

In order to construct finite fields in AbstractAlgebra.jl, one must first construct the
field itself. This is accomplished with the following constructors.

```@docs
GF(p::T) where T <: Integer
```

Here are some examples of creating a finite field and making use of the resulting
parent object to coerce various elements into the field.

**Examples**

```jldoctest
julia> F = GF(13)

julia> g = F(3)

julia> h = F(g)

```

## Basic field functionality

The finite field module in AbstractAlgebra.jl implements the full Field interface.

We give some examples of such functionality.

**Examples**

```jldoctest
julia> F = GF(13)

julia> h = zero(F)

julia> k = one(F)

julia> isone(k) == true

julia> iszero(f) == false

julia> U = base_ring(F)

julia> V = base_ring(h)

julia> T = parent(h)

julia> h == deepcopy(h)

julia> h = h + 2

julia> m = inv(k)

```

## Basic manipulation of fields and elements

```@docs
gen{T <: Integer}(F::GFField{T})
```

```@docs
order(F::GFField)
```

```@docs
degree(F::GFField)
```

**Examples**

```jldoctest
julia> F = GF(13)

julia> d = degree(F)

julia> n = order(F)

julia> g = gen(F)

```


