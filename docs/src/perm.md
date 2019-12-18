```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

# Permutations and Permutation groups

AbstractAlgebra.jl provides rudimentary native support for permutation groups (implemented in `src/generic/PermGroups.jl`). All functionality of permutations is accesible in the `Generic` submodule.

Permutations are represented internally via vector of integers, wrapped in type `Perm{T}`, where `T<:Integer` carries the information on the type of elements of a permutation. Permutation groups are singleton parent objects of type `SymmetricGroup{T}` and are used mostly to store the length of a permutation, since it is not included in the permutation type.

Permutation groups are created using the `SymmetricGroup` (inner) constructor.

Both `SymmetricGroup` and `Perm` and can be parametrized by any type `T<:Integer` .
By default the parameter is the `Int`-type native to the systems architecture.
However, if you are sure that your permutations are small enough to fit into smaller integer type (such as `Int32`, `UInt16`, or even `Int8`), you may choose to change the parametrizing type accordingly.
In practice this may result in decreased memory footprint (when storing multiple permutations) and noticable faster performance, if your workload is heavy in operations on permutations, which e.g. does not fit into cache of your cpu.

All the permutation group types belong to the `Group` abstract type and the corresponding permutation element types belong to the `GroupElem` abstract type.

```@docs
Generic.setpermstyle
```

## Permutations constructors

There are several methods to construct permutations in AbstractAlgebra.jl.

* The easiest way is to directly call to the `Perm` (inner) constructor:

```@docs
Generic.Perm
```

  Since the parent object can be reconstructed from the permutation itself, you can work with permutations without explicitly constructing the parent object.

* The other way is to first construct the permutation group they belong to.
  This is accomplished with the inner constructor `SymmetricGroup(n::Integer)` which
  constructs the permutation group on $n$ symbols and returns the parent object
  representing the group.

```@docs
Generic.SymmetricGroup
```

  A vector of integers can be then coerced to a permutation by calling a parent permutation group on it.
  The advantage is that the vector is automatically converted to the integer type fixed at the creation of the parent object.

  **Examples:**

```jldoctest
julia> G = SymmetricGroup(BigInt(5)); p = G([2,3,1,5,4])
(1,2,3)(4,5)

julia> typeof(p)
Perm{BigInt}

julia> H = SymmetricGroup(UInt16(5)); r = H([2,3,1,5,4])
(1,2,3)(4,5)

julia> typeof(r)
Perm{UInt16}

julia> H()
()
```

  By default the coercion checks for non-unique values in the vector, but this can be switched off with `G([2,3,1,5,4], false)`.

* Finally there is a `perm"..."` string macro to construct a permutation from a string input.

```@docs
@perm_str
```

## Permutation interface

The following basic functionality is provided by the default permutation group
implementation in AbstractAlgebra.jl, to support construction of other generic
constructions over permutation groups.
Any custom permutation group implementation in AbstractAlgebra.jl should provide these functions along with the usual group element arithmetic and comparison.

```@docs
parent(::Perm)
elem_type(::SymmetricGroup)
parent_type(::Perm)
```

A custom implementation also needs to implement `hash(::Perm, ::UInt)` and (possibly) `deepcopy_internal(::Perm, ::ObjectIdDict)`.

!!! note

    Permutation group elements are mutable and so returning shallow copies is not sufficient.

```julia
getindex(a::Perm, n::Integer)
```

Allow access to entry $n$ of the given permutation via the syntax `a[n]`.
Note that entries are $1$-indexed.

```julia
setindex!(a::Perm, d::Integer, n::Integer)
```

Set the $n$-th entry of the given permutation to $d$.
This allows Julia to provide the syntax `a[n] = d` for setting entries of a permutation. Entries are $1$-indexed.

!!! note
    Using `setindex!` invalidates the cycle decomposition cached in a permutation, which will be computed the next time it is needed.

Given the parent object `G` for a permutation group, the following coercion functions are provided to coerce various arguments into the permutation group.
Developers provide these by overloading the permutation group parent objects.

```julia
G()
```

Return the identity permutation.

```julia
G(A::Vector{<:Integer})
```

Return the permutation whose entries are given by the elements of the supplied
vector.

```julia
G(p::Perm)
```

Take a permutation that is already in the permutation group and simply return
it. A copy of the original is not made if not necessary.

## Basic manipulation

Numerous functions are provided to manipulate permutation group elements.

```@docs
cycles(::Perm)
```

Cycle structure is cached in a permutation, since once available, it provides a convenient shortcut in many other algorithms.

```@docs
parity(::Perm)
sign(::Perm)
permtype(::Perm)
order(::Perm)
order(::Generic.SymmetricGroup)
```

Note that even an `Int64` can be easily overflowed when computing with permutation groups.
Thus, by default, `order` returns (always correct) `BigInt`s.
If you are sure that the computation will not overflow, you may use `order(::Type{T}, ...)` to perform computations with machine integers.
Julia's standard promotion rules apply for the returned value.

Since `SymmetricGroup` implements the iterator protocol, you may iterate over all permutations via a simple loop:

```
for p in SymmetricGroup(n)
   ...
end
```
Iteration over all permutations in reasonable time, (i.e. in terms of minutes) is possible when $n ≤ 13$.

You may also use the non-allocating `Generic.elements!` function for $n ≤ 14$ (or even $15$ if you are patient enough), which is an order of magnitude faster.

```@docs
Generic.elements!(::Generic.SymmetricGroup)
```

However, since all permutations yielded by `elements!` are aliased (modified "in-place"), `collect(Generic.elements!(SymmetricGroup(n)))` returns a vector of identical permutations.

!!! note
    If you intend to use or store elements yielded by `elements!` you need to **deepcopy** them explicitly.

## Arithmetic operators

```@docs
*(::Perm{T}, ::Perm{T}) where T
^(::Perm, n::Integer)
inv(::Perm)
```

Permutations parametrized by different types can be multiplied, and follow the standard julia integer promotion rules:

```jldoctest
g = rand(SymmetricGroup(Int8(5)));
h = rand(SymmetricGroup(UInt32(5)));
typeof(g*h)

# output
Perm{UInt32}
```

## Coercion

The following coercions are available for `G::SymmetricGroup` parent objects.
Each of the methods perform basic sanity checks on the input which can be switched off by the second argument.

**Examples**

```jldoctest
julia> SymmetricGroup(4)()
()
```
> Return the identity element of `G`.


```julia
(G::SymmetricGroup)(::AbstractVector{<:Integer}[, check=true])
```
> Turn a vector of integers into a permutation (performing conversion, if necessary).


```julia
(G::SymmetricGroup)(::Perm[, check=true])
```
> Coerce a permutation `p` into group `G` (performing the conversion, if necessary).
> If `p` is already an element of `G` no copy is performed.

```julia
(G::SymmetricGroup)(::String[, check=true])
```
> Parse the string input e.g. copied from the output of GAP.
> The method uses the same logic as the `perm"..."` macro.
> The string is sanitized and checked for disjoint cycles.
> Both `string(p::Perm)` (if `setpermstyle(:cycles)`) and `string(cycles(p::Perm))` are valid input for this method.

```julia
(G::SymmetricGroup{T})(::CycleDec{T}[, check=true]) where T
```
> Turn a cycle decomposition object into a permutation.

## Comparison

```@docs
==(::Perm, ::Perm)
==(::Generic.SymmetricGroup, ::Generic.SymmetricGroup)
```

## Misc
```@docs
rand(::Generic.SymmetricGroup)
Generic.matrix_repr(::Perm)
Generic.emb(::Generic.SymmetricGroup, ::Vector{Int}, ::Bool)
Generic.emb!(::Perm, ::Perm, V)
```
