```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
```

AbstractAlgebra.jl provides rudimentary native support for permutation groups (implemented in `src/generic/PermGroups.jl`) and their characters through Young Tableaux (`src/generic/YoungTabs.jl`).

# Permutations and Permutation groups

Permutations in AbstractAlgebra.jl are represented internally via vector (list) of integers, wrapped in type `perm{T}`, where `T<:Integer` carries the information on the type of elements of a permutation. Permutation groups are singleton parent object of type `PermGroup{T}` and are used mostly to store the length of a permutation, since it is not included in the permutation type.

Permutation groups are created using the `PermGroup` (inner) constructor.
However, for convenience we define
```
PermutationGroup = PermGroup
```
so that permutation groups can be created using `PermutationGroup` instead of `PermGroup`.

The types of permutations in AbstractAlgebra.jl are given by the following table, along with
the libraries that provide them and the associated types of the parent objects.

 Library | Group                          | Element type  | Parent type
---------|--------------------------------|---------------|---------------------
Native   | $S_n$                          | `perm`        | `PermGroup`

Both `PermGroup` and `perm` and can be parametrized by any type `T<:Integer` .
By default the parameter is the `Int`-type native to the systems architecture.
However, if You are sure that Your permutations are small enough to fit into smaller integer type (such as `Int32`, `Uint16`, or even `Int8`), You may choose to change the paramtrizing type accordingly.
In practice this may result in decreased memory footprint (when storing multiple permutations) and noticable faster performance, if Your workload is heavy in operations on permutations, which e.g. does not fit into cache of Your cpu.

All the permutation group types belong to the `Group` abstract type and the corresponding permutation element types belong to the `GroupElem` abstract type.

```@docs
setpermstyle
```

## Permutations constructors

There are several methods to to construct permutations in AbstractAlgebra.jl.

* Since the parent object can be reconstructed from the permutation object, the easiest way is to directly call to `perm` constructor:

```@docs
Generic.perm
```

* The other way is to first construct the permutation group they belong to.
  This is accomplished with the inner constructor `PermGroup(n::Integer)` which
  constructs the permutation group on $n$ symbols and returns the parent object
  representing the group.

```@docs
Generic.PermGroup
```

  A vector of integers can be then coerced to a permutation via call to parent. The advantage is that the vector is automatically converted to the integer type fixed at the creation of the parent object.

  **Examples:**

```jldoctest
julia> G = PermutationGroup(BigInt(5)); p = G([2,3,1,5,4])
(1,2,3)(4,5)

julia> typeof(p)
AbstractAlgebra.Generic.perm{BigInt}

julia> H = PermutationGroup(UInt16(5)); r = H([2,3,1,5,4])
(1,2,3)(4,5)

julia> typeof(r)
AbstractAlgebra.Generic.perm{UInt16}

julia> H()
()
```

  By default the coercion checks for non-unique values in the vector, but this can be switched off with `G([2,3,1,5,4], false)`.

* Finally there is a `perm"..."` string macro to construct permutation from string input.

```@docs
@perm_str
```

## Permutation interface

The following basic functionality is provided by the default permutation group
implementation in AbstractAlgebra.jl, to support construction of other generic
constructions over permutation groups.
Any custom permutation group implementation in AbstractAlgebra.jl should provide these functions along with the usual group element arithmetic and comparison.

```@docs
parent(::perm)
elem_type(::PermGroup)
parent_type(::perm)
```

A custom implementation also needs to implement `hash(::perm, ::UInt)` and (possibly) `deepcopy_internal(::perm, ::ObjectIdDict)`.
Note that permutation group elements are mutable and so returning shallow copies is not sufficient.

```julia
getindex(a::perm, n::Int)
```

Allows access to entry $n$ of the given permutation via the syntax `a[n]`.
Note that entries are $1$-indexed.

```julia
setindex!(a::perm, d::Int, n::Int)
```

Set the $n$-th entry of the given permutation to $d$.
This allows Julia to provide the syntax `a[n] = d` for setting entries of a permuation.
Note that entries are $1$-indexed.

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
G(p::perm)
```

Take a permutation that is already in the permutation group and simply return
it. A copy of the original is not made if not necessary.

## Basic manipulation

Numerous functions are provided to manipulate permutation group elements.

```@docs
cycles(::perm)
```

Cycle structure is stored in a permutation, since once available, it provides a convenient shortcut in many other algorithms.

**WARNING:** You should never modify a permutation (e.g. via `setindex`) for which cycle structure has been already computed.
Cycle decomposition is computed once for the lifetime of a permutation.

```@docs
parity(::perm)
sign(::perm)
permtype(::perm)
order(::perm)
order(::PermGroup)
```

Note that even an `Int64` can be easily overflowd when computing with permutation groups.
Thus, by default, `order` returns (always correct) `BigInt`s.
If You are sure that the computation will not overflow, You may use `order(::Type{T}, ...)` to perform computations with machine integers.
Julias standard promotion rules apply for the returned value.

Iteration over all permutations in the permutation group $S_n$ can be achieved with

```@docs
elements(::PermGroup)
```
Iteration in reasonable time (i.e. in terms of minutes) is possible for $S_n$ when $n ≤ 13$.
You may also use the non-allocating `Generic.elements!(::PermGroup)` for $n = 14$ (or even $15$ if You are patient enough), which is an order of mangitude faster, since all permutations are generated "in-place":

```jldoctest
julia> collect(elements(PermGroup(3)))
6-element Array{AbstractAlgebra.Generic.perm{Int64},1}:
 ()
 (1,2)
 (1,3,2)
 (2,3)
 (1,2,3)
 (1,3)

julia> for p in Generic.elements!(PermGroup(3)); println(p.d); end
[1, 2, 3]
[2, 1, 3]
[3, 1, 2]
[1, 3, 2]
[2, 3, 1]
[3, 2, 1]
```

However, `collect(Generic.elements!(PermGroup(n)))` returns a vector of identical permutations:

```jldoctest
julia> A = collect(Generic.elements!(PermGroup(3))); length(A)
6

julia> unique(A)
1-element Array{AbstractAlgebra.Generic.perm{Int64},1}:
 (1,3)
```
If You intend to use or store elements yielded by `elements!` You need to **copy** them explicitely.
Note be careful not to use `cycles` and `elements!`, as the cycle structue will be computed only for the first element.
```julia
for p in Generic.elements!(PermGroup(3))
   cycles(p)
end
```
will compute cycle decomposition exactly **once**.

## Arithmetic operators

```@docs
*(::perm, ::perm)
^(::perm, n::Integer)
inv(::perm)
```

Permutations parametrized by different types can be multiplied, and follow the standard julia integer promotion rules:

```jldoctest
g = rand(PermGroup(Int8(5)));
h = rand(PermGroup(UInt32(5)));
typeof(g*h)

# output
AbstractAlgebra.Generic.perm{Int64}
```

## Coercion

The following coercions are available for `G::PermGroup` parent objects.
Each of the method performs basic sanity checks on the input which can be switched off by the second argument.

```julia
(G::PermGroup)()
```
returns the identity element of `G`.


```julia
(G::PermGrup)(::Vector{<:Integer}[, check=true])
```
turns a vector od integers into a permutation (performing conversion, if necessary).


```julia
(G::PermGroup)(::perm{<:Integer}[, check=true])
```
coerces a permutation `p` into group $G$ (performing the conversion, if necessary).
If `p` is already an element of `G` no copy is performed.

```julia
(G::PermGroup)(::String[, check=true])
```
parses the string input e.g. copied from the output of GAP.
The method is the same as `perm"..."` macro.
The string is sanitized and checked for disjoint cycles.
Both `string(p)` (if `setpermstyle(:cycles)`) and `string(cycles(p))` are valid input for this method.

```julia
(G::PermGroup{T})(::CycleDec{T}[, check=true]) where T
```
turns a cycle decomposition object into a permutation.

## Comparison

```@docs
==(::perm, ::perm)
==(::PermGroup, ::PermGroup)
```

## Misc
```@docs
rand(::PermGroup)
Generic.matrix_repr(::perm)
Generic.emb(::PermGroup, ::Vector{Int}, ::Bool)
Generic.emb!(::perm, ::perm, V)
```

# Characters
