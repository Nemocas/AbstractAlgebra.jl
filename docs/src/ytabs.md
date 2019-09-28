```@meta
CurrentModule = AbstractAlgebra
DocTestSetup = quote
    using AbstractAlgebra
end
DocTestFilters = r"[0-9\.]+ seconds \(.*\)"
```

# Partitions and Young tableaux

AbstractAlgebra.jl provides basic support for computations with Young tableaux, skew diagrams and the characters of permutation groups (implemented `src/generic/YoungTabs.jl`).
All functionality of permutations is accesible in the `Generic` submodule.

## Partitions

The basic underlying object for those concepts is `Partition` of a number $n$, i.e. a sequence of positive integers $n_1, \ldots, n_k$ which sum to $n$.
Partitions in AbstractAlgebra.jl are represented internally by non-increasing `Vector`s of `Int`s.
Partitions are printed using the standard notation, i.e. $9 = 4 + 2 + 1 + 1 + 1$ is shown as $4_1 2_1 1_3$ with the subscript indicating the count of a summand in the partition.

```@docs
Generic.Partition
```

### Array interface

`Partition` is a concrete subtype of `AbstractVector{Int}` and implements the following standard Array interface:

```@docs
size(::Generic.Partition)
getindex(::Generic.Partition, i::Integer)
setindex!(::Generic.Partition, v::Integer, i::Integer)
```
These functions work on the level of `p.part` vector.
Additionally `setindex!` will try to prevent uses which result in non-valid (i.e. non-decreasing) partition vectors.

One can easily iterate over all partitions of $n$ using the `AllParts` type:

```@docs
Generic.AllParts
```

The number of all partitions can be computed by the hidden function `_numpart`.
Much faster implementation is available in [Nemo.jl](http://nemocas.github.io/Nemo.jl/latest/arb.html#Nemo.numpart-Tuple{Int64,ArbField}).

```@docs
Generic._numpart
```

Since `Partition` is a subtype of `AbstractVector` generic functions which operate on vectors should work in general.
However the meaning of `conj` has been changed to agree with the traditional understanding of conjugation of `Partitions`:
```@docs
conj(::Generic.Partition)
conj(::Generic.Partition, v::Vector)
```

## Young Diagrams and Young Tableaux

Mathematicaly speaking Young diagram is a diagram which consists of rows of square boxes such that the number of boxes in each row is no less than the number of boxes in the previous row.
For example partition $4_1 3_2 1$ represents the following diagram.
```
┌───┬───┬───┬───┐
│   │   │   │   │
├───┼───┼───┼───┘
│   │   │   │
├───┼───┼───┤
│   │   │   │
├───┼───┴───┘
│   │
└───┘
```
Young Tableau is formally a bijection between the set of boxes of a Young Diagram and the set $\{1, \ldots, n\}$.
If a bijection is increasing along rows and columns of the diagram it is referred to as **standard**.
For example
```
┌───┬───┬───┬───┐
│ 1 │ 2 │ 3 │ 4 │
├───┼───┼───┼───┘
│ 5 │ 6 │ 7 │
├───┼───┼───┤
│ 8 │ 9 │10 │
├───┼───┴───┘
│11 │
└───┘
```
is a standard Young tableau of $4_1 3_2 1$ where the bijection assigns consecutive natural numbers to consecutive (row-major) cells.

### Constructors

In AbstractAlgebra.jl Young tableau are implemented as essentially row-major sparse matrices, i.e. `YoungTableau <: AbstractArray{Int,2}` but only the defining `Partition` and the (row-major) fill-vector is stored.

```@docs
Generic.YoungTableau
```
For convenience there exists an alternative constructor of `YoungTableau`, which accepts a vector of integers and constructs `Partition` internally.
```
YoungTableau(p::Vector{Integer}[, fill=collect(1:sum(p))])
```

### Array interface

To make `YoungTableaux` array-like we implement the following functions:

```@docs
size(::Generic.YoungTableau)
getindex(::Generic.YoungTableau, n::Integer)
```
Also the double-indexing corresponds to `(row, column)` access to an abstract array.
```jldoctest
julia> y = YoungTableau([4,3,1])
┌───┬───┬───┬───┐
│ 1 │ 2 │ 3 │ 4 │
├───┼───┼───┼───┘
│ 5 │ 6 │ 7 │
├───┼───┴───┘
│ 8 │
└───┘

julia> y[1,2]
2

julia> y[2,3]
7

julia> y[3,2]
0
```

Functions defined for `AbstractArray` type based on those (e.g. `length`) should work.
Again, as in the case of `Partition` the meaning of `conj` is altered to
reflect the usual meaning for Young tableaux:
```@docs
conj(::Generic.YoungTableau)
```

### Pretty-printing

Similarly to permutations we have two methods of displaying Young Diagrams:

```@docs
Generic.setyoungtabstyle
```

### Ulitility functions
```@docs
matrix_repr(::Generic.YoungTableau)
fill!(::Generic.YoungTableau, ::AbstractVector{<:Integer})
```

## Characters of permutation groups

Irreducible characters (at least over field of characteristic $0$) of the full group of permutations $S_n$ correspond via [Specht modules](https://en.wikipedia.org/wiki/Specht_module) to partitions of $n$.

```@docs
character(::Generic.Partition)
character(lambda::Generic.Partition, p::Generic.Perm)
character(lambda::Generic.Partition, mu::Generic.Partition)
```
The values computed by characters are cached in an internal dictionary `Dict{Tuple{BitVector,Vector{Int}}, BigInt}`.
Note that all of the above functions return `BigInts`.
If you are sure that the computations do not overflow, variants of the last two functions using `Int` are available:
```
character(::Type{Int}, lambda::Partition, p::Perm[, check::Bool=true])
character(::Type{Int}, lambda::Partition, mu::Partition[, check::Bool=true])
```

The dimension $\dim \lambda$ of the irreducible module corresponding to partition $\lambda$ can be computed using [Hook length formula](https://en.wikipedia.org/wiki/Hook_length_formula)

```@docs
Generic.rowlength
Generic.collength
hooklength
dim(::Generic.YoungTableau)
```

The character associated with `Y.part` can also be used to compute the dimension, but as it is expected the Murnaghan-Nakayama is much slower even though (due to caching) consecutive calls are fast:

```jldoctest
julia> λ = Partition(collect(12:-1:1))
12₁11₁10₁9₁8₁7₁6₁5₁4₁3₁2₁1₁

julia> @time dim(YoungTableau(λ))
  0.224430 seconds (155.77 k allocations: 7.990 MiB)
9079590132732747656880081324531330222983622187548672000

julia> @time dim(YoungTableau(λ))
  0.000038 seconds (335 allocations: 10.734 KiB)
9079590132732747656880081324531330222983622187548672000

julia> G = PermutationGroup(sum(λ))
Permutation group over 78 elements

julia> @time character(λ, G())
 24.154105 seconds (58.13 M allocations: 3.909 GiB, 42.84% gc time)
9079590132732747656880081324531330222983622187548672000

julia> @time character(λ, G())
  0.001439 seconds (195 allocations: 24.453 KiB)
9079590132732747656880081324531330222983622187548672000
```

### Low-level functions and characters

As mentioned above `character` functions use the Murnaghan-Nakayama rule for evaluation.
The implementation follows
> Dan Bernstein,
> The computational complexity of rules for the character table of $S_n$
> *Journal of Symbolic Computation*, **37** (6), 2004, p. 727-748,
implementing the following functions. For precise definitions and meaning please consult the paper cited.

```@docs
Generic.partitionseq
isrimhook(::BitVector, ::Int, ::Int)
Generic.MN1inner
```

## Skew Diagrams

Skew diagrams are formally differences of two Young diagrams. Given $\lambda$ and $\mu$, two partitions of $n+m$ and $m$ (respectively). Suppose that each of cells of $\mu$ is a cell of $\lambda$ (i.e. parts of $\mu$ are no greater than the corresponding parts of $\lambda$). Then the skew diagram denoted by $\lambda/\mu$ is the set theoretic difference the of sets of boxes, i.e. is a diagram with exactly $n$ boxes:

```@docs
Generic.SkewDiagram
```

`SkewDiagram` implements array interface with the following functions:

```@docs
size(xi::Generic.SkewDiagram)
in(t::Tuple{T,T}, xi::Generic.SkewDiagram) where T<:Integer
getindex(xi::Generic.SkewDiagram, n::Integer)
```

The support for skew diagrams is very rudimentary. The following functions are available:

```@docs
isrimhook(::Generic.SkewDiagram)
leglength
matrix_repr(::Generic.SkewDiagram)
```
