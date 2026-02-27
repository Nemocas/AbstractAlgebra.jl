###############################################################################
#
#   Cycle Decomposition
#
###############################################################################

@doc raw"""
    CycleDec{T}(ccycles, cptrs, n) where T

Cycle decomposition of a permutation.
* `ccycles`: an array of consecutive entries of cycles;
* `cptrs`: an array of pointers to the locations where cycles begin: ```ccycles[cptrs[i], cptrs[i+1]-1]` contains the i-th cycle;
* `n`: the number of cycles;
"""
struct CycleDec{T<:Integer}
   ccycles::Vector{T}
   cptrs::Vector{T}
   n::T
end

@doc raw"""
    Perm{T<:Integer}

The type of permutations.
Fieldnames:
* `d::Vector{T}` - vector representing the permutation
* `modified::Bool` - bit to check the validity of cycle decomposition
* `cycles::CycleDec{T}` - (cached) cycle decomposition

A permutation $p$ consists of a vector (`p.d`) of $n$ integers from $1$ to $n$.
If the $i$-th entry of the vector is $j$, this corresponds to $p$ sending $i \to j$.
The cycle decomposition (`p.cycles`) is computed on demand and should never be
accessed directly. Use [`cycles(p)`](@ref) instead.

There are two inner constructors of `Perm`:

* `Perm(n::T)` constructs the trivial `Perm{T}`-permutation of length $n$.
* `Perm(v::AbstractVector{<:Integer} [,check=true])` constructs a permutation
  represented by `v`. By default `Perm` constructor checks if the vector
  constitutes a valid permutation. To skip the check call `Perm(v, false)`.

# Examples
```jldoctest
julia> Perm([1,2,3])
()

julia> g = Perm(Int32[2,3,1])
(1,2,3)

julia> typeof(g)
Perm{Int32}
```
"""
mutable struct Perm{T<:Integer} <: AbstractPerm
   d::Vector{T}
   modified::Bool
   cycles::CycleDec{T}

   function Perm(n::T) where T<:Integer
      @assert isconcretetype(T)
      return new{T}(collect(Base.OneTo(n)), false)
   end

   function Perm(v::AbstractVector{T}, check::Bool=true) where T<:Integer
      @assert isconcretetype(T)
      if check
         Base.Set(v) != Base.Set(1:length(v)) && error("Unable to coerce to permutation:
         non-unique elements in array")
      end
      return new{T}(v, false)
   end
end

