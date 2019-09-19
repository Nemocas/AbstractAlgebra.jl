###############################################################################
#
#   GenericTypes.jl : Generic Types
#
###############################################################################

###############################################################################
#
#   Cycle Decomposition
#
###############################################################################

@doc Markdown.doc"""
    CycleDec{T}(ccycles, cptrs, n) where T
> Cycle decomposition of a permutation.
> * `ccycles`: an array of consecutive entries of cycles;
> * `cptrs`: an array of pointers to the locations where cycles begin: ```ccycles[cptrs[i], cptrs[i+1]-1]` contains the i-th cycle;
> * `n`: the number of cycles;
"""
struct CycleDec{T}
   ccycles::Vector{T}
   cptrs::Vector{Int}
   n::Int
end

###############################################################################
#
#   PermGroup / perm
#
###############################################################################

@doc Markdown.doc"""
    PermGroup{T<:Integer}
> The permutation group singleton type.
> `PermGroup(n)` constructs the permutation group $S_n$ on $n$-symbols. The type of elements of the group is inferred from the type of `n`.

# Examples:
```jldoctest; setup = :(using AbstractAlgebra)
julia> G = PermGroup(5)
Permutation group over 5 elements

julia> elem_type(G)
perm{Int64}

julia> H = PermGroup(UInt16(5))
Permutation group over 5 elements

julia> elem_type(H)
perm{UInt16}
```
"""
struct PermGroup{T<:Integer} <: AbstractAlgebra.Group
   n::T
end

@doc Markdown.doc"""
    perm{T<:Integer}
> The type of permutations.
> Fieldnames:
> * `d::Vector{T}` - vector representing the permutation
> * `modified::Bool` - bit to check the validity of cycle decomposition
> * `cycles::CycleDec{T}` - (cached) cycle decomposition
>
> Permutation $p$ consists of a vector (`p.d`) of $n$ integers from $1$ to $n$.
> If the $i$-th entry of the vector is $j$, this corresponds to $p$ sending $i \to j$.
> The cycle decomposition (`p.cycles`) is computed on demand and should never be
> accessed directly. Use [`cycles(p)`](@ref) instead.
>
> There are two inner constructors of `perm`:
>
> * `perm(n::T)` constructs the trivial `perm{T}`-permutation of length $n$.
> * `perm(v::Vector{T<:Integer}[,check=true])` constructs a permutation
> represented by `v`. By default `perm` constructor checks if the vector
> constitutes a valid permutation. To skip the check call `perm(v, false)`.

# Examples:
```jldoctest; setup = :(using AbstractAlgebra)
julia> perm([1,2,3])
()

julia> g = perm(Int32[2,3,1])
(1,2,3)

julia> typeof(g)
perm{Int32}
```
"""
mutable struct perm{T<:Integer} <: AbstractAlgebra.GroupElem
   d::Array{T, 1}
   modified::Bool
   cycles::CycleDec{T}

   function perm(n::T) where T<:Integer
      return new{T}(collect(T, 1:n), false)
   end

   function perm(v::AbstractVector{T}, check::Bool=true) where T<:Integer
      if check
         Set(v) != Set(1:length(v)) && error("Unable to coerce to permutation:
         non-unique elements in array")
      end
      return new{T}(v, false)
   end
end

@doc Markdown.doc"""
    AllPerms(n::T) where T
> Return an iterator over arrays representing all permutations of `1:n`.
> Similar to `Combinatorics.permutations(1:n)`
"""
struct AllPerms{T<:Integer}
   all::Int
   c::Vector{Int}
   elts::perm{T}

   function AllPerms(n::T) where T
      new{T}(Int(factorial(n)), ones(Int, n), perm(collect(T, 1:n), false))
   end
end

###############################################################################
#
#   Partition
#
###############################################################################

@doc Markdown.doc"""
    Partition(part::Vector{<:Integer}[, check::Bool=true]) <: AbstractVector{Int}
> Represent integer partition in the non-increasing order.
>
> `part` will be sorted, if necessary. Checks for validity of input can be skipped by calling the (inner) constructor with `false` as the second argument.
>
> Functionally `Partition` is a thin wrapper over `Vector{Int}`.

> Fieldnames:
>  * `n::Int` - the partitioned number
>  * `part::Vector{Int}` - a non-increasing sequence of summands of `n`.

# Examples:
```jldoctest; setup = :(using AbstractAlgebra)
julia> p = Partition([4,2,1,1,1])
4₁2₁1₃

julia> p.n == sum(p.part)
true
```
"""
mutable struct Partition <: AbstractVector{Int}
   n::Int
   part::Vector{Int}

   function Partition(part::AbstractVector{T}, check::Bool=true) where T<:Integer
      if check
         all(diff(part) .<= 0) || sort!(part, rev=true)
         if length(part) > 0
            part[end] >= 1 || throw(ArgumentError("Found non-positive entry in partition!"))
         end
      end
      return new(sum(part), part)
   end
end

@doc Markdown.doc"""
    AllParts(n::Integer)
> Return an iterator over all integer `Partition`s of `n`.
> Partitions are produced in ascending order according to RuleAsc (Algorithm 3.1) from
>
> > Jerome Kelleher and Barry O’Sullivan,
> > *Generating All Partitions: A Comparison Of Two Encodings*
> > ArXiv:0909.2331
>
> See also `Combinatorics.partitions(1:n)`.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> ap = AllParts(5);


julia> collect(ap)
7-element Array{AbstractAlgebra.Generic.Partition,1}:
 1₅
 2₁1₃
 3₁1₂
 2₂1₁
 4₁1₁
 3₁2₁
 5₁
```
"""
struct AllParts
    n::Int
    part::Vector{Int}
    AllParts(n::Integer) = new(n, zeros(Int,n))
end

###############################################################################
#
#   SkewDiagram
#
###############################################################################

@doc Markdown.doc"""
    SkewDiagram(lambda::Partition, mu::Partition) <: AbstractArray{Int, 2}
> Implements a skew diagram, i.e. a difference of two Young diagrams
> represented by partitions `lambda` and `mu`.
> (below dots symbolise the removed entries)

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> l = Partition([4,3,2])
4₁3₁2₁

julia> m = Partition([3,1,1])
3₁1₂

julia> xi = SkewDiagram(l,m)
3×4 AbstractAlgebra.Generic.SkewDiagram:
 ⋅  ⋅  ⋅  1
 ⋅  1  1
 ⋅  1

```
"""
struct SkewDiagram <: AbstractArray{Int, 2}
   lam::Partition
   mu::Partition

   function SkewDiagram(lambda, mu)
      lambda.n >= mu.n || throw("Can't create SkewDiagram: $mu is partition of  $(mu.n) > $(lambda.n).")
      length(lambda) >= length(mu) || throw("Can't create SkewDiagram: $mu is longer than $(lambda)!")
      for (l, m) in zip(lambda, mu)
         l >= m || throw("a row of $mu is longer than a row of $lambda")
      end
      return new(lambda, mu)
   end
end

###############################################################################
#
#   Young Tableau
#
###############################################################################

@doc Markdown.doc"""
    YoungTableau(part::Partition[, fill::Vector{Int}=collect(1:sum(part))])  <: AbstractArray{Int, 2}
> Return the Young tableaux of partition `part`, filled linearly
> by `fill` vector. Note that `fill` vector is in **row-major** format.
>
> Fields:
> * `part` - the partition defining Young diagram
> * `fill` - the row-major fill vector: the entries of the diagram.

# Examples:
```jldoctest; setup = :(using AbstractAlgebra)
julia> p = Partition([4,3,1]); y = YoungTableau(p)
┌───┬───┬───┬───┐
│ 1 │ 2 │ 3 │ 4 │
├───┼───┼───┼───┘
│ 5 │ 6 │ 7 │
├───┼───┴───┘
│ 8 │
└───┘

julia> y.part
4₁3₁1₁

julia> y.fill
8-element Array{Int64,1}:
 1
 2
 3
 4
 5
 6
 7
 8
```
"""
struct YoungTableau <: AbstractArray{Int, 2}
   part::Partition
   fill::Vector{Int}

   function YoungTableau(part::Partition, fill::Vector{T}=collect(1:sum(part))) where T<:Integer
      sum(part) == length(fill) || throw(ArgumentError("Can't fill Young digaram of $part with $fill: different number of elemnets."))

      # _, fill = conj(part, fill)

      return new(part, fill)
   end
end

###############################################################################
#
#   PolyRing / Poly
#
###############################################################################

mutable struct PolyRing{T <: RingElement} <: AbstractAlgebra.PolyRing{T}
   base_ring::Ring
   S::Symbol

   function PolyRing{T}(R::Ring, s::Symbol, cached::Bool = true) where T <: RingElement
      if cached && haskey(PolyID, (R, s))
         return PolyID[R, s]::PolyRing{T}
      else
         z = new{T}(R, s)
         if cached
           PolyID[R, s] = z
         end
         return z
      end
   end
end

const PolyID = Dict{Tuple{Ring, Symbol}, Ring}()

mutable struct Poly{T <: RingElement} <: AbstractAlgebra.PolyElem{T}
   coeffs::Array{T, 1}
   length::Int
   parent::PolyRing{T}

   Poly{T}() where T <: RingElement = new{T}(Array{T}(undef, 0), 0)

   function Poly{T}(b::Array{T, 1}) where T <: RingElement
      z = new{T}(b)
      z.length = normalise(z, length(b))
      return z
   end

   Poly{T}(a::T) where T <: RingElement = iszero(a) ? new{T}(Array{T}(undef, 0), 0) : new{T}([a], 1)
end

###############################################################################
#
#   NCPolyRing / NCPoly
#
###############################################################################

mutable struct NCPolyRing{T <: NCRingElem} <: AbstractAlgebra.NCPolyRing{T}
   base_ring::NCRing
   S::Symbol

   function NCPolyRing{T}(R::NCRing, s::Symbol, cached::Bool = true) where T <: NCRingElem
      if cached && haskey(NCPolyID, (R, s))
         return NCPolyID[R, s]::NCPolyRing{T}
      else
         z = new{T}(R, s)
         if cached
           NCPolyID[R, s] = z
         end
         return z
      end
   end
end

const NCPolyID = Dict{Tuple{NCRing, Symbol}, NCRing}()

mutable struct NCPoly{T <: NCRingElem} <: AbstractAlgebra.NCPolyElem{T}
   coeffs::Array{T, 1}
   length::Int
   parent::NCPolyRing{T}

   NCPoly{T}() where T <: NCRingElem = new{T}(Array{T}(undef, 0), 0)

   function NCPoly{T}(b::Array{T, 1}) where T <: NCRingElem
      z = new{T}(b)
      z.length = normalise(z, length(b))
      return z
   end

   NCPoly{T}(a::T) where T <: NCRingElem = iszero(a) ? new{T}(Array{T}(undef, 0), 0) : new{T}([a], 1)
end

const PolynomialElem{T} = Union{AbstractAlgebra.PolyElem{T}, AbstractAlgebra.NCPolyElem{T}}

###############################################################################
#
#   MPolyRing / MPoly
#
###############################################################################

# S is a Symbol which can take the values:
# :lex
# :revlex
# :deglex
# :degrevlex
#
# T is an Int which is the number of variables
# (plus one if ordered by total degree)

mutable struct MPolyRing{T <: RingElement} <: AbstractAlgebra.MPolyRing{T}
   base_ring::Ring
   S::Array{Symbol, 1}
   ord::Symbol
   num_vars::Int
   N::Int

   function MPolyRing{T}(R::Ring, s::Array{Symbol, 1}, ord::Symbol, N::Int,
                         cached::Bool = true) where T <: RingElement
      if cached && haskey(MPolyID, (R, s, ord, N))
         return MPolyID[R, s, ord, N]::MPolyRing{T}
      else
         z = new{T}(R, s, ord, length(s), N)
         if cached
           MPolyID[R, s, ord, N] = z
         end
         return z
      end
   end
end

const MPolyID = Dict{Tuple{Ring, Array{Symbol, 1}, Symbol, Int}, Ring}()

mutable struct MPoly{T <: RingElement} <: AbstractAlgebra.MPolyElem{T}
   coeffs::Array{T, 1}
   exps::Array{UInt, 2}
   length::Int
   parent::MPolyRing{T}

   function MPoly{T}(R::MPolyRing) where T <: RingElement
      N = R.N
      return new{T}(Array{T}(undef, 0), Array{UInt}(undef, N, 0), 0, R)
   end

   MPoly{T}(R::MPolyRing, a::Array{T, 1}, b::Array{UInt, 2}) where T <: RingElement = new{T}(a, b, length(a), R)

   function MPoly{T}(R::MPolyRing, a::T) where T <: RingElement
      N = R.N
      return iszero(a) ? new{T}(Array{T}(undef, 0), Array{UInt}(undef, N, 0), 0, R) :
                                          new{T}([a], zeros(UInt, N, 1), 1, R)
   end
end

# Iterators

struct MPolyCoeffs{T <: AbstractAlgebra.MPolyElem}
   poly::T
end

struct MPolyExponentVectors{T <: AbstractAlgebra.MPolyElem}
   poly::T
end

struct MPolyTerms{T <: AbstractAlgebra.MPolyElem}
   poly::T
end

struct MPolyMonomials{T <: AbstractAlgebra.MPolyElem}
   poly::T
end

mutable struct MPolyBuildCtx{T <: AbstractAlgebra.MPolyElem, S}
  poly::T
  state::S

  function MPolyBuildCtx(R::T, s::S) where {S, T <: AbstractAlgebra.MPolyRing}
    return new{elem_type(T), S}(R())
  end
end

###############################################################################
#
#   SparsePolyRing / SparsePoly
#
###############################################################################

mutable struct SparsePolyRing{T <: RingElement} <: AbstractAlgebra.Ring
   base_ring::Ring
   S::Symbol
   num_vars::Int

   function SparsePolyRing{T}(R::Ring, s::Symbol, cached::Bool = true) where T <: RingElement
      if cached && haskey(SparsePolyID, (R, s))
         return SparsePolyID[R, s]::SparsePolyRing{T}
      else
         z = new{T}(R, s)
         if cached
           SparsePolyID[R, s] = z
         end
         return z
      end
   end
end

const SparsePolyID = Dict{Tuple{Ring, Symbol}, SparsePolyRing}()

mutable struct SparsePoly{T <: RingElement} <: AbstractAlgebra.RingElem
   coeffs::Array{T, 1}
   exps::Array{UInt}
   length::Int
   parent::SparsePolyRing{T}

   SparsePoly{T}() where T <: RingElement = new{T}(Array{T}(undef, 0), Array{UInt}(undef, 0), 0)

   SparsePoly{T}(a::Array{T, 1}, b::Array{UInt, 1}) where T <: RingElement = new{T}(a, b, length(a))

   SparsePoly{T}(a::T) where T <: RingElement = iszero(a) ? new{T}(Array{T}(undef, 0), Array{UInt}(undef, 0), 0) :
                                               new{T}([a], [UInt(0)], 1)
end

###############################################################################
#
#   ResRing / Res
#
###############################################################################

mutable struct ResRing{T <: RingElement} <: AbstractAlgebra.ResRing{T}
   base_ring::Ring
   modulus::T

   function ResRing{T}(modulus::T, cached::Bool = true) where T <: RingElement
      c = canonical_unit(modulus)
      if !isone(c)
        modulus = divexact(modulus, c)
      end
      if cached && haskey(ModulusDict, (parent(modulus), modulus))
         return ModulusDict[parent(modulus), modulus]::ResRing{T}
      else
         z = new{T}(parent(modulus), modulus)
         if cached
            ModulusDict[parent(modulus), modulus] = z
         end
         return z
      end
   end
end

const ModulusDict = Dict{Tuple{Ring, RingElement}, Ring}()

mutable struct Res{T <: RingElement} <: AbstractAlgebra.ResElem{T}
   data::T
   parent::ResRing{T}

   Res{T}(a::T) where T <: RingElement = new{T}(a)
end

###############################################################################
#
#   ResField / ResFieldElem
#
###############################################################################

mutable struct ResField{T <: RingElement} <: AbstractAlgebra.ResField{T}
   base_ring::Ring
   modulus::T

   function ResField{T}(modulus::T, cached::Bool = true) where T <: RingElement
      c = canonical_unit(modulus)
      if !isone(c)
        modulus = divexact(modulus, c)
      end
      if cached && haskey(ModulusFieldDict, (parent(modulus), modulus))
         return ModulusFieldDict[parent(modulus), modulus]::ResField{T}
      else
         z = new{T}(parent(modulus), modulus)
         if cached
            ModulusFieldDict[parent(modulus), modulus] = z
         end
         return z
      end
   end
end

const ModulusFieldDict = Dict{Tuple{Ring, RingElement}, Field}()

mutable struct ResF{T <: RingElement} <: AbstractAlgebra.ResFieldElem{T}
   data::T
   parent::ResField{T}

   ResF{T}(a::T) where T <: RingElement = new{T}(a)
end

###############################################################################
#
#   RelSeriesRing / RelSeries
#
###############################################################################

mutable struct RelSeriesRing{T <: RingElement} <: AbstractAlgebra.SeriesRing{T}
   base_ring::Ring
   prec_max::Int
   S::Symbol

   function RelSeriesRing{T}(R::Ring, prec::Int, s::Symbol, cached::Bool = true) where T <: RingElement
      if cached && haskey(RelSeriesID, (R, prec, s))
         return RelSeriesID[R, prec, s]::RelSeriesRing{T}
      else
         z = new{T}(R, prec, s)
         if cached
            RelSeriesID[R, prec, s] = z
         end
         return z
      end
   end
end

const RelSeriesID = Dict{Tuple{Ring, Int, Symbol}, Ring}()

mutable struct RelSeries{T <: RingElement} <: AbstractAlgebra.RelSeriesElem{T}
   coeffs::Array{T, 1}
   length::Int
   prec::Int
   val::Int
   parent::RelSeriesRing{T}

   function RelSeries{T}(a::Array{T, 1}, length::Int, prec::Int, val::Int) where T <: RingElement
      new{T}(a, length, prec, val)
   end

   RelSeries{T}(a::RelSeries{T}) where T <: RingElement = a
end

###############################################################################
#
#   AbsSeriesRing / AbsSeries
#
###############################################################################

mutable struct AbsSeriesRing{T <: RingElement} <: AbstractAlgebra.SeriesRing{T}
   base_ring::Ring
   prec_max::Int
   S::Symbol

   function AbsSeriesRing{T}(R::Ring, prec::Int, s::Symbol, cached::Bool = true) where T <: RingElement
      if cached && haskey(AbsSeriesID, (R, prec, s))
         return AbsSeriesID[R, prec, s]::AbsSeriesRing{T}
      else
         z = new{T}(R, prec, s)
         if cached
            AbsSeriesID[R, prec, s] = z
         end
         return z
      end
   end
end

const AbsSeriesID = Dict{Tuple{Ring, Int, Symbol}, Ring}()

mutable struct AbsSeries{T <: RingElement} <: AbstractAlgebra.AbsSeriesElem{T}
   coeffs::Array{T, 1}
   length::Int
   prec::Int
   parent::AbsSeriesRing{T}

   AbsSeries{T}(a::Array{T, 1}, length::Int, prec::Int) where T <: RingElement = new{T}(a, length, prec)
   AbsSeries{T}(a::AbsSeries{T}) where T <: RingElement = a
end

###############################################################################
#
#   LaurentSeriesRing / LarentSeriesRingElem
#
###############################################################################

mutable struct LaurentSeriesRing{T <: RingElement} <: AbstractAlgebra.Ring
   base_ring::Ring
   prec_max::Int
   S::Symbol

   function LaurentSeriesRing{T}(R::Ring, prec::Int, s::Symbol, cached::Bool = true) where T <: RingElement
      if cached && haskey(LaurentSeriesID, (R, prec, s))
         return LaurentSeriesID[R, prec, s]::LaurentSeriesRing{T}
      else
         z = new{T}(R, prec, s)
         if cached
            LaurentSeriesID[R, prec, s] = z
         end
         return z
      end
   end
end

const LaurentSeriesID = Dict{Tuple{Ring, Int, Symbol}, Ring}()

mutable struct LaurentSeriesRingElem{T <: RingElement} <: AbstractAlgebra.RingElem
   coeffs::Array{T, 1}
   length::Int
   prec::Int
   val::Int
   scale::Int
   parent::LaurentSeriesRing{T}

   function LaurentSeriesRingElem{T}(a::Array{T, 1}, length::Int, prec::Int, val::Int, scale::Int) where T <: RingElement
      new{T}(a, length, prec, val, scale)
   end
end

###############################################################################
#
#   LaurentSeriesField / LarentSeriesFieldElem
#
###############################################################################

mutable struct LaurentSeriesField{T <: FieldElement} <: AbstractAlgebra.Field
   base_ring::Field
   prec_max::Int
   S::Symbol

   function LaurentSeriesField{T}(R::Field, prec::Int, s::Symbol, cached::Bool = true) where T <: FieldElement
      if cached && haskey(LaurentSeriesID, (R, prec, s))
         return LaurentSeriesID[R, prec, s]::LaurentSeriesField{T}
      else
         z = new{T}(R, prec, s)
         if cached
            LaurentSeriesID[R, prec, s] = z
         end
         return z
      end
   end
end

mutable struct LaurentSeriesFieldElem{T <: FieldElement} <: AbstractAlgebra.FieldElem
   coeffs::Array{T, 1}
   length::Int
   prec::Int
   val::Int
   scale::Int
   parent::LaurentSeriesField{T}

   function LaurentSeriesFieldElem{T}(a::Array{T, 1}, length::Int, prec::Int, val::Int, scale::Int) where T <: FieldElement
      new{T}(a, length, prec, val, scale)
   end
end

const LaurentSeriesElem{T} = Union{LaurentSeriesRingElem{T}, LaurentSeriesFieldElem{T}} where T <: RingElement

###############################################################################
#
#   PuiseuxSeriesRing / PuiseuxSeriesRingElem
#
###############################################################################

mutable struct PuiseuxSeriesRing{T <: RingElement} <: AbstractAlgebra.Ring
   laurent_ring::Ring

   function PuiseuxSeriesRing{T}(R::LaurentSeriesRing{T}, cached::Bool = true) where T <: RingElement
      if cached && haskey(PuiseuxSeriesID, R)
         return PuiseuxSeriesID[R]::PuiseuxSeriesRing{T}
      else
         z = new{T}(R)
         if cached
            PuiseuxSeriesID[R] = z
         end
         return z
      end
   end
end

const PuiseuxSeriesID = Dict{Ring, Ring}()

mutable struct PuiseuxSeriesRingElem{T <: RingElement} <: AbstractAlgebra.RingElem
   data::LaurentSeriesRingElem{T}
   scale::Int
   parent::PuiseuxSeriesRing{T}

   function PuiseuxSeriesRingElem{T}(d::LaurentSeriesRingElem{T}, scale::Int) where T <: RingElement
      new{T}(d, scale)
   end
end

###############################################################################
#
#   PuiseuxSeriesField / PuiseuxSeriesFieldElem
#
###############################################################################

mutable struct PuiseuxSeriesField{T <: FieldElement} <: AbstractAlgebra.Field
   laurent_ring::Field

   function PuiseuxSeriesField{T}(R::LaurentSeriesField{T}, cached::Bool = true) where T <: FieldElement
      if cached && haskey(PuiseuxSeriesID, R)
         return PuiseuxSeriesID[R]::PuiseuxSeriesField{T}
      else
         z = new{T}(R)
         if cached
            PuiseuxSeriesID[R] = z
         end
         return z
      end
   end
end

mutable struct PuiseuxSeriesFieldElem{T <: FieldElement} <: AbstractAlgebra.FieldElem
   data::LaurentSeriesFieldElem{T}
   scale::Int
   parent::PuiseuxSeriesField{T}

   function PuiseuxSeriesFieldElem{T}(d::LaurentSeriesFieldElem{T}, scale::Int) where T <: FieldElement
      new{T}(d, scale)
   end
end

const PuiseuxSeriesElem{T} = Union{PuiseuxSeriesRingElem{T}, PuiseuxSeriesFieldElem{T}} where T <: RingElement

###############################################################################
#
#   FracField / Frac
#
###############################################################################

mutable struct FracField{T <: RingElem} <: AbstractAlgebra.FracField{T}
   base_ring::Ring

   function FracField{T}(R::Ring, cached::Bool = true) where T <: RingElem
      if cached && haskey(FracDict, R)
         return FracDict[R]::FracField{T}
      else
         z = new{T}(R)
         if cached
            FracDict[R] = z
         end
         return z
      end
   end
end

const FracDict = Dict{Ring, Ring}()

mutable struct Frac{T <: RingElem} <: AbstractAlgebra.FracElem{T}
   num::T
   den::T
   parent::FracField{T}

   Frac{T}(num::T, den::T) where T <: RingElem = new{T}(num, den)
end

###############################################################################
#
#   MatSpace / Mat
#
###############################################################################

# All MatSpaceElem's and views thereof
abstract type Mat{T} <: MatElem{T} end

# not really a mathematical ring
mutable struct MatSpace{T <: RingElement} <: AbstractAlgebra.MatSpace{T}
   nrows::Int
   ncols::Int
   base_ring::Ring

   function MatSpace{T}(R::Ring, r::Int, c::Int, cached::Bool = true) where T <: RingElement
      if cached && haskey(MatDict, (R, r, c))
         return MatDict[R, r, c]::MatSpace{T}
      else
         z = new{T}(r, c, R)
         if cached
            MatDict[R, r, c] = z
         end
         return z
      end
   end
end

const MatDict = Dict{Tuple{Ring, Int, Int}, MatSpace}()

mutable struct MatSpaceElem{T <: RingElement} <: Mat{T}
   entries::Array{T, 2}
   base_ring::Ring

   function MatSpaceElem{T}(A::Array{T, 2}) where T <: RingElement
      return new{T}(A)
    end

   function MatSpaceElem{T}(A::AbstractArray{T, 2}) where T <: RingElement
      return new{T}(Array(A))
   end

   function MatSpaceElem{T}(r::Int, c::Int, A::Array{T, 1}) where T <: RingElement
      t = Array{T}(undef, r, c)
      for i = 1:r
         for j = 1:c
            t[i, j] = A[(i - 1) * c + j]
         end
      end
      return new{T}(t)
   end
end

mutable struct MatSpaceView{T <: RingElement, V, W} <: Mat{T}
   entries::SubArray{T, 2, Array{T, 2}, V, W}
   base_ring::Ring
end

###############################################################################
#
#   MatAlgebra / MatAlgElem
#
###############################################################################

mutable struct MatAlgebra{T <: RingElement} <: AbstractAlgebra.MatAlgebra{T}
   n::Int
   base_ring::Ring

   function MatAlgebra{T}(R::Ring, n::Int, cached::Bool = true) where T <: RingElement
      if haskey(MatAlgDict, (R, n))
         return MatAlgDict[R, n]::MatAlgebra{T}
      else
         z = new{T}(n, R)
         if cached
            MatAlgDict[R, n] = z
         end
         return z
      end
   end
end

const MatAlgDict = Dict{Tuple{Ring, Int}, NCRing}()

mutable struct MatAlgElem{T <: RingElement} <: AbstractAlgebra.MatAlgElem{T}
   entries::Array{T, 2}
   base_ring::Ring

   function MatAlgElem{T}(A::Array{T, 2}) where T <: RingElement
      return new{T}(A)
   end

   function MatAlgElem{T}(n::Int, A::Array{T, 1}) where T <: RingElement
      t = Array{T}(undef, n, n)
      for i = 1:n
         for j = 1:n
            t[i, j] = A[(i - 1) * c + j]
         end
      end
      return new{T}(t)
   end
end

const MatrixElem{T} = Union{AbstractAlgebra.MatElem{T}, AbstractAlgebra.MatAlgElem{T}}

###############################################################################
#
#   CompositeMap
#
###############################################################################

Map(::Type{T}) where T <: AbstractAlgebra.Map = supertype(T)
Map(::Type{S}) where S <: AbstractAlgebra.SetMap = Map{D, C, <:S, T} where {D, C, T}

mutable struct CompositeMap{D, C} <: AbstractAlgebra.Map{D, C, AbstractAlgebra.SetMap, CompositeMap}
   domain::D
   codomain::C
   map1::AbstractAlgebra.Map
   map2::AbstractAlgebra.Map

   function CompositeMap(map1::AbstractAlgebra.Map{D, U}, map2::AbstractAlgebra.Map{U, C}) where {D, U, C}
      return new{D, C}(domain(map1), codomain(map2), map1, map2)
   end
end

###############################################################################
#
#   FunctionalMap
#
###############################################################################

mutable struct FunctionalMap{D, C} <: AbstractAlgebra.Map{D, C, AbstractAlgebra.FunctionalMap, FunctionalMap}
    domain::D
    codomain::C
    image_fn::Function
end

###############################################################################
#
#   IdentityMap
#
###############################################################################

struct IdentityMap{D} <: AbstractAlgebra.Map{D, D, AbstractAlgebra.IdentityMap, IdentityMap}
   domain::D
end

###############################################################################
#
#   FunctionalCompositeMap
#
###############################################################################

mutable struct FunctionalCompositeMap{D, C} <: AbstractAlgebra.Map{D, C, AbstractAlgebra.FunctionalMap, FunctionalCompositeMap}
   domain::D
   codomain::C
   map1::AbstractAlgebra.Map
   map2::AbstractAlgebra.Map
   fn_cache::Function

   function FunctionalCompositeMap(map1::Map(AbstractAlgebra.FunctionalMap){D, U}, map2::Map(AbstractAlgebra.FunctionalMap){U, C}) where {D, U, C}
      return new{D, C}(domain(map1), codomain(map2), map1, map2)
   end
end

###############################################################################
#
#   MapWithSection
#
###############################################################################

mutable struct MapWithSection{D, C} <: AbstractAlgebra.Map{D, C, AbstractAlgebra.SetMap, MapWithSection}
   map::AbstractAlgebra.Map
   section::AbstractAlgebra.Map

   function MapWithSection(map::AbstractAlgebra.Map{D, C}, section::AbstractAlgebra.Map{C, D}) where {D, C}
      (domain(map) != codomain(section) || codomain(map) != domain(section)) &&
error("Maps not compatible")
      return new{D, C}(map, section)
   end

   function MapWithSection(map::AbstractAlgebra.Map{D, C}) where {D, C}
      return new{D, C}(map)
   end
end

###############################################################################
#
#   MapWithRetraction
#
###############################################################################

mutable struct MapWithRetraction{D, C} <: AbstractAlgebra.Map{D, C, AbstractAlgebra.SetMap, MapWithRetraction}
   map::AbstractAlgebra.Map
   retraction::AbstractAlgebra.Map

   function MapWithRetraction(map::AbstractAlgebra.Map{D, C}, retraction::AbstractAlgebra.Map{C, D}) where {D, C}
      (domain(map) != codomain(retraction) || codomain(map) != domain(retraction)) &&
error("Maps not compatible")
      return new{D, C}(map, retraction)
   end

   function MapWithRetraction(map::AbstractAlgebra.Map{D, C}) where {D, C}
      return new{D, C}(map)
   end
end

###############################################################################
#
#   MapCache
#
###############################################################################

mutable struct MapCache{D, C, S, T, De, Ce} <: AbstractAlgebra.Map{D, C, S, T}
   map::AbstractAlgebra.Map{D, C}
   limit::Int
   enabled::Bool
   image_cache::Dict{De, Ce}

   function MapCache(f::AbstractAlgebra.Map{D, C, S, T}, limit::Int, enabled::Bool) where {D, C, S, T}
      De = elem_type(D)
      Ce = elem_type(C)
      r = new{D, C, S, T, De, Ce}(f, limit)
      if enabled
         r.image_cache = Dict{De, Ce}()
      end
      r.enabled = enabled
      return r
   end
end

###############################################################################
#
#   FreeModule/free_module_elem
#
###############################################################################

mutable struct FreeModule{T <: Union{RingElement, NCRingElem}} <: AbstractAlgebra.FPModule{T}
   rank::Int
   base_ring::NCRing

   function FreeModule{T}(R::NCRing, rank::Int, cached::Bool = true) where T <: Union{RingElement, NCRingElem}
      if cached && haskey(FreeModuleDict, (R, rank))
         return FreeModuleDict[R, rank]::FreeModule{T}
      else
         z = new{T}(rank, R)
         if cached
            FreeModuleDict[R, rank] = z
         end
         return z
      end
   end
end

const FreeModuleDict = Dict{Tuple{NCRing, Int}, FreeModule}()

mutable struct free_module_elem{T <: Union{RingElement, NCRingElem}} <: AbstractAlgebra.FPModuleElem{T}
    v::AbstractAlgebra.MatElem{T}
    parent::FreeModule{T}

    function free_module_elem{T}(m::AbstractAlgebra.FPModule{T}, v::AbstractAlgebra.MatElem{T}) where T <: Union{RingElement, NCRingElem}
       z = new{T}(v, m)
    end
end

###############################################################################
#
#   ModuleHomomorphism
#
###############################################################################

mutable struct ModuleHomomorphism{T <: RingElement} <: AbstractAlgebra.Map{AbstractAlgebra.FPModule{T}, AbstractAlgebra.FPModule{T}, AbstractAlgebra.FPModuleHomomorphism, ModuleHomomorphism}

   domain::AbstractAlgebra.FPModule{T}
   codomain::AbstractAlgebra.FPModule{T}
   matrix::AbstractAlgebra.MatElem{T}
   image_fn::Function

   function ModuleHomomorphism{T}(D::AbstractAlgebra.FPModule{T}, C::AbstractAlgebra.FPModule{T}, m::AbstractAlgebra.MatElem{T}) where T <: RingElement
      z = new(D, C, m, x::AbstractAlgebra.FPModuleElem{T} -> C(x.v*m))
   end
end

mutable struct ModuleIsomorphism{T <: RingElement} <: AbstractAlgebra.Map{AbstractAlgebra.FPModule{T}, AbstractAlgebra.FPModule{T}, AbstractAlgebra.FPModuleHomomorphism, ModuleIsomorphism}

   domain::AbstractAlgebra.FPModule{T}
   codomain::AbstractAlgebra.FPModule{T}
   matrix::AbstractAlgebra.MatElem{T}
   inverse_matrix::AbstractAlgebra.MatElem{T}
   image_fn::Function
   inverse_image_fn::Function

   function ModuleIsomorphism{T}(D::AbstractAlgebra.FPModule{T}, C::AbstractAlgebra.FPModule{T}, m::AbstractAlgebra.MatElem{T}, minv::AbstractAlgebra.MatElem{T}) where T <: RingElement
      z = new(D, C, m, minv, x::AbstractAlgebra.FPModuleElem{T} -> C(x.v*m), y::AbstractAlgebra.FPModuleElem{T} -> D(y.v*minv))
   end
end

###############################################################################
#
#   Submodule/submodule_elem
#
###############################################################################

mutable struct Submodule{T <: RingElement} <: AbstractAlgebra.FPModule{T}
   m::AbstractAlgebra.FPModule{T}
   gens::Vector{<:AbstractAlgebra.FPModuleElem{T}}
   rels::Vector{<:AbstractAlgebra.MatElem{T}}
   gen_cols::Vector{Int}
   pivots::Vector{Int}
   base_ring::Ring
   map::ModuleHomomorphism{T}

   function Submodule{T}(M::AbstractAlgebra.FPModule{T}, gens::Vector{<:AbstractAlgebra.FPModuleElem{T}}, rels::Vector{<:AbstractAlgebra.MatElem{T}}, gen_cols::Vector{Int}, pivots::Vector{Int}) where T <: RingElement
      z = new{T}(M, gens, rels, gen_cols, pivots, base_ring(M))
   end
end

mutable struct submodule_elem{T <: RingElement} <: AbstractAlgebra.FPModuleElem{T}
   v::AbstractAlgebra.MatElem{T}
   parent::Submodule{T}

   function submodule_elem{T}(m::AbstractAlgebra.FPModule{T}, v::AbstractAlgebra.MatElem{T}) where T <: RingElement
      z = new{T}(v, m)
   end
end

###############################################################################
#
#   QuotientModule/quotient_module_elem
#
###############################################################################

mutable struct QuotientModule{T <: RingElement} <: AbstractAlgebra.FPModule{T}
   m::AbstractAlgebra.FPModule{T}
   rels::Vector{<:AbstractAlgebra.MatElem{T}}
   gen_cols::Vector{Int} # which original columns correspond to gens of quotient
   pivots::Vector{Int} # pivot column of each culled relation in new rels matrix
   base_ring::Ring
   map::ModuleHomomorphism{T}

   function QuotientModule{T}(M::AbstractAlgebra.FPModule{T}, combined_rels::AbstractAlgebra.MatElem{T}) where T <: RingElement
      # concatenate relations in M and new rels
      R = base_ring(M)
      # remove zero rows and all rows/cols corresponding to unit pivots
      gen_cols, culled, pivots = cull_matrix(combined_rels)
      # put all the culled relations into new relations
      new_rels = [matrix(R, 1, length(gen_cols),
                    [combined_rels[culled[i], gen_cols[j]]
                       for j in 1:length(gen_cols)]) for i = 1:length(culled)]
      # create quotient module
      z = new{T}(M, new_rels, gen_cols, pivots, base_ring(M))
      return z
   end
end

mutable struct quotient_module_elem{T <: RingElement} <: AbstractAlgebra.FPModuleElem{T}
   v::AbstractAlgebra.MatElem{T}
   parent::QuotientModule{T}

   function quotient_module_elem{T}(m::AbstractAlgebra.FPModule{T}, v::AbstractAlgebra.MatElem{T}) where T <: RingElement
      z = new{T}(v, m)
   end
end

###############################################################################
#
#   SNFModule/snf_module_elem
#
###############################################################################

mutable struct SNFModule{T <: RingElement} <: AbstractAlgebra.FPModule{T}
   m::AbstractAlgebra.FPModule{T}
   gens::Vector{<:AbstractAlgebra.FPModuleElem{T}}
   invariant_factors::Vector{T}
   base_ring::Ring
   map::ModuleIsomorphism{T}

   function SNFModule{T}(M::AbstractAlgebra.FPModule{T}, gens::Vector{<:AbstractAlgebra.FPModuleElem{T}}, invariant_factors::Vector{T}) where T <: RingElement
      return new{T}(M, gens, invariant_factors, base_ring(M))
   end
end

mutable struct snf_module_elem{T <: RingElement} <: AbstractAlgebra.FPModuleElem{T}
   v::AbstractAlgebra.MatElem{T}
   parent::SNFModule{T}

   function snf_module_elem{T}(m::AbstractAlgebra.FPModule{T}, v::AbstractAlgebra.MatElem{T}) where T <: RingElement
      z = new{T}(v, m)
   end
end

###############################################################################
#
#   DirectSumModule/direct_sum_module_elem
#
###############################################################################

mutable struct DirectSumModule{T <: RingElement} <: AbstractAlgebra.FPModule{T}
   m::Vector{<:AbstractAlgebra.FPModule{T}}
   rels::Vector{<:AbstractAlgebra.MatElem{T}}
   inj::Vector{<:ModuleHomomorphism{T}}
   pro::Vector{<:ModuleHomomorphism{T}}

   function DirectSumModule{T}(m::Vector{<:AbstractAlgebra.FPModule{T}}, rels::Vector{<:AbstractAlgebra.MatElem{T}}) where T <: RingElement
      return new{T}(m, rels)
   end
end

mutable struct direct_sum_module_elem{T <: RingElement} <: AbstractAlgebra.FPModuleElem{T}
   v::AbstractAlgebra.MatElem{T}
   parent::DirectSumModule{T}

   function direct_sum_module_elem{T}(m::AbstractAlgebra.FPModule{T}, v::AbstractAlgebra.MatElem{T}) where T <: RingElement
      z = new{T}(v, m)
   end
end
