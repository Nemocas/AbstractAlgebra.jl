###############################################################################
#
#   GenericTypes.jl : Generic Types
#
###############################################################################

###############################################################################
#
#   PermGroup / perm
#
###############################################################################

const PermID = ObjectIdDict()

mutable struct PermGroup{T<:Integer} <: AbstractAlgebra.Group
  n::T

  function PermGroup(n::T, cached=true) where T<:Integer
     if haskey(PermID, n)
        return PermID[n]::PermGroup{T}
     else
        z = new{T}(n)
        if cached
           PermID[n] = z
        end
        return z
     end
  end
end

mutable struct perm{T<:Integer} <: AbstractAlgebra.GroupElem
  d::Array{T, 1}
  cycles::Vector{Vector{T}}
  parent::PermGroup{T}

  function perm(n::T) where T<:Integer
     return new{T}(collect(1:n))
  end

  function perm(a::Array{T, 1}) where T<:Integer
     return new{T}(a)
  end
end

doc"""
    AllPerms(n::Int)
> Returns an iterator over arrays representing all permutations of `1:n`.
> Similar to `Combinatorics.permutations(1:n)`
"""
struct AllPerms{T}
  n::T
  all::Int

  function AllPerms(n::T) where {T<:Integer}
     return new{T}(n, factorial(n))
  end
end

###############################################################################
#
#   Partition
#
###############################################################################

doc"""
    Partition(part::Vector{Int}, check::Bool=true)
> Partition represents integer partition into numbers in non-increasing order.
> It is a thin wrapper over `Vector{Int}`
"""
mutable struct Partition{T} <: AbstractVector{Integer}
   n::T
   part::Vector{T}

   function Partition(part::Vector{T}, check::Bool=true) where T
      if check
         all(diff(part) .<= zero(T)) || throw("Partition must be decreasing!")
         if length(part) > zero(T)
            part[end] >= one(T) || throw("Found non-positive entry in partition!")
         end
      end
      return new{T}(sum(part), part)
   end
end

doc"""
    AllParts(n::Int)
> Returns an iterator over all integer `Partition`s of `n`. They come in
> ascending order. See also `Combinatorics.partitions(n)`.
"""
struct AllParts{T}
    n::T
end

###############################################################################
#
#   SkewDiagram
#
###############################################################################

doc"""
    SkewDiagram(lambda::Partition, mu::Partition)
> Implements a skew diagram, i.e. a difference of two Young diagrams
> represented by partitions `lambda` and `mu`.
"""
struct SkewDiagram
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

doc"""
    YoungTableau(part::Partition, tab::Array{Int, 2})
> Returns the Young tableaux of partition `part` of `n`, filled linearly
> (row-major) by `fill` vector.
"""
struct YoungTableau <: AbstractArray{Int, 2}
   part::Partition
   tab::Array{Int,2}
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
      if haskey(PolyID, (R, s))
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

   Poly{T}() where T <: RingElement = new{T}(Array{T}(0), 0)

   function Poly{T}(b::Array{T, 1}) where T <: RingElement
      z = new{T}(b)
      z.length = normalise(z, length(b))
      return z
   end

   Poly{T}(a::T) where T <: RingElement = iszero(a) ? new{T}(Array{T}(0), 0) : new{T}([a], 1)
end

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
      if haskey(MPolyID, (R, s, ord, N))
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
      return new{T}(Array{T}(0), Array{UInt}(N, 0), 0, R)
   end

   MPoly{T}(R::MPolyRing, a::Array{T, 1}, b::Array{UInt, 2}) where T <: RingElement = new{T}(a, b, length(a), R)

   function MPoly{T}(R::MPolyRing, a::T) where T <: RingElement
      N = R.N
      return iszero(a) ? new{T}(Array{T}(0), Array{UInt}(N, 0), 0, R) :
                                          new{T}([a], zeros(UInt, N, 1), 1, R)
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
      if haskey(SparsePolyID, (R, s))
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

   SparsePoly{T}() where T <: RingElement = new{T}(Array{T}(0), Array{UInt}(0), 0)

   SparsePoly{T}(a::Array{T, 1}, b::Array{UInt, 1}) where T <: RingElement = new{T}(a, b, length(a))

   SparsePoly{T}(a::T) where T <: RingElement = iszero(a) ? new{T}(Array{T}(0), Array{UInt}(0), 0) :
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
      if haskey(ModulusDict, (parent(modulus), modulus))
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
      if haskey(ModulusFieldDict, (parent(modulus), modulus))
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
      if haskey(RelSeriesID, (R, prec, s))
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
      if haskey(AbsSeriesID, (R, prec, s))
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
      if haskey(LaurentSeriesID, (R, prec, s))
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
      if haskey(LaurentSeriesID, (R, prec, s))
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
      if haskey(PuiseuxSeriesID, R)
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
      if haskey(PuiseuxSeriesID, R)
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
      if haskey(FracDict, R)
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

# not really a mathematical ring
mutable struct MatSpace{T <: RingElement} <: AbstractAlgebra.MatSpace{T}
   rows::Int
   cols::Int
   base_ring::Ring

   function MatSpace{T}(R::Ring, r::Int, c::Int, cached::Bool = true) where T <: RingElement
      if haskey(MatDict, (R, r, c))
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

const MatDict = Dict{Tuple{Ring, Int, Int}, Ring}()

mutable struct Mat{T <: RingElement} <: AbstractAlgebra.MatElem{T}
   entries::Array{T, 2}
   base_ring::Ring

   function Mat{T}(A::Array{T, 2}) where T <: RingElement
      return new{T}(A)
   end

   function Mat{T}(r::Int, c::Int, A::Array{T, 1}) where T <: RingElement
      t = Array{T}(r, c)
      for i = 1:r
         for j = 1:c
            t[i, j] = A[(i - 1) * c + j]
         end
      end
      return new{T}(t)
   end
end

###############################################################################
#
#   CompositeMap
#
###############################################################################

mutable struct CompositeMap{D, C} <: AbstractAlgebra.Map{D, C, CompositeMap}
   domain::D
   codomain::C
   map1::AbstractAlgebra.Map
   map2::AbstractAlgebra.Map

   function CompositeMap(map1::AbstractAlgebra.Map{U, C}, map2::AbstractAlgebra.Map{D, U}) where {D, U, C}
      return new{D, C}(domain(map2), codomain(map1), map1, map2)
   end
end

###############################################################################
#
#   FunctionalMap
#
###############################################################################

mutable struct FunctionalMap{D, C} <: AbstractAlgebra.Map{D, C, AbstractAlgebra.FunctionalMap}
    image_fn::Function
    domain::D
    codomain::C
end

###############################################################################
#
#   IdentityMap
#
###############################################################################

struct IdentityMap{D} <: AbstractAlgebra.Map{D, D, AbstractAlgebra.IdentityMap}
   domain::D
end

###############################################################################
#
#   FunctionalCompositeMap
#
###############################################################################

mutable struct FunctionalCompositeMap{D, C} <: AbstractAlgebra.Map{D, C, AbstractAlgebra.FunctionalMap}
   domain::D
   codomain::C
   map1::AbstractAlgebra.Map
   map2::AbstractAlgebra.Map
   fn_cache::Function

   function FunctionalCompositeMap(map1::AbstractAlgebra.Map{U, C, <:AbstractAlgebra.FunctionalMap}, map2::AbstractAlgebra.Map{D, U, <:AbstractAlgebra.FunctionalMap}) where {D, U, C}
      return new{D, C}(domain(map2), codomain(map1), map1, map2)
   end
end

###############################################################################
#
#   MapWithSection
#
###############################################################################

mutable struct MapWithSection{D, C} <: AbstractAlgebra.Map{D, C, MapWithSection}
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

mutable struct MapWithRetraction{D, C} <: AbstractAlgebra.Map{D, C, MapWithRetraction}
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

mutable struct MapCache{D, C, T, De, Ce} <: AbstractAlgebra.Map{D, C, T}
   map::AbstractAlgebra.Map{D, C}
   limit::Int
   image_cache::Dict{De, Ce}

   function MapCache(f::AbstractAlgebra.Map{D, C, T}, limit::Int, enable::Bool) where {D, C, T}
      De = elem_type(D)
      Ce = elem_type(C)
      r = new{D, C, T, De, Ce}(f, limit)
      if enable
         r.image_cache = Dict{De, Ce}()
      end
      return r
   end
end

