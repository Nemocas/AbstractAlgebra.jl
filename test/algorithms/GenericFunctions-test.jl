module GenericTest

# the following code should be essentially copied and pasted from
# docs/src/ring_interface.md to make sure that it actually works

using AbstractAlgebra, Test

using Random: Random, SamplerTrivial, GLOBAL_RNG
using RandomExtensions: RandomExtensions, Make2, AbstractRNG

import AbstractAlgebra: parent_type, elem_type, base_ring, parent, isdomain_type,
       isexact_type, canonical_unit, isequal, divexact, zero!, mul!, add!, addeq!,
       get_cached!, isunit, characteristic, Ring, RingElem, expressify

import Base: show, +, -, *, ^, ==, inv, isone, iszero, one, zero, rand,
             deepcopy_internal, hash

import Main: test_Ring_interface, test_EuclideanRing_interface, test_elem

mutable struct ConstPolyRing{T <: RingElement} <: Ring
   base_ring::Ring

   function ConstPolyRing{T}(R::Ring, cached::Bool) where T <: RingElement
      return get_cached!(ConstPolyID, R, cached) do
         new{T}(R)
      end::ConstPolyRing{T}
   end
end

const ConstPolyID = AbstractAlgebra.CacheDictType{Ring, ConstPolyRing}()
   
mutable struct ConstPoly{T <: RingElement} <: RingElem
   c::T
   parent::ConstPolyRing{T}

   function ConstPoly{T}(c::T) where T <: RingElement
      return new(c)
   end
end

# Data type and parent object methods

parent_type(::Type{ConstPoly{T}}) where T <: RingElement = ConstPolyRing{T}

elem_type(::Type{ConstPolyRing{T}}) where T <: RingElement = ConstPoly{T}

base_ring(R::ConstPolyRing) = R.base_ring

parent(f::ConstPoly) = f.parent

isdomain_type(::Type{ConstPoly{T}}) where T <: RingElement = isdomain_type(T)

isexact_type(::Type{ConstPoly{T}}) where T <: RingElement = isexact_type(T)

function hash(f::ConstPoly, h::UInt)
   r = 0x65125ab8e0cd44ca
   return xor(r, hash(f.c, h))
end

function deepcopy_internal(f::ConstPoly{T}, d::IdDict) where T <: RingElement
   r = ConstPoly{T}(deepcopy_internal(f.c, d))
   r.parent = f.parent # parent should not be deepcopied
   return r
end

# Basic manipulation

zero(R::ConstPolyRing) = R()

one(R::ConstPolyRing) = R(1)

iszero(f::ConstPoly) = iszero(f.c)

isone(f::ConstPoly) = isone(f.c)

isunit(f::ConstPoly) = isunit(f.c)

characteristic(R::ConstPolyRing) = characteristic(base_ring(R))

# Canonical unit

canonical_unit(f::ConstPoly) = canonical_unit(f.c)

# String I/O

function show(io::IO, R::ConstPolyRing)
   print(io, "Constant polynomials over ")
   show(io, base_ring(R))
end

function show(io::IO, f::ConstPoly)
   print(io, f.c)
end

# Expressification (optional)

function expressify(R::ConstPolyRing; context = nothing)
   return Expr(:sequence, Expr(:text, "Constant polynomials over "),
                          expressify(base_ring(R), context = context))
end

function expressify(f::ConstPoly; context = nothing)
   return expressify(f.c, context = context)
end

# Unary operations

function -(f::ConstPoly)
   R = parent(f)
   return R(-f.c)
end

# Binary operations

function +(f::ConstPoly{T}, g::ConstPoly{T}) where T <: RingElement
   parent(f) != parent(g) && error("Incompatible rings")
   R = parent(f)
   return R(f.c + g.c)
end

function -(f::ConstPoly{T}, g::ConstPoly{T}) where T <: RingElement
   parent(f) != parent(g) && error("Incompatible rings")
   R = parent(f)
   return R(f.c - g.c)
end

function *(f::ConstPoly{T}, g::ConstPoly{T}) where T <: RingElement
   parent(f) != parent(g) && error("Incompatible rings")
   R = parent(f)
   return R(f.c*g.c)
end

# Comparison

function ==(f::ConstPoly{T}, g::ConstPoly{T}) where T <: RingElement
   parent(f) != parent(g) && error("Incompatible rings")
   return f.c == g.c
end

function isequal(f::ConstPoly{T}, g::ConstPoly{T}) where T <: RingElement
   parent(f) != parent(g) && error("Incompatible rings")
   return isequal(f.c, g.c)
end

# Powering need not be implemented if * is

# Exact division

function divexact(f::ConstPoly{T}, g::ConstPoly{T}; check::Bool = true) where T <: RingElement
   parent(f) != parent(g) && error("Incompatible rings")
   R = parent(f)
   return R(divexact(f.c, g.c, check = check))
end

# Inverse

function inv(f::ConstPoly)
   R = parent(f)
   return R(AbstractAlgebra.inv(f.c))
end

# Unsafe operators

function zero!(f::ConstPoly)
   f.c = zero(base_ring(parent(f)))
   return f
end

function mul!(f::ConstPoly{T}, g::ConstPoly{T}, h::ConstPoly{T}) where T <: RingElement
   f.c = g.c*h.c
   return f
end

function add!(f::ConstPoly{T}, g::ConstPoly{T}, h::ConstPoly{T}) where T <: RingElement
   f.c = g.c + h.c
   return f
end

function addeq!(f::ConstPoly{T}, g::ConstPoly{T}) where T <: RingElement
   f.c += g.c
   return f
end

# Random generation

RandomExtensions.maketype(R::ConstPolyRing, _) = elem_type(R)

rand(rng::AbstractRNG, sp::SamplerTrivial{<:Make2{ConstPoly,ConstPolyRing}}) =
        sp[][1](rand(rng, sp[][2]))

rand(rng::AbstractRNG, R::ConstPolyRing, n::UnitRange{Int}) = R(rand(rng, n))

rand(R::ConstPolyRing, n::UnitRange{Int}) = rand(Random.GLOBAL_RNG, R, n)

# Promotion rules

promote_rule(::Type{ConstPoly{T}}, ::Type{ConstPoly{T}}) where T <: RingElement = ConstPoly{T}

function promote_rule(::Type{ConstPoly{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? ConstPoly{T} : Union{}
end

# Constructors

function (R::ConstPolyRing{T})() where T <: RingElement
   r = ConstPoly{T}(base_ring(R)(0))
   r.parent = R
   return r
end

function (R::ConstPolyRing{T})(c::Integer) where T <: RingElement
   r = ConstPoly{T}(base_ring(R)(c))
   r.parent = R
   return r
end

# Needed to prevent ambiguity
function (R::ConstPolyRing{T})(c::T) where T <: Integer
   r = ConstPoly{T}(base_ring(R)(c))
   r.parent = R
   return r
end

function (R::ConstPolyRing{T})(c::T) where T <: RingElement
   base_ring(R) != parent(c) && error("Unable to coerce element")
   r = ConstPoly{T}(c)
   r.parent = R
   return r
end

function (R::ConstPolyRing{T})(f::ConstPoly{T}) where T <: RingElement
   R != parent(f) && error("Unable to coerce element")
   return f
end

# Parent constructor

function ConstantPolynomialRing(R::Ring, cached::Bool=true)
   T = elem_type(R)
   return ConstPolyRing{T}(R, cached)
end

# we need only divrem to satsify the
# Euclidean interface

function Base.divrem(a::ConstPoly{elem_type(ZZ)}, b::ConstPoly{elem_type(ZZ)})
   parent(a) != parent(b) && error("Incompatible rings")
   q, r = AbstractAlgebra.divrem(a.c, b.c)
   return parent(a)(q), parent(a)(r)
end

####

function test_elem(R::ConstPolyRing{elem_type(ZZ)})
   n = rand(1:999)
   return R(rand(-n:n))
end

@testset "GenericFunctions.Ring_interface" begin
   test_Ring_interface(ConstantPolynomialRing(ZZ))
end

@testset "GenericFunctions.EuclideanRing_interface" begin
   test_EuclideanRing_interface(ConstantPolynomialRing(ZZ))
end

end
