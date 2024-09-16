module GenericTest

# the following code should be essentially copied and pasted from
# docs/src/ring_interface.md to make sure that it actually works

using AbstractAlgebra, Test

using ..Random: Random, SamplerTrivial, GLOBAL_RNG
using ..RandomExtensions: RandomExtensions, Make2, AbstractRNG

import AbstractAlgebra: Ring
import AbstractAlgebra: RingElem
import AbstractAlgebra: add!
import AbstractAlgebra: base_ring
import AbstractAlgebra: base_ring_type
import AbstractAlgebra: canonical_unit
import AbstractAlgebra: characteristic
import AbstractAlgebra: divexact
import AbstractAlgebra: elem_type
import AbstractAlgebra: expressify
import AbstractAlgebra: get_cached!
import AbstractAlgebra: is_domain_type
import AbstractAlgebra: is_exact_type
import AbstractAlgebra: is_unit
import AbstractAlgebra: isequal
import AbstractAlgebra: mul!
import AbstractAlgebra: parent
import AbstractAlgebra: parent_type
import AbstractAlgebra: zero!

import Base: *
import Base: +
import Base: -
import Base: ==
import Base: ^
import Base: deepcopy_internal
import Base: hash
import Base: inv
import Base: isone
import Base: iszero
import Base: one
import Base: rand
import Base: show
import Base: zero

import ..test_Ring_interface
import ..test_EuclideanRing_interface
import ..test_elem

@attributes mutable struct ConstPolyRing{T <: RingElement} <: Ring
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

base_ring_type(::Type{ConstPolyRing{T}}) where T <: RingElement = parent_type(T)

base_ring(R::ConstPolyRing) = R.base_ring::base_ring_type(R)

parent(f::ConstPoly) = f.parent

is_domain_type(::Type{ConstPoly{T}}) where T <: RingElement = is_domain_type(T)

is_exact_type(::Type{ConstPoly{T}}) where T <: RingElement = is_exact_type(T)

function hash(f::ConstPoly, h::UInt)
   r = 0x65125ab8e0cd44ca
   return xor(r, hash(f.c, h))
end

function deepcopy_internal(f::ConstPoly{T}, dict::IdDict) where T <: RingElement
   r = ConstPoly{T}(deepcopy_internal(f.c, dict))
   r.parent = f.parent # parent should not be deepcopied
   return r
end

# Basic manipulation

zero(R::ConstPolyRing) = R()

one(R::ConstPolyRing) = R(1)

iszero(f::ConstPoly) = iszero(f.c)

isone(f::ConstPoly) = isone(f.c)

is_unit(f::ConstPoly) = is_unit(f.c)

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

# Random generation

RandomExtensions.maketype(R::ConstPolyRing, _) = elem_type(R)

rand(rng::AbstractRNG, sp::SamplerTrivial{<:Make2{ConstPoly,ConstPolyRing}}) =
        sp[][1](rand(rng, sp[][2]))

rand(rng::AbstractRNG, R::ConstPolyRing, n::AbstractUnitRange{Int}) = R(rand(rng, n))

rand(R::ConstPolyRing, n::AbstractUnitRange{Int}) = rand(Random.GLOBAL_RNG, R, n)

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
