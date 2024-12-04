###############################################################################
#
#   Residue.jl : generic residue rings (modulo a principal ideal)
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{EuclideanRingResidueRingElem{T}}) where T <: RingElement = EuclideanRingResidueRing{T}

elem_type(::Type{EuclideanRingResidueRing{T}}) where {T <: RingElement} = EuclideanRingResidueRingElem{T}

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{EuclideanRingResidueRingElem{T}}, ::Type{EuclideanRingResidueRingElem{T}}) where T <: RingElement = EuclideanRingResidueRingElem{T}

function promote_rule(::Type{EuclideanRingResidueRingElem{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? EuclideanRingResidueRingElem{T} : Union{}
end

###############################################################################
#
#   Parent object call overloading
#
###############################################################################

function (a::EuclideanRingResidueRing{T})(b::RingElement) where {T <: RingElement}
   return a(base_ring(a)(b))
end

function (a::EuclideanRingResidueRing{T})() where {T <: RingElement}
   z = EuclideanRingResidueRingElem{T}(zero(base_ring(a)))
   z.parent = a
   return z
end

function (a::EuclideanRingResidueRing{T})(b::Integer) where {T <: RingElement}
   z = EuclideanRingResidueRingElem{T}(mod(base_ring(a)(b), modulus(a)))
   z.parent = a
   return z
end

function (a::EuclideanRingResidueRing{T})(b::T) where {T <: RingElem}
   base_ring(a) !== parent(b) && error("Operation on incompatible objects")
   z = EuclideanRingResidueRingElem{T}(mod(b, modulus(a)))
   z.parent = a
   return z
end

function (a::EuclideanRingResidueRing{T})(b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}
   a !== parent(b) && error("Operation on incompatible objects")
   return b
end

################################################################################
#
#  Map
#
################################################################################

domain(f::EuclideanRingResidueMap) = f.domain

codomain(f::EuclideanRingResidueMap) = f.codomain

function image(f::EuclideanRingResidueMap, a)
  parent(a) !== domain(f) && error("Not an element of the domain")
  return codomain(f)(a)
end

(f::EuclideanRingResidueMap)(a) = image(f, a)

function preimage(f::EuclideanRingResidueMap, a)
  parent(a) != codomain(f) && error("Not an element of the codomain")
  return lift(a)
end

###############################################################################
#
#   Some functions for residue rings of polynomial rings
#
###############################################################################

function gen(R::Union{EuclideanRingResidueRing{T}, EuclideanRingResidueField{T}}) where {T<:PolyRingElem}
   return R(gen(base_ring(R)))
end

# TODO: the names `gen` (algebra generator) and `gens` (module generators) are
# very unfortunate
function gens(R::Union{EuclideanRingResidueRing{T}, EuclideanRingResidueField{T}}) where {T<:PolyRingElem} ## probably needs more cases
   ## as the other residue functions
   g = gen(R)
   r = Vector{typeof(g)}()
   push!(r, one(R))
   if degree(modulus(R)) == 1
      return r
   end
   push!(r, g)
   for i = 2:degree(modulus(R))-1
      push!(r, r[end] * g)
   end
   return r
end

function characteristic(R::Union{EuclideanRingResidueRing{T}, EuclideanRingResidueField{T}}) where {T<:PolyRingElem}
   return characteristic(base_ring(base_ring(R)))
end

function size(R::Union{EuclideanRingResidueRing{T}, EuclideanRingResidueField{T}}) where {T<:PolyRingElem}
   return size(base_ring(base_ring(R)))^degree(modulus(R))
end

function rand(R::Union{EuclideanRingResidueRing{T}, EuclideanRingResidueField{T}}) where {T<:PolyRingElem}
   r = rand(base_ring(base_ring(R)))
   g = gen(R)
   for i = 1:degree(modulus(R))
      r = r * g + rand(base_ring(base_ring(R)))
   end
   return r
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(a::T) where {T <: EuclideanRingResidueRingElem}
   a.data = zero!(a.data)
   return a
end

function mul!(c::T, a::T, b::T) where {T <: EuclideanRingResidueRingElem}
   c.data = mod(data(a)*data(b), modulus(a))
   return c
end

function add!(c::T, a::T, b::T) where {T <: EuclideanRingResidueRingElem}
   c.data = mod(data(a) + data(b), modulus(a))
   return c
end

euclid(n::EuclideanRingResidueRingElem) = degree(gcd(data(n), modulus(n)))

#horrible - and copied from fmpz_mod
#don't know how to seriously simplify it
#maybe a durect gcdx should be added as well
function Base.divrem(n::T, m::T) where {T <: EuclideanRingResidueRingElem}
  @assert !iszero(m)
  R = parent(n)
  e = euclid(m)
  if iszero(e)
    fl, q = divides(n, m)
    @assert fl
    return q, zero(R)
  end

  S = typeof(data(n))

  cp = coprime_base(S[data(n), data(m), modulus(m)])::Vector{S}

  q = Vector{Tuple{S, S}}()
  inf = -1
  for i=1:length(cp)
    is_unit(cp[i]) && continue
    v = valuation(modulus(R), cp[i])::Int
    if v != 0
      pk = cp[i]^v
      nv = iszero(data(n) % pk) ? inf : valuation(data(n) % pk, cp[i])
      mv = iszero(data(m) % pk) ? inf : valuation(data(m) % pk, cp[i])
      if nv < mv
        push!(q, (pk, zero(data(n))))
      else
        if nv === inf
          push!(q, (pk, one(data(n))))
        else
          push!(q, (pk, divexact(data(n) % pk, cp[i]^nv)))
        end
      end
    end
  end
  qq = R(crt([x[2] for x = q], [x[1] for x = q])::S)::T
  #need to adjust the leading term of qq so it cancelles:
  # x * lc(qq)*lc(m) = lc(n)
  # so x = lc(n)/lc(qq)/lc(m)
  if !is_zero(qq)
    qq *= leading_coefficient(data(n)) //leading_coefficient(data(qq)) // leading_coefficient(data(m))
  end
  rr = n-qq*m
  @assert n == qq*m+rr
  @assert rr == 0 || euclid(rr) < e
  return (qq,rr)::Tuple{T, T}
end

