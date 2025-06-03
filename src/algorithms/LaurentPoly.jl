###############################################################################
#
#   LaurentPoly.jl : Generic algorithms for abstract Laurent polynomials
#
###############################################################################

# required methods without default implementation:
# coeff, setcoeff!, map_coefficients, gen

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

is_domain_type(::Type{<:LaurentPolyRingElem{T}}) where {T} = is_domain_type(T)

is_exact_type(::Type{<:LaurentPolyRingElem{T}}) where {T} = is_exact_type(T)

###############################################################################
#
#   Basic manipulation
#
###############################################################################

const hashp_seed = UInt === UInt64 ? 0xcdaf0e0b5ade239b : 0x5ade239b

function Base.hash(p::LaurentPolyRingElem, h::UInt)
   for i in terms_degrees(p)
      c = coeff(p, i)
      if !iszero(c)
         h = hash(i, h)
         h = hash(c, h)
      end
   end
   hash(hashp_seed, h)
end

# required implementation
"""
    terms_degrees(p::LaurentPolyRingElem) -> AbstractVector{<:Integer}

Return a vector containing at least all the degrees of the non-null
terms of `p`.
"""
function terms_degrees end

# return a terms_degrees vector valid for both polys
function terms_degrees(p::LaurentPolyRingElem, q::LaurentPolyRingElem)
   mdp = terms_degrees(p)
   mdq = terms_degrees(q)
   T = promote_type(eltype(mdp), eltype(mdq))
   @assert T <: Signed # -one(T) must not wrap around
   if isempty(mdq)
      mdp, mdq = mdq, mdp
   end
   # if only one is empty, it's now mdp
   if isempty(mdp) # extrema requires non-empty
      if isempty(mdq)
         zero(T):-one(T)
      else
         UnitRange{T}(extrema(mdq)...) # mdq not returned for type-stability
      end
   else # none is empty
      minp, maxp = extrema(mdp)
      minq, maxq = extrema(mdq)
      min(minp, minq):max(maxp, maxq)
   end
end

# like terms_degrees, but return a "strict" range, whose min/max
# corresponds to non-zero coeffs, or an empty range for p == 0
function degrees_range(p::LaurentPolyRingElem)
   mds = terms_degrees(p)
   T = eltype(mds)
   minp, maxp = isempty(mds) ? (T(0), T(0)) : extrema(mds)
   while iszero(coeff(p, minp))
      minp += 1
      minp > maxp && return minp:maxp
   end
   while iszero(coeff(p, maxp))
      maxp -= 1
   end
   minp:maxp
end


# return the degree of the unique term, throw if not a term
function term_degree(p::LaurentPolyRingElem)
   isnull = true
   local deg
   mds = terms_degrees(p)
   isempty(mds) && throw(DomainError(p, "not a term"))
   for d in mds
      if !iszero(coeff(p, d))
         !isnull && throw(DomainError(p, "not a term"))
         deg = d
         isnull = false
      end
   end
   deg
end

function leading_coefficient(p::LaurentPolyRingElem)
   dr = degrees_range(p)
   isempty(dr) ? zero(base_ring(p)) : coeff(p, last(dr))
end

function trailing_coefficient(p::LaurentPolyRingElem)
   dr = degrees_range(p)
   isempty(dr) ? zero(base_ring(p)) : coeff(p, first(dr))
end

gens(R::LaurentPolyRing) = [gen(R)]

is_gen(p::LaurentPolyRingElem) = p == gen(parent(p))

# whether p is a (monic) monomial of degree i, non-recursively:
# return true iff `p` has only one non-null coefficient
# (of degree i) at the outer layer, and this coefficient is one
# (as an element of the base ring)
function is_monomial(p::LaurentPolyRingElem, i::Integer)
   dr = degrees_range(p)
   length(dr) == 1 || return false
   dr[] == i || return false
   isone(coeff(p, i))
end

function is_monomial(p::LaurentPolyRingElem)
   dr = degrees_range(p)
   length(dr) == 1 || return false
   isone(coeff(p, dr[]))
end

function is_monomial_recursive(p::LaurentPolyRingElem)
   dr = degrees_range(p)
   length(dr) == 1 || return false
   is_monomial_recursive(coeff(p, dr[]))
end

function is_unit(f::T) where {T <: LaurentPolyRingElem}
  # **NOTE**  f.poly is not normalized so that the degree 0 coeff is non-zero
  is_trivial(parent(f)) && return true  # coeffs in zero ring
  unit_seen = false
  for i in 0:degree(f.poly)
    if is_nilpotent(coeff(f.poly, i))
      continue
    end
    if unit_seen || !is_unit(coeff(f.poly, i))
      return false
    end
    unit_seen = true
  end
  return unit_seen
end


function is_nilpotent(f::T) where {T <: LaurentPolyRingElem}
  return is_nilpotent(f.poly);
end

ConformanceTests._implements(::Type{LaurentPolyRingElem{T}}, ::typeof(is_unit)) where T = _implements(T, is_unit) && _implements(T, is_nilpotent)

ConformanceTests._implements(::Type{LaurentPolyRingElem{T}}, ::typeof(is_nilpotent)) where T = _implements(T, is_nilpotent)


###############################################################################
#
#   Comparisons
#
###############################################################################

==(p::LaurentPolyRingElem, q::LaurentPolyRingElem) =
   all(i -> coeff(p, i) == coeff(q, i), terms_degrees(p, q))


###############################################################################
#
#   Approximation
#
###############################################################################

function Base.isapprox(p::LaurentPolyRingElem, q::LaurentPolyRingElem; atol::Real=sqrt(eps()))
   check_parent(p, q)
   all(terms_degrees(p, q)) do d
      isapprox(coeff(p, d), coeff(q, d); atol=atol)
   end
end

Base.isapprox(p::LaurentPolyRingElem{T}, q::T; atol::Real=sqrt(eps())) where {T} =
   isapprox(p, parent(p)(q); atol=atol)

Base.isapprox(q::T, p::LaurentPolyRingElem{T}; atol::Real=sqrt(eps())) where {T} =
   isapprox(p, q; atol=atol)


################################################################################
#
#  Change base ring
#
################################################################################

change_base_ring(R::Ring, p::LaurentPolyRingElem) = map_coefficients(R, p)

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(x::LaurentPolyRingElem) = canonical_unit(leading_coefficient(x))

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, mime::MIME"text/plain", p::LaurentPolyRing)
  @show_name(io, p)
  @show_special(io, mime, p)
  print(io, "Univariate Laurent polynomial ring in ", var(p))
  println(io)
  io = pretty(io)
  print(io, Indent(), "over ", Lowercase(), base_ring(p))
  print(io, Dedent())
end

function show(io::IO, p::LaurentPolyRing)
  @show_name(io, p)
  @show_special(io, p)
  if is_terse(io)
    print(io, "Univariate Laurent polynomial ring")
  else
    io = pretty(io)
    print(io, "Univariate Laurent polynomial ring in ", var(p), " over ")
    print(terse(io), Lowercase(), base_ring(p))
  end
end

