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

base_ring(p::LaurentPolyElem) = base_ring(parent(p))

isdomain_type(::Type{<:LaurentPolyElem{T}}) where {T} = isdomain_type(T)

isexact_type(::Type{<:LaurentPolyElem{T}}) where {T} = isexact_type(T)

function check_parent(a::LaurentPolyElem, b::LaurentPolyElem, throw::Bool = true)
   c = parent(a) == parent(b)
   c || !throw || error("incompatible Laurent polynomial rings")
   return c
end


###############################################################################
#
#   Basic manipulation
#
###############################################################################

const hashp_seed = UInt === UInt64 ? 0xcdaf0e0b5ade239b : 0x5ade239b

function Base.hash(p::LaurentPolyElem, h::UInt)
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
    terms_degrees(p::LaurentPolyElem) -> AbstractVector{<:Integer}

Return a vector containing at least all the degrees of the non-null
terms of `p`.
"""
function terms_degrees end

# return a terms_degrees vector valid for both polys
function terms_degrees(p::LaurentPolyElem, q::LaurentPolyElem)
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
function degrees_range(p::LaurentPolyElem)
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
function term_degree(p::LaurentPolyElem)
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

function leading_coefficient(p::LaurentPolyElem)
   dr = degrees_range(p)
   isempty(dr) ? zero(base_ring(p)) : coeff(p, last(dr))
end

function trailing_coefficient(p::LaurentPolyElem)
   dr = degrees_range(p)
   isempty(dr) ? zero(base_ring(p)) : coeff(p, first(dr))
end

gens(R::LaurentPolynomialRing) = [gen(R)]

isgen(p::LaurentPolyElem) = p == gen(parent(p))

# whether p is a (monic) monomial of degree i, non-recursively:
# return true iff `p` has only one non-null coefficient
# (of degree i) at the outer layer, and this coefficient is one
# (as an element of the base ring)
function ismonomial(p::LaurentPolyElem, i::Integer)
   dr = degrees_range(p)
   length(dr) == 1 || return false
   dr[] == i || return false
   isone(coeff(p, i))
end

function ismonomial(p::LaurentPolyElem)
   dr = degrees_range(p)
   length(dr) == 1 || return false
   isone(coeff(p, dr[]))
end

function ismonomial_recursive(p::LaurentPolyElem)
   dr = degrees_range(p)
   length(dr) == 1 || return false
   ismonomial_recursive(coeff(p, dr[]))
end

function isunit(p::LaurentPolyElem)
   dr = degrees_range(p)
   length(dr) == 1 || return false
   isunit(coeff(p, dr[]))
end

###############################################################################
#
#   Comparisons
#
###############################################################################

==(p::LaurentPolyElem, q::LaurentPolyElem) =
   all(i -> coeff(p, i) == coeff(q, i), terms_degrees(p, q))


###############################################################################
#
#   Approximation
#
###############################################################################

function Base.isapprox(p::LaurentPolyElem, q::LaurentPolyElem; atol::Real=sqrt(eps()))
   check_parent(p, q)
   all(terms_degrees(p, q)) do d
      isapprox(coeff(p, d), coeff(q, d); atol=atol)
   end
end

Base.isapprox(p::LaurentPolyElem{T}, q::T; atol::Real=sqrt(eps())) where {T} =
   isapprox(p, parent(p)(q); atol=atol)

Base.isapprox(q::T, p::LaurentPolyElem{T}; atol::Real=sqrt(eps())) where {T} =
   isapprox(p, q; atol=atol)


################################################################################
#
#  Change base ring
#
################################################################################

change_base_ring(R::Ring, p::LaurentPolyElem) = map_coefficients(R, p)

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(x::LaurentPolyElem) = canonical_unit(leading_coefficient(x))

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, p::LaurentPolynomialRing)
   print(io, "Univariate Laurent Polynomial Ring in ")
   print(io, string(var(p)))
   print(io, " over ")
   print(IOContext(io, :compact => true), base_ring(p))
end

