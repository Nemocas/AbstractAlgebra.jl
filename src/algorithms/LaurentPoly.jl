###############################################################################
#
#   LaurentPoly.jl : Generic algorithms for abstract Laurent polynomials
#
###############################################################################

# required methods without default implementation:
# coeff, setcoeff!, gen

base_ring(p::LaurentPolyElem) = base_ring(parent(p))

###############################################################################
#
#   Basic manipulation
#
###############################################################################

const hashp_seed = UInt === UInt64 ? 0xcdaf0e0b5ade239b : 0x5ade239b

function Base.hash(p::LaurentPolyElem, h::UInt)
   for i in monomials_degrees(p)
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
    monomials_degrees(p::LaurentPolyElem) -> AbstractVector{<:Integer}

> Return a vector containing at least all the degrees of the non-null
> monomials of `p`.
"""
function monomials_degrees end

# return a monomials_degrees vector valid for both polys
function monomials_degrees(p::LaurentPolyElem, q::LaurentPolyElem)
   mdp = monomials_degrees(p)
   mdq = monomials_degrees(q)
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
         UnitRange{T}(extrema(mdq)...)
      end
   else # none is empty
      minp, maxp = extrema(mdp)
      minq, maxq = extrema(mdq)
      min(minp, minq):max(maxp, maxq)
   end
end

# like monomials_degrees, but return a "strict" range, whose min/max
# corresponds to non-zero coeffs
function degrees_range(p::LaurentPolyElem)
   mds = monomials_degrees(p)
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


# return the degree of the unique monomial, throw if not a monomial
function monomial_degree(p::LaurentPolyElem)
   isnull = true
   local deg
   for d in monomials_degrees(p)
      if !iszero(coeff(p, d))
         !isnull && throw(DomainError(p, "not a monomial"))
         deg = d
         isnull = false
      end
   end
   deg
end

function lead(p::LaurentPolyElem)
   dr = degrees_range(p)
   isempty(dr) ? zero(base_ring(p)) : coeff(p, last(dr))
end

function trail(p::LaurentPolyElem)
   dr = degrees_range(p)
   isempty(dr) ? zero(base_ring(p)) : coeff(p, first(dr))
end

gens(R::LaurentPolynomialRing) = [gen(R)]

isgen(p::LaurentPolyElem) = p == gen(parent(p))

# whether p is a (monic) monomial of degree i
function ismonomial(p::LaurentPolyElem, i::Integer; rec::Bool=true)
   rec && error("not implemented")
   dr = degrees_range(p)
   length(dr) == 1 || return false
   dr[] == i || return false
   isone(coeff(p, i))
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
   all(i -> coeff(p, i) == coeff(q, i), monomials_degrees(p, q))


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
