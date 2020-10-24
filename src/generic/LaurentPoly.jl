###############################################################################
#
#   LaurentPoly.jl : Generic Laurent polynomials over rings
#
###############################################################################

import AbstractAlgebra: terms_degrees
using AbstractAlgebra: term_degree, degrees_range

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{LaurentPolyWrap{T, PE}}) where {T, PE} =
   LaurentPolyWrapRing{T, parent_type(PE)}

elem_type(::Type{LaurentPolyWrapRing{T, PR}}) where {T, PR} =
   LaurentPolyWrap{T, elem_type(PR)}

parent(p::LaurentPolyWrap) = LaurentPolyWrapRing(parent(p.poly))

base_ring(R::LaurentPolyWrapRing) = base_ring(R.polyring)

var(R::LaurentPolyWrapRing) = var(R.polyring)

symbols(R::LaurentPolyWrapRing) = symbols(R.polyring)

nvars(R::LaurentPolyWrapRing) = nvars(R.polyring)

characteristic(R::LaurentPolyWrapRing) = characteristic(R.polyring)


###############################################################################
#
#   Basic manipulation
#
###############################################################################

terms_degrees(p::LaurentPolyWrap) = p.mindeg .+ (0:degree(p.poly))

"""
    trail_degree(p::LaurentPolyElem)

Return the degree of the term with lowest degree in `p`.
The result is undefined when `p` is null.
"""
function trail_degree(p::LaurentPolyWrap)
   # TODO: implement in terms of trail_degree for polynomials
   first(degrees_range(p))
end

"""
    lead_degree(p::LaurentPolyElem)

Return the degree of the term with highest degree in `p`.
The result is undefined when `p` is null.
"""
lead_degree(p::LaurentPolyWrap) = p.mindeg + degree(p.poly)

coeff(p::LaurentPolyWrap, i::Int) =
   i < p.mindeg ? zero(base_ring(p)) : coeff(p.poly, i - p.mindeg)

function _enable_deg!(p::LaurentPolyWrap, i::Int)
   diff = p.mindeg - i
   if diff > 0
      p.mindeg = i
      p.poly = shift_left(p.poly, diff)
   end
   nothing
end

# the underlying storage is adjusted (increased) to allow setting the coeff
function setcoeff!(p::LaurentPolyWrap, i::Int, a)
   _enable_deg!(p, i)
   setcoeff!(p.poly, i - p.mindeg, a)
   p
end

iszero(p::LaurentPolyWrap) = iszero(p.poly)

isone(p::LaurentPolyWrap) = ismonomial(p, 0)

zero(R::LaurentPolyWrapRing) = LaurentPolyWrap(zero(R.polyring))
one(R::LaurentPolyWrapRing) = LaurentPolyWrap(one(R.polyring))

gen(R::LaurentPolyWrapRing) = LaurentPolyWrap(gen(R.polyring))

isgen(p::LaurentPolyWrap) = ismonomial(p, 1)

# only an optimization over the default Base implementation (maybe 1.4 speed-up)
deepcopy_internal(p::LaurentPolyWrap, dict::IdDict) =
   LaurentPolyWrap(deepcopy_internal(p.poly, dict), p.mindeg)


###############################################################################
#
#   Unary and Binary operations
#
###############################################################################

-(p::LaurentPolyWrap) = LaurentPolyWrap(-p.poly, p.mindeg)

function +(p::LaurentPolyWrap, q::LaurentPolyWrap)
   if p.mindeg > q.mindeg
      p, q = q, p
   end
   p_, q_ = p.poly, q.poly
   if p.mindeg < q.mindeg
      q_ = shift_left(q_, q.mindeg - p.mindeg)
   end
   LaurentPolyWrap(p_ + q_, p.mindeg)
end

-(p::LaurentPolyWrap, q::LaurentPolyWrap) = p + (-q) # TODO: optimize

*(p::LaurentPolyWrap{T}, q::LaurentPolyWrap{T}) where {T} = LaurentPolyWrap(p.poly * q.poly, p.mindeg + q.mindeg)

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

*(p::LaurentPolyWrap{T}, a::T) where {T<:RingElem} = LaurentPolyWrap(p.poly * a, p.mindeg)
*(a::T, p::LaurentPolyWrap{T}) where {T<:RingElem} = p * a

*(p::LaurentPolyWrap, a::Union{Integer,Rational,AbstractFloat}) = LaurentPolyWrap(p.poly * a, p.mindeg)
*(a::Union{Integer,Rational,AbstractFloat}, p::LaurentPolyWrap) = p * a

+(p::LaurentPolyWrap, a::RingElement) = p + LaurentPolyWrap(one(p.poly) * a)
+(a::RingElement, p::LaurentPolyWrap) = p + a

###############################################################################
#
#   Powering
#
###############################################################################

function ^(p::LaurentPolyWrap, e::Integer)
   if e >= 0
      LaurentPolyWrap(p.poly^e, p.mindeg * e)
   else
      # p must be a term, whose coeff is invertible
      deg = term_degree(p)
      c = coeff(p, deg)
      # the following is to allow x^-3 even if 1^-3 is failing
      c = isone(c) ? c : c^e
      LaurentPolyWrap(c * one(p.poly), deg * e)
   end
end


###############################################################################
#
#   Evaluation
#
###############################################################################

function evaluate(p::LaurentPolyWrap, b::RingElement)
   z = evaluate(p.poly, b)
   s = b^p.mindeg
   s * z
end


###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(p::LaurentPolyWrap)
   q = zero!(p.poly)
   if q !== p.poly
      LaurentPolyWrap(q, 0)
   else
      p.mindeg = 0
      p
   end
end

function mul!(c::LaurentPolyWrap{T}, a::LaurentPolyWrap{T}, b::LaurentPolyWrap{T}) where T
   am = a.mindeg
   bm = b.mindeg
   d = mul!(c.poly, a.poly, b.poly)
   if d === c.poly
      c.mindeg = am + bm
      c
   else
      LaurentPolyWrap(d, am + bm)
   end
end

function addeq!(c::LaurentPolyWrap{T}, a::LaurentPolyWrap{T}) where T
   # TODO: optimize (together with +)
   d = c + a
   c.poly = d.poly
   c.mindeg = d.mindeg
   c
end

function add!(c::LaurentPolyWrap{T}, a::LaurentPolyWrap{T}, b::LaurentPolyWrap{T}) where T
   # TODO: optimize
   d = a + b
   c.poly = d.poly
   c.mindeg = d.mindeg
   c
end

###############################################################################
#
#   Shifting
#
###############################################################################

# return a copy of `f` whose underlying poly has a constant term
# (this maximizes the .mindeg field)
function canonicalize(f::LaurentPolyWrap)
   td = trail_degree(f)
   tdp = td - f.mindeg # trail degree for f.poly
   LaurentPolyWrap(shift_right(f.poly, tdp), td)
end

function shift_left(f::LaurentPolyWrap, n::Integer)
   n < 0 && throw(DomainError(n, "n must be >= 0"))
   f = canonicalize(f) # this ensures the underlying polynomial is copied
   LaurentPolyWrap(f.poly, f.mindeg + n)
end

function shift_right(f::LaurentPolyWrap, n::Integer)
   n < 0 && throw(DomainError(n, "n must be >= 0"))
   f = canonicalize(f) # this ensures the underlying polynomial is copied
   LaurentPolyWrap(f.poly, f.mindeg - n)
end

###############################################################################
#
#   Random elements
#
###############################################################################

RandomExtensions.maketype(S::LaurentPolyWrapRing, _, _) = elem_type(S)

function RandomExtensions.make(S::LaurentPolyWrapRing, v1, vs...)
   R = S.polyring
   if length(vs) == 1 && vs[1] isa Integer && elem_type(R) == Random.gentype(v1)
     Make(S, v1, vs[1]) # forward to default Make constructor
   else
      degrees_range = v1
      m = minimum(degrees_range)
      degrees_range = degrees_range .- m
      make(S, make(R, degrees_range, vs...), m)
   end
end

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make3{<:LaurentPolyWrap, <:LaurentPolyWrapRing}})
   v, m = sp[][2:end]
   LaurentPolyWrap(rand(rng, v), m)
end

rand(rng::AbstractRNG, S::LaurentPolyWrapRing, degrees_range, v...) =
   rand(rng, make(S, degrees_range, v...))

rand(S::LaurentPolyWrapRing, degrees_range, v...) =
   rand(GLOBAL_RNG, S, degrees_range, v...)


###############################################################################
#
#   Promotion rules
#
###############################################################################

# TODO: add tests

promote_rule(::Type{L}, ::Type{L}) where {L <: LaurentPolyWrap} = L

function promote_rule(::Type{L}, ::Type{U}) where
       {T<:RingElement, L <: LaurentPolyWrap{T}, U <: RingElement}
   promote_rule(T, U) == T ? L : Union{}
end


################################################################################
#
#  map_coeffs
#
################################################################################

map_coeffs(f, p::LaurentPolyWrap) = LaurentPolyWrap(map_coeffs(f, p.poly), p.mindeg)

###############################################################################
#
#   Parent object call overload
#
###############################################################################

(R::LaurentPolyWrapRing)(b::RingElement) = LaurentPolyWrap(R.polyring(b))

(R::LaurentPolyWrapRing)() = LaurentPolyWrap(R.polyring())

function (R::LaurentPolyWrapRing)(p::LaurentPolyWrap)
   parent(p) == R ? p :
                    LaurentPolyWrap(R.polyring(p.poly), p.mindeg)
end

###############################################################################
#
#   LaurentPolynomialRing constructor
#
###############################################################################

@doc doc"""
    LaurentPolynomialRing(R::AbstractAlgebra.Ring, s::AbstractString)

Given a base ring `R` and string `s` specifying how the generator (variable)
should be printed, return a tuple `S, x` representing the new Laurent polynomial
ring $S = R[x, 1/x]$ and the generator $x$ of the ring. The parent
object `S` will depend only on `R` and `x`.
"""
function LaurentPolynomialRing(R::AbstractAlgebra.Ring, s::AbstractString)
   P, x = AbstractAlgebra.PolynomialRing(R, s)
   LaurentPolyWrapRing(P), LaurentPolyWrap(x)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function AbstractAlgebra.expressify(y::LaurentPolyWrap, S = var(parent(y));
                                                        context = nothing)
   x = y.poly
   mindeg = y.mindeg
   len = length(x)
   sum = Expr(:call, :+)
   for i in 1:len
      c = coeff(x, len - i)
      k = len - i + mindeg
      if !iszero(c)
         if k == 0
            xk = 1
         elseif k == 1
            xk = S
         else
            xk = Expr(:call, :^, S, k)
         end
         if isone(c)
            push!(sum.args, Expr(:call, :*, xk))
         else
            push!(sum.args, Expr(:call, :*, expressify(c, context = context), xk))
         end
      end
   end
   return sum
end

function Base.show(io::IO, ::MIME"text/plain", a::LaurentPolyWrap)
  print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

function Base.show(io::IO, a::LaurentPolyWrap)
  print(io, AbstractAlgebra.obj_to_string(a, context = io))
end
