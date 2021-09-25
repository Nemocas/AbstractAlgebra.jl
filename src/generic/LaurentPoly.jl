###############################################################################
#
#   LaurentPoly.jl : Generic Laurent polynomials over rings
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{LaurentPolyWrap{T, PE, LR}}) where {T, PE, LR} = LR

elem_type(::Type{LaurentPolyWrapRing{T, PR}}) where {T, PR} =
                  LaurentPolyWrap{T, elem_type(PR), LaurentPolyWrapRing{T, PR}}

parent(p::LaurentPolyWrap) = p.parent

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
function set_coefficient!(p::LaurentPolyWrap, i::Int, a)
   _enable_deg!(p, i)
   p.poly = set_coefficient!(p.poly, i - p.mindeg, a)
   return p
end

iszero(p::LaurentPolyWrap) = iszero(p.poly)

isone(p::LaurentPolyWrap) = ismonomial(p, 0)

zero(R::LaurentPolyWrapRing) = LaurentPolyWrap(R, zero(R.polyring))
one(R::LaurentPolyWrapRing) = LaurentPolyWrap(R, one(R.polyring))

gen(R::LaurentPolyWrapRing) = LaurentPolyWrap(R, gen(R.polyring))

isgen(p::LaurentPolyWrap) = ismonomial(p, 1)

function deepcopy_internal(p::LaurentPolyWrap, dict::IdDict)
   return LaurentPolyWrap(p.parent, deepcopy_internal(p.poly, dict), p.mindeg)
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

###############################################################################
#
#   Unary operations
#
###############################################################################

-(p::LaurentPolyWrap) = LaurentPolyWrap(parent(p), -p.poly, p.mindeg)

###############################################################################
#
#   Binary operations
#
###############################################################################

function +(p::LaurentPolyWrap, q::LaurentPolyWrap)
   if p.mindeg > q.mindeg
      p, q = q, p
   end
   p_, q_ = p.poly, q.poly
   if p.mindeg < q.mindeg
      q_ = shift_left(q_, q.mindeg - p.mindeg)
   end
   LaurentPolyWrap(parent(p), p_ + q_, p.mindeg)
end

-(p::LaurentPolyWrap, q::LaurentPolyWrap) = p + (-q) # TODO: optimize

*(p::LaurentPolyWrap{T}, q::LaurentPolyWrap{T}) where {T} = LaurentPolyWrap(parent(p), p.poly * q.poly, p.mindeg + q.mindeg)

function _remove_gen(a::LaurentPolyWrap)
   for i in 0:3
      if !iszero(coeff(a.poly, i))
         return i, (i == 0 ? a.poly : shift_right(a.poly, i))
      end
   end
   return remove(a.poly, gen(parent(a.poly)))
end

function divexact(a::LaurentPolyWrap{T}, b::LaurentPolyWrap{T}; check::Bool = true) where T
   vb, ub = _remove_gen(b)
   f = divexact(a.poly, ub, check = check)
   return LaurentPolyWrap(parent(a), f, a.mindeg - b.mindeg - vb)
end

function divides(a::LaurentPolyWrap{T}, b::LaurentPolyWrap{T}) where T
   vb, ub = _remove_gen(b)
   ok, f = divides(a.poly, ub)
   return ok, LaurentPolyWrap(parent(a), f, a.mindeg - b.mindeg - vb)
end

function isdivisible_by(a::LaurentPolyWrap{T}, b::LaurentPolyWrap{T}) where T
   iszero(b) && return iszero(a)
   vb, ub = _remove_gen(b)
   # should use isdivisible_by here, but it throws on ZZ[x]
   return divides(a.poly, ub)[1]
end

function Base.inv(p::LaurentPolyWrap)
   isunit(p) || throw(NotInvertibleError(p))
   v, g = _remove_gen(p)
   return LaurentPolyWrap(parent(p), inv(g), -p.mindeg-v)
end

function isunit(p::LaurentPolyWrap)
   iszero(p) && return false
   v, g = _remove_gen(p)
   return length(g) < 2
end

function Base.divrem(p::LaurentPolyWrap{T}, q::LaurentPolyWrap{T}) where T
   iszero(q) && error(DivideError())
   iszero(p) && return one(parent(p)), p
   #euc structure: write p (and q) as unit * poly, so remove "x" from p.poly
   # the degree is then the euc function
   vp, up = _remove_gen(p)
   vq, uq = _remove_gen(q)
   qq, rr = divrem(up, uq)
   return LaurentPolyWrap(parent(p), qq, p.mindeg+vp-q.mindeg-vq),
          LaurentPolyWrap(parent(p), rr, p.mindeg+vp)
end

function canonical_unit(p::LaurentPolyWrap)
   iszero(p) && return one(parent(p))
   R = parent(p.poly)
   v, _ = _remove_gen(p)
   return LaurentPolyWrap(parent(p), R(canonical_unit(p.poly)), p.mindeg + v)
end

function gcd(p::LaurentPolyWrap{T}, q::LaurentPolyWrap{T}) where T
   parent(p) == parent(q) || error("Incompatible parents")
   if iszero(p)
      return divexact(q, canonical_unit(q))
   elseif iszero(q)
      return divexact(p, canonical_unit(p))
   end
   vp, up = _remove_gen(p)
   vq, uq = _remove_gen(q)
   return LaurentPolyWrap(parent(p), gcd(up, uq), 0)
end

function gcdx(a::LaurentPolyWrap{T}, b::LaurentPolyWrap{T}) where T
   parent(a) == parent(b) || error("Incompatible parents")
   R = parent(a)
   if iszero(a)
      if iszero(b)
         return zero(R), zero(R), zero(R)
      else
         t = canonical_unit(b)
         return divexact(b, t), zero(R), inv(t)
      end
   elseif iszero(b)
      t = canonical_unit(a)
      return divexact(a, t), inv(t), zero(R)
   end
   va, ua = _remove_gen(a)
   vb, ub = _remove_gen(b)
   g, s, t = gcdx(ua, ub)
   return LaurentPolyWrap(R, g, 0), LaurentPolyWrap(R, s, -a.mindeg - va),
                                    LaurentPolyWrap(R, t, -b.mindeg - vb)
end

function lcm(p::LaurentPolyWrap{T}, q::LaurentPolyWrap{T}) where T
   return LaurentPolyWrap(parent(p), lcm(p.poly, q.poly), 0)
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

*(p::LaurentPolyWrap{T}, a::T) where {T<:RingElem} = LaurentPolyWrap(parent(p), p.poly * a, p.mindeg)
*(a::T, p::LaurentPolyWrap{T}) where {T<:RingElem} = p * a

*(p::LaurentPolyWrap, a::Union{Integer,Rational,AbstractFloat}) = LaurentPolyWrap(parent(p), p.poly * a, p.mindeg)
*(a::Union{Integer,Rational,AbstractFloat}, p::LaurentPolyWrap) = p * a

+(p::LaurentPolyWrap, a::RingElement) = p + LaurentPolyWrap(parent(p), one(p.poly) * a)
+(a::RingElement, p::LaurentPolyWrap) = p + a

###############################################################################
#
#   Powering
#
###############################################################################

function ^(p::LaurentPolyWrap, e::Integer)
   if e >= 0
      LaurentPolyWrap(parent(p), p.poly^e, p.mindeg * e)
   else
      # p must be a term, whose coeff is invertible
      deg = term_degree(p)
      c = coeff(p, deg)
      # the following is to allow x^-3 even if 1^-3 is failing
      c = isone(c) ? c : c^e
      LaurentPolyWrap(parent(p), c * one(p.poly), deg * e)
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
#   Derivative
#
###############################################################################

function derivative(a::LaurentPolyWrap)
   dp = derivative(a.poly)
   n = a.mindeg
   if n != 0
      dp = n*a.poly + shift_left(dp, 1)
      n -= 1
   end
   return LaurentPolyWrap(parent(a), dp, n)
end


###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(z::LaurentPolyWrap)
   z.poly = zero!(z.poly)
   z.mindeg = 0
   return z
end

function mul!(z::LaurentPolyWrap{T}, a::LaurentPolyWrap{T}, b::LaurentPolyWrap{T}) where T
   z.mindeg = a.mindeg + b.mindeg
   z.poly = mul!(z.poly, a.poly, b.poly)
   return z
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
   return LaurentPolyWrap(parent(f), shift_right(f.poly, tdp), td)
end

function shift_left(f::LaurentPolyWrap, n::Integer)
   n < 0 && throw(DomainError(n, "n must be >= 0"))
   f = canonicalize(f) # this ensures the underlying polynomial is copied
   f.mindeg += n
   return f
end

function shift_right(f::LaurentPolyWrap, n::Integer)
   n < 0 && throw(DomainError(n, "n must be >= 0"))
   f = canonicalize(f) # this ensures the underlying polynomial is copied
   f.mindeg -= n
   return f
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
   R, v, m = sp[][1:end]
   LaurentPolyWrap(R, rand(rng, v), m)
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

function promote_rule(::Type{LaurentPolyWrap{S, T, V}}, ::Type{U}) where {S, T, U, V}
   promote_rule(T, U) == T ? LaurentPolyWrap{S, T, V} : Union{}
end

################################################################################
#
#  Change base ring / map_coefficients
#
################################################################################

function AbstractAlgebra._map(g, p::LaurentPolyWrap, Rx::LaurentPolyWrapRing)
   return LaurentPolyWrap(Rx,
                        AbstractAlgebra._map(g, p.poly, Rx.polyring), p.mindeg)
end

function change_base_ring(R::Ring, p::LaurentPolyWrap; cached::Bool = true,
                 parent::LaurentPolyWrapRing = LaurentPolyWrapRing(
         AbstractAlgebra._change_poly_ring(R, parent(p.poly), cached), cached))

   return AbstractAlgebra._map(R, p, parent)
end

function map_coefficients(g, p::LaurentPolyWrap; cached::Bool = true,
                       parent::LaurentPolyWrapRing = LaurentPolyWrapRing(
                      AbstractAlgebra._make_parent(g, p.poly, cached), cached))
   return AbstractAlgebra._map(g, p, parent)
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

(R::LaurentPolyWrapRing)(b::RingElement) = LaurentPolyWrap(R, R.polyring(b))

(R::LaurentPolyWrapRing)(b::RingElement, d::Int) = LaurentPolyWrap(R, R.polyring(b), d)

(R::LaurentPolyWrapRing)() = LaurentPolyWrap(R, R.polyring())

function (R::LaurentPolyWrapRing)(p::LaurentPolyWrap)
   parent(p) == R ? p : LaurentPolyWrap(R, R.polyring(p.poly), p.mindeg)
end

###############################################################################
#
#   LaurentPolynomialRing constructor
#
###############################################################################

function LaurentPolynomialRing(R::AbstractAlgebra.Ring, s::Symbol; cached::Bool = true)
   P, x = AbstractAlgebra.PolynomialRing(R, s, cached = cached)
   R = LaurentPolyWrapRing(P, cached)
   R, LaurentPolyWrap(R, x)
end

