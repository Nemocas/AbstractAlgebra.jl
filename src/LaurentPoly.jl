###############################################################################
#
#   LaurentPoly.jl : Laurent polynomials over rings
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{LaurentPolyRingElem{T, PE, LR}}) where {T, PE, LR} = LR

elem_type(::Type{LaurentPolyRing{T, PR}}) where {T, PR} =
                  LaurentPolyRingElem{T, elem_type(PR), LaurentPolyRing{T, PR}}

parent(p::LaurentPolyRingElem) = p.parent

base_ring_type(::Type{<:LaurentPolyRing{T}}) where {T} = poly_ring_type(T)
base_ring(R::LaurentPolyRing) = R.polyring::base_ring_type(R)

coefficient_ring_type(::Type{LaurentPolyRing{T, PR}}) where {T, PR} = coefficient_ring_type(PR)
coefficient_ring(R::LaurentPolyRing) = coefficient_ring(base_ring(R))

var(R::LaurentPolyRing) = var(base_ring(R))

symbols(R::LaurentPolyRing) = symbols(base_ring(R))

number_of_variables(R::LaurentPolyRing) = 1
number_of_generators(R::LaurentPolyRing) = number_of_variables(R)

characteristic(R::LaurentPolyRing) = characteristic(coefficient_ring(R))
is_known(::typeof(characteristic), R::LaurentPolyRing) = is_known(characteristic, coefficient_ring(R))

is_finite(R::LaurentPolyRing) = is_trivial(R)
is_known(::typeof(is_finite), R::LaurentPolyRing) = is_known(is_trivial, R)

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

"""
    terms_degrees(p::LaurentPolyRingElem) -> AbstractVector{<:Integer}

Return a vector containing at least all the degrees of the non-null
terms of `p`.
"""
terms_degrees(p::LaurentPolyRingElem) = p.mindeg .+ (0:degree(p.poly))

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
   deg = -1
   mds = terms_degrees(p)
   isempty(mds) && throw(DomainError(p, "not a term"))
   for d in mds
      if !iszero(coeff(p, d))
         !isnull && throw(DomainError(p, "not a term"))
         deg = d
         isnull = false
      end
   end
   return deg
end

function leading_coefficient(p::LaurentPolyRingElem)
   dr = degrees_range(p)
   isempty(dr) ? zero(coefficient_ring(p)) : coeff(p, last(dr))
end

function trailing_coefficient(p::LaurentPolyRingElem)
   dr = degrees_range(p)
   isempty(dr) ? zero(coefficient_ring(p)) : coeff(p, first(dr))
end

constant_coefficient(p::LaurentPolyRingElem) = coeff(p, 0)

function is_constant(p::LaurentPolyRingElem)
   dr = degrees_range(p)
   return isempty(dr) || dr == 0:0
end

gens(R::LaurentPolyRing) = [gen(R)]

length(p::LaurentPolyRingElem) = length(degrees_range(p))

is_term(p::LaurentPolyRingElem) = isone(length(p))

# whether p is a (monic) monomial of degree i, non-recursively:
# return true iff `p` has only one non-null coefficient
# (of degree i) at the outer layer, and this coefficient is one
# (as an element of the base ring)
is_monomial(p::LaurentPolyRingElem, i::Integer) = is_term(p) && isone(coeff(p, i))

function is_monomial(p::LaurentPolyRingElem)
   dr = degrees_range(p)
   return isone(length(dr)) && isone(coeff(p, dr[]))
end

function is_monomial_recursive(p::LaurentPolyRingElem)
   dr = degrees_range(p)
   return isone(length(dr)) && is_monomial_recursive(coeff(p, dr[]))
end

function is_unit(f::T) where {T <: LaurentPolyRingElem}
  # **NOTE**  f.poly is not normalized so that the degree 0 coeff is non-zero
  is_trivial(parent(f)) && return true  # coeffs in zero ring
  is_zero(f) && return false
  dr = degrees_range(f)
  is_one(length(dr)) && return is_unit(coeff(f, dr[]))
  is_domain_type(coefficient_ring(parent(f))) && return false
  # For coefficient rings with zero divisors, things are more complex;
  # if someone needs it we can implement some more, just ask for it
  # see also <https://github.com/Nemocas/AbstractAlgebra.jl/issues/2359>
  throw(NotImplementedError(is_unit, f))
end


function is_nilpotent(f::T) where {T <: LaurentPolyRingElem}
  return is_nilpotent(f.poly)
end

"""
    trail_degree(p::LaurentPolyRingElem)

Return the degree of the term with lowest degree in `p`.
The result is undefined when `p` is null.
"""
function trail_degree(p::LaurentPolyRingElem)
   # TODO: implement in terms of trail_degree for polynomials
   first(degrees_range(p))
end

"""
    lead_degree(p::LaurentPolyRingElem)

Return the degree of the term with highest degree in `p`.
The result is undefined when `p` is null.
"""
lead_degree(p::LaurentPolyRingElem) = p.mindeg + degree(p.poly)

coeff(p::LaurentPolyRingElem, i::Int) =
   i < p.mindeg ? zero(coefficient_ring(p)) : coeff(p.poly, i - p.mindeg)

function _enable_deg!(p::LaurentPolyRingElem, i::Int)
   diff = p.mindeg - i
   if diff > 0
      p.mindeg = i
      p.poly = shift_left(p.poly, diff)
   end
   nothing
end

# the underlying storage is adjusted (increased) to allow setting the coeff
function set_coefficient!(p::LaurentPolyRingElem, i::Int, a)
   _enable_deg!(p, i)
   p.poly = set_coefficient!(p.poly, i - p.mindeg, a)
   return p
end

iszero(p::LaurentPolyRingElem) = iszero(p.poly)

isone(p::LaurentPolyRingElem) = is_monomial(p, 0)

zero(R::LaurentPolyRing) = LaurentPolyRingElem(R, zero(base_ring(R)))
one(R::LaurentPolyRing) = LaurentPolyRingElem(R, one(base_ring(R)))

gen(R::LaurentPolyRing) = LaurentPolyRingElem(R, one(base_ring(R)), 1)

is_gen(p::LaurentPolyRingElem) = is_monomial(p, 1)

function deepcopy_internal(p::LaurentPolyRingElem, dict::IdDict)
   return LaurentPolyRingElem(p.parent, deepcopy_internal(p.poly, dict), p.mindeg)
end

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
  print(io, Indent(), "over ", Lowercase(), coefficient_ring(p))
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
    print(terse(io), Lowercase(), coefficient_ring(p))
  end
end

function expressify(y::LaurentPolyRingElem, S = var(parent(y));
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

function Base.show(io::IO, ::MIME"text/plain", a::LaurentPolyRingElem)
   print(io, obj_to_string(a, context = io))
end

function Base.show(io::IO, a::LaurentPolyRingElem)
   print(io, obj_to_string(a, context = io))
end

###############################################################################
#
#   Unary operations
#
###############################################################################

-(p::LaurentPolyRingElem) = LaurentPolyRingElem(parent(p), -p.poly, p.mindeg)

###############################################################################
#
#   Binary operations
#
###############################################################################

function +(p::LaurentPolyRingElem, q::LaurentPolyRingElem)
   if p.mindeg > q.mindeg
      p, q = q, p
   end
   p_, q_ = p.poly, q.poly
   if p.mindeg < q.mindeg
      q_ = shift_left(q_, q.mindeg - p.mindeg)
   end
   LaurentPolyRingElem(parent(p), p_ + q_, p.mindeg)
end

-(p::LaurentPolyRingElem, q::LaurentPolyRingElem) = p + (-q) # TODO: optimize

*(p::LaurentPolyRingElem{T}, q::LaurentPolyRingElem{T}) where {T} = LaurentPolyRingElem(parent(p), p.poly * q.poly, p.mindeg + q.mindeg)

function _remove_gen(a::LaurentPolyRingElem)
   for i in 0:3
      if !iszero(coeff(a.poly, i))
         return i, (i == 0 ? a.poly : shift_right(a.poly, i))
      end
   end
   # we need something finite on zero input
   if iszero(a.poly)
      return (0, a.poly)
   end
   return remove(a.poly, gen(parent(a.poly)))
end

function divexact(a::LaurentPolyRingElem{T}, b::LaurentPolyRingElem{T}; check::Bool = true) where T
   vb, ub = _remove_gen(b)
   f = divexact(a.poly, ub, check = check)
   return LaurentPolyRingElem(parent(a), f, a.mindeg - b.mindeg - vb)
end

function divides(a::LaurentPolyRingElem{T}, b::LaurentPolyRingElem{T}) where T
   vb, ub = _remove_gen(b)
   ok, f = divides(a.poly, ub)
   return ok, LaurentPolyRingElem(parent(a), f, a.mindeg - b.mindeg - vb)
end

function is_divisible_by(a::LaurentPolyRingElem{T}, b::LaurentPolyRingElem{T}) where T
   iszero(b) && return iszero(a)
   vb, ub = _remove_gen(b)
   # should use is_divisible_by here, but it throws on ZZ[x]
   return divides(a.poly, ub)[1]
end

# may throw even when p is invertible
function Base.inv(p::LaurentPolyRingElem)
   v, g = _remove_gen(p)
   return LaurentPolyRingElem(parent(p), inv(g), -p.mindeg-v)
end

is_zero_divisor(p::LaurentPolyRingElem) = is_zero_divisor(p.poly)

function is_zero_divisor_with_annihilator(p::LaurentPolyRingElem)
   f, b = is_zero_divisor_with_annihilator(p.poly)
   return f, LaurentPolyRingElem(parent(p), b)
end

function Base.divrem(p::LaurentPolyRingElem{T}, q::LaurentPolyRingElem{T}) where T
   iszero(q) && error(DivideError())
   iszero(p) && return one(parent(p)), p
   #euc structure: write p (and q) as unit * poly, so remove "x" from p.poly
   # the degree is then the euc function
   vp, up = _remove_gen(p)
   vq, uq = _remove_gen(q)
   qq, rr = divrem(up, uq)
   return LaurentPolyRingElem(parent(p), qq, p.mindeg+vp-q.mindeg-vq),
          LaurentPolyRingElem(parent(p), rr, p.mindeg+vp)
end

function canonical_unit(p::LaurentPolyRingElem)
   iszero(p) && return one(parent(p))
   R = parent(p.poly)
   v, _ = _remove_gen(p)
   return LaurentPolyRingElem(parent(p), R(canonical_unit(p.poly)), p.mindeg + v)
end

function gcd(p::LaurentPolyRingElem{T}, q::LaurentPolyRingElem{T}) where T
   check_parent(p, q)
   if iszero(p)
      return divexact(q, canonical_unit(q))
   elseif iszero(q)
      return divexact(p, canonical_unit(p))
   end
   vp, up = _remove_gen(p)
   vq, uq = _remove_gen(q)
   return LaurentPolyRingElem(parent(p), gcd(up, uq), 0)
end

function gcdx(a::LaurentPolyRingElem{T}, b::LaurentPolyRingElem{T}) where T
   check_parent(a, b)
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
   return LaurentPolyRingElem(R, g, 0), LaurentPolyRingElem(R, s, -a.mindeg - va),
                                    LaurentPolyRingElem(R, t, -b.mindeg - vb)
end

function lcm(p::LaurentPolyRingElem{T}, q::LaurentPolyRingElem{T}) where T
   return LaurentPolyRingElem(parent(p), lcm(p.poly, q.poly), 0)
end

function factor(a::LaurentPolyRingElem)
   R = parent(a)
   va, ua = _remove_gen(a)
   f = factor(ua)
   d = Dict{typeof(a), Int}()
   for (p, e) in f
      d[LaurentPolyRingElem(R, p, 0)] = e
   end
   return Fac(LaurentPolyRingElem(R, unit(f)::typeof(a.poly), a.mindeg + va), d)
end

function factor_squarefree(a::LaurentPolyRingElem)
   R = parent(a)
   va, ua = _remove_gen(a)
   f = factor_squarefree(ua)
   d = Dict{typeof(a), Int}()
   for (p, e) in f
      d[LaurentPolyRingElem(R, p, 0)] = e
   end
   return Fac(LaurentPolyRingElem(R, unit(f)::typeof(a.poly), a.mindeg + va), d)
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

*(p::LaurentPolyRingElem{T}, a::T) where {T<:RingElem} = LaurentPolyRingElem(parent(p), p.poly * a, p.mindeg)
*(a::T, p::LaurentPolyRingElem{T}) where {T<:RingElem} = p * a

*(p::LaurentPolyRingElem, a::JuliaRingElement) = LaurentPolyRingElem(parent(p), p.poly * a, p.mindeg)
*(a::JuliaRingElement, p::LaurentPolyRingElem) = p * a

+(p::LaurentPolyRingElem, a::RingElement) = p + LaurentPolyRingElem(parent(p), one(p.poly) * a)
+(a::RingElement, p::LaurentPolyRingElem) = p + a

###############################################################################
#
#   Powering
#
###############################################################################

function ^(p::LaurentPolyRingElem, e::Integer)
   if e >= 0
      LaurentPolyRingElem(parent(p), p.poly^e, p.mindeg * e)
   else
      # p must be a term, whose coeff is invertible
      deg = term_degree(p)
      c = coeff(p, deg)
      # the following is to allow x^-3 even if 1^-3 is failing
      c = isone(c) ? c : c^e
      LaurentPolyRingElem(parent(p), c * one(p.poly), deg * e)
   end
end


###############################################################################
#
#   Evaluation
#
###############################################################################

function evaluate(p::LaurentPolyRingElem, b::RingElement)
   z = evaluate(p.poly, b)
   s = b^p.mindeg
   s * z
end


###############################################################################
#
#   Derivative
#
###############################################################################

function derivative(a::LaurentPolyRingElem)
   dp = derivative(a.poly)
   n = a.mindeg
   if n != 0
      dp = n*a.poly + shift_left(dp, 1)
      n -= 1
   end
   return LaurentPolyRingElem(parent(a), dp, n)
end


###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(z::LaurentPolyRingElem)
   z.poly = zero!(z.poly)
   z.mindeg = 0
   return z
end

function mul!(z::LaurentPolyRingElem{T}, a::LaurentPolyRingElem{T}, b::LaurentPolyRingElem{T}) where T
   z.mindeg = a.mindeg + b.mindeg
   z.poly = mul!(z.poly, a.poly, b.poly)
   return z
end

function add!(c::LaurentPolyRingElem{T}, a::LaurentPolyRingElem{T}, b::LaurentPolyRingElem{T}) where T
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
function _canonicalize(f::LaurentPolyRingElem)
   td = trail_degree(f)
   tdp = td - f.mindeg # trail degree for f.poly
   return LaurentPolyRingElem(parent(f), shift_right(f.poly, tdp), td)
end

function shift_left(f::LaurentPolyRingElem, n::Integer)
   n < 0 && throw(DomainError(n, "n must be >= 0"))
   f = _canonicalize(f) # this ensures the underlying polynomial is copied
   f.mindeg += n
   return f
end

function shift_right(f::LaurentPolyRingElem, n::Integer)
   n < 0 && throw(DomainError(n, "n must be >= 0"))
   f = _canonicalize(f) # this ensures the underlying polynomial is copied
   f.mindeg -= n
   return f
end

###############################################################################
#
#   Inflation
#
###############################################################################

function inflate(f::LaurentPolyRingElem, shift::Int64, defl::Int64)
   return LaurentPolyRingElem(parent(f), inflate(f.poly, defl), f.mindeg*defl+shift)
end

inflate(f::LaurentPolyRingElem, defl::Int64) = inflate(f, 0, defl)

###############################################################################
#
#   Random elements
#
###############################################################################

RandomExtensions.maketype(S::LaurentPolyRing, _, _) = elem_type(S)

function RandomExtensions.make(S::LaurentPolyRing, v1, vs...)
   R = base_ring(S)
   if length(vs) == 1 && vs[1] isa Integer && elem_type(R) == Random.gentype(v1)
      Make(S, v1, vs[1]) # forward to default Make constructor
   else
      degrees_range = v1
      m = minimum(degrees_range)
      degrees_range = degrees_range .- m
      Make(S, make(R, degrees_range, vs...), m)
   end
end

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make3{<:LaurentPolyRingElem, <:LaurentPolyRing}})
   R, v, m = sp[][1:end]
   LaurentPolyRingElem(R, rand(rng, v), m)
end

rand(rng::AbstractRNG, S::LaurentPolyRing,
     degrees_range::AbstractUnitRange{Int}, v...) =
   rand(rng, make(S, degrees_range, v...))

rand(S::LaurentPolyRing, degrees_range::AbstractUnitRange{Int}, v...) =
   rand(Random.default_rng(), S, degrees_range, v...)

###############################################################################
#
#   Conformance test element generation
#
###############################################################################

function ConformanceTests.generate_element(R::LaurentPolyRing)
   n = rand(0:10)
   if n == 0
      return zero(R)
   else
      m = rand(0:5)
      rand(R, -m:n-m, -99:99)
   end
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

# TODO: add tests

promote_rule(::Type{LaurentPolyRingElem{S, T, V}},
             ::Type{LaurentPolyRingElem{S, T, V}}) where {S, T, V} = LaurentPolyRingElem{S, T, V}

function promote_rule(::Type{LaurentPolyRingElem{S, T, V}}, ::Type{U}) where {S, T, U, V}
   promote_rule(T, U) == T ? LaurentPolyRingElem{S, T, V} : Union{}
end

################################################################################
#
#  Change base ring / map_coefficients
#
################################################################################

function _map(g::T, p::LaurentPolyRingElem, Rx::LaurentPolyRing) where T
   return LaurentPolyRingElem(Rx,
                        _map(g, p.poly, base_ring(Rx)), p.mindeg)
end

function change_base_ring(R::Ring, p::LaurentPolyRingElem; cached::Bool = true,
                 parent::LaurentPolyRing = LaurentPolyRing(
         _change_poly_ring(R, parent(p.poly), cached), cached))

   return _map(R, p, parent)
end

function map_coefficients(g::T, p::LaurentPolyRingElem; cached::Bool = true,
                       parent::LaurentPolyRing = LaurentPolyRing(
                      _make_parent(g, p.poly, cached), cached)) where T
   return _map(g, p, parent)
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

(R::LaurentPolyRing)(b::RingElement) = LaurentPolyRingElem(R, base_ring(R)(b))

(R::LaurentPolyRing)(b::RingElement, d::Int) = LaurentPolyRingElem(R, base_ring(R)(b), d)

(R::LaurentPolyRing)() = LaurentPolyRingElem(R, base_ring(R)())

function (R::LaurentPolyRing)(p::LaurentPolyRingElem)
   parent(p) == R ? p : LaurentPolyRingElem(R, base_ring(R)(p.poly), p.mindeg)
end

###############################################################################
#
#   laurent_polynomial_ring constructor
#
###############################################################################

@doc raw"""
    laurent_polynomial_ring(R::Ring, s::VarName)

Given a base ring `R` and string `s` specifying how the generator (variable)
should be printed, return a tuple `S, x` representing the new Laurent polynomial
ring $S = R[x, 1/x]$ and the generator $x$ of the ring.

# Examples
```julia
julia> R, x = laurent_polynomial_ring(ZZ, :x)
(Univariate Laurent Polynomial Ring in x over Integers, x)

julia> 2x^-3 + x^2
x^2 + 2*x^-3

julia> rand(R, -3:3, -9:9)
-3*x^2 - 8*x + 4 + 3*x^-1 - 6*x^-2 + 9*x^-3
```
"""
function laurent_polynomial_ring(R::Ring, s::VarName; cached::Bool = true)
   @req !is_trivial(R) "Zero rings are currently not supported as coefficient ring."
   P = poly_ring(R, Symbol(s); cached)
   L = LaurentPolyRing(P, cached)
   return L, gen(L)
end
