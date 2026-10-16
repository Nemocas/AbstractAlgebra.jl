###############################################################################
#
#   LaurentMPoly.jl: Multivariate Laurent polynomials
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{LaurentMPolyRingElem{T, PE, LR}}) where {T, PE, LR} = LR

elem_type(::Type{LaurentMPolyRing{T, PR}}) where {T, PR} =
                LaurentMPolyRingElem{T, elem_type(PR), LaurentMPolyRing{T, PR}}

parent(p::LaurentMPolyRingElem) = p.parent

base_ring_type(::Type{<:LaurentMPolyRing{T}}) where {T} = mpoly_ring_type(T)
base_ring(R::LaurentMPolyRing) = R.mpolyring::base_ring_type(R)

coefficient_ring_type(::Type{LaurentMPolyRing{T, PR}}) where {T, PR} = coefficient_ring_type(PR)
coefficient_ring(R::LaurentMPolyRing) = coefficient_ring(base_ring(R))

symbols(R::LaurentMPolyRing) = symbols(base_ring(R))

number_of_variables(R::LaurentMPolyRing) = number_of_variables(base_ring(R))
number_of_generators(R::LaurentMPolyRing) = number_of_variables(base_ring(R))

is_univariate(R::LaurentMPolyRing) = is_univariate(base_ring(R))

laurent_mpoly_type(::Type{T}) where T<:RingElement = LaurentMPolyRingElem{T, mpoly_type(T), LaurentMPolyRing{T, mpoly_ring_type(T)}}
laurent_mpoly_type(::Type{S}) where S<:Ring = laurent_mpoly_type(elem_type(S))
laurent_mpoly_type(x) = laurent_mpoly_type(typeof(x)) # to stop this method from eternally recursing on itself, we better add ...
laurent_mpoly_type(::Type{T}) where T = throw(ArgumentError("Type `$T` must be subtype of `RingElement`."))
laurent_mpoly_type(T::Type{Union{}}) = throw(MethodError(laurent_mpoly_type, (T,)))

laurent_mpoly_ring_type(x) = parent_type(laurent_mpoly_type(x))

characteristic(R::LaurentMPolyRing) = characteristic(coefficient_ring(R))
is_known(::typeof(characteristic), R::LaurentMPolyRing) = is_known(characteristic, coefficient_ring(R))

is_finite(R::LaurentMPolyRing) = is_trivial(coefficient_ring(R)) || (nvars(R) == 0 && is_finite(coefficient_ring(R)))
is_known(::typeof(is_finite), R::LaurentMPolyRing) =
  is_known(is_trivial, coefficient_ring(R)) &&
    (is_trivial(coefficient_ring(R)) || nvars(R) > 0 || is_known(is_finite, coefficient_ring(R)))

is_domain_type(::Type{<:LaurentMPolyRingElem{T}}) where {T} = is_domain_type(T)
is_exact_type(::Type{<:LaurentMPolyRingElem{T}}) where {T} = is_exact_type(T)

###############################################################################
#
#   String I/O
#
###############################################################################

function expressify(a::LaurentMPolyRingElem, x = symbols(parent(a)); context = nothing)
    sum = Expr(:call, :+)
    n = nvars(parent(a))
    for (c, v) in zip(coefficients(a), exponent_vectors(a))
        prod = Expr(:call, :*)
        if !isone(c)
            push!(prod.args, expressify(c, context = context))
        end
        for i in 1:n
            if isone(v[i])
                push!(prod.args, x[i])
            elseif !iszero(v[i])
                push!(prod.args, Expr(:call, :^, x[i], v[i]))
            end
        end
        push!(sum.args, prod)
    end
    return sum
end

@enable_all_show_via_expressify LaurentMPolyRingElem

function show(io::IO, mime::MIME"text/plain", p::LaurentMPolyRing)
  @show_name(io, p)
  @show_special(io, mime, p)

  print(io, "Multivariate Laurent polynomial ring")
  print(io, " in ", ItemQuantity(nvars(p), "variable"), " ")
  print_abbreviated_list(io, symbols(p))
  println(io)
  io = pretty(io)
  print(io, Indent(), "over ", Lowercase(), coefficient_ring(p))
  print(io, Dedent())
end

function show(io::IO, p::LaurentMPolyRing)
  @show_name(io, p)
  @show_special(io, p)
  if is_terse(io)
    print(io, "Multivariate Laurent polynomial ring")
  else
    io = pretty(io)
    print(io, "Multivariate Laurent polynomial ring in ", ItemQuantity(nvars(p), "variable"))
    print(terse(io), " over ", Lowercase(), coefficient_ring(p))
  end
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function deepcopy_internal(a::LaurentMPolyRingElem, dict::IdDict)
   return LaurentMPolyRingElem(a.parent, deepcopy_internal(a.mpoly, dict),
                                     deepcopy_internal(a.mindegs, dict))
end

function Base.hash(a::LaurentMPolyRingElem, h::UInt)
    (ap, ad) = _normalize(a)
    h = hash(ap, h)
    h = xor(h, 0x9c64b62806a3d51d%UInt)
    return hash(ad, h)
end

function length(a::LaurentMPolyRingElem)
    return length(a.mpoly)
end

function is_term(a::LaurentMPolyRingElem)
    return is_term(a.mpoly)
end

function is_monomial(a::LaurentMPolyRingElem)
    return is_monomial(a.mpoly)
end

function zero(R::LaurentMPolyRing)
    return LaurentMPolyRingElem(R, zero(base_ring(R)))
end

function one(R::LaurentMPolyRing)
    return LaurentMPolyRingElem(R, one(base_ring(R)))
end

function gen(R::LaurentMPolyRing, i::Int)
    mindegs = zeros(Int, nvars(R))
    mindegs[i] = 1
    return LaurentMPolyRingElem(R, one(base_ring(R)), mindegs)
end

function iszero(a::LaurentMPolyRingElem)
    return iszero(a.mpoly)
end

function isone(a::LaurentMPolyRingElem)
    isone(length(a.mpoly)) || return false
    isone(leading_coefficient(a.mpoly)) || return false
    e = leading_exponent_vector(a.mpoly)
    for i in 1:length(e)
        e[i] == -a.mindegs[i] || return false
    end
    return true
end

function _var_index(a::LaurentMPolyRingElem)
    isone(length(a)) || return 0
    isone(leading_coefficient(a)) || return 0
    e = leading_exponent_vector(a.mpoly)
    found = 0
    for i in 1:length(e)
        s = e[i] + a.mindegs[i]
        if isone(s)
            found == 0 || return 0
            found = i
        elseif !iszero(s)
            return 0
        end
    end
    return found
end

function var_index(a::LaurentMPolyRingElem)
    z = _var_index(a)
    iszero(z) && error("Not a variable in var_index")
    return z
end

function is_gen(a::LaurentMPolyRingElem)
    return !iszero(_var_index(a))
end

function var_indices(a::LaurentMPolyRingElem)
    iszero(a) && return Int[]
    lo, _ = _mindegs(a.mpoly)
    hi = degrees(a.mpoly)
    # variable i is absent iff it has exponent 0 in every term
    return [i for i in 1:nvars(parent(a)) if lo[i] != hi[i] || lo[i] != -a.mindegs[i]]
end

function gens(R::LaurentMPolyRing)
    return [gen(R, i) for i in 1:nvars(R)]
end

function is_constant(a::LaurentMPolyRingElem)
    length(a) > 1 && return false
    return iszero(a) || iszero(first(exponent_vectors(a)))
end

vars(p::LaurentMPolyRingElem) = [gen(parent(p), i) for i in var_indices(p)]

is_univariate(p::LaurentMPolyRingElem) = length(var_indices(p)) <= 1

degree(f::T, x::T) where T <: LaurentMPolyRingElem = degree(f, var_index(x))

###############################################################################
#
#   Degrees
#
###############################################################################

# Unlike for polynomials, -1 is a valid degree, so zero has no degree.

function degree(a::LaurentMPolyRingElem, i::Int)
    @req !iszero(a) "Zero polynomial does not have a degree"
    return degree(a.mpoly, i) + a.mindegs[i]
end

function degrees(a::LaurentMPolyRingElem)
    @req !iszero(a) "Zero polynomial does not have a degree"
    return degrees(a.mpoly) .+ a.mindegs
end

function total_degree(a::LaurentMPolyRingElem)
    @req !iszero(a) "Zero polynomial does not have a degree"
    return total_degree(a.mpoly) + sum(a.mindegs)
end

function Base.inv(a::LaurentMPolyRingElem)
    (ap, ad) = _normalize(a)
    return LaurentMPolyRingElem(parent(a), inv(ap), neg!(ad, ad))
end

function is_unit(f::T) where {T <: LaurentMPolyRingElem}
  # **NOTE** f.mpoly is not normalized in any way
  is_trivial(parent(f)) && return true  # coeffs in zero ring
  is_zero(f) && return false
  is_one(length(f.mpoly)) && return is_unit(coeff(f.mpoly,1))
  is_domain_type(coefficient_ring(parent(f))) && return false
  # For coefficient rings with zero divisors, things are more complex;
  # if someone needs it we can implement some more, just ask for it
  # see also <https://github.com/Nemocas/jl/issues/2359>
  throw(NotImplementedError(is_unit, f))
end


function is_nilpotent(f::T) where {T <: LaurentMPolyRingElem}
  return is_nilpotent(f.mpoly)
end


is_zero_divisor(p::LaurentMPolyRingElem) = is_zero_divisor(p.mpoly)

function is_zero_divisor_with_annihilator(p::LaurentMPolyRingElem)
   f, b = is_zero_divisor_with_annihilator(p.mpoly)
   return f, LaurentMPolyRingElem(parent(p), b)
end

###############################################################################
#
#   Arithmetic
#
###############################################################################

function ==(a::LaurentMPolyRingElem, b::LaurentMPolyRingElem)
    check_parent(a, b, false) || return false
    if a.mindegs == b.mindegs
        return a.mpoly == b.mpoly
    end
    g, x, y = _gcdhelper(a, b)
    return x == y
end

function +(a::LaurentMPolyRingElem, b::LaurentMPolyRingElem)
    check_parent(a, b)
    if a.mindegs == b.mindegs
        return LaurentMPolyRingElem(parent(a), a.mpoly + b.mpoly, a.mindegs)
    end
    g, x, y = _gcdhelper(a, b)
    z, d = _normalize(x + y)
    return LaurentMPolyRingElem(parent(a), z, add!(d, d, g))
end

function -(a::LaurentMPolyRingElem, b::LaurentMPolyRingElem)
    check_parent(a, b)
    if a.mindegs == b.mindegs
        return LaurentMPolyRingElem(parent(a), a.mpoly - b.mpoly, a.mindegs)
    end
    g, x, y = _gcdhelper(a, b)
    z, d = _normalize(x - y)
    return LaurentMPolyRingElem(parent(a), z, add!(d, d, g))
end

function -(a::LaurentMPolyRingElem)
    LaurentMPolyRingElem(parent(a), -a.mpoly, a.mindegs)
end

# For ^ and *, if the inputs are normalized and the base ring is a domain,
# then the output will be normalized with this code
function ^(a::LaurentMPolyRingElem, b::Integer)
    # possible promotion of the vector and then conversion back to Vector{Int}
    if b >= 0
        return LaurentMPolyRingElem(parent(a), a.mpoly^b, a.mindegs*b)
    else
        ap, ad = _normalize(a)
        return LaurentMPolyRingElem(parent(a), inv(ap)^-b, ad*b)
    end
end

function *(a::LaurentMPolyRingElem{T}, b::LaurentMPolyRingElem{T}) where {T}
    check_parent(a, b)
    return LaurentMPolyRingElem(parent(a), a.mpoly*b.mpoly, a.mindegs + b.mindegs)
end

function divides(a::LaurentMPolyRingElem, b::LaurentMPolyRingElem)
    check_parent(a, b)
    (bp, bd) = _normalize(b)
    flag, q = divides(a.mpoly, bp)
    return flag, LaurentMPolyRingElem(parent(a), q, sub!(bd, a.mindegs, bd))
end

function divexact(a::LaurentMPolyRingElem, b::LaurentMPolyRingElem; check::Bool=true)
    check_parent(a, b)
    (bp, bd) = _normalize(b)
    q = divexact(a.mpoly, bp, check=check)
    return LaurentMPolyRingElem(parent(a), q, sub!(bd, a.mindegs, bd))
end

function gcd(a::LaurentMPolyRingElem, b::LaurentMPolyRingElem)
    check_parent(a, b)
    ap, ad = _normalize(a)
    bp, bd = _normalize(b)
    return LaurentMPolyRingElem(parent(a), gcd(ap, bp), zero!(ad))
end

function Base.divrem(a::LaurentMPolyRingElem, b::LaurentMPolyRingElem)
    check_parent(a, b)
    error("divrem not implemented for LaurentMPoly")
end

function factor(a::LaurentMPolyRingElem)
   R = parent(a)
   ap, ad = _normalize(a)
   f = factor(ap)
   d = Dict{typeof(a), Int}()
   for (p, e) in f
      d[LaurentMPolyRingElem(R, p)] = e
   end
   return Fac(LaurentMPolyRingElem(R, unit(f)::typeof(a.mpoly), ad), d)
end

function factor_squarefree(a::LaurentMPolyRingElem)
   R = parent(a)
   ap, ad = _normalize(a)
   f = factor_squarefree(ap)
   d = Dict{typeof(a), Int}()
   for (p, e) in f
      d[LaurentMPolyRingElem(R, p)] = e
   end
   return Fac(LaurentMPolyRingElem(R, unit(f)::typeof(a.mpoly), ad), d)
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

function canonical_unit(a::LaurentMPolyRingElem)
    amin, aiszero = _mindegs(a.mpoly)
    aiszero && return one(parent(a))
    return LaurentMPolyRingElem(parent(a), parent(a.mpoly)(canonical_unit(a.mpoly)),
                                       add!(amin, amin, a.mindegs))
end

###############################################################################
#
#   Evaluation
#
###############################################################################

function evaluate(a::LaurentMPolyRingElem, b::Vector)
    length(b) == nvars(parent(a)) || error("Number of variables does not match number of values")
    (ap, ad) = _normalize(a)
    z = evaluate(ap, b)
    for i in 1:nvars(parent(a))
        if !iszero(ad[i])
            z *= b[i]^ad[i]
        end
    end
    return z
end

# The values are coerced into the parent of `a`; a variable occurring with a
# negative exponent needs a unit value.
function evaluate(a::LaurentMPolyRingElem, vars::Vector{Int}, vals::Vector{<:RingElement})
   R = parent(a)
   @req allunique(vars) "Variables not unique"
   @req length(vars) == length(vals) "Number of variables does not match number of values"
   @req all(i -> 1 <= i <= nvars(R), vars) "Variable index not in range"

   subs = gens(R)
   for (i, v) in zip(vars, vals)
      subs[i] = R(v)
   end
   return evaluate(a, subs)
end

function evaluate(a::T, vars::Vector{T}, vals::Vector{<:RingElement}) where T <: LaurentMPolyRingElem
   return evaluate(a, [var_index(x) for x in vars], vals)
end

###############################################################################
#
#   Derivative
#
###############################################################################

# this has a chance of generating non-normalized output from normalized input
function derivative(a::LaurentMPolyRingElem, j::Int)
    z = derivative(a.mpoly, j)
    e = copy(a.mindegs)
    if !iszero(e[j])
        z = gen(parent(a.mpoly), j)*z + e[j]*a.mpoly
        e[j] -= 1
    end
    return LaurentMPolyRingElem(parent(a), z, e)
end

function derivative(a::LaurentMPolyRingElem{T}, x::LaurentMPolyRingElem{T}) where T <: RingElement
   check_parent(a, x)
   return derivative(a, var_index(x))
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(z::LaurentMPolyRingElem)
   z.mpoly = zero(parent(z.mpoly))
   z.mindegs = zeros(Int, nvars(parent(z)))
   return z
end

###############################################################################
#
#   Internal helpers
#
###############################################################################

function zero!(z::Vector{Int})
    for i in 1:length(z)
        z[i] = 0
    end
    return z
end

function min_broadcast!(z::Vector{Int}, a::Vector{Int}, b::Vector{Int})
    for i in 1:length(z)
        z[i] = min(a[i], b[i])
    end
    return z
end

function add!(z::Vector{Int}, a::Vector{Int}, b::Vector{Int})
    for i in 1:length(z)
        z[i] = a[i] + b[i]
    end
    return z
end

function sub!(z::Vector{Int}, a::Vector{Int}, b::Vector{Int})
    for i in 1:length(z)
        z[i] = a[i] - b[i]
    end
    return z
end

function neg!(z::Vector{Int}, a::Vector{Int})
    for i in 1:length(z)
        z[i] = -a[i]
    end
    return z
end

# min broadcasted over all all exponent vectors
# a return of (0, true) indicates the min is infinite
# TODO specialize for implementations
function _mindegs(a::MPolyRingElem)
    d = zeros(Int, nvars(parent(a)))
    first = true
    for e in exponent_vectors(a)
        if first
            d = e
        else
            min_broadcast!(d, d, e)
        end
        first = false
    end
    return d, first
end

# TODO specialize for implementations
function _divexact_by_exponent_vector(a::MPolyRingElem, d::Vector)
    for i in 1:length(d)
        if d[i] > 0
            a = divexact(a, gen(parent(a), i)^d[i])
        elseif d[i] < 0
            a = a * gen(parent(a), i)^-d[i]
        end
    end
    return a
end

# return equivalent members (mpoly, mindegs) where the mpoly is normalized
function _normalize(a::LaurentMPolyRingElem)
    ap, ad = _normalize(a.mpoly)
    return ap, add!(ad, ad, a.mindegs)
end

# A normalized mpoly is either zero or not divisible by any gen
function _normalize(a::MPolyRingElem)
    d, isinf = _mindegs(a)
    return _divexact_by_exponent_vector(a, d), d
end

# pull out a monomial common factor g, leaving mpolys
function _gcdhelper(a::LaurentMPolyRingElem, b::LaurentMPolyRingElem)
    amin, aiszero = _mindegs(a.mpoly)
    amin = add!(amin, amin, a.mindegs)
    if aiszero
        bp, bd = _normalize(b)
        return (bd, a.mpoly, bp)
    end
    bmin, biszero = _mindegs(b.mpoly)
    bmin = add!(bmin, bmin, b.mindegs)
    if biszero
        ap, ad = _normalize(a)
        return (ad, ap, b.mpoly)
    end
    g = min.(amin, bmin)
    amin = sub!(amin, g, a.mindegs)
    bmin = sub!(bmin, g, b.mindegs)
    return (g, _divexact_by_exponent_vector(a.mpoly, amin),
               _divexact_by_exponent_vector(b.mpoly, bmin))
end

###############################################################################
#
#   Iterators
#
###############################################################################

#### coefficients

function leading_coefficient(a::LaurentMPolyRingElem)
    return leading_coefficient(a.mpoly)
end

function coefficients(a::LaurentMPolyRingElem)
    return coefficients(a.mpoly)
end

function trailing_coefficient(a::LaurentMPolyRingElem)
    return trailing_coefficient(a.mpoly)
end

content(a::LaurentMPolyRingElem) = content(a.mpoly)

# the monomial ordering is invariant under multiplication by a monomial, so
# the terms of a.mpoly are in the same order as those of a
tail(a::LaurentMPolyRingElem) = LaurentMPolyRingElem(parent(a), tail(a.mpoly), copy(a.mindegs))

is_homogeneous(a::LaurentMPolyRingElem) = is_homogeneous(a.mpoly)

function constant_coefficient(a::LaurentMPolyRingElem)
    e = -a.mindegs
    any(x -> x < 0, e) && return zero(coefficient_ring(parent(a)))
    return coeff(a.mpoly, e)
end

function coeff(a::LaurentMPolyRingElem, e::Vector{Int})
  return coeff(a.mpoly, e - a.mindegs)
end

#### exponent vectors

function leading_exponent_vector(a::LaurentMPolyRingElem)
    e = leading_exponent_vector(a.mpoly)
    return add!(e, e, a.mindegs)
end

struct LaurentMPolyExponentVectors{T, S}
    poly::T
    it::S
end

function exponent_vectors(a::LaurentMPolyRingElem)
    t = exponent_vectors(a.mpoly)
    return LaurentMPolyExponentVectors{typeof(a), typeof(t)}(a, t)
end

function Base.iterate(a::LaurentMPolyExponentVectors)
    t = Base.iterate(a.it)
    return isnothing(t) ? t : (add!(t[1], t[1], a.poly.mindegs), t[2])
end

function Base.iterate(a::LaurentMPolyExponentVectors, state)
    t = Base.iterate(a.it, state)
    return isnothing(t) ? t : (add!(t[1], t[1], a.poly.mindegs), t[2])
end

function Base.eltype(::Type{LaurentMPolyExponentVectors{T, S}}) where {T, S}
    return Vector{Int}
end

function Base.length(a::LaurentMPolyExponentVectors)
    return length(a.it)
end

#### monomials

function leading_monomial(a::LaurentMPolyRingElem)
    return LaurentMPolyRingElem(parent(a), leading_monomial(a.mpoly), a.mindegs)
end

struct LaurentMPolyMonomials{T, S}
    poly::T
    it::S
end

function monomials(a::LaurentMPolyRingElem)
    t = monomials(a.mpoly)
    return LaurentMPolyMonomials{typeof(a), typeof(t)}(a, t)
end

function Base.iterate(a::LaurentMPolyMonomials)
    t = Base.iterate(a.it)
    return isnothing(t) ? t :
                (LaurentMPolyRingElem(parent(a.poly), t[1], a.poly.mindegs), t[2])
end

function Base.iterate(a::LaurentMPolyMonomials, state)
    t = Base.iterate(a.it, state)
    return isnothing(t) ? t :
                (LaurentMPolyRingElem(parent(a.poly), t[1], a.poly.mindegs), t[2])
end

function Base.eltype(::Type{LaurentMPolyMonomials{T, S}}) where {T, S}
    return T
end

function Base.length(a::LaurentMPolyMonomials)
    return length(a.it)
end

#### terms

function leading_term(a::LaurentMPolyRingElem)
    return LaurentMPolyRingElem(parent(a), leading_term(a.mpoly), a.mindegs)
end

struct LaurentMPolyTerms{T, S}
    poly::T
    it::S
end

function terms(a::LaurentMPolyRingElem)
    t = terms(a.mpoly)
    return LaurentMPolyTerms{typeof(a), typeof(t)}(a, t)
end

function Base.iterate(a::LaurentMPolyTerms)
    t = Base.iterate(a.it)
    return isnothing(t) ? t :
                 (LaurentMPolyRingElem(parent(a.poly), t[1], a.poly.mindegs), t[2])
end

function Base.iterate(a::LaurentMPolyTerms, state)
    t = Base.iterate(a.it, state)
    return isnothing(t) ? t :
                 (LaurentMPolyRingElem(parent(a.poly), t[1], a.poly.mindegs), t[2])
end

function Base.eltype(::Type{LaurentMPolyTerms{T, S}}) where {T, S}
    return T
end

function Base.length(a::LaurentMPolyTerms)
    return length(a.it)
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::LaurentMPolyRing{T})() where T <: RingElement
    return zero(a)
end

function (a::LaurentMPolyRing{T})(b::RingElement) where T <: RingElement
   return LaurentMPolyRingElem(a, base_ring(a)(b))
end

function (a::LaurentMPolyRing{T})(b::MPolyRingElem{T}) where T <: RingElement
   parent(b) == base_ring(a) || error("Unable to coerce polynomial")
   return LaurentMPolyRingElem(a, b)
end

function (a::LaurentMPolyRing{T})(b::LaurentMPolyRingElem{T}) where T <: RingElement
   parent(b) == a || error("Unable to coerce polynomial")
   return b
end

function (a::LaurentMPolyRing{T})(b::Vector{T}, e::Vector{Vector{Int}}) where T <: RingElement
    if isempty(e)
        return zero(a)
    end
    n = nvars(a)
    m = copy(e[1])
    for i in 1:length(e)
        length(e[i]) == n || error("Exponent vector $i has length $(length(m[i])) (expected $(n))")
        min_broadcast!(m, m, e[i])
    end
    return LaurentMPolyRingElem(a, base_ring(a)(b, map(x -> x - m, e)), m)
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

# rule for R -> R[...]
function promote_rule(::Type{LaurentMPolyRingElem{T, PE, LR}}, ::Type{T}) where {T <: NCRingElement, PE, LR}
  return LaurentMPolyRingElem{T, PE, LR}
end

# If U can be promoted to R[x,y], then U can be promoted to R[x,1/x,y,1/y].
# Handles promotion from R[x,y] to R[x,1/x,y,1/y] as well.
function promote_rule(::Type{LaurentMPolyRingElem{T, PE, LR}}, ::Type{U}) where {T, PE, LR, U <: NCRingElement}
   promote_rule(PE, U) == PE ? LaurentMPolyRingElem{T, PE, LR} : Union{}
end

###############################################################################
#
#  Change base ring / map_coefficients
#
###############################################################################

function _map(g::T, p::LaurentMPolyRingElem, R::LaurentMPolyRing) where T
   return LaurentMPolyRingElem(R, _map(g, p.mpoly, base_ring(R)),
                              p.mindegs)
end

function change_base_ring(
    R::Ring,
    p::LaurentMPolyRingElem;
    cached::Bool = true,
    parent::LaurentMPolyRing = LaurentMPolyRing(
       _change_mpoly_ring(R, parent(p.mpoly), cached), cached))
   return _map(R, p, parent)
end

function map_coefficients(g::T, p::LaurentMPolyRingElem; cached::Bool = true,
                       parent::LaurentMPolyRing = LaurentMPolyRing(
                        _change_mpoly_ring(AbstractAlgebra.parent(g(zero(base_ring(p.mpoly)))), AbstractAlgebra.parent(p.mpoly), cached),
                            cached)) where T
   return _map(g, p, parent)
end

change_coefficient_ring(R::Ring, p::LaurentMPolyRingElem; kwargs...) =
   change_base_ring(R, p; kwargs...)

###############################################################################
#
#   Random elements
#
###############################################################################

function RandomExtensions.maketype(S::LaurentMPolyRing, _, _, _)
    return elem_type(S)
end

function RandomExtensions.make(S::LaurentMPolyRing,
                               term_range::AbstractUnitRange{Int},
                               exp_bound::AbstractUnitRange{Int},
                               vs...)
   R = coefficient_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, term_range, exp_bound, vs[1])
   else
      Make(S, term_range, exp_bound, make(R, vs...))
   end
end

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make4{<:RingElement,
                                         <:LaurentMPolyRing,
                                         <:AbstractUnitRange{Int},
                                         <:AbstractUnitRange{Int}}})
   S, term_range, exp_bound, v = sp[][1:end]
   f = zero(S)
   g = gens(S)
   R = coefficient_ring(S)
   for i = 1:rand(rng, term_range)
      term = one(S)
      for j = 1:length(g)
         term *= g[j]^rand(rng, exp_bound)
      end
      term *= rand(rng, v)
      f += term
   end
   return f
end

function rand(rng::AbstractRNG, S::LaurentMPolyRing,
              term_range::AbstractUnitRange{Int}, exp_bound::AbstractUnitRange{Int}, v...)
   rand(rng, make(S, term_range, exp_bound, v...))
end

function rand(S::LaurentMPolyRing,
              term_range::AbstractUnitRange{Int},
              exp_bound::AbstractUnitRange{Int}, v...)
   rand(Random.default_rng(), S, term_range, exp_bound, v...)
end

###############################################################################
#
#   Conformance test element generation
#
###############################################################################

function ConformanceTests.generate_element(R::LaurentMPolyRing{BigInt})
  n = rand(1:10)
  # R: length between 1 and 9
  # R: exponents between -n and n
  # ZZ: coeffs between -99 and 99
  rand(R, 1:9, -n:n, -99:99)
end

function ConformanceTests.generate_element(R::LaurentMPolyRing{<:ResElem{BigInt}})
  n = rand(1:5)
  # R: length between 1 and 9
  # R: exponents between -n and n
  # ZZ/6ZZ: coeffs between ??? <- TODO
  rand(R, 1:4, -n:n, 1:10)
end

###############################################################################
#
#   Inflation
#
###############################################################################

function inflate(f::LaurentMPolyRingElem, shift::Vector{Int}, defl::Vector{Int})
   return LaurentMPolyRingElem(parent(f), inflate(f.mpoly, defl), f.mindegs.*defl.+shift)
end

function inflate(f::LaurentMPolyRingElem, defl::Vector{Int})
   return LaurentMPolyRingElem(parent(f), inflate(f.mpoly, defl), f.mindegs.*defl)
end

function inflate(f::LaurentMPolyRingElem, vars::Vector{Int}, shifts::Vector{Int}, defls::Vector{Int})
   mindegs = copy(f.mindegs)
   for (var, shift, defl) in zip(vars, shifts, defls)
      mindegs[var] *= defl
      mindegs[var] += shift
   end
   return LaurentMPolyRingElem(parent(f), inflate(f.mpoly, vars, zeros(Int, length(shifts)), defls), mindegs)
end

function deflation(f::LaurentMPolyRingElem)
   shift, defl = deflation(f.mpoly)
   iszero(f) && return shift, defl
   return shift .+ f.mindegs, defl
end

function deflate(f::LaurentMPolyRingElem, shift::Vector{Int}, defl::Vector{Int})
   defl = [iszero(d) ? 1 : d for d in defl]
   B = MPolyBuildCtx(parent(f))
   for (c, e) in zip(coefficients(f), exponent_vectors(f))
      push_term!(B, c, div.(e .- shift, defl))
   end
   return finish(B)
end

deflate(f::LaurentMPolyRingElem, defl::Vector{Int}) = deflate(f, zeros(Int, length(defl)), defl)

###############################################################################
#
#   Square root
#
###############################################################################

# Over a domain, a square s^2 with s = m*r for a monomial m and a polynomial r
# not divisible by any variable normalizes to m^2 and r^2. Over a non-domain
# this fails, e.g. (3 + x)^2 = 6x + x^2 over ZZ/9.
function is_square_with_sqrt(a::LaurentMPolyRingElem)
   R = parent(a)
   is_domain_type(coefficient_ring(R)) || throw(NotImplementedError(:is_square_with_sqrt, a))
   iszero(a) && return true, zero(R)
   ap, ad = _normalize(a)
   all(iseven, ad) || return false, zero(R)
   flag, r = is_square_with_sqrt(ap)
   flag || return false, zero(R)
   return true, LaurentMPolyRingElem(R, r, div.(ad, 2))
end

is_square(a::LaurentMPolyRingElem) = is_square_with_sqrt(a)[1]

function Base.sqrt(a::LaurentMPolyRingElem; check::Bool = true)
   flag, q = is_square_with_sqrt(a)
   check && !flag && error("Not a square in square root")
   return q
end

###############################################################################
#
#   Comparison
#
###############################################################################

function Base.isless(a::LaurentMPolyRingElem, b::LaurentMPolyRingElem)
   check_parent(a, b)
   (!is_monomial(a) || !is_monomial(b)) && error("Not monomials in comparison")
   # shift both to polynomials; the monomial ordering is invariant under this
   m = min.(a.mindegs, b.mindegs)
   return isless(_divexact_by_exponent_vector(a.mpoly, m .- a.mindegs),
                 _divexact_by_exponent_vector(b.mpoly, m .- b.mindegs))
end

###############################################################################
#
#  Universal laurent polynomial ring methods
#
###############################################################################

@doc raw"""
    _upgrade(p::LaurentMPolyRingElem{T}, R::LaurentMPolyRing{T}) where {T}

Return an element of `R` which is obtained from `p` by mapping the $i$-th variable
of `parent(p)` to the $i$-th variable of `R`.
For this to work, `R` needs to have at least as many variables as `parent(p)`.
"""
function _upgrade(p::LaurentMPolyRingElem{T}, R::LaurentMPolyRing{T}) where {T}
   n = nvars(R) - nvars(parent(p))
   n < 0 && error("Too few variables")
   return LaurentMPolyRingElem(R, _upgrade(p.mpoly, R.mpolyring), vcat(p.mindegs, zeros(Int, n)))
end

@doc raw"""
    _add_gens(R::LaurentMPolyRing, varnames::Vector{Symbol})

Return a new uncached multivariate Laurent polynomial ring which has the same properties
as `R` but `varnames` as additional generators.
"""
function _add_gens(R::LaurentMPolyRing, varnames::Vector{Symbol})
   return LaurentMPolyRing(_add_gens(R.mpolyring, varnames), false)
end

###############################################################################
#
#   laurent_polynomial_ring constructor
#
###############################################################################

@doc raw"""
    laurent_polynomial_ring(R::Ring, varnames...; cached::Bool = true)

Given a base ring `R` and variable names `varnames...`, say `:x, :y, :z`, return
a tuple `S, x, y, z` representing the new ring $S = R[x, 1/x, y, 1/y, z, 1/z]$
and the generators $x, y, z$ of the ring.

By default (`cached=true`), the output `S` will be cached, i.e. if
`laurent_polynomial_ring ` is invoked again with the same arguments, the same
(*identical*) ring is returned. Setting `cached` to `false` ensures a distinct
new ring is returned, and will also prevent it from being cached.

For information about the many ways to specify `varnames...` refer to [`polynomial_ring`](@ref) or the
specification in [`AbstractAlgebra.@varnames_interface`](@ref).

# Examples

```jldoctest
julia> S, (x, y) = laurent_polynomial_ring(ZZ, [:x, :y])
(Multivariate Laurent polynomial ring in 2 variables over integers, LaurentMPolyRingElem{BigInt, AbstractAlgebra.Generic.MPoly{BigInt}, LaurentMPolyRing{BigInt, AbstractAlgebra.Generic.MPolyRing{BigInt}}}[x, y])

julia> (x + y)*x^-1
1 + x^-1*y
```
"""
function laurent_polynomial_ring(R::Ring, s::Vector{Symbol}; cached::Bool = true)
   @req !is_trivial(R) "Zero rings are currently not supported as coefficient ring."
   P = poly_ring(R, s; cached)
   L = LaurentMPolyRing(P, cached)
   return L, gens(L)
end

@varnames_interface laurent_polynomial_ring(R::Ring, s)

"""
    @laurent_polynomial_ring(R::Ring, varnames...; cached=true)

Return the ring from [`laurent_polynomial_ring`](@ref laurent_polynomial_ring(::Ring, ::Vector{Symbol}))
and introduce the generators into the current scope.

# Examples

```jldoctest
julia> S = @laurent_polynomial_ring(ZZ, [:u, :v])
Multivariate Laurent polynomial ring in 2 variables u, v
  over integers

julia> u*v^-1
u*v^-1
```
"""
:(@laurent_polynomial_ring)
