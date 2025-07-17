###############################################################################
#
#   LaurentMPoly.jl : Generic multivariate Laurent polynomials over rings
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{LaurentMPolyWrap{T, PE, LR}}) where {T, PE, LR} = LR

elem_type(::Type{LaurentMPolyWrapRing{T, PR}}) where {T, PR} =
                LaurentMPolyWrap{T, elem_type(PR), LaurentMPolyWrapRing{T, PR}}

parent(p::LaurentMPolyWrap) = p.parent

base_ring_type(::Type{<:LaurentMPolyWrapRing{T}}) where {T} = parent_type(T)
base_ring(R::LaurentMPolyWrapRing) = base_ring(R.mpolyring)::base_ring_type(R)

coefficient_ring_type(::Type{LaurentMPolyWrapRing{T, PR}}) where {T, PR} = coefficient_ring_type(PR)
coefficient_ring(R::LaurentMPolyWrapRing) = coefficient_ring(R.mpolyring)

symbols(R::LaurentMPolyWrapRing) = symbols(R.mpolyring)

number_of_variables(R::LaurentMPolyWrapRing) = number_of_variables(R.mpolyring)
number_of_generators(R::LaurentMPolyWrapRing) = number_of_variables(R.mpolyring)

characteristic(R::LaurentMPolyWrapRing) = characteristic(R.mpolyring)

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function deepcopy_internal(a::LaurentMPolyWrap, dict::IdDict)
   return LaurentMPolyWrap(a.parent, deepcopy_internal(a.mpoly, dict),
                                     deepcopy_internal(a.mindegs, dict))
end

function Base.hash(a::LaurentMPolyWrap, h::UInt)
    (ap, ad) = _normalize(a)
    h = hash(ap, h)
    h ^= 0x9c64b62806a3d51d%UInt
    return hash(ad, h)
end

function length(a::LaurentMPolyWrap)
    return length(a.mpoly)
end

function zero(R::LaurentMPolyWrapRing)
    return LaurentMPolyWrap(R, zero(R.mpolyring))
end

function one(R::LaurentMPolyWrapRing)
    return LaurentMPolyWrap(R, one(R.mpolyring))
end

function gen(R::LaurentMPolyWrapRing, i::Int)
    return LaurentMPolyWrap(R, gen(R.mpolyring, i))
end

function iszero(a::LaurentMPolyWrap)
    return iszero(a.mpoly)
end

function isone(a::LaurentMPolyWrap)
    isone(length(a.mpoly)) || return false
    isone(leading_coefficient(a.mpoly)) || return false
    e = leading_exponent_vector(a.mpoly)
    for i in 1:length(e)
        e[i] == -a.mindegs[i] || return false
    end
    return true
end

function _var_index(a::LaurentMPolyWrap)
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

function var_index(a::LaurentMPolyWrap)
    z = _var_index(a)
    iszero(z) && error("Not a variable in var_index")
    return z
end

function is_gen(a::LaurentMPolyWrap)
    return !iszero(_var_index(a))
end

function Base.inv(a::LaurentMPolyWrap)
    (ap, ad) = _normalize(a)
    return LaurentMPolyWrap(parent(a), inv(ap), neg!(ad, ad))
end

function is_unit(f::T) where {T <: LaurentMPolyRingElem}
  # **NOTE** f.mpoly is not normalized in any way
  is_trivial(parent(f)) && return true  # coeffs in zero ring
  unit_seen = false
  for i in 1:length(f.mpoly)
    if is_nilpotent(coeff(f.mpoly, i))
      continue
    end
    if unit_seen || !is_unit(coeff(f.mpoly, i))
      return false
    end
    unit_seen = true
  end
  return unit_seen
end

ConformanceTests._implements(::Type{LaurentMPolyRingElem{T}}, ::typeof(is_unit)) where T = _implements(T, is_unit) && _implements(T, is_nilpotent)

function is_nilpotent(f::T) where {T <: LaurentMPolyRingElem}
  return is_nilpotent(f.mpoly);
end

ConformanceTests._implements(::Type{LaurentMPolyRingElem{T}}, ::typeof(is_nilpotent)) where T = _implements(T, is_nilpotent)

is_zero_divisor(p::LaurentMPolyWrap) = is_zero_divisor(p.mpoly)

ConformanceTests._implements(::Type{LaurentMPolyRingElem{T}}, ::typeof(is_zero_divisor)) where T = _implements(T, is_zero_divisor)

function is_zero_divisor_with_annihilator(p::LaurentMPolyWrap)
   f, b = is_zero_divisor_with_annihilator(p.mpoly)
   return f, LaurentMPolyWrap(parent(p), b)
end

###############################################################################
#
#   Arithmetic
#
###############################################################################

function ==(a::LaurentMPolyWrap, b::LaurentMPolyWrap)
    check_parent(a, b, false) || return false
    if a.mindegs == b.mindegs
        return a.mpoly == b.mpoly
    end
    g, x, y = _gcdhelper(a, b)
    return x == y
end

function +(a::LaurentMPolyWrap, b::LaurentMPolyWrap)
    check_parent(a, b)
    if a.mindegs == b.mindegs
        return LaurentMPolyWrap(parent(a), a.mpoly + b.mpoly, a.mindegs)
    end
    g, x, y = _gcdhelper(a, b)
    z, d = _normalize(x + y)
    return LaurentMPolyWrap(parent(a), z, add!(d, d, g))
end

function -(a::LaurentMPolyWrap, b::LaurentMPolyWrap)
    check_parent(a, b)
    if a.mindegs == b.mindegs
        return LaurentMPolyWrap(parent(a), a.mpoly - b.mpoly, a.mindegs)
    end
    g, x, y = _gcdhelper(a, b)
    z, d = _normalize(x - y)
    return LaurentMPolyWrap(parent(a), z, add!(d, d, g))
end

function -(a::LaurentMPolyWrap)
    LaurentMPolyWrap(parent(a), -a.mpoly, a.mindegs)
end

# For ^ and *, if the inputs are normalized and the base ring is a domain,
# then the output will be normalized with this code
function ^(a::LaurentMPolyWrap, b::Integer)
    # possible promotion of the vector and then conversion back to Vector{Int}
    if b >= 0
        return LaurentMPolyWrap(parent(a), a.mpoly^b, a.mindegs*b)
    else
        ap, ad = _normalize(a)
        return LaurentMPolyWrap(parent(a), inv(ap)^-b, ad*b)
    end
end

function *(a::LaurentMPolyWrap, b::LaurentMPolyWrap)
    check_parent(a, b)
    return LaurentMPolyWrap(parent(a), a.mpoly*b.mpoly, a.mindegs + b.mindegs)
end

function divides(a::LaurentMPolyWrap, b::LaurentMPolyWrap)
    check_parent(a, b)
    (bp, bd) = _normalize(b)
    flag, q = divides(a.mpoly, bp)
    return flag, LaurentMPolyWrap(parent(a), q, sub!(bd, a.mindegs, bd))
end

function divexact(a::LaurentMPolyWrap, b::LaurentMPolyWrap; check::Bool=true)
    check_parent(a, b)
    (bp, bd) = _normalize(b)
    q = divexact(a.mpoly, bp, check=check)
    return LaurentMPolyWrap(parent(a), q, sub!(bd, a.mindegs, bd))
end

function gcd(a::LaurentMPolyWrap, b::LaurentMPolyWrap)
    check_parent(a, b)
    ap, ad = _normalize(a)
    bp, bd = _normalize(b)
    return LaurentMPolyWrap(parent(a), gcd(ap, bp), zero!(ad))
end

function Base.divrem(a::LaurentMPolyWrap, b::LaurentMPolyWrap)
    check_parent(a, b)
    error("divrem not implemented for LaurentMPoly")
end

function factor(a::LaurentMPolyWrap)
   R = parent(a)
   ap, ad = _normalize(a)
   f = factor(ap)
   d = Dict{typeof(a), Int}()
   for (p, e) in f
      d[LaurentMPolyWrap(R, p)] = e
   end
   return Fac(LaurentMPolyWrap(R, unit(f), ad), d)
end

function factor_squarefree(a::LaurentMPolyWrap)
   R = parent(a)
   ap, ad = _normalize(a)
   f = factor_squarefree(ap)
   d = Dict{typeof(a), Int}()
   for (p, e) in f
      d[LaurentMPolyWrap(R, p)] = e
   end
   return Fac(LaurentMPolyWrap(R, unit(f), ad), d)
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

function canonical_unit(a::LaurentMPolyWrap)
    amin, aiszero = _mindegs(a.mpoly)
    aiszero && return one(parent(a))
    return LaurentMPolyWrap(parent(a), parent(a.mpoly)(canonical_unit(a.mpoly)),
                                       add!(amin, amin, a.mindegs))
end

###############################################################################
#
#   Evaluation
#
###############################################################################

function evaluate(a::LaurentMPolyWrap, b::Vector)
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

###############################################################################
#
#   Derivative
#
###############################################################################

# this has a chance of generating non-normalized output from normalized input
function derivative(a::LaurentMPolyWrap, j::Int)
    z = derivative(a.mpoly, j)
    e = copy(a.mindegs)
    if !iszero(e[j])
        z = gen(parent(a.mpoly), j)*z + e[j]*a.mpoly
        e[j] -= 1
    end
    return LaurentMPolyWrap(parent(a), z, e)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(z::LaurentMPolyWrap)
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
function _normalize(a::LaurentMPolyWrap)
    ap, ad = _normalize(a.mpoly)
    return ap, add!(ad, ad, a.mindegs)
end

# A normalized mpoly is either zero or not divisible by any gen
function _normalize(a::MPolyRingElem)
    d, isinf = _mindegs(a)
    return _divexact_by_exponent_vector(a, d), d
end

# pull out a monomial common factor g, leaving mpolys
function _gcdhelper(a::LaurentMPolyWrap, b::LaurentMPolyWrap)
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

function leading_coefficient(a::LaurentMPolyWrap)
    return leading_coefficient(a.mpoly)
end

function coefficients(a::LaurentMPolyWrap)
    return coefficients(a.mpoly)
end

function constant_coefficient(a::LaurentMPolyWrap)
    e = -a.mindegs
    any(x -> x < 0, e) && return zero(coefficient_ring(parent(a)))
    return coeff(a.mpoly, e)
end

#### exponent vectors

function leading_exponent_vector(a::LaurentMPolyWrap)
    e = leading_exponent_vector(a.mpoly)
    return add!(e, e, a.mindegs)
end

struct LaurentMPolyWrapExponentVectors{T, S}
    poly::T
    it::S
end

function exponent_vectors(a::LaurentMPolyWrap)
    t = exponent_vectors(a.mpoly)
    return LaurentMPolyWrapExponentVectors{typeof(a), typeof(t)}(a, t)
end

function Base.iterate(a::LaurentMPolyWrapExponentVectors)
    t = Base.iterate(a.it)
    return isnothing(t) ? t : (add!(t[1], t[1], a.poly.mindegs), t[2])
end

function Base.iterate(a::LaurentMPolyWrapExponentVectors, state)
    t = Base.iterate(a.it, state)
    return isnothing(t) ? t : (add!(t[1], t[1], a.poly.mindegs), t[2])
end

function Base.eltype(::Type{LaurentMPolyWrapExponentVectors{T, S}}) where {T, S}
    return Vector{Int}
end

function Base.length(a::LaurentMPolyWrapExponentVectors)
    return length(a.it)
end

#### monomials

function leading_monomial(a::LaurentMPolyWrap)
    return LaurentMPolyWrap(parent(a), leading_monomial(a.mpoly), a.mindegs)
end

struct LaurentMPolyWrapMonomials{T, S}
    poly::T
    it::S
end

function monomials(a::LaurentMPolyWrap)
    t = monomials(a.mpoly)
    return LaurentMPolyWrapMonomials{typeof(a), typeof(t)}(a, t)
end

function Base.iterate(a::LaurentMPolyWrapMonomials)
    t = Base.iterate(a.it)
    return isnothing(t) ? t :
                (LaurentMPolyWrap(parent(a.poly), t[1], a.poly.mindegs), t[2])
end

function Base.iterate(a::LaurentMPolyWrapMonomials, state)
    t = Base.iterate(a.it, state)
    return isnothing(t) ? t :
                (LaurentMPolyWrap(parent(a.poly), t[1], a.poly.mindegs), t[2])
end

function Base.eltype(::Type{LaurentMPolyWrapMonomials{T, S}}) where {T, S}
    return T
end

function Base.length(a::LaurentMPolyWrapMonomials)
    return length(a.it)
end

#### terms

function leading_term(a::LaurentMPolyWrap)
    return LaurentMPolyWrap(parent(a), leading_term(a.mpoly), a.mindegs)
end

struct LaurentMPolyWrapTerms{T, S}
    poly::T
    it::S
end

function terms(a::LaurentMPolyWrap)
    t = terms(a.mpoly)
    return LaurentMPolyWrapTerms{typeof(a), typeof(t)}(a, t)
end

function Base.iterate(a::LaurentMPolyWrapTerms)
    t = Base.iterate(a.it)
    return isnothing(t) ? t :
                 (LaurentMPolyWrap(parent(a.poly), t[1], a.poly.mindegs), t[2])
end

function Base.iterate(a::LaurentMPolyWrapTerms, state)
    t = Base.iterate(a.it, state)
    return isnothing(t) ? t :
                 (LaurentMPolyWrap(parent(a.poly), t[1], a.poly.mindegs), t[2])
end

function Base.eltype(::Type{LaurentMPolyWrapTerms{T, S}}) where {T, S}
    return T
end

function Base.length(a::LaurentMPolyWrapTerms)
    return length(a.it)
end

###############################################################################
#
#   Build Context
#
###############################################################################

mutable struct LaurentMPolyBuildCtx{T, S}
    coeffs::Vector{T}
    exps::Vector{Vector{Int}}
    parent::S
end

function MPolyBuildCtx(R::AbstractAlgebra.LaurentMPolyRing{T}) where T
    return LaurentMPolyBuildCtx{T, typeof(R)}(T[], Vector{Int}[], R)
end

function push_term!(B::LaurentMPolyBuildCtx{T, S}, c::U, expv::Vector{Int}) where {S, T, U}
    length(expv) == nvars(B.parent) || error("length of exponent vector should match the number of variables")
    push!(B.coeffs, c)
    push!(B.exps, expv)
    return B
end

function finish(B::LaurentMPolyBuildCtx{T, S}) where {T, S}
    res = B.parent(B.coeffs, B.exps)
    B.coeffs = T[]
    B.exps = Vector{Int}[]
    return res
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::LaurentMPolyWrapRing{T})() where T <: RingElement
    return zero(a)
end

function (a::LaurentMPolyWrapRing{T})(b::RingElement) where T <: RingElement
   return LaurentMPolyWrap(a, a.mpolyring(b))
end

function (a::LaurentMPolyWrapRing{T})(b::MPolyRingElem{T}) where T <: RingElement
   parent(b) == a.mpolyring || error("Unable to coerce polynomial")
   return LaurentMPolyWrap(a, b)
end

function (a::LaurentMPolyWrapRing{T})(b::LaurentMPolyWrap{T}) where T <: RingElement
   parent(b) == a || error("Unable to coerce polynomial")
   return b
end

function (a::LaurentMPolyWrapRing{T})(b::Vector{T}, e::Vector{Vector{Int}}) where T <: RingElement
    if isempty(e)
        return zero(a)
    end
    n = nvars(a)
    m = copy(e[1])
    for i in 1:length(e)
        length(e[i]) == n || error("Exponent vector $i has length $(length(m[i])) (expected $(n))")
        min_broadcast!(m, m, e[i])
    end
    return LaurentMPolyWrap(a, a.mpolyring(b, map(x -> x - m, e)), m)
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

# If U can be promoted to R[x,y], then U can be promoted to R[x,1/x,y,1/y].
# Handles promotion from R[x,y] to R[x,1/x,y,1/y] as well.
function promote_rule(::Type{LaurentMPolyWrap{T, PE, LR}}, ::Type{U}) where {T, PE, LR, U}
   promote_rule(PE, U) == PE ? LaurentMPolyWrap{T, PE, LR} : Union{}
end

################################################################################
#
#  Change base ring / map_coefficients
#
################################################################################

function AbstractAlgebra._map(g::T, p::LaurentMPolyWrap, R::LaurentMPolyWrapRing) where T
   return LaurentMPolyWrap(R, AbstractAlgebra._map(g, p.mpoly, R.mpolyring),
                              p.mindegs)
end

function change_base_ring(
    R::Ring,
    p::LaurentMPolyWrap;
    cached::Bool = true,
    parent::LaurentMPolyWrapRing = LaurentMPolyWrapRing(
       AbstractAlgebra._change_mpoly_ring(R, parent(p.mpoly), cached), cached))
   return AbstractAlgebra._map(R, p, parent)
end

function map_coefficients(g::T, p::LaurentMPolyWrap; cached::Bool = true,
                       parent::LaurentMPolyWrapRing = LaurentMPolyWrapRing(
                        AbstractAlgebra._change_mpoly_ring(AbstractAlgebra.parent(g(zero(base_ring(p.mpoly)))), AbstractAlgebra.parent(p.mpoly), cached),
                            cached)) where T
   return AbstractAlgebra._map(g, p, parent)
end

###############################################################################
#
#   laurent_polynomial_ring constructor
#
###############################################################################

function laurent_polynomial_ring(R::AbstractAlgebra.Ring, s::Vector{Symbol}; cached::Bool = true)
   @req !is_trivial(R) "Zero rings are currently not supported as coefficient ring."
   P, x = AbstractAlgebra.polynomial_ring(R, s, cached = cached)
   R = LaurentMPolyWrapRing(P, cached)
   R, map(p -> LaurentMPolyWrap(R, p), x)
end

