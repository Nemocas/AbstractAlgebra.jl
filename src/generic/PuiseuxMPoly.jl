#################################################################################
#
# Puiseux polynomials
# ===================
#
# A multivariate Puiseux polynomial ring is just a wrapper around a normal
# multivariate Laurent polynomial ring for the normal poly below.
#
# A Puiseux polynomial can be written of the form
#
#   t^{k//d} * (c_0 * t^(0//d) + ... + c_r * t^(r//d)).
#
# We store it as
# - a normal Laurent poly: c_0 * t^(k_0) + ... + c_r * t^(k_r)
# - a scale: d
#
# The representation above is normalized if
# - gcd(d, exponents of poly) = 1
# - c_0 != 0
# Every Puiseux polynomial has a unique normalized representation.
#
# With the exception of normalized! and rescale, all function inputs are assumed
# to be normalized and all function outputs will be normalized.
#
#################################################################################

function puiseux_polynomial_ring_elem(
    Kt::PuiseuxMPolyRing,
    f::LaurentMPolyRingElem,
    d::Int = Int(1);
    skip_normalization::Bool = false
    )
    pf = PuiseuxMPolyRingElem(Kt, f, d)
    if !skip_normalization
        normalize!(pf)
    end
    return pf
end

#################################################################################
#
# Properties
#
#################################################################################

elem_type(::Type{PuiseuxMPolyRing{T}}) where T <: RingElement = PuiseuxMPolyRingElem{T}
parent_type(::Type{PuiseuxMPolyRingElem{T}}) where T <: RingElement = PuiseuxMPolyRing{T}

base_ring_type(::Type{PuiseuxMPolyRing{T}}) where T <: RingElement = Generic.LaurentMPolyWrapRing{T, mpoly_ring_type(T)}
coefficient_ring_type(::Type{PuiseuxMPolyRing{T}}) where T = parent_type(T)

characteristic(R::PuiseuxMPolyRing) = characteristic(base_ring(R))

symbols(R::PuiseuxMPolyRing) = symbols(base_ring(R))

function promote_rule(::Type{PuiseuxMPolyRingElem{S}}, ::Type{T}) where {S <: RingElement, T <: RingElement}
    if PuiseuxMPolyRingElem{S} === T
        return T
    end
    SS = elem_type(coefficient_ring_type(parent_type(S)))
    return promote_rule(SS, T)
end

#################################################################################
#
# Getters
#
#################################################################################

base_ring(R::PuiseuxMPolyRing{T}) where T = R.baseRing::laurent_mpoly_ring_type(T)
coefficient_ring(R::PuiseuxMPolyRing) = coefficient_ring(base_ring(R))

Base.parent(f::PuiseuxMPolyRingElem) = f.parent
poly(f::PuiseuxMPolyRingElem{T}) where {T} = f.poly::laurent_mpoly_type(T)
scale(f::PuiseuxMPolyRingElem) = f.scale

#################################################################################
#
# Setters
#
#################################################################################

# WARNING: input is not assumed to be normalized
function normalize!(f::PuiseuxMPolyRingElem)
    if iszero(f)
        if !isone(scale(f))
            f.scale = one(ZZ)
            return true
        else
            return false
        end
    end

    # make sure scale is correct,
    # i.e., gcd of numerators (= exponents of poly + shift) and denominatos (= scale) is 1
    gcdExponents = gcd(vcat(scale(f), reduce(vcat, [e for e in exponent_vectors(poly(f))])))
    if gcdExponents > 1
        # TODO: use delflate for laurent polys when it is implemented

        vars = gens(base_ring(parent(f)))
        f.poly = sum(c*prod(vars[i].^(Int.(div.(e, gcdExponents))[i]) for i in 1:nvars(parent(f)))
            for (c, e) in zip(coefficients(poly(f)), exponent_vectors(poly(f)))
            )

        f.scale = div(f.scale, gcdExponents)
    end

    return gcdExponents > 1
end

# WARNING: output may not be normalized
function rescale(f::PuiseuxMPolyRingElem, newScale::Int)
    @req newScale > 0 "new scale must be positive"

    # we assume f is normalized
    if newScale == scale(f)
        return f
    end

    newScaleMultipleOfCurrentScale, scaleQuotient = divides(newScale, scale(f))
    @req newScaleMultipleOfCurrentScale "new scale must be a multiple of the current scale"

    newPoly = evaluate(poly(f), gens(base_ring(parent(f))).^(scaleQuotient))

    # when updating to the latest version of AbstractAlgebra, use the new inflate function as below
    # newPoly = inflate(poly(f), [Int(scaleQuotient) for i in 1:nvars(parent(f))])
    return PuiseuxMPolyRingElem(parent(f), newPoly, newScale)
end

#################################################################################
#
# Conversions
#
#################################################################################

# The next function is required but not tested in AbstractAlgebra.
# The following code errors without it:
#   QQt,(t,) = puiseux_polynomial_ring(QQ,["t"]);
#   QQtx, x = polynomial_ring(QQt,3);
#   prod(x .^ rand(1:9,3)) # calls mul_johnson in AA/src/generic/MPoly,jl
function (R::PuiseuxMPolyRing)()
    return zero(R)
end

function (Kt::PuiseuxMPolyRing)(c::RingElement)
    return PuiseuxMPolyRingElem(Kt, base_ring(Kt)(c))
end

function (Kt::PuiseuxMPolyRing{T})(ct::PuiseuxMPolyRingElem{T}) where T <: RingElement
    @req parent(ct) === Kt "parents must be equal" 
    return ct
end

# The next function is required but not tested in AbstractAlgebra.
# The following code errors without it:
# K = algebraic_closure(QQ);
# Kz, z = polynomial_ring(K, "z");
# C = roots(rand(Int8)*z^2+rand(Int8)*z+rand(Int8))
# for _ in 1:99
#     C = vcat(C,roots(rand(Int8)*z^2+rand(Int8)*z+rand(Int8)))
# end
# Kt,(t,) = puiseux_polynomial_ring(K,["t"]);
# Ct = [ Kt(c) * t^rand(Int8) for c in C ]
# Ktx,x = polynomial_ring(Kt,3);
# f = sum([ rand(Ct) * prod(x .^ rand(1:9,3))  for _ in 1:9])
# evaluate(f, [Ktx(1), Ktx(1), Ktx(1)+x[3]]) # runs ^(::MPoly,Int) in AA/src/generic/MPoly,jl

function Base.hash(f::PuiseuxMPolyRingElem, h::UInt)
    normalize!(f)
    return hash((parent(f), poly(f), scale(f)), h)
end

gens(R::PuiseuxMPolyRing) = puiseux_polynomial_ring_elem.(Ref(R), gens(base_ring(R)))
ngens(R::PuiseuxMPolyRing) = ngens(base_ring(R))
nvars(R::PuiseuxMPolyRing) = nvars(base_ring(R))
zero(R::PuiseuxMPolyRing) = puiseux_polynomial_ring_elem(R, zero(base_ring(R)); skip_normalization=true)
one(R::PuiseuxMPolyRing) = puiseux_polynomial_ring_elem(R, one(base_ring(R)); skip_normalization=true)
iszero(f::PuiseuxMPolyRingElem) = iszero(poly(f))
isone(f::PuiseuxMPolyRingElem) = isone(poly(f)) && scale(f) == 1

function Base.:(==)(f::PuiseuxMPolyRingElem, g::PuiseuxMPolyRingElem)
    check_parent(f, g)
    return poly(f) == poly(g) && scale(f) == scale(g)
end

function Base.deepcopy_internal(f::PuiseuxMPolyRingElem, dict::IdDict)
    return puiseux_polynomial_ring_elem(parent(f), deepcopy_internal(poly(f), dict), deepcopy_internal(scale(f), dict); skip_normalization=true)
end

coefficients(f::PuiseuxMPolyRingElem) = coefficients(poly(f))
exponent_vectors(f::PuiseuxMPolyRingElem) = [ e .// scale(f) for e in exponent_vectors(poly(f)) ]
monomials(f::PuiseuxMPolyRingElem) = puiseux_polynomial_ring_elem.(Ref(parent(f)), monomials(poly(f)), Ref(scale(f)))

Base.length(f::PuiseuxMPolyRingElem) = length(poly(f))

function valuation(f::PuiseuxMPolyRingElem)
    @req nvars(parent(f)) == 1 "valuation is only defined for univariate Puiseux polynomials"
    if iszero(f)
        return PosInf()
    end
    return minimum(e[1] for e in exponent_vectors(f))
end

is_univariate(R::PuiseuxMPolyRing) = is_univariate(base_ring(R))
is_gen(f::PuiseuxMPolyRingElem) = is_gen(poly(f)) && scale(f) == 1
is_term(f::PuiseuxMPolyRingElem) = is_term(poly(f))
is_monomial(f::PuiseuxMPolyRingElem) = is_monomial(poly(f).mpoly)
is_unit(f::PuiseuxMPolyRingElem) = is_monomial(f) && is_unit(leading_coefficient(poly(f)))

is_nilpotent(f::PuiseuxMPolyRingElem) = is_nilpotent(poly(f))

#################################################################################
#
# Printing
#
#################################################################################

function expressify(a::PuiseuxMPolyRingElem, x = symbols(parent(a)); context = nothing)
   sum = Expr(:call, :+)
   n = nvars(parent(a))
   for (c, v) in zip(coefficients(a), exponent_vectors(a))
      prod = Expr(:call, :*)
      if !isone(c)
         push!(prod.args, expressify(c, context = context))
      end
      for i in 1:n
         if v[i] != 1
            push!(prod.args, Expr(:call, :^, x[i], expressify(v[i]; context = context)))
         elseif v[i] == 1
            push!(prod.args, x[i])
         end
      end
      push!(sum.args, prod)
   end
   return sum
end

@enable_all_show_via_expressify PuiseuxMPolyRingElem

function show(io::IO, mime::MIME"text/plain", p::PuiseuxMPolyRing)
  @show_name(io, p)
  @show_special(io, mime, p)

  max_vars = 5 # largest number of variables to print
  n = nvars(p)
  print(io, "Puiseux polynomial ring")
  print(io, " in ", ItemQuantity(nvars(p), "variable"), " ")
  if n > max_vars
    join(io, symbols(p)[1:max_vars - 1], ", ")
    println(io, ", ..., ", symbols(p)[n])
  else
    join(io, symbols(p), ", ")
    println(io)
  end
  io = pretty(io)
  print(io, Indent(), "over ", Lowercase(), coefficient_ring(p))
  print(io, Dedent())
end

function show(io::IO, p::PuiseuxMPolyRing)
  @show_name(io, p)
  @show_special(io, p)
  if is_terse(io)
    print(io, "Multivariate polynomial ring")
  else
    io = pretty(io)
    print(io, "Multivariate polynomial ring in ", ItemQuantity(nvars(p), "variable"))
    print(terse(io), " over ", Lowercase(), base_ring(p))
  end
end

#################################################################################
#
# Operations
#
#################################################################################

function Base.:+(f::PuiseuxMPolyRingElem, g::PuiseuxMPolyRingElem)
    check_parent(f, g)
    if iszero(f)
        return g
    elseif iszero(g)
        return f
    end

    # rescale both to the lcm of their scales
    newScale = lcm(scale(f), scale(g))
    frescaled = rescale(f, newScale)
    grescaled = rescale(g, newScale)

    newPoly = poly(frescaled) + poly(grescaled)

    # normalize output, in case of cancellations
    fplusg = PuiseuxMPolyRingElem(
        parent(f),
        newPoly,
        newScale)
    normalize!(fplusg)
    return fplusg
end

function Base.:-(f::PuiseuxMPolyRingElem)
    return PuiseuxMPolyRingElem(parent(f), -poly(f), scale(f))
end

function Base.:-(f::PuiseuxMPolyRingElem, g::PuiseuxMPolyRingElem)
    check_parent(f, g)
    return f + (-g)
end

function Base.:*(f::PuiseuxMPolyRingElem, g::PuiseuxMPolyRingElem)
    check_parent(f, g)
    if iszero(f) || iszero(g)
        return zero(parent(f))
    elseif isone(f)
        return g
    elseif isone(g)
        return f
    end

    # multiply scales and polys
    newScale = scale(f)*scale(g)

    newPoly = evaluate(poly(f), gens(base_ring(parent(f))).^scale(g)) * evaluate(poly(g), gens(base_ring(parent(f))).^scale(f))

    # when updating to the latest version of AbstractAlgebra, use the new inflate function as below
    # newPoly = inflate(poly(f), [Int(scale(g)) for i in 1:nvars(parent(f))]) * inflate(poly(g), [Int(scale(f)) for i in 1:nvars(parent(f))])

    return puiseux_polynomial_ring_elem(parent(f), newPoly, newScale)
end

function Base.:^(f::PuiseuxMPolyRingElem, a::Rational)

    if denominator(a) == 1
        return f^numerator(a)
    end

    @req length(f) == 1 "only monomials can be exponentiated to rational powers"
    @req isone(first(coefficients(f))) "only monomials with coefficient 1 can be exponentiated to rational powers"

    return puiseux_polynomial_ring_elem(
        parent(f),
        poly(f)^numerator(a),
        scale(f)*denominator(a)
    )
end

#
# The next function converts the exponent to a Rational{BigInt}, which
# is not compatible with exponentiating LaurentMPolyRingElems

# function Base.:^(f::PuiseuxMPolyRingElem, a::Rational{Int})
#     return f^(QQ(a))
# end

function Base.:^(f::PuiseuxMPolyRingElem, a::Int)
    if a == 0
        return one(parent(f))
    end
    if a == 1
        return f
    end

    @req a >= 0 || length(f) == 1 "only monomials can be exponentiated to negative powers"
    return puiseux_polynomial_ring_elem(
        parent(f),
        poly(f)^a,
        scale(f)
    )
end

function Base.:^(f::PuiseuxMPolyRingElem, a::Integer)
    return f^(Int(a))
end


function divexact(f::PuiseuxMPolyRingElem{T}, a::T; check::Bool = true) where {T <: RingElement}
    @req !iszero(a) "division by zero"
    @req parent(a) === coefficient_ring(f) "coefficient rings must agree"
    return puiseux_polynomial_ring_elem(parent(f), poly(f)*1//a, scale(f); skip_normalization=true)
end

function divexact(f::PuiseuxMPolyRingElem, a::Integer; check::Bool = true)
    return divexact(f, coefficient_ring(f)(a); check = check)
end

# Note that this method should also be able to divide a product by one of its factors
# e.g., divexact(a*b, a) should return b
function divexact(f::PuiseuxMPolyRingElem, g::PuiseuxMPolyRingElem; check::Bool = true)
    check_parent(f, g)
    @req !iszero(g) "division by zero"

    if iszero(f)
        return zero(parent(f))
    end
    if isone(g)
        return f
    end

    # multiply scales and divide poly(f) by coefficient of poly(g)
    newScale = scale(f)*scale(g)

    vars = gens(base_ring(parent(f)))
    newPoly = divexact(evaluate(poly(f), vars.^scale(g)), evaluate(poly(g), vars.^scale(f)); check = check)

    # when updating to the latest version of AbstractAlgebra, use the new inflate function as below
    # newPoly = divexact(inflate(poly(f), [Int(scale(g)) for i in 1:nvars(parent(f))]), inflate(poly(g), [Int(scale(f)) for i in 1:nvars(parent(f))]))

    return puiseux_polynomial_ring_elem(parent(f), newPoly, newScale)
end

# The following function is required for running the Conformance Tests
function ConformanceTests.generate_element(R::PuiseuxMPolyRing)
    f = ConformanceTests.generate_element(base_ring(R).mpolyring)
    f_laurent = base_ring(R)(f)
    scale = Int(rand(ZZ, 1:10))
    return puiseux_polynomial_ring_elem(R, f_laurent, scale)
end
