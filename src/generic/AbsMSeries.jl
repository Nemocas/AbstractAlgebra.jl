###############################################################################
#
#   AbsMSeries.jl : Multivariate power series over rings, capped absolute
#                   precision
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

function O(R::AbsMSeriesRing{T}, prec::Int) where T <: RingElement
    prec < 0 && error("Precision must be nonnegative")
    return R(poly_ring(R)(), fill(prec, nvars(R)))
end

function parent_type(::Type{AbsMSeries{T, S}}) where {T <: RingElement, S}
    return AbsMSeriesRing{T, S}
end

function elem_type(::Type{AbsMSeriesRing{T, S}}) where {T <: RingElement, S}
    return AbsMSeries{T, S}
end

function check_parent(a::AbsMSeries, b::AbsMSeries, throw::Bool = true)
    c = parent(a) != parent(b)
    c && throw &&
            error("Incompatible multivariate series rings in series operation")
    return !c
 end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

poly(a::AbsMSeries{T, S}) where {T <: RingElement, S} = a.poly::S

function poly_ring(R::AbsMSeriesRing{T, S}) where {T <: RingElement, S}
   return R.poly_ring::parent_type(S) 
end

@doc Markdown.doc"""
    length(a::AbsMSeries)

Return the number of nonzero terms in the series $a$.
"""
length(a::AbsMSeries) = length(poly(a))

@doc Markdown.doc"""
    nvars(R::AbsMSeriesRing)

Return the number of variables in the series ring.
"""
nvars(R::AbsMSeriesRing) = nvars(poly_ring(R))

@doc Markdown.doc"""
    precision(a::AbsMSeries)

Return a vector of precisions, one for each variable in the series ring.
"""
precision(a::AbsMSeries) = a.prec

@doc Markdown.doc"""
    set_precision!(a::AbsMSeries, prec::Vector{Int})

Set the precisions of the variables in the given series to the values in the
vector `prec`. The precisions must be non-negative. The series will be
truncated to the new precisions. The mutated series is returned.
"""
function set_precision!(a::AbsMSeries, prec::Vector{Int})
    length(prec) != length(precision(a)) &&
                         error("Array length not equal to number of variables")
    if !exponents_lt(precision(a), prec)
        a.poly = truncate_poly(a.poly, prec)
    end
    a.prec = prec
    return a
end

@doc Markdown.doc"""
    max_precision(R::AbsMSeriesRing)

Return a vector of precision caps, one for each variable in the ring.
Arithmetic operations will be performed to precisions not exceeding these
values.
"""
max_precision(R::AbsMSeriesRing) = R.prec_max

@doc Markdown.doc"""
    valuation(a::AbsMSeries)

Return the valuation of $a$ as a vector of integers, one for each variable.
"""
function valuation(a::AbsMSeries)
    p = poly(a)
    prec = precision(a)
    val = prec
    for v in exponent_vectors(p)
        if iszero(v)
            return v
        end
        v = exponents_clamp_zero_to_prec(v, prec)
        val = min.(v, val)
    end
    return val
end

@doc Markdown.doc"""
    coeff(a::AbsMSeries, n::Int)

Return the coefficient of the $n$-th nonzero term of the series (or zero if
there are fewer than $n$ nonzero terms). Terms are numbered from the least
significant term, i.e. the first term displayed when the series is printed.
"""
function coeff(a::AbsMSeries, n::Int)
    return coeff(poly(a), length(a) - n + 1)
end

iszero(a::AbsMSeries) = length(poly(a)) == 0

function isone(a::AbsMSeries)
    return isone(poly(a)) || iszero(precision(a))
end

zero(R::AbsMSeriesRing) = R(0)

one(R::AbsMSeriesRing) = R(1)

@doc Markdown.doc"""
    isunit(a::AbsMSeries)

Return `true` if the series is a unit in its series ring, i.e. if its constant
term is a unit in the base ring.
"""
isunit(a::AbsMSeries) = isunit(constant_coefficient(poly(a)))

@doc Markdown.doc"""
    gen(R::AbsMSeriesRing, i::Int)

Return the $i$-th generator (variable) of the series ring $R$. Numbering starts
from $1$ for the most significant variable.
"""
function gen(R::AbsMSeriesRing, i::Int)
    S = poly_ring(R)
    prec = [R.prec_max[ind] for ind in 1:nvars(R)]
    x = R.prec_max[i] > 1 ? gen(S, i) : S()
    return R(x, prec)
end

@doc Markdown.doc"""
    gens(R::AbsMSeriesRing)

Return a vector of the generators (variables) of the series ring $R$, starting
with the most significant.
"""
gens(R::AbsMSeriesRing) = [gen(R, i) for i in 1:nvars(R)]


@doc Markdown.doc"""
    isgen(a::AbsMSeries)

Return true if the series $a$ is a generator of its parent series ring.
"""
function isgen(a::AbsMSeries)
    R = parent(a)
    p = poly(a)
    v = vars(p)
    return length(v) == 1 && length(p) == 1 &&
          isone(leading_coefficient(p)) && sum(first(exponent_vectors(p))) == 1
end

function deepcopy_internal(a::AbsMSeries, dict::IdDict)
    return parent(a)(deepcopy_internal(poly(a), dict), precision(a))
end

function Base.hash(a::AbsMSeries, h::UInt)
    b = 0xf7f073b6c9e1d560
    return xor(b, hash(poly(a), h))
end

###############################################################################
#
#   Iterators
#
###############################################################################

@doc Markdown.doc"""
    coefficients(a::AbsMSeries)

Return an array of the nonzero coefficients of the series, in the order they
would be displayed, i.e. least significant term first.
"""
function coefficients(a::AbsMSeries)
    return reverse!(collect(coefficients(poly(a))))
end

@doc Markdown.doc"""
    exponent_vectors(a::AbsMSeries)

Return an array of the exponent vectors of the nonzero terms of the series, in
the order they would be displayed, i.e. least significant term first.
"""
function exponent_vectors(a::AbsMSeries)
    return reverse!(collect(exponent_vectors(poly(a))))
end

###############################################################################
#
#   Coefficients, terms and exponent vectors
#
###############################################################################

# set the exponent vector of the underlying polynomial (used internally)
function set_exponent_vector!(a::AbsMSeries, i::Int64, v::Vector{Int64})
    a.poly = set_exponent_vector!(poly(a), i, v)
    return a
end

# set the coefficient of the underlying polynomial (used internally)
function setcoeff!(a::AbsMSeries{T}, i::Int64, c::T) where T <: RingElement
    a.poly = setcoeff!(poly(a), i, c)
    return a
end

# used by MPolyBuildCtx by evaluation
function sort_terms!(a::AbsMSeries)
    a.poly = sort_terms!(poly(a))
    return a
end

# used by MPolyBuildCtx by evaluation
function combine_like_terms!(a::AbsMSeries)
    a.poly = combine_like_terms!(poly(a))
    return a
end

###############################################################################
#
#   Exponent helper functions (not exported)
#
###############################################################################

function exponents_clamp_zero_to_prec(a::Vector{Int}, prec::Vector{Int})
    return [a[i] == 0 ? prec[i] : a[i] for i in 1:length(a)]
end

function exponents_lt(v::Vector{Int}, p::Vector{Int})
    return all(((x, y),) -> x < y, zip(v, p))
end

###############################################################################
#
#   Truncation
#
###############################################################################

function truncate_poly(a::MPolyElem, prec::Vector{Int})
    R = parent(a)
    ctx = MPolyBuildCtx(R)
    for (c, v) in zip(coefficients(a), exponent_vectors(a))
        if exponents_lt(v, prec)
            push_term!(ctx, c, v)
        end
    end
    return finish(ctx)
end

@doc Markdown.doc"""
    truncate(a::AbstractAlgebra.AbsMSeries, prec::Vector{Int})

Return $a$ truncated to (absolute) precisions given by the vector `prec`.
"""
function truncate(a::AbsMSeries, prec::Vector{Int})
    R = parent(a)
    length(prec) != nvars(R) &&
             error("Array length not equal to number of variables in truncate")
    trunc_needed = false
    p = precision(a)
    for i = 1:nvars(R)
        if prec[i] < p[i]
            trunc_needed = true
            break
        end
    end
    if !trunc_needed
        return a
    end
    prec = min.(prec, p)
    q = truncate_poly(poly(a), prec)
    return R(q, prec)
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::AbsMSeries)
    R = parent(a)
    return R(-poly(a), precision(a))
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::AbsMSeries, b::AbsMSeries)
    check_parent(a, b)
    R = parent(a)
    prec = min.(precision(a), precision(b))
    z = truncate_poly(poly(a) + poly(b), prec)
    return R(z, prec)
end

function -(a::AbsMSeries, b::AbsMSeries)
    check_parent(a, b)
    R = parent(a)
    prec = min.(precision(a), precision(b))
    z = truncate_poly(poly(a) - poly(b), prec)
    return R(z, prec)
end

function *(a::AbsMSeries, b::AbsMSeries)
    check_parent(a, b)
    R = parent(a)
    prec = min.(precision(a) .+ valuation(b), precision(b) .+ valuation(a))
    prec = min.(prec, max_precision(R))
    z = truncate_poly(poly(a)*poly(b), prec)
    return R(z, prec)
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(a::T, b::AbsMSeries{T}) where {T <: RingElem}
    R = parent(b)
    return R(a*poly(b), precision(b)) 
end

function *(a::Union{Integer, Rational, AbstractFloat}, b::AbsMSeries)
    R = parent(b)
    return R(a*poly(b), precision(b)) 
end

*(a::AbsMSeries{T}, b::T) where T <: RingElem = b*a
 
*(a::AbsMSeries, b::Union{Integer, Rational, AbstractFloat}) = b*a

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::AbsMSeries, b::Int) where T <: RingElement
    b < 0 && throw(DomainError(b, "Can't take negative power"))
    R = parent(a)
    prec = precision(a)
    if b == 0
        p = one(poly_ring(R))
        p = truncate_poly(p, prec)
        return R(p, prec)
    elseif isconstant(poly(a))
        return R(poly(a)^b, precision(a))
    elseif b == 1
        return deepcopy(a)
    end
    bit = ~((~UInt(0)) >> 1)
    while (UInt(bit) & b) == 0
        bit >>= 1
    end
    z = a
    bit >>= 1
    while bit !=0
        z = z*z
        if (UInt(bit) & b) != 0
            z *= a
        end
        bit >>= 1
    end
    return z
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(x::AbsMSeries{T}, y::AbsMSeries{T}) where T <: RingElement
    check_parent(x, y)
    prec = min.(precision(x), precision(y))
    p1 = truncate_poly(poly(x), prec)
    p2 = truncate_poly(poly(y), prec)
    return p1 == p2
end

function isequal(x::AbsMSeries{T}, y::AbsMSeries{T}) where T <: RingElement
    check_parent(x, y)
    prec = precision(x)
    prec == precision(y) || return false
    return truncate_poly(poly(x), prec) == truncate_poly(poly(y), prec)
end

###############################################################################
#
#   Inverse
#
###############################################################################

@doc Markdown.doc"""
    Base.inv(x::AbsMSeries)

Return the inverse of the series $x$. An exception is raised if the series is
not a unit.
"""
function Base.inv(x::AbsMSeries)
    !isunit(x) && error("Not a unit")
    R = parent(x)
    prec = [1 for n in 1:nvars(R)]
    cinv = inv(coeff(x, 1))
    xinv = R(poly_ring(R)(cinv), prec)
    two = R(poly_ring(R)(2), prec)
    # lift each variable in turn
    for var = nvars(R):-1:1
        nvar = precision(x)[var]
        var_prec = [nvar]
        while nvar != 1
            nvar = div(nvar + 1, 2)
            push!(var_prec, nvar)
        end
        # list var quadratically
        for i = length(var_prec) - 1:-1:1
            prec[var] = var_prec[i]
            two = set_precision!(two, prec)
            xinv = set_precision!(xinv, prec)
            xinv = (two - x*xinv)*xinv
        end
    end
    return xinv
end

###############################################################################
#
#   Exact division
#
###############################################################################

@doc Markdown.doc"""
    divexact(x::AbsMSeries{T}, y::AbsMSeries{T}) where T <: RingElement

Return the exact quotient of the series $x$ by the series $y$. This function
currently assumes $y$ is an invertible series.
"""
function divexact(x::AbsMSeries{T}, y::AbsMSeries{T}) where T <: RingElement
    check_parent(x, y)
    return x*inv(y)
end

###############################################################################
#
#   Evaluation
#
###############################################################################

@doc Markdown.doc"""
    evaluate(a::U, vars::Vector{Int}, vals::Vector{U}) where {T <: RingElement, U <: AbsMSeries{T}}

Evaluate the series expression by substituting in the supplied values in
the array `vals` for the corresponding variables with indices given by the
array `vars`. The values must be in the same ring as $a$.
"""
function evaluate(a::U, vars::Vector{Int}, vals::Vector{U}) where
                                         {T <: RingElement, U <: AbsMSeries{T}}
    R = parent(a)
    unique(vars) != vars && error("Variables not unique")
    length(vars) != length(vals) &&
        error("Number of variables does not match number of values")
    for i = 1:length(vars)
        if vars[i] < 1 || vars[i] > nvars(parent(a))
            error("Variable index not in range")
        end
        parent(vals[i]) !== R && error("Element not in series ring")
    end
 
    if length(vars) == 0
        return a
    end
 
    S = parent(a)
    R = base_ring(a)
    return AbstractAlgebra._evaluate(a, S, R, vars, vals)
end

@doc Markdown.doc"""
    evaluate(a::U, vars::Vector{U}, vals::Vector{U}) where {T <: RingElement, U <: AbsMSeries{T}}

Evaluate the series expression by substituting in the supplied values in
the array `vals` for the corresponding variables given by the array `vars`.
The values must be in the same ring as $a$.
"""
function evaluate(a::U, vars::Vector{U}, vals::Vector{U}) where
                                        {T <: RingElement, U <: AbsMSeries{T}}
    varidx = Int[var_index(poly(x)) for x in vars]
    return evaluate(a, varidx, vals)
end

@doc Markdown.doc"""
    evaluate(a::U, vals::Vector{U}) where {T <: RingElement, U <: AbsMSeries{T}}

Evaluate the series expression by substituting in the supplied values in
the array `vals` for the variables the series ring to which $a$ belongs. The
values must be in the same ring as $a$.
"""
function evaluate(a::U, vals::Vector{U}) where
                                         {T <: RingElement, U <: AbsMSeries{T}}
    R = parent(a)
    return evaluate(a, [i for i in 1:nvars(R)], vals)
end

###############################################################################
#
#   Unsafe operators
#
###############################################################################

function addeq!(a::AbsMSeries{T}, b::AbsMSeries{T}) where T <: RingElement
    R = parent(a)
    prec = min.(precision(a), precision(b))
    a.poly = addeq!(a.poly, b.poly)
    a.poly = truncate_poly(a.poly, prec)
    a.prec = prec
    return a
end

function mul!(c::AbsMSeries{T}, a::AbsMSeries{T}, b::AbsMSeries{T}) where
                                                            T <: RingElement
    R = parent(a)
    prec = min.(precision(a) .+ valuation(b), precision(b) .+ valuation(a))
    prec = min.(prec, max_precision(R))
    c.poly = mul!(c.poly, a.poly, b.poly)
    c.poly = truncate_poly(c.poly, prec)
    c.prec = prec
    return c
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{AbsMSeries{T}}, ::Type{AbsMSeries{T}}) where
                                               T <: RingElement = AbsMSeries{T}

function promote_rule(::Type{AbsMSeries{T}}, ::Type{U}) where
                                           {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? AbsMSeries{T} : Union{}
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (R::AbsMSeriesRing{T, S})(x::S, prec::Vector{Int}) where
                          {T <: RingElement, S <: AbstractAlgebra.MPolyElem{T}}
    for v in prec
        v < 0 && error("Precision must be non-negative")
    end
    s = AbsMSeries{T, S}(R, x, prec)
    return s
end

function (R::AbsMSeriesRing)()
    return R(poly_ring(R)(), max_precision(R))
end

function (R::AbsMSeriesRing{T})(x::T) where T <: RingElem
    return R(poly_ring(R)(x), max_precision(R))
end

function (R::AbsMSeriesRing{T})(x::AbsMSeries{T}) where T <: RingElement
    parent(x) != R && error("Unable to coerce")
    return x
end

function (R::AbsMSeriesRing)(b::Union{Integer, Rational, AbstractFloat})
    return R(poly_ring(R)(b), max_precision(R))
end

###############################################################################
#
#   PowerSeriesRing constructor
#
###############################################################################

function PowerSeriesRing(R::AbstractAlgebra.Ring, prec::Vector{Int},
                  s::Vector{T}; cached=true, model=:capped_absolute) where
                                                                    T <: Symbol
    str = [String(a) for a in s]
    U = elem_type(R)
 
    S, _ = AbstractAlgebra.PolynomialRing(R, str)
    V = elem_type(S)

    if model == :capped_absolute
       parent_obj = AbsMSeriesRing{U, V}(S, prec, s, cached)
    else
       error("Unknown model")
    end
 
    return tuple(parent_obj, gens(parent_obj))
end

function PowerSeriesRing(R::AbstractAlgebra.Ring, prec::Int,
                  s::Vector{T}; cached=true, model=:capped_absolute) where
                                                                    T <: Symbol
    prec_vec = [prec for v in s]
    return PowerSeriesRing(R, prec_vec, s; cached=cached, model=model)
end
