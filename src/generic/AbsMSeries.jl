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
    R.weighted_prec != -1 && error("Operation not possible in weighted rings")
    prec < 0 && error("Precision must be non-negative")
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

@doc raw"""
    weights(R::AbsMSeriesRing)

Return a vector of weights which the variables are weighted with.
"""
function weights(R::AbsMSeriesRing)
   R.weighted_prec == -1 && error("Not a weighted ring")
   return R.prec_max # prec doubles as weights in weighted mode
end

@doc raw"""
    length(a::AbsMSeries)

Return the number of nonzero terms in the series $a$.
"""
length(a::AbsMSeries) = length(poly(a))

@doc raw"""
    number_of_variables(R::AbsMSeriesRing)

Return the number of variables in the series ring.
"""
number_of_variables(R::AbsMSeriesRing) = number_of_variables(poly_ring(R))

number_of_generators(R::AbsMSeriesRing) = number_of_generators(poly_ring(R))

@doc raw"""
    precision(a::AbsMSeries)

Return a vector of precisions, one for each variable in the series ring.
If the ring is weighted the weighted precision is returned instead.
"""
function precision(a::AbsMSeries)
   S = parent(a)
   if S.weighted_prec == -1
      return a.prec
   else
      return S.weighted_prec
   end
end

@doc raw"""
    set_precision!(a::AbsMSeries, prec::Vector{Int})

Set the precisions of the variables in the given series to the values in the
vector `prec`. The precisions must be non-negative. The series will be
truncated to the new precisions. The mutated series is returned.
"""
function set_precision!(a::AbsMSeries, prec::Vector{Int})
    parent(a).weighted_prec != -1 && error("Operation not possible in weighted rings")
    length(prec) != length(a.prec) &&
                         error("Array length not equal to number of variables")
    if !exponents_lt(a.prec, prec)
        a.poly = truncate_poly(a.poly, prec)
    end
    a.prec = prec
    return a
end

@doc raw"""
    max_precision(R::AbsMSeriesRing)

Return a vector of precision caps, one for each variable in the ring.
Arithmetic operations will be performed to precisions not exceeding these
values.
"""
function max_precision(R::AbsMSeriesRing)
   R.weighted_prec != -1 && error("Operation not possible in weighted rings")
   return R.prec_max
end

@doc raw"""
    valuation(a::AbsMSeries)

Return the valuation of $a$ as a vector of integers, one for each variable.
"""
function valuation(a::AbsMSeries)
   parent(a).weighted_prec != -1 && error("Operation not possible in weighted rings")
   p = poly(a)
    prec = a.prec
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

@doc raw"""
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

@doc raw"""
    is_unit(a::AbsMSeries)

Return `true` if the series is a unit in its series ring, i.e. if its constant
term is a unit in the base ring.
"""
is_unit(a::AbsMSeries) = is_unit(constant_coefficient(poly(a)))

@doc raw"""
    gen(R::AbsMSeriesRing, i::Int)

Return the $i$-th generator (variable) of the series ring $R$. Numbering starts
from $1$ for the most significant variable.
"""
function gen(R::AbsMSeriesRing, i::Int)
   @boundscheck 1 <= i <= nvars(R) || throw(ArgumentError("variable index out of range"))
   S = poly_ring(R)
   if R.weighted_prec == -1
      prec = [R.prec_max[ind] for ind in 1:nvars(R)]
      x = R.prec_max[i] > 1 ? gen(S, i) : S()
   else
      w = weights(R)
      prec = [0 for ind in 1:nvars(R)]
      x = R.weighted_prec > w[i] ? gen(S, i) : S()
   end
   return R(x, prec)
end

@doc raw"""
    gens(R::AbsMSeriesRing)

Return a vector of the generators (variables) of the series ring $R$, starting
with the most significant.
"""
gens(R::AbsMSeriesRing) = [gen(R, i) for i in 1:nvars(R)]


@doc raw"""
    is_gen(a::AbsMSeries)

Return true if the series $a$ is a generator of its parent series ring.
"""
function is_gen(a::AbsMSeries)
    R = parent(a)
    p = poly(a)
    v = vars(p)
    return length(v) == 1 && length(p) == 1 &&
          isone(leading_coefficient(p)) && sum(first(exponent_vectors(p))) == 1
end

function deepcopy_internal(a::AbsMSeries, dict::IdDict)
    return parent(a)(deepcopy_internal(poly(a), dict), a.prec)
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

@doc raw"""
    coefficients(a::AbsMSeries)

Return an array of the nonzero coefficients of the series, in the order they
would be displayed, i.e. least significant term first.
"""
function coefficients(a::AbsMSeries)
    return reverse!(collect(coefficients(poly(a))))
end

@doc raw"""
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

function exponents_lt(v::Vector{Int}, w::Vector{Int}, p::Int)
   return sum(v .* w) < p
end

###############################################################################
#
#   Truncation
#
###############################################################################

function truncate_poly(a::MPolyRingElem, prec::Vector{Int}, weighted_prec::Int=-1)
    R = parent(a)
    ctx = MPolyBuildCtx(R)
    for (c, v) in zip(coefficients(a), exponent_vectors(a))
        if weighted_prec == -1
            if exponents_lt(v, prec)
                push_term!(ctx, c, v)
            end
        else
            if exponents_lt(v, prec, weighted_prec)
                push_term!(ctx, c, v)
            end
        end  
    end
    return finish(ctx)
end

@doc raw"""
    truncate(a::AbstractAlgebra.AbsMSeries, prec::Vector{Int})

Return $a$ truncated to (absolute) precisions given by the vector `prec`.
"""
function truncate(a::AbsMSeries, prec::Vector{Int})
    R = parent(a)
    R.weighted_prec != -1 && error("Operation not permitted")
    length(prec) != nvars(R) &&
             error("Array length not equal to number of variables in truncate")
    p = a.prec
    prec = min.(prec, p)
    if prec == p
        # no truncation needed
        return a
    else
        return R(truncate_poly(poly(a), prec), prec)
    end
end

@doc raw"""
    truncate(a::AbstractAlgebra.AbsMSeries, prec::Int)

Return $a$ truncated to precision `prec`. This either truncates by weight in
the weighted cases or truncates each variable to precision `prec` in the
unweighted case.
"""
function truncate(a::AbsMSeries, prec::Int)
    R = parent(a)
    if R.weighted_prec == -1
        return truncate(a, [prec for i in 1:nvars(R)])
    else
        return R(truncate_poly(poly(a), weights(R), prec),
                 [0 for i in 1:nvars(R)]) #??
    end
end


###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::AbsMSeries)
    R = parent(a)
    return R(-poly(a), a.prec)
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::AbsMSeries, b::AbsMSeries)
    check_parent(a, b)
    R = parent(a)
    if R.weighted_prec == -1
        prec = min.(a.prec, b.prec)
        z = truncate_poly(poly(a) + poly(b), prec)
    else
        z = poly(a) + poly(b)
        prec = a.prec
    end
    return R(z, prec)
end

function -(a::AbsMSeries, b::AbsMSeries)
    check_parent(a, b)
    R = parent(a)
    if R.weighted_prec == -1
        prec = min.(a.prec, b.prec)
        z = truncate_poly(poly(a) - poly(b), prec)
    else
        z = poly(a) - poly(b)
        prec = a.prec
    end
    return R(z, prec)
end

function *(a::AbsMSeries, b::AbsMSeries)
    check_parent(a, b)
    R = parent(a)
    if R.weighted_prec == -1
        prec = min.(a.prec .+ valuation(b), b.prec .+ valuation(a))
        prec = min.(prec, max_precision(R))
        z = truncate_poly(poly(a)*poly(b), prec)
    else
        z = truncate_poly(poly(a)*poly(b), weights(R), R.weighted_prec)
        prec = a.prec
    end
    return R(z, prec)
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(a::T, b::AbsMSeries{T}) where {T <: RingElem}
    R = parent(b)
    return R(a*poly(b), b.prec) 
end

function *(a::Union{Integer, Rational, AbstractFloat}, b::AbsMSeries)
    R = parent(b)
    return R(a*poly(b), b.prec) 
end

*(a::AbsMSeries{T}, b::T) where T <: RingElem = b*a
 
*(a::AbsMSeries, b::Union{Integer, Rational, AbstractFloat}) = b*a

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::AbsMSeries, b::Int)
    b < 0 && throw(DomainError(b, "Can't take negative power"))
    R = parent(a)
    if b == 0
        p = one(poly_ring(R))
        if R.weighted_prec == -1
            p = truncate_poly(p, a.prec)
        end
        return R(p, a.prec)
    elseif is_constant(poly(a))
        return R(poly(a)^b, a.prec)
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
    R = parent(x)
    if R.weighted_prec == -1
        prec = min.(x.prec, y.prec)
        p1 = truncate_poly(poly(x), prec)
        p2 = truncate_poly(poly(y), prec)
    else
        p1 = poly(x)
        p2 = poly(y)
    end
    return p1 == p2
end

function isequal(x::AbsMSeries{T}, y::AbsMSeries{T}) where T <: RingElement
    check_parent(x, y)
    R = parent(x)
    if R.weighted_prec == -1
        prec = x.prec
        prec == y.prec || return false
        return truncate_poly(poly(x), prec) == truncate_poly(poly(y), prec)
    else
        return x == y
    end
end

###############################################################################
#
#   Inverse
#
###############################################################################


function Base.inv(a::AbsMSeries)
    R = parent(a)
    ainv = R(inv(constant_coefficient(poly(a))))
    if R.weighted_prec == -1
        # use the precision stored in the polynomial
        # arithmetic uses precision stored in polynomials
        max_n = sum(a.prec)
        cur_n = 1
        # 1-a*ainv = er where each monomial in er has total degree >= cur_n
        # Furthermore, we only care about the terms in er where the exponent
        # on variable i is restricted to the range [0, a.prec[i]).
        # Therefore, with max_n = sum(a.prec), we are done if cur_n >= max_n
        while true
            cur_n *= 2
            trunc = [min(a.prec[i], cur_n) for i in 1:nvars(R)]
            set_precision!(ainv, trunc)
            e = 2 - truncate(a, trunc)*ainv
            (trunc == a.prec && isone(e)) && break
            ainv = e*ainv
            cur_n >= max_n && break
        end
    else
        # use the precision stored in the parent
        # arithmetic uses precision stored in parent
        max_n = R.weighted_prec
        cur_n = minimum(R.prec_max)
        @assert cur_n > 0
        # 1-a*ainv = er where each monomial in er has weight >= cur_n
        while true
            cur_n *= 2
            trunc = min(max_n, cur_n)
            e = 2 - truncate(a, trunc)*ainv
            (trunc == max_n && isone(e)) && break
            ainv = e*ainv
            cur_n >= max_n && break
        end
    end
    return ainv
end

###############################################################################
#
#   Exact division
#
###############################################################################

@doc raw"""
    divexact(x::AbsMSeries{T}, y::AbsMSeries{T}; check::Bool=true) where T <: RingElement

Return the exact quotient of the series $x$ by the series $y$. This function
currently assumes $y$ is an invertible series.
"""
function divexact(x::AbsMSeries{T}, y::AbsMSeries{T}; check::Bool=true) where T <: RingElement
    check_parent(x, y)
    return x*inv(y)
end

###############################################################################
#
#   Evaluation
#
###############################################################################

@doc raw"""
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

@doc raw"""
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

@doc raw"""
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
    if R.weighted_prec == -1
        a.poly = truncate_poly(a.poly, prec)
        a.prec = prec
    end
    return a
end

function mul!(c::AbsMSeries{T}, a::AbsMSeries{T}, b::AbsMSeries{T}) where
                                                            T <: RingElement
    R = parent(a)
    if R.weighted_prec == -1
        prec = min.(a.prec .+ valuation(b), b.prec .+ valuation(a))
        prec = min.(prec, max_precision(R))
        c.poly = mul!(c.poly, a.poly, b.poly)
        c.poly = truncate_poly(c.poly, prec)
        c.prec = prec
    else
        c.poly = mul!(c.poly, a.poly, b.poly)
        c.poly = truncate_poly(c.poly, weights(R), R.weighted_prec)
    end 
    return c
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

function promote_rule(::Type{AbsMSeries{T, V}}, ::Type{U}) where
                                        {V, T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? AbsMSeries{T, V} : Union{}
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (R::AbsMSeriesRing{T, S})(x::S, prec::Vector{Int}) where
                          {T <: RingElement, S <: AbstractAlgebra.MPolyRingElem{T}}
    for v in prec
        v < 0 && error("Precision must be non-negative")
    end
    s = AbsMSeries{T, S}(R, x, prec)
    return s
end

function (R::AbsMSeriesRing)()
    if R.weighted_prec == -1
        return R(poly_ring(R)(), max_precision(R))
    else
        return R(poly_ring(R)(), [0 for i in 1:nvars(R)])
    end
end

function (R::AbsMSeriesRing{T})(x::T) where T <: RingElem
    if R.weighted_prec == -1
        return R(poly_ring(R)(x), max_precision(R))
    else
        return R(poly_ring(R)(x), [0 for i in 1:nvars(R)])
    end
end

function (R::AbsMSeriesRing{T})(x::AbsMSeries{T}) where T <: RingElement
    parent(x) != R && error("Unable to coerce")
    return x
end

function (R::AbsMSeriesRing)(b::Union{Integer, Rational, AbstractFloat})
    if R.weighted_prec == -1
        return R(poly_ring(R)(b), max_precision(R))
    else
        return R(poly_ring(R)(b),  [0 for i in 1:nvars(R)])
    end
end

###############################################################################
#
#   power_series_ring constructor
#
###############################################################################

function power_series_ring(R::AbstractAlgebra.Ring, prec::Vector{Int},
                  s::Vector{Symbol}; cached::Bool=true, model=:capped_absolute)
    U = elem_type(R)

    S, _ = AbstractAlgebra.polynomial_ring(R, s)
    V = elem_type(S)

    model === :capped_absolute || error("Unknown model")

    parent_obj = AbsMSeriesRing{U, V}(S, prec, s, cached)

    return parent_obj, gens(parent_obj)
end

function power_series_ring(R::AbstractAlgebra.Ring, prec::Int,
        s::Vector{Symbol}; weights::Union{Vector{Int}, Nothing}=nothing,
        cached::Bool=true, model=:capped_absolute)
    U = elem_type(R)

    S, _ = AbstractAlgebra.polynomial_ring(R, s)
    V = elem_type(S)

    model === :capped_absolute || error("Unknown model")

    if weights === nothing
        parent_obj = AbsMSeriesRing{U, V}(S, [prec for _ in s], s, cached)
    else
        parent_obj = AbsMSeriesRing{U, V}(S, weights, prec, s, cached)
    end

    return parent_obj, gens(parent_obj)
end

power_series_ring(R::AbstractAlgebra.Ring, weights::Vector{Int},
    prec::Int, s::Vector{Symbol}; kw...
    ) = power_series_ring(R, prec, s; weights, kw...)
