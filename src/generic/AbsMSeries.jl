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

function O(a::AbstractAlgebra.AbsMSeriesElem{T}) where T <: RingElement
    if iszero(a)
       return deepcopy(a)    # 0 + O(x^n)
    end
    R = parent(a)
    p = poly(a)
    v = vars(p)
    exps = first(exponent_vector(p, 1))
    (length(v) != 1 || length(p) != 1 || !isone(lc(p))) &&
                                               error("Not a pure power in O()")
    ind = var_index(v[1])
    prec = [i == ind ? exps[i] : R.prec_max[i] for i in 1:length(exps)]
    return R(poly(a), prec)
end

parent_type(::Type{AbsMSeries{T}}) where {T <: RingElement} = AbsMSeriesRing{T}

elem_type(::Type{AbsMSeriesRing{T}}) where {T <: RingElement} = AbsMSeries{T}

###############################################################################
#
#   Basic manipulation
#
###############################################################################

poly(a::AbsMSeries) = a.poly

length(a::AbsMSeries) = length(poly(a))

nvars(R::AbsMSeriesRing) = nvars(R.poly_ring)

precision(a::AbsMSeries) = a.prec

max_precision(R::AbsMSeriesRing) = R.prec_max

function coeff(a::AbsMSeries, n::Int)
    return coeff(poly(a), n)
end

function gen(R::AbsMSeriesRing, i::Int)
    S = R.poly_ring
    prec = [i == ind ? 1 : R.prec_max[i] for ind in 1:nvars(R)]
    x = gen(S, i)
    return R(x, prec)
end

gens(R::AbsMSeriesRing) = [gen(R, i) for i in 1:nvars(R)]

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (R::AbsMSeriesRing{T})(x::MPoly{T}, prec::Vector{Int}) where T <: RingElement
    s = AbsMSeries{T}(x, prec)
    s.parent = R
end

###############################################################################
#
#   PowerSeriesRing constructor
#
###############################################################################

function PowerSeriesRing(R::AbstractAlgebra.Ring, prec::Vector{Int}, s::Vector{T}; cached=true, model=:capped_absolute) where T <: AbstractString
    sym = [Symbol(a) for a in s]
    U = elem_type(R)
 
    S, _ = AbstractAlgebra.PolynomialRing(R, s)

    if model == :capped_absolute
       parent_obj = AbsMSeriesRing{U}(R, S, prec, sym, cached)
    else
       error("Unknown model")
    end
 
    return tuple(parent_obj, gens(parent_obj))
 end