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
    prec = [R.prec_max[ind] for ind in 1:nvars(R)]
    x = R.prec_max[i] > 1 ? gen(S, i) : S()
    return R(x, prec)
end

gens(R::AbsMSeriesRing) = [gen(R, i) for i in 1:nvars(R)]

vars(R::AbstractAlgebra.MSeriesRing) = R.sym

parent(a::AbstractAlgebra.MSeriesElem) = a.parent

function base_ring(R::AbstractAlgebra.MSeriesRing{T}) where T <: RingElement
    return R.base_ring::parent_type(T)
end

base_ring(a::AbstractAlgebra.MSeriesElem) = base_ring(parent(a))

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function AbstractAlgebra.expressify(a::AbstractAlgebra.AbsMSeriesElem,
                                        x = vars(parent(a)); context = nothing)
    sum = Expr(:call, :+)

    push!(sum.args, expressify(a.poly, context = context))

    for i in 1:nvars(parent(a))
        push!(sum.args, Expr(:call, :O, Expr(:call, :^, x[i], a.prec[i])))
    end

    return sum
end

function Base.show(io::IO, a::AbstractAlgebra.MSeriesElem)
    print(io, AbstractAlgebra.obj_to_string(a, context = io))
end
  
function Base.show(io::IO, ::MIME"text/plain", a::AbstractAlgebra.MSeriesElem)
    print(io, AbstractAlgebra.obj_to_string(a, context = io))
end
  
function show(io::IO, a::AbstractAlgebra.MSeriesRing)
    v = join([String(s) for s in vars(a)], ", ")
    print(io, "Multivariate power series ring in ", v, " over ")
    print(IOContext(io, :compact => true), base_ring(a))
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (R::AbsMSeriesRing{T})(x::MPoly{T}, prec::Vector{Int}) where T <: RingElement
    s = AbsMSeries{T}(x, prec)
    s.parent = R
    return s
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