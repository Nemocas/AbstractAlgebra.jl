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
       return deepcopy(a)
    end
    R = parent(a)
    p = poly(a)
    v = vars(p)
    (length(v) != 1 || length(p) != 1 || !isone(lc(p))) &&
                                               error("Not a pure power in O()")
    ind = var_index(v[1])
    exps = first(exponent_vectors(p))
    prec = [i == ind ? exps[i] : R.prec_max[i] for i in 1:length(exps)]
    return R(parent(p)(), prec)
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

function set_precision!(a::AbsMSeries, prec::Vector{Int})
    length(prec) != length(precision(a)) &&
                         error("Array length not equal to number of variables")
    if !exponents_lt(precision(a), prec)
        a.poly = truncate_poly(a.poly, prec)
    end
    a.prec = prec
    return a
end

max_precision(R::AbsMSeriesRing) = R.prec_max

@doc Markdown.doc"""
    valuation(a::AbsMSeries)

Return the valuation of $a$ as a vector of integers.
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

function coeff(a::AbsMSeries, n::Int)
    return coeff(poly(a), length(a) - n + 1)
end

iszero(a::AbsMSeries) = length(poly(a)) == 0

function isone(a::AbsMSeries)
    return isone(poly(a)) || iszero(precision(a))
end

zero(R::AbsMSeriesRing) = R(0)

one(R::AbsMSeriesRing) = R(1)

isunit(a::AbsMSeries) = isunit(constant_coefficient(poly(a)))

function gen(R::AbsMSeriesRing, i::Int)
    S = R.poly_ring
    prec = [R.prec_max[ind] for ind in 1:nvars(R)]
    x = R.prec_max[i] > 1 ? gen(S, i) : S()
    return R(x, prec)
end

gens(R::AbsMSeriesRing) = [gen(R, i) for i in 1:nvars(R)]

function isgen(a::AbsMSeries)
    R = parent(a)
    p = poly(a)
    v = vars(p)
    return length(v) == 1 && length(p) == 1 &&
                          isone(lc(p)) && sum(first(exponent_vectors(p))) == 1
end

vars(R::AbstractAlgebra.MSeriesRing) = R.sym

parent(a::AbstractAlgebra.MSeriesElem) = a.parent

function base_ring(R::AbstractAlgebra.MSeriesRing{T}) where T <: RingElement
    return R.base_ring::parent_type(T)
end

base_ring(a::AbstractAlgebra.MSeriesElem) = base_ring(parent(a))

function deepcopy_internal(a::AbsMSeries, dict::IdDict)
    return parent(a)(deepcopy(poly(a)), precision(a))
end

function Base.hash(a::AbsMSeries, h::UInt)
    b = 0xf7f073b6c9e1d560
    return xor(b, hash(poly(a), h))
end

function characteristic(a::AbstractAlgebra.MSeriesRing)
    return characteristic(base_ring(a))
 end

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function AbstractAlgebra.expressify(a::AbstractAlgebra.AbsMSeriesElem,
                                        x = vars(parent(a)); context = nothing)
    apoly = poly(a)

    poly_sum = Expr(:call, :+)
    n = nvars(parent(apoly))
    
    iter = zip(coeffs(apoly), exponent_vectors(apoly))
    cv = reverse(collect(iter))
    
    for (c, v) in cv
        prod = Expr(:call, :*)
        if !isone(c)
            push!(prod.args, expressify(c, context = context))
        end
        for i in n:-1:1
            if v[i] > 1
                push!(prod.args, Expr(:call, :^, x[i], v[i]))
            elseif v[i] == 1
                push!(prod.args, x[i])
            end
        end
        push!(poly_sum.args, prod)
    end
    
    sum = Expr(:call, :+)

    push!(sum.args, poly_sum)

    for i in nvars(parent(a)):-1:1
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
#   Exponent helper functions (not exported)
#
###############################################################################

function exponents_clamp_zero_to_prec(a::Vector{Int}, prec::Vector{Int})
    return [a[i] == 0 ? prec[i] : a[i] for i in 1:length(a)]
end

function exponents_lt(v::Vector{Int}, p::Vector{Int})
    for i = 1:length(v)
        if v[i] >= p[i]
            return false
        end
    end
    return true
end

###############################################################################
#
#   Truncation
#
###############################################################################

function truncate_poly(a::MPolyElem, prec::Vector{Int})
    R = parent(a)
    ctx = MPolyBuildCtx(R)
    for (c, v) in zip(coeffs(a), exponent_vectors(a))
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
    R = parent(a)
    prec = min.(precision(a), precision(b))
    z = truncate_poly(poly(a) + poly(b), prec)
    return R(z, prec)
end

function -(a::AbsMSeries, b::AbsMSeries)
    R = parent(a)
    prec = min.(precision(a), precision(b))
    z = truncate_poly(poly(a) - poly(b), prec)
    return R(z, prec)
end

function *(a::AbsMSeries, b::AbsMSeries)
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
        p = one(R.poly_ring)
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
    prec = min.(precision(x), precision(y))
    p1 = truncate_poly(poly(x), prec)
    p2 = truncate_poly(poly(y), prec)
    return p1 == p2
end

function isequal(x::AbsMSeries{T}, y::AbsMSeries{T}) where T <: RingElement
    return precision(x) == precision(y) && poly(x) == poly(y)
end

###############################################################################
#
#   Inverse
#
###############################################################################

function Base.inv(x::AbsMSeries)
    !isunit(x) && error("Not a unit")
    R = parent(x)
    prec = [1 for n in 1:nvars(R)]
    cinv = inv(coeff(x, 1))
    xinv = R(R.poly_ring(cinv), prec)
    two = R(R.poly_ring(2), prec)
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

function divexact(x::AbsMSeries{T}, y::AbsMSeries{T}) where T <: RingElement
    return x*inv(y)
end

###############################################################################
#
#   Random elements
#
###############################################################################

RandomExtensions.maketype(S::AbstractAlgebra.MSeriesRing, _, _) = elem_type(S)

function RandomExtensions.make(S::AbstractAlgebra.MSeriesRing,
                                      term_range::UnitRange{Int}, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, term_range, vs[1])
   else
      make(S, term_range, make(R, vs...))
   end
end

function rand(rng::AbstractRNG, sp::SamplerTrivial{<:Make3{
                 <:RingElement,<:AbstractAlgebra.MSeriesRing,UnitRange{Int}}})
   S, term_range, v = sp[][1:end]
   f = S()
   g = gens(S)
   R = base_ring(S)
   prec = max_precision(S)
   for i = 1:rand(rng, term_range)
      term = S(1)
      for j = 1:length(g)
         term *= g[j]^rand(rng, 0:prec[j])
      end
      term *= rand(rng, v)
      f += term
   end
   return f
end

function rand(rng::AbstractRNG, S::AbstractAlgebra.MSeriesRing,
                                             term_range::UnitRange{Int}, v...)
   rand(rng, make(S, term_range, v...))
end

function rand(S::AbstractAlgebra.MSeriesRing, term_range, v...)
   rand(GLOBAL_RNG, S, term_range, v...)
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{AbsMSeries{T}}, ::Type{AbsMSeries{T}}) where T <: RingElement = AbsMSeries{T}

function promote_rule(::Type{AbsMSeries{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? AbsMSeries{T} : Union{}
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (R::AbsMSeriesRing{T})(x::MPoly{T}, prec::Vector{Int}) where T <: RingElement
    for v in prec
        v < 0 && error("Precision must be non-negative")
    end
    s = AbsMSeries{T}(x, prec)
    s.parent = R
    return s
end

function (R::AbsMSeriesRing)()
    return R(R.poly_ring(), max_precision(R))
end

function (R::AbsMSeriesRing{T})(x::T) where T <: RingElem
    return R(R.poly_ring(x), max_precision(R))
end

function (R::AbsMSeriesRing)(b::Union{Integer, Rational, AbstractFloat})
    return R(R.poly_ring(b), max_precision(R))
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