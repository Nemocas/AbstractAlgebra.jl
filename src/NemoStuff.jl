export neg!

function neg!(w::Vector{Int})
    w .*= -1
end

sub!(z::T, x::T, y::T) where {T} = x - y

sub!(z::Rational{Int}, x::Rational{Int}, y::Int) = x - y

neg!(z::Rational{Int}, x::Rational{Int}) = -x

add!(z::Rational{Int}, x::Rational{Int}, y::Int) = x + y

mul!(z::Rational{Int}, x::Rational{Int}, y::Int) = x * y

is_negative(n::T) where T<:Real = n < zero(T)
is_positive(n::T) where T<:Real = n > zero(T)

# TODO (CF):
# should be Bernstein'ed: this is slow for large valuations
# returns the maximal v s.th. z mod p^v == 0 and z div p^v
#   also useful if p is not prime....
#
# TODO: what happens to z = 0???

function remove(z::Rational{T}, p::T) where {T<:Integer}
    z == 0 && return (0, z)
    v, d = remove(denominator(z), p)
    w, n = remove(numerator(z), p)
    return w - v, n // d
end

function valuation(z::Rational{T}, p::T) where {T<:Integer}
    z == 0 && error("Not yet implemented")
    v = valuation(denominator(z), p)
    w = valuation(numerator(z), p)
    return w - v
end

base_ring(::Vector{Int}) = Int

################################################################################
#
#  Zero matrix constructors
#
################################################################################

function zero_matrix(::Type{MatElem}, R::Ring, n::Int, m::Int)
    return zero_matrix(R, n, m)
end

function identity_matrix(::Type{MatElem}, R::Ring, n::Int)
    return identity_matrix(R, n)
end

################################################################################
#
#  Unsafe functions for generic matrices
#
################################################################################

#function zero!(a::MatElem)
#  for i in 1:nrows(a)
#    for j in 1:ncols(a)
#      a[i, j] = zero!(a[i, j])
#    end
#  end
#  return a
#end

function mul!(c::MatElem, a::MatElem, b::MatElem)
    ncols(a) != nrows(b) && error("Incompatible matrix dimensions")
    nrows(c) != nrows(a) && error("Incompatible matrix dimensions")
    ncols(c) != ncols(b) && error("Incompatible matrix dimensions")

    if c === a || c === b
        d = parent(a)()
        return mul!(d, a, b)
    end

    t = base_ring(a)()
    for i = 1:nrows(a)
        for j = 1:ncols(b)
            c[i, j] = zero!(c[i, j])
            for k = 1:ncols(a)
                c[i, j] = addmul_delayed_reduction!(c[i, j], a[i, k], b[k, j], t)
            end
            c[i, j] = reduce!(c[i, j])
        end
    end
    return c
end

function mul!(c::MatElem, a::MatElem, b::RingElement)
    nrows(c) != nrows(a) && error("Incompatible matrix dimensions")

    if c === a || c === b
        d = parent(a)()
        return mul!(d, a, b)
    end

    t = base_ring(a)()
    for i = 1:nrows(a)
        for j = 1:ncols(a)
            c[i, j] = mul!(c[i, j], a[i, j], b)
        end
    end
    return c
end

function add!(c::MatElem, a::MatElem, b::MatElem)
    parent(a) != parent(b) && error("Parents don't match.")
    parent(c) != parent(b) && error("Parents don't match.")
    for i = 1:nrows(c)
        for j = 1:ncols(c)
            c[i, j] = add!(c[i, j], a[i, j], b[i, j])
        end
    end
    return c
end

function addmul!(z::T, x::T, y::T) where {T<:RingElement}
    zz = parent(z)()
    zz = mul!(zz, x, y)
    return addeq!(z, zz)
end

#TODO: should be done in Nemo/AbstractAlgebra s.w.
#      needed by ^ (the generic power in Base using square and multiply)
Base.copy(f::Generic.MPoly) = deepcopy(f)
Base.copy(f::Generic.Poly) = deepcopy(f)
Base.copy(a::PolyRingElem) = deepcopy(a)
Base.copy(a::SeriesElem) = deepcopy(a)

###############################################################################
#
#  Sub
#
###############################################################################

function sub(M::MatElem, rows::Vector{Int}, cols::Vector{Int})
    N = zero_matrix(base_ring(M), length(rows), length(cols))
    for i = 1:length(rows)
        for j = 1:length(cols)
            N[i, j] = M[rows[i], cols[j]]
        end
    end
    return N
end

function sub(M::MatElem{T}, r::AbstractUnitRange{<:Integer}, c::AbstractUnitRange{<:Integer}) where {T}
    z = similar(M, length(r), length(c))
    for i in 1:length(r)
        for j in 1:length(c)
            z[i, j] = M[r[i], c[j]]
        end
    end
    return z
end

function sub(M::Generic.Mat, rows::AbstractUnitRange{Int}, cols::AbstractUnitRange{Int})
    @assert step(rows) == 1 && step(cols) == 1
    z = zero_matrix(base_ring(M), length(rows), length(cols))
    for i in rows
        for j in cols
            z[i-first(rows)+1, j-first(cols)+1] = M[i, j]
        end
    end
    return z
end

gens(L::SimpleNumField{T}) where {T} = [gen(L)]

function gen(L::SimpleNumField{T}, i::Int) where {T}
    i == 1 || error("index must be 1")
    return gen(L)
end

function Base.getindex(L::SimpleNumField{T}, i::Int) where {T}
    if i == 0
        return one(L)
    elseif i == 1
        return gen(L)
    else
        error("index has to be 0 or 1")
    end
end

number_of_generators(L::SimpleNumField{T}) where {T} = 1

is_unit(a::NumFieldElem) = !iszero(a)

canonical_unit(a::NumFieldElem) = a

"""
The coefficients of `f` when viewed as a univariate polynomial in the `i`-th
variable.
"""
function coefficients(f::MPolyRingElem, i::Int)
    d = degree(f, i)
    cf = [MPolyBuildCtx(parent(f)) for j = 0:d]
    for (c, e) = zip(coefficients(f), exponent_vectors(f))
        a = e[i]
        e[i] = 0
        push_term!(cf[a+1], c, e)
    end
    return map(finish, cf)
end

function change_base_ring(p::MPolyRingElem{T}, g, new_polynomial_ring) where {T<:RingElement}
    cvzip = zip(coefficients(p), exponent_vectors(p))
    M = MPolyBuildCtx(new_polynomial_ring)
    for (c, v) in cvzip
        res = g(c)
        if !iszero(res)
            push_term!(M, g(c), v)
        end
    end
    return finish(M)::elem_type(new_polynomial_ring)
end

#check with Nemo/ Dan if there are better solutions
#the block is also not used here I think
#functionality to view mpoly as upoly in variable `i`, so the
#coefficients are mpoly's without variable `i`.
function leading_coefficient(f::MPolyRingElem, i::Int)
    g = MPolyBuildCtx(parent(f))
    d = degree(f, i)
    for (c, e) = zip(coefficients(f), exponent_vectors(f))
        if e[i] == d
            e[i] = 0
            push_term!(g, c, e)
        end
    end
    return finish(g)
end

"""
`content` as a polynomial in the variable `i`, i.e. the gcd of all the
coefficients when viewed as univariate polynomial in `i`.
"""
function content(f::MPolyRingElem, i::Int)
    return reduce(gcd, coefficients(f, i))
end

function content(a::PolyRingElem{<:FieldElem})
    return one(base_ring(a))
end

function canonical_unit(a::SeriesElem)
    iszero(a) && return one(parent(a))
    v = valuation(a)
    v == 0 && return a
    v > 0 && return shift_right(a, v)
    return shift_left(a, -v)
end

# should be Nemo/AA
# TODO: symbols vs strings
#       lift(PolyRing, Series)
#       lift(FracField, Series)
#       (to be in line with lift(ZZ, padic) and lift(QQ, padic)
#TODO: some of this would only work for Abs, not Rel, however, this should be fine here
function map_coefficients(f::T, a::RelPowerSeriesRingElem; parent::SeriesRing) where T
    c = typeof(f(coeff(a, 0)))[]
    for i = 0:pol_length(a)-1
        push!(c, f(polcoeff(a, i)))
    end
    b = parent(c, length(c), precision(a), valuation(a))
    return b
end

#=
function map_coefficients(f, a::RelPowerSeriesRingElem)
  d = f(coeff(a, 0))
  T = parent(a)
  if parent(d) == base_ring(T)
    S = T
  else
    S = power_series_ring(parent(d), max_precision(T), string(var(T)), cached = false)[1]
  end
  c = typeof(d)[d]
  for i=1:pol_length(a)-1
    push!(c, f(polcoeff(a, i)))
  end
  b = S(c, length(c), precision(a), valuation(a))
  return b
end
=#
function lift(R::PolyRing{S}, s::SeriesElem{S}) where {S}
    t = R()
    for x = 0:pol_length(s)
        setcoeff!(t, x, polcoeff(s, x))
    end
    return shift_left(t, valuation(s))
end

#TODO: this is for rings, not for fields, maybe different types?
function Base.gcd(a::T, b::T) where {T<:SeriesElem}
    iszero(a) && iszero(b) && return a
    iszero(a) && return gen(parent(a))^valuation(b)
    iszero(b) && return gen(parent(a))^valuation(a)
    return gen(parent(a))^min(valuation(a), valuation(b))
end

function Base.lcm(a::T, b::T) where {T<:SeriesElem}
    iszero(a) && iszero(b) && return a
    iszero(a) && return a
    iszero(b) && return b
    return gen(parent(a))^max(valuation(a), valuation(b))
end

function gen(R::Union{Generic.EuclideanRingResidueRing{T},Generic.EuclideanRingResidueField{T}}) where {T<:PolyRingElem}
    return R(gen(base_ring(R)))
end

function characteristic(R::Union{Generic.EuclideanRingResidueRing{T},Generic.EuclideanRingResidueField{T}}) where {T<:PolyRingElem}
    return characteristic(base_ring(base_ring(R)))
end

function size(R::Union{Generic.EuclideanRingResidueRing{T},Generic.EuclideanRingResidueField{T}}) where {T<:ResElem}
    return size(base_ring(base_ring(R)))^degree(modulus(R))
end

function size(R::Union{Generic.EuclideanRingResidueRing{T},Generic.EuclideanRingResidueField{T}}) where {T<:PolyRingElem}
    return size(base_ring(base_ring(R)))^degree(R.modulus)
end

function rand(R::Union{Generic.EuclideanRingResidueRing{T},Generic.EuclideanRingResidueField{T}}) where {T<:PolyRingElem}
    r = rand(base_ring(base_ring(R)))
    g = gen(R)
    for i = 1:degree(R.modulus)
        r = r * g + rand(base_ring(base_ring(R)))
    end
    return r
end

function gens(R::Union{Generic.EuclideanRingResidueRing{T},Generic.EuclideanRingResidueField{T}}) where {T<:PolyRingElem} ## probably needs more cases
    ## as the other residue functions
    g = gen(R)
    r = Vector{typeof(g)}()
    push!(r, one(R))
    if degree(R.modulus) == 1
        return r
    end
    push!(r, g)
    for i = 2:degree(R.modulus)-1
        push!(r, r[end] * g)
    end
    return r
end

promote_rule(::Type{LocElem{T}}, ::Type{T}) where {T} = LocElem{T}

promote_rule(::Type{T}, ::Type{S}) where {S<:NumFieldElem,T<:Integer} = S

promote_rule(::Type{S}, ::Type{T}) where {S<:NumFieldElem,T<:Integer} = S

function mulmod(a::S, b::S, mod::Vector{S}) where {S<:MPolyRingElem{T}} where {T<:RingElem}
    return Base.divrem(a * b, mod)[2]
end

Base.:\(f::Map, x) = preimage(f, x)

export set_precision, set_precision!

function set_precision(f::PolyRingElem{T}, n::Int) where {T<:SeriesElem}
    g = parent(f)()
    for i = 0:length(f)
        setcoeff!(g, i, set_precision(coeff(f, i), n))
    end
    return g
end

function set_precision!(f::PolyRingElem{T}, n::Int) where {T<:SeriesElem}
    for i = 0:length(f)
        setcoeff!(f, i, set_precision!(coeff(f, i), n))
    end
    return f
end

function Base.minimum(::typeof(precision), a::Vector{<:SeriesElem})
    return minimum(map(precision, a))
end

function Base.maximum(::typeof(precision), a::Vector{<:SeriesElem})
    return maximum(map(precision, a))
end

#TODO: in Nemo, rename to setprecision
#      fix/report series add for different length
function set_precision(a::SeriesElem, i::Int)
    b = deepcopy(a)
    set_precision!(b, i)
    return b
end

Random.gentype(::Type{T}) where {T<:FinField} = elem_type(T)

import LinearAlgebra
LinearAlgebra.dot(a::NCRingElem, b::NCRingElem) = a * b

transpose!(A::MatrixElem) = transpose(A)

function Base.div(f::PolyRingElem, g::PolyRingElem)
    q, r = divrem(f, g)
    return q
end

function Base.rem(f::PolyRingElem, g::PolyRingElem)
    return mod(f, g)
end

###############################################################################
#
#   Random functions
#
###############################################################################

Random.Sampler(::Type{RNG}, K::FinField, n::Random.Repetition) where {RNG<:AbstractRNG} =
    Random.SamplerSimple(K, Random.Sampler(RNG, BigInt(0):BigInt(characteristic(K) - 1), n))

function rand(rng::AbstractRNG, Ksp::Random.SamplerSimple{<:FinField})
    K = Ksp[]
    r = degree(K)
    alpha = gen(K)
    res = zero(K)
    for i = 0:(r-1)
        c = rand(rng, Ksp.data)
        res += c * alpha^i
    end
    return res
end

function (R::Generic.PolyRing{T})(x::Generic.RationalFunctionFieldElem{T,U}) where {T<:RingElem,U}
    @assert isone(denominator(x))
    @assert parent(numerator(x)) === R
    return numerator(x)
end
function (R::PolyRing{T})(x::Generic.RationalFunctionFieldElem{T,U}) where {T<:RingElem,U}
    @assert isone(denominator(x))
    @assert parent(numerator(x)) === R
    return numerator(x)
end
