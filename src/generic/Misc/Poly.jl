##############################################################################
#
#  Basic manipulation
#
##############################################################################

function ispower(a::PolyElem, n::Int)
    # not the best algorithm... but it works generically
    # probably a equal-degree-factorisation would be good + some more gcd's
    # implement some Newton-type algo?
    degree(a) % n == 0 || return false, a
    fl, x = ispower(lead(a), n)
    fl || return false, a
    f = factor(a)
    all(i -> i % n == 0, values(f.fac)) || return false, a
    return true, x*prod(p^div(k, n) for (p, k) = f.fac)
end

@doc Markdown.doc"""
    valence(f::PolyElem) -> RingElem

 The last non-zero coefficient of $f$.
"""
function valence(f::PolyElem)
    for i = 0:degree(f)
        c = coeff(f, i)
        if !iszero(c)
            return c
        end
    end
    return c
end

################################################################################
#
#  Iterators
#
################################################################################

struct PolyCoeffs{T <: RingElem}
    f::T
end
  
function coefficients(f::PolyElem)
    return PolyCoeffs(f)
end
  
function Base.iterate(PC::PolyCoeffs{<:PolyElem}, st::Int = -1)
    st += 1
    if st > degree(PC.f)
        return nothing
    else
        return coeff(PC.f, st), st
    end
end
  
Base.IteratorEltype(M::PolyElem) = Base.HasEltype()

Base.eltype(M::PolyElem{T}) where {T} = T
  
Base.IteratorSize(M::PolyCoeffs{<:PolyElem}) = Base.HasLength()

Base.length(M::PolyCoeffs{<:PolyElem}) = degree(M.f) + 1
  
function Base.lastindex(a::PolyCoeffs{<:PolyElem})
    return degree(a.f)
end
  
function Base.getindex(a::PolyCoeffs{<:PolyElem}, i::Int)
    return coeff(a.f, i)
end

##############################################################################
#
#  Binary arithmetic
#
##############################################################################

# computes the top n coeffs of the product only
function mulhigh_n(a::PolyElem{T}, b::PolyElem{T}, n::Int) where {T}
    # sum a_i t^i and sum b_j t^j
    # want (i, j) s.th. i + j >= deg a + deg b - n
    r = parent(a)()
    for i=max(degree(a)-n, 0):degree(a)
        for j = max(degree(a) + degree(b) - n - i, 0):degree(b)
            setcoeff!(r, i+j, coeff(r, i+j) + coeff(a, i)*coeff(b, j))
        end
    end
    return r
end

# assuming b divides a, compute the last n coeffs of the quotient
# will produce garbage otherwise
# div(a, b) mod x^n
function divexact_low(a::PolyElem{T}, b::PolyElem{T}, n::Int) where {T}
    r = parent(a)()
    a = truncate(a, n)
    b = truncate(b, n)
    for i = 0:n - 1
        q = divexact(trail(a), trail(b))
        setcoeff!(r, i, q)
        a = shift_right(a - q*b, 1)
        b = truncate(b, n - i - 1)
        # truncate both a and b to n - i - 1 (for generic polys one could just change the length)
    end
    return r
end

#computes the top coeffs starting with x^n
function divhigh(a::PolyElem{T}, b::PolyElem{T}, n0::Int) where {T}
    r = parent(a)()
    n = degree(a) - degree(b) - n0
    fit!(r, degree(a) - degree(b))
    a = deepcopy(a)
    k = degree(a) - n0
    da = degree(a)
    for i=0:n
        if degree(a) < degree(b)
            break
        end
        q = divexact(coeff(a, da), lead(b))
        setcoeff!(r, da - degree(b), q)
        for j=da:-1:max(k, da - degree(b))
            setcoeff!(a, j, coeff(a, j)-q*coeff(b, j-da+degree(b)))
        end
        da -= 1
    end
    if iszero(r)
        return r
    end
    set_length!(r, normalise(r, length(r) - 1))
    return r
end

################################################################################
#
#  Deflate and inflate
#
################################################################################

@doc Markdown.doc"""
    inflate(f::PolyElem, n::Int64) -> PolyElem

Given a polynomial $f$ in $x$, return $f(x^n)$, i.e. multiply
all exponents by $n$.
"""
function inflate(x::PolyElem, n::Int64)
    y = parent(x)()
    for i = 0:degree(x)
        setcoeff!(y, n*i, coeff(x, i))
    end
    return y
end

@doc Markdown.doc"""
    deflate(f::PolyElem, n::Int64) -> PolyElem

Given a polynomial $f$ in $x^n$, write it as a polynomial in $x$, i.e. divide
all exponents by $n$.
"""
function deflate(x::PolyElem, n::Int64)
    y = parent(x)()
    for i = 0:div(degree(x), n)
        setcoeff!(y, i, coeff(x, n*i))
    end
    return y
end

@doc Markdown.doc"""
    deflate(x::PolyElem) -> PolyElem

Deflate the polynomial $f$ maximally, i.e. find the largest $n$ s.th.
$f$ can be deflated by $n$, i.e. $f$ is actually a polynomial in $x^n$.
"""
function deflate(x::PolyElem)
    g = 0
    for i = 0:degree(x)
        if coeff(x, i) != 0
            g = gcd(g, i)
            if g==1
                return x, 1
            end
        end
    end
    return deflate(x, g), g
end

################################################################################
#
#  Power sums
#
################################################################################

@doc Markdown.doc"""
    polynomial_to_power_sums(f::PolyElem{T}, n::Int=degree(f)) -> Array{T, 1}

Uses Newton (or Newton-Girard) formulas to compute the first $n$
power sums from the coefficients of $f$.
"""
function polynomial_to_power_sums(f::PolyElem{T}, n::Int=degree(f)) where T <: FieldElem
    d = degree(f)
    R = base_ring(f)
    S = AbstractAlgebra.PowerSeriesRing(R, n + 1, "gen(S)", cached = false, model =:capped_absolute)[1]
    # Beware: converting to power series and derivative do not commute
    A = S()
    B = S()
    fit!(A, d)
    fit!(B, d + 1)
    for i=1:d
        c = coeff(f, i)
        setcoeff!(A, d - i, i*c)
        setcoeff!(B, d - i, c)
    end
    setcoeff!(B, d, coeff(f, 0))
    A.prec = n + 1
    B.prec = n + 1
    L = A*inv(B)
    s = T[coeff(L, i) for i = 1:n]
    return s
end

# plain vanilla recursion
function polynomial_to_power_sums(f::PolyElem{T}, n::Int=degree(f)) where T
    if n == 0
        return elem_type(base_ring(f))[]
    end
    d = degree(f)
    R = base_ring(f)
    E = T[(-1)^i*coeff(f, d-i) for i=0:min(d, n)] #should be the elementary symm.
    while length(E) <= n
        push!(E, R(0))
    end
    P = T[]
    push!(P, E[1+1])
    for k = 2:n
        push!(P, (-1)^(k-1)*k*E[k+1] + sum((-1)^(k-1+i)*E[k-i+1]*P[i] for i=1:k-1))
    end
    return P
end

@doc Markdown.doc"""
    power_sums_to_polynomial(P::Array{T, 1}) -> PolyElem{T}

Uses the Newton (or Newton-Girard) identities to obtain the polynomial
coefficients (the elementary symmetric functions) from the power sums.
"""
function power_sums_to_polynomial(P::Array{T, 1}) where T <: FieldElem
    d = length(P)
    R = parent(P[1])
    S = AbstractAlgebra.PowerSeriesRing(R, d, "gen(S)")[1] # capped_absolute
    s = S(P, length(P), d, 0)
    r = -integral(s)
    r1 = exp(r)
    @assert iszero(valuation(r1))
    Rx, x = AbstractAlgebra.PolynomialRing(R, "x", cached = false)
    return Rx([polcoeff(r1, d - i) for i = 0:d])
end

function power_sums_to_polynomial(P::Array{T, 1}) where T
    E = T[one(parent(P[1]))]
    R = parent(P[1])
    last_non_zero = 0
    for k = 1:length(P)
        push!(E, divexact(sum((-1)^(i-1)*E[k-i+1]*P[i] for i=1:k), R(k)))
        if E[end] != 0
            last_non_zero = k
        end
    end
    E = E[1:last_non_zero+1]
    d = length(E) # the length of the resulting poly...
    for i = 1:div(d, 2)
        E[i], E[d - i + 1] = (-1)^(d - i)*E[d - i + 1], (-1)^(i - 1)*E[i]
    end
    if isodd(d)
        E[div(d + 1, 2)] *= (-1)^div(d, 2)
    end
    return AbstractAlgebra.PolynomialRing(R, "x", cached = false)[1](E)
end

################################################################################
#
#  Factorisation
#
################################################################################

function factor(f::PolyElem, R::Field)
    Rt = AbstractAlgebra.PolynomialRing(R, "t", cached = false)[1]
    f1 = change_base_ring(R, f, parent = Rt)
    return factor(f1)
end

function factor(f::FracElem, R::Ring)
    fn = factor(R(numerator(f)))
    fd = factor(R(denominator(f)))
    fn.unit = divexact(fn.unit, fd.unit)
    for (k,v) = fd.fac
        if Base.haskey(fn.fac, k)
            fn.fac[k] -= v
        else
            fn.fac[k] = -v
        end
    end
    return fn
end

################################################################################
#
#  Roots
#
################################################################################

function roots(f::PolyElem)
    lf = factor(f)
    rts = Vector{elem_type(base_ring(f))}()
    for (p, e) in lf
        if degree(p) == 1
            push!(rts, -divexact(trail(p), lead(p)))
        end
    end
    return rts
end

function roots(f::PolyElem, R::Field)
    Rt = AbstractAlgebra.PolynomialRing(R, "t", cached = false)[1]
    f1 = change_base_ring(R, f, parent = Rt)
    return roots(f1)
end

function sturm_sequence(f::PolyElem{<:FieldElem})
    g = f
    h = derivative(g)
    seq = typeof(f)[g,h]
    while true
        r = rem(g, h)
        if r != 0
            push!(seq, -r)
            g, h = h, -r
        else 
            break
        end
    end
    return seq
end
  
################################################################################
#
#  Random polynomial
#
################################################################################

@doc Markdown.doc"""
    Base.rand(Rt::PolyRing{T}, n::Int) where T <: ResElem{S} where S <: Integer -> PolyElem{T}

Find a random polynomial of degree=$n$.
"""
function Base.rand(Rt::PolyRing{T}, n::Int) where T <: ResElem{S} where S <: Integer
    f = Rt()
    R = base_ring(Rt)
    for i = 0:n
        setcoeff!(f, i, rand(R))
    end
    return f
end
