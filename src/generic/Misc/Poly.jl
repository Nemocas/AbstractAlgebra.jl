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
    fl, x = ispower(leading_coefficient(a), n)
    fl || return false, a
    f = factor(a)
    all(i -> i % n == 0, values(f.fac)) || return false, a
    return true, x*prod(p^div(k, n) for (p, k) = f.fac)
end

################################################################################
#
#  Power sums
#
################################################################################

function polynomial_to_power_sums(f::PolyElem{T}, n::Int=degree(f)) where T <: FieldElement
    degree(f) < 1 && error("Polynomial has no roots")
    !ismonic(f) && error("Requires monic polynomial")
    iszero(constant_coefficient(f)) && error("Requires nonzero constant coefficient")
    n < 0 && throw(DomainError(n, "number of terms must be nonnegative"))
    d = degree(f)
    R = base_ring(f)
    # Beware: converting to power series and derivative do not commute
    dfc = collect(Base.Iterators.take(reverse!(collect(coefficients(derivative(f)))), n + 1))
    A = abs_series(R, dfc, length(dfc), n + 1; cached=false)
    S = parent(A)
    fc = collect(Base.Iterators.take(reverse!(collect(coefficients(f))), n + 1))
    B = S(fc, length(fc), n + 1)
    L = A*inv(B)
    s = T[coeff(L, i) for i = 1:n]
    return s
end

@doc Markdown.doc"""
    polynomial_to_power_sums(f::PolyElem{T}, n::Int=degree(f)) where T <: RingElement -> Array{T, 1}

Uses Newton (or Newton-Girard) formulas to compute the first $n$
sums of powers of the roots of $f$ from the coefficients of $f$, starting
with the sum of (first powers of) the roots. The input polynomial must be
monic, at least degree $1$ and have nonzero constant coefficient.
"""
function polynomial_to_power_sums(f::PolyElem{T}, n::Int=degree(f)) where T <: RingElement
    # plain vanilla recursion
    degree(f) < 1 && error("Polynomial has no roots")
    !ismonic(f) && error("Requires monic polynomial")
    iszero(constant_coefficient(f)) && error("Requires nonzero constant coefficient")
    n < 0 && throw(DomainError(n, "number of terms must be nonnegative"))
    d = degree(f)
    R = base_ring(f)
    if n == 0
       return elem_type(R)[]
    end
    if n == 1
       return [-coeff(f, d - 1)]
    end
    E = T[(-1)^i*coeff(f, d - i) for i = 0:min(d, n)] # elementary symm. polys
    while length(E) <= n
        push!(E, R())
    end
    P = T[]
    push!(P, E[1 + 1])
    for k = 2:n
        push!(P, (-1)^(k - 1)*k*E[k + 1] +
		 sum((-1)^(k - 1 + i)*E[k - i + 1]*P[i] for i = 1:k - 1))
    end
    return P
end

@doc Markdown.doc"""
    power_sums_to_polynomial(P::Array{T, 1}; parent::AbstractAlgebra.PolyRing{T}=
           AbstractAlgebra.PolynomialRing(parent(P[1])) where T <: RingElement -> PolyElem{T}

Uses the Newton (or Newton-Girard) identities to obtain the polynomial
with given sums of powers of roots. The list must be nonempty and contain
`degree(f)` entries where $f$ is the polynomial to be recovered. The list
must start with the sum of first powers of the roots.
"""
function power_sums_to_polynomial(P::Array{T, 1}; 
                          parent::AbstractAlgebra.PolyRing{T}=
	    AbstractAlgebra.PolynomialRing(parent(P[1]))) where T <: RingElement
   return power_sums_to_polynomial(P, parent)
end

function power_sums_to_polynomial(P::Array{T, 1}, Rx::AbstractAlgebra.PolyRing{T}) where T <: FieldElement
    d = length(P)
    R = base_ring(Rx)
    s = rel_series(R, P, d, d, 0)
    r = -integral(s)
    r1 = exp(r)
    @assert iszero(valuation(r1))
    return Rx([polcoeff(r1, d - i) for i = 0:d])
end

function power_sums_to_polynomial(P::Array{T, 1}, Rx::AbstractAlgebra.PolyRing{T}) where T <: RingElement
    E = T[one(parent(P[1]))]
    R = parent(P[1])
    last_non_zero = 0
    for k = 1:length(P)
        push!(E, divexact(sum((-1)^(i - 1)*E[k - i + 1]*P[i] for i = 1:k), R(k)))
        if E[end] != 0
            last_non_zero = k
        end
    end
    E = E[1:last_non_zero + 1]
    d = length(E) # the length of the resulting polynomial
    for i = 1:div(d, 2)
        E[i], E[d - i + 1] = (-1)^(d - i)*E[d - i + 1], (-1)^(i - 1)*E[i]
    end
    if isodd(d)
        E[div(d + 1, 2)] *= (-1)^div(d, 2)
    end
    return Rx(E)
end

################################################################################
#
#  Factorisation
#
################################################################################

function factor(f::PolyElem, R::Field)
    Rt = AbstractAlgebra.PolynomialRing(R)
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

@doc doc"""
    roots(f::PolyElem)

Returns the roots of the polynomial `f` in the base ring of `f` as an array.
"""
function roots(f::PolyElem)
    lf = factor(f)
    rts = Vector{elem_type(base_ring(f))}()
    for (p, e) in lf
        if degree(p) == 1
            push!(rts, -divexact(constant_coefficient(p), leading_coefficient(p)))
        end
    end
    return rts
end

@doc doc"""
    roots(f::PolyElem, R::Field)

Returns the roots of the polynomial `f` in the field `R` as an array.
"""
function roots(f::PolyElem, R::Field)
    Rt = AbstractAlgebra.PolynomialRing(R)
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
 
