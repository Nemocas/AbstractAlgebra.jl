##############################################################################
#
#  Basic manipulation
#
##############################################################################

function is_power(a::PolyRingElem, n::Int)
    # not the best algorithm... but it works generically
    # probably a equal-degree-factorisation would be good + some more gcd's
    # implement some Newton-type algo?
    degree(a) % n == 0 || return false, a
    f = factor(a)
    all(i % n == 0 for (_, i) in f) || return false, a
    fl, x = is_power(constant_coefficient(f.unit), n)
    fl || return false, a
    return true, x*prod(p^div(k, n) for (p, k) = f)
end

################################################################################
#
#  Factorisation
#
################################################################################

function factor(R::Field, f::PolyRingElem)
    Rt = AbstractAlgebra.poly_ring(R)
    f1 = change_base_ring(R, f; parent = Rt)
    return factor(f1)
end

function factor(R::Ring, f::Union{FracElem, Rational})
    fn = factor(R(numerator(f)))
    fd = factor(R(denominator(f)))
    un = divexact(unit(fn), unit(fd))
    fndc = Dict{elem_type(R), Int}(p => e for (p, e) in fn)
    for (k, v) = fd
        if Base.haskey(fndc, k)
            fndc[k] -= v
        else
            fndc[k] = -v
        end
    end
    return Fac(un, sort!(collect(fndc); by = x -> pretty_sort(x[1])))
end

################################################################################
#
#  Roots
#
################################################################################

@doc raw"""
    roots(f::PolyRingElem)

Returns the roots of the polynomial `f` in the base ring of `f` as an array.
"""
function roots(f::PolyRingElem)
    lf = factor(f)
    rts = Vector{elem_type(base_ring(f))}()
    for (p, e) in lf
        if degree(p) == 1
            fl, q = divides(constant_coefficient(p), leading_coefficient(p))
            fl && push!(rts, -q)
        end
    end
    return rts
end

@doc raw"""
    roots(R::Field, f::PolyRingElem)

Returns the roots of the polynomial `f` in the field `R` as an array.
"""
function roots(R::Field, f::PolyRingElem)
    Rt = AbstractAlgebra.poly_ring(R)
    f1 = change_base_ring(R, f, parent = Rt)
    return roots(f1)
end

@doc raw"""
    roots_with_multiplicities(f::PolyRingElem) -> Vector{Tuple{T, Int}}

Return the roots of `f` in its base ring, whose element type is `T`, each with
its multiplicity: the largest $e$ such that $(x - r)^e$ divides `f`. The order
of the roots is not specified. The base ring must be exact, and `f` non-zero.
"""
function roots_with_multiplicities(f::PolyRingElem)
    @req !iszero(f) "Argument must be non-zero"
    @req is_exact_type(elem_type(base_ring(f))) "Multiplicities can only be computed over exact rings"
    x = gen(parent(f))
    # `roots` does not promise to list each root once
    return [(r, valuation(f, x - r)) for r in unique(roots(f))]
end

function roots(a::FinFieldElem, i::Int)
  _, x = polynomial_ring(parent(a), cached = false)
  return roots(x^i-a)
end

function sturm_sequence(f::PolyRingElem{<:FieldElem})
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

