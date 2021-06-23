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
#  Factorisation
#
################################################################################

function factor(f::PolyElem, R::Field)
    Rt = AbstractAlgebra.PolyRing(R)
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
    Rt = AbstractAlgebra.PolyRing(R)
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
 
