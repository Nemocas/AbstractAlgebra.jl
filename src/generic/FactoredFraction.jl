###############################################################################
#
#   FactoredFraction.jl : fraction fields in factored form
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

function parent(a::FactoredFracFieldElem{T}) where T <: RingElement
    return a.parent
end

function parent_type(::Type{FactoredFracFieldElem{T}}) where {T <: RingElement}
    return FactoredFracField{T}
end

function elem_type(::Type{FactoredFracField{T}}) where {T <: RingElement}
    return FactoredFracFieldElem{T}
end

base_ring_type(::Type{FactoredFracField{T}}) where T <: NCRingElement = parent_type(T)

function base_ring(F::FactoredFracField{T}) where T <: RingElement
    return F.base_ring::parent_type(T)
end

function characteristic(F::FactoredFracField{T}) where T <: RingElement
   return characteristic(base_ring(F))
end

###############################################################################
#
#   Constructors and Parent Call Overloads
#
###############################################################################

function (F::FactoredFracField{T})() where T
    return _make_base_elem(F, zero(base_ring(F)))
end

function (F::FactoredFracField{T})(a::T) where T <: RingElem
    parent(a) == base_ring(F) || error("Could not coerce into $F")
    return _make_base_elem(F, a)
end

function (F::FactoredFracField{T})(a::Integer) where T <: RingElem
    return _make_base_elem(F, base_ring(F)(a))
end

function (F::FactoredFracField{T})(a::Integer) where T <: Integer
    return _make_base_elem(F, base_ring(F)(a))
end

function (F::FactoredFracField{T})(a::Rational) where T <: Integer
    return F(numerator(a), denominator(a))
end

function (F::FactoredFracField{T})(a::Rational) where T <: RingElem
    return F(numerator(a), denominator(a))
end

function (F::FactoredFracField{T})(a::AbstractAlgebra.Generic.FracFieldElem{T}) where T <: RingElement
    base_ring(F) == base_ring(a) || error("Could not coerce into $F")
    return _append_pow!(_make_base_elem(F, numerator(a)), denominator(a), -1)
end

function (F::FactoredFracField{T})(a::FactoredFracFieldElem{T}) where T
    F == parent(a) || error("Could not coerce into $F")
    return a
end

# construction of the fraction a/b
function (F::FactoredFracField{T})(a::RingElement, b::RingElement) where T <: RingElement
    a = base_ring(F)(a)
    b = base_ring(F)(b)
    return _append_pow!(_make_base_elem(F, a), b, -1)
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.deepcopy_internal(a::FactoredFracFieldElem{T}, dict::IdDict) where T <: RingElement
   return FactoredFracFieldElem{T}(deepcopy_internal(a.unit, dict),
                          deepcopy_internal(a.terms, dict),
                          a.parent)
end

# the non-expanding hash function would have to normalise the bases, and then
# either sort the bases or combine the hashes of the bases in a commutative way

function Base.numerator(a::FactoredFracFieldElem, canonicalise::Bool=true)
    z = unit(a)
    for (b, e) in a
        if canonicalise
            u = canonical_unit(b)
            z *= _pow(u, e)
            b = divexact(b, u)
        end
        if e > 0
            z *= b^e
        end
    end
    return z
end

function Base.denominator(a::FactoredFracFieldElem, canonicalise::Bool=true)
    z = one(base_ring(a))
    for (b, e) in a
        if e < 0
            e = Base.checked_neg(e)
            if canonicalise
                b = divexact(b, canonical_unit(b))
            end
            z *= b^e
        end
    end
    return z
end

function one(F::FactoredFracField{T}) where T
    FactoredFracFieldElem{T}(one(base_ring(F)), FactoredFracTerm{T}[], F)
end

function zero(F::FactoredFracField{T}) where T
    FactoredFracFieldElem{T}(zero(base_ring(F)), FactoredFracTerm{T}[], F)
end

function is_unit(a::FactoredFracFieldElem{T}) where T
    return !iszero(a)
end

function iszero(a::FactoredFracFieldElem{T}) where T
    return iszero(a.unit)
end

function isone(a::FactoredFracFieldElem{T}) where T
    iszero(a.unit) && false
    for i in a.terms
        # if some base appears to a non-zero power and is not a unit and is
        # relatively prime to all other bases, then a cannot be one.
        ok = !iszero(i.exp) && !is_unit(i.base)
        for j in a.terms
            if j !== i
                ok = ok && is_unit(gcd(i.base, j.base))
            end
        end
        ok && return false
    end
    z = normalise(a)
    return isempty(z.terms) && isone(z.unit)
end

###############################################################################
#
#   Fac Interface
#
###############################################################################

function unit(a::FactoredFracFieldElem)
    return a.unit
end

function Base.iterate(a::FactoredFracFieldElem)
    t = Base.iterate(a.terms)
    return isnothing(t) ? t : ((t[1].base, t[1].exp), t[2])
end

function Base.iterate(a::FactoredFracFieldElem, b)
    t = Base.iterate(a.terms, b)
    return isnothing(t) ? t : ((t[1].base, t[1].exp), t[2])
end

function Base.eltype(::Type{FactoredFracFieldElem{T}}) where T
    return Tuple{T, Int}
end

function Base.length(a::FactoredFracFieldElem)
    Base.length(a.terms)
end

function push_term!(a::FactoredFracFieldElem{T}, b::T, e::Int) where T <: RingElement
    parent(b) == base_ring(a) || error("Incompatible parents")
    push!(a.terms, FactoredFracTerm{T}(b, e))
    return a
end

function push_term!(a::FactoredFracFieldElem{T}, b, e::Int) where T <: RingElement
    return push_term!(a, base_ring(a)(b), e)
end

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function expressify(a::FactoredFracFieldElem; context = nothing)
    n = Expr(:call, :*, expressify(a.unit, context = context))
    d = Expr(:call, :*)
    for t in a.terms
        b = expressify(t.base; context = context)
        e = Base.checked_abs(t.exp)
        push!((e == t.exp ? n : d).args, isone(e) ? b : Expr(:call, :^, b, e))
    end
    return length(d.args) < 2 ? n : Expr(:call, :/, n, d)
end

@enable_all_show_via_expressify FactoredFracFieldElem

function expressify(a::FactoredFracField; context = nothing)
    return Expr(:sequence, Expr(:text, "Factored fraction field of "),
                           expressify(base_ring(a); context = context))
end

@enable_all_show_via_expressify FactoredFracField

###############################################################################
#
#   Comparison, Addition, Subtraction, GCD
#
###############################################################################

function ==(a::FactoredFracFieldElem{T}, b::FactoredFracFieldElem{T}) where T
    (g, x, y) = _gcdhelper(a, b)
    return x == y
end

function +(a::FactoredFracFieldElem{T}, b::Rational) where T <: RingElem
   return a + parent(a)(b)
end

function +(a::FactoredFracFieldElem{T}, b::T) where T <: RingElem
   return a + parent(a)(b)
end

function +(b::Rational, a::FactoredFracFieldElem{T}) where T <: RingElem
   return parent(a)(b) + a
end

function +(b::T, a::FactoredFracFieldElem{T}) where T <: RingElem
   return parent(a)(b) + a
end

function +(a::FactoredFracFieldElem{T}, b::FactoredFracFieldElem{T}) where T
    (g, x, y) = _gcdhelper(a, b)
    return mul_by_base_elem(g, x + y)
end

function -(a::FactoredFracFieldElem{T}, b::Rational) where T <: RingElem
   return a - parent(a)(b)
end

function -(a::FactoredFracFieldElem{T}, b::T) where T <: RingElem
   return a - parent(a)(b)
end

function -(b::Rational, a::FactoredFracFieldElem{T}) where T <: RingElem
   return parent(a)(b) - a
end

function -(b::T, a::FactoredFracFieldElem{T}) where T <: RingElem
   return parent(a)(b) - a
end

function -(a::FactoredFracFieldElem{T}, b::FactoredFracFieldElem{T}) where T
    (g, x, y) = _gcdhelper(a, b)
    return mul_by_base_elem(g, x - y)
end

function -(a::FactoredFracFieldElem{T}) where T
    return FactoredFracFieldElem{T}(-a.unit, a.terms, a.parent)
end

function gcd(a::FactoredFracFieldElem{T}, b::FactoredFracFieldElem{T}) where T
    (g, x, y) = _gcdhelper(a, b)
    return mul_by_base_elem(g, gcd(x, y))
end

###############################################################################
#
#   Multiplication
#
###############################################################################

function mul_by_base_elem(a::FactoredFracFieldElem{T}, b::T) where T <: RingElement
    F = parent(a)
    if iszero(b)
        return zero(F)
    else
        z = FactoredFracFieldElem{T}(a.unit, map(copy, a.terms), F)
        _append_pow_normalise!(z, b, 1, 1)
        return z
    end
end

function *(a::FactoredFracFieldElem{T}, b::T) where T <: RingElem
    parent(b) == base_ring(a) || error("Incompatible rings")
    return mul_by_base_elem(a, b)
end

function *(a::FactoredFracFieldElem{T}, b::Integer) where T <: RingElem
    return mul_by_base_elem(a, base_ring(a)(b))
end

function *(b::T, a::FactoredFracFieldElem{T}) where T <: RingElem
   return mul_by_base_elem(a, b)
end

function *(b::Integer, a::FactoredFracFieldElem{T}) where T <: RingElem
   return mul_by_base_elem(a, base_ring(a)(b))
end

function *(b::FactoredFracFieldElem{T}, c::FactoredFracFieldElem{T}) where T <: RingElement
    parent(b) == parent(c) || error("Incompatible rings")
    input_is_good = _bases_are_coprime(b) && _bases_are_coprime(b)
    z = FactoredFracFieldElem{T}(b.unit*c.unit, FactoredFracTerm{T}[], parent(b))
    if iszero(z.unit)
        return z
    end
    b = map(copy, b.terms)
    c = map(copy, c.terms)
    i = 1
    while i <= length(b)
        j = 1
        while j <= length(c)
            (g, b[i].base, c[j].base) = gcd_with_cofactors(b[i].base, c[j].base)
            if is_unit(g)
                z.unit *= _pow(g, Base.checked_add(b[i].exp, c[j].exp))
            else
                b[i].base = _append_coprimefac!(z, b[i].base, b[i].exp, g)
                c[j].base = _append_coprimefac!(z, c[j].base, c[j].exp, g)
            end
            if is_unit(c[j].base)
                z.unit *= _pow(c[j].base, c[j].exp)
                c[j] = c[end]
                pop!(c)
            else
                j += 1
            end
        end
        if is_unit(b[i].base)
            z.unit *= _pow(b[i].base, b[i].exp)
            b[i] = b[end]
            pop!(b)
        else
            i += 1
        end
    end
    _append_pow!(z, b, 1)
    _append_pow!(z, c, 1)
    @assert !input_is_good || _bases_are_coprime(z)
    return z
end

###############################################################################
#
#   Division
#
###############################################################################

function Base.inv(a::FactoredFracFieldElem{T}) where T
    z = FactoredFracTerm{T}[]
    for i in a.terms
        push!(z, FactoredFracTerm{T}(i.base, Base.checked_neg(i.exp)))
    end
    return FactoredFracFieldElem{T}(inv(a.unit), z, parent(a))
end

function divexact(a::FactoredFracFieldElem{T}, b::FactoredFracFieldElem{T}; check::Bool = true) where T
    return a*inv(b)
end

function divexact(a::Integer, b::FactoredFracFieldElem{T}; check::Bool = true) where T
    return divexact(parent(b)(a), b, check = check)
end

function divexact(a::T, b::FactoredFracFieldElem{T}; check::Bool = true) where T <: RingElem
    return divexact(parent(b)(a), b, check = check)
end

#=
function divexact(a::FactoredFracFieldElem{T}, b::Rational; check::Bool = true) where T <: RingElem
   return divexact(a, parent(a)(b), check = check)
end

function divexact(a::Rational, b::FactoredFracFieldElem{T}; check::Bool = true) where T <: RingElem
   return divexact(parent(b)(a), b, check = check)
end
=#

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::FactoredFracFieldElem{T}, b::Int) where T
    z = FactoredFracTerm{T}[]
    if !iszero(b)
        for t in a.terms
            push!(z, FactoredFracTerm(t.base, Base.checked_mul(b, t.exp)))
        end
    end
    return FactoredFracFieldElem{T}(_pow(a.unit, b), z, parent(a))
end

##############################################################################
#
#  Evaluation
#
##############################################################################

function evaluate(f::FactoredFracFieldElem{T}, v::Vector{U}) where {T <: RingElement, U <: RingElement}
    z = evaluate(unit(f), v)
    for (b, e) in f
        z *= evaluate(b, v)^e
    end
    return z
end

function evaluate(f::FactoredFracFieldElem{T}, v::U) where {T <: RingElement, U <: RingElement}
    z = evaluate(unit(f), v)
    for (b, e) in f
        z *= evaluate(b, v)^e
    end
    return z
end

##############################################################################
#
#  Derivative
#
##############################################################################

# Return the derivative with respect to the `i`-th variable.
function derivative(a::FactoredFracFieldElem{T}, i::Int) where {T <: MPolyRingElem}
    z = FactoredFracFieldElem{T}(one(base_ring(a)), FactoredFracTerm{T}[], parent(a))
    p = unit(a)
    for (b, e) in a
        p *= b
        if !isone(e)
            push!(z.terms, FactoredFracTerm{T}(b, Base.checked_sub(e, 1)))
        end
    end
    s = divexact(p, unit(a))*derivative(unit(a), i)
    for (b, e) in a
        s += divexact(p, b)*e*derivative(b, i)
    end
    return _append_pow_normalise!(z, s::T, 1, 1)
end

################################################################################
#
#   Remove and valuation
#
################################################################################

function remove(a::FactoredFracFieldElem{T}, p::T) where T <: RingElement
    z = FactoredFracFieldElem{T}(unit(a), FactoredFracTerm{T}[], parent(a))
    v = 0
    for (b, e) in a
        v1, b1 = remove(b, p)
        v += e*v1
        _append_pow!(z, b1, e)
    end
    return (v, z)
end

###############################################################################
#
#   Unsafe operators and functions: These are defined generically on FracElem
#   in terms of the numerator/denominator interface, so we need some overrides.
#
###############################################################################

function zero!(c::FactoredFracFieldElem)
    c.unit = zero!(c.unit)
    empty!(c.terms)
    return c
end

function mul!(c::FactoredFracFieldElem{T}, a::FactoredFracFieldElem{T}, b::FactoredFracFieldElem{T}) where T <: RingElement
    return a*b
end

function addeq!(a::FactoredFracFieldElem{T}, b::FactoredFracFieldElem{T}) where T <: RingElement
    return a + b
end

function add!(c::FactoredFracFieldElem{T}, a::FactoredFracFieldElem{T}, b::FactoredFracFieldElem{T}) where T <: RingElement
    return a + b
end

###############################################################################
#
#   Internal helpers
#
###############################################################################

function _make_base_elem(F::FactoredFracField{T}, a::T) where T <: RingElement
    if iszero(a) || is_unit(a)
        return FactoredFracFieldElem{T}(a, FactoredFracTerm{T}[], F)
    else
        return FactoredFracFieldElem{T}(one(base_ring(F)), [FactoredFracTerm{T}(a, 1)], F)
    end
end

# return a version of a with coprime bases
function normalise(a::FactoredFracFieldElem{T}) where T
    z = FactoredFracFieldElem{T}(a.unit, FactoredFracTerm{T}[], parent(a))
    if !iszero(z.unit)
        for i in a.terms
            _append_pow_normalise!(z, i.base, i.exp, 1)
        end
    end
    return z
end

# this power tries to invert for negative exponents (and not just throw)
function _pow(a::T, e::Int) where T <: RingElement
    if e < 0
        return divexact(one(parent(a)), a)^Base.checked_neg(e)
    else
        return a^e
    end
end

function copy(a::FactoredFracTerm{T}) where T
    return FactoredFracTerm{T}(a.base, a.exp)
end

function _bases_are_coprime(a::FactoredFracFieldElem{T}) where T
    a = a.terms
    i = 1
    while i <= length(a)
        j = i + 1
        while j <= length(a)
            if !is_unit(gcd(a[i].base, a[j].base))
                return false
            end
            j += 1
        end
        i += 1
    end
    return true
end

# bases are coprime and none are units
function _bases_are_nice(a::FactoredFracFieldElem{T}) where T
    if !_bases_are_coprime(a)
        return false
    end
    for i in a.terms
        if is_unit(i.base)
            return false
        end
    end
    return true
end

# z *= f^e with no extra normalization
function _append_pow!(z::FactoredFracFieldElem{T}, f::Vector{FactoredFracTerm{T}}, e::Int) where T
    for i in f
        ie = Base.checked_mul(i.exp, e)
        if !iszero(ie)
            push!(z.terms, FactoredFracTerm{T}(i.base, ie))
        end
    end
    return z
end

function _append_pow!(z::FactoredFracFieldElem{T}, f::T, e::Int) where T
    if !iszero(e)
        if is_unit(f)
            z.unit *= _pow(f, e)
        else
            push!(z.terms, FactoredFracTerm{T}(f, e))
        end
    end
    return z
end

# z *= a^e with normalization
function _append_pow_normalise!(z::FactoredFracFieldElem{T}, a::T, e::Int, i::Int) where T
    iszero(e) && return z
    if iszero(a)
        z.unit = a
        empty!(z.terms)
        return z
    end
    input_is_good = _bases_are_coprime(z)
    l = z.terms
    while i <= length(l) && !is_unit(a)
        (g, lbar, abar) = gcd_with_cofactors(l[i].base, a)
        # (g*lbar)^l[i].exp * (g*abar)^e
        # (lbar)^l[i].exp * (g)^(l[i].exp+e) * (abar)^e
        if is_unit(g)
            i += 1
        elseif is_unit(lbar)
            a = abar
            z.unit *= _pow(lbar, l[i].exp)
            l[i].base = g
            l[i].exp = Base.checked_add(l[i].exp, e)
            if iszero(l[i].exp)
                l[i] = l[end]
                pop!(l)
            end
        elseif is_unit(abar)
            z.unit *= _pow(abar, e)
            l[i].base = lbar
            a = g
            e = Base.checked_add(l[i].exp, e)
            if iszero(e)
                return
            end
        else
            a = abar
            _append_pow_normalise!(z, g, e, i)
        end
    end
    if is_unit(a)
        z.unit *= _pow(a, e)
    else
        push!(l, FactoredFracTerm{T}(a, e))
    end
    @assert !input_is_good || _bases_are_coprime(z)
    return z
end

# multiply l by p^e and return a such that
# p*a = b*g and gcd(a, g) == 1
function _append_coprimefac!(l::FactoredFracFieldElem{T}, b::T, e::Int, g::T) where T
    @assert !iszero(b)
    _append_pow_normalise!(l, g, e, 1)
    a = b
    r = gcd(a, g)
    while !is_unit(r)
        @assert !iszero(r)
        _append_pow_normalise!(l, r, e, 1)
        a = divexact(a, r)
        r = gcd(r, a)
    end
    return a
end

# pull out some common factor b = z*bbar, c = z*cbar
# where bbar and cbar are elements of the base ring.
function _gcdhelper(b::FactoredFracFieldElem{T}, c::FactoredFracFieldElem{T}) where T
    F = parent(b)
    z = FactoredFracFieldElem(one(base_ring(F)), FactoredFracTerm{T}[], F)   # gcd(b,c)
    bbar = b.unit                       # expanded form of eventual b/gcd(b,c)
    cbar = c.unit                       # expanded form of eventual c/gcd(b,c)
    b = map(copy, b.terms)
    c = map(copy, c.terms)
    i = 1
    while i <= length(b)
        j = 1
        while j <= length(c)
            (g, b[i].base, c[j].base) = gcd_with_cofactors(b[i].base, c[j].base)
            if !is_unit(g)
                e = Base.checked_sub(b[i].exp, c[j].exp)
                if e >= 0
                    _append_pow_normalise!(z, g, c[j].exp, 1)
                    if e > 0
                        push!(b, FactoredFracTerm{T}(g, e))
                    end
                else
                    _append_pow_normalise!(z, g, b[i].exp, 1)
                    push!(c, FactoredFracTerm{T}(g, Base.checked_neg(e)))
                end
            end
            if is_unit(c[j].base)
                cbar *= _pow(c[j].base, c[j].exp)
                c[j] = c[end]
                pop!(c)
            else
                j += 1
            end
        end
        if is_unit(b[i].base)
            bbar *= _pow(b[i].base, b[i].exp)
            b[i] = b[end]
            pop!(b)
        else
            i += 1
        end
    end
    i = 1
    while i <= length(b)
        if b[i].exp < 0
            _append_pow_normalise!(z, b[i].base, b[i].exp, 1)
            cbar *= b[i].base^Base.checked_neg(b[i].exp)
        else
            bbar *= b[i].base^b[i].exp
        end
        i += 1
    end
    j = 1
    while j <= length(c)
        if c[j].exp < 0
            _append_pow_normalise!(z, c[j].base, c[j].exp, 1)
            bbar *= c[j].base^Base.checked_neg(c[j].exp)
        else
            cbar *= c[j].base^c[j].exp
        end
        j += 1
    end
    return (z, bbar, cbar)
end

###############################################################################
#
#   FactoredFractionField constructor
#
###############################################################################

function FactoredFractionField(R::AbstractAlgebra.Ring; cached::Bool=true)
   return FactoredFracField{AbstractAlgebra.elem_type(R)}(R, cached)
end 

