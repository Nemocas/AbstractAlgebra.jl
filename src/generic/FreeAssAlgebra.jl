###############################################################################
#
#   FreeAssAlgebra.jl : free associative algebra R<x1,...,xn>
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

function parent_type(::Type{FreeAssAlgElem{T}}) where T <: RingElement
    return FreeAssAlgebra{T}
end

function elem_type(::Type{FreeAssAlgebra{T}}) where T <: RingElement
    return FreeAssAlgElem{T}
end

function parent(a::FreeAssAlgElem)
    return a.parent
end

base_ring_type(::Type{FreeAssAlgebra{T}}) where T <: RingElement = parent_type(T)

base_ring(a::FreeAssAlgebra{T}) where T <: RingElement = a.base_ring::parent_type(T)

function symbols(a::FreeAssAlgebra)
    return a.S
end

function number_of_variables(a::FreeAssAlgebra)
    return length(a.S)
end

function length(a::FreeAssAlgElem)
    return a.length
end

function check_parent(
    a::FreeAssAlgElem{T},
    b::FreeAssAlgElem{T},
    throw::Bool = true,
) where T <: RingElement
    b = parent(a) != parent(b)
    b & throw && error("Incompatible rings in operation")
    return !b
end

###############################################################################
#
# Basic manipulation
#
###############################################################################

function Base.deepcopy_internal(a::FreeAssAlgElem{T}, dict::IdDict) where T <: RingElement
    return FreeAssAlgElem{T}(
        a.parent,
        deepcopy_internal(a.coeffs, dict),
        deepcopy_internal(a.exps, dict),
        a.length,
    )
end

function zero(a::FreeAssAlgebra{T}) where T
    return FreeAssAlgElem{T}(a, T[], Vector{Int}[], 0)
end

function one(a::FreeAssAlgebra{T}) where T
    c = one(base_ring(a))
    !iszero(c) || return zero(a)
    return FreeAssAlgElem{T}(a, [c], [Int[]], 1)
end

function iszero(a::FreeAssAlgElem{T}) where T
    return length(a) == 0
end

function isone(a::FreeAssAlgElem{T}) where T
    if length(a) < 1
        return isone(zero(base_ring(a)))
    else
        return a.length == 1 && isone(a.coeffs[1]) && isempty(a.exps[1])
    end
end

function number_of_generators(a::FreeAssAlgebra{T}) where T
    return number_of_variables(a)
end

function gen(a::FreeAssAlgebra{T}, i::Int) where T
    @boundscheck 1 <= i <= ngens(a) || throw(ArgumentError("variable index out of range"))
    c = one(base_ring(a))
    iszero(c) && return zero(a)
    return FreeAssAlgElem{T}(a, T[c], [Int[i]], 1)
end

function gens(a::FreeAssAlgebra{T}) where T <: RingElement
    return [gen(a, i) for i in 1:ngens(a)]
end

function is_gen(a::FreeAssAlgElem{T}) where T
    if length(a) < 1
        return iszero(one(base_ring(a)))
    else
        return a.length == 1 && isone(a.coeffs[1]) && length(a.exps[1]) == 1
    end
end

function is_constant(a::FreeAssAlgElem{T}) where T
    return length(a) == 0 || (length(a) == 1 && isempty(a.exps[1]))
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{FreeAssAlgElem{T}}, ::Type{FreeAssAlgElem{T}}) where T <: RingElement =
    FreeAssAlgElem{T}

function promote_rule(
    ::Type{FreeAssAlgElem{T}},
    ::Type{U},
) where {T <: RingElement, U <: RingElement}
    promote_rule(T, U) == T ? FreeAssAlgElem{T} : Union{}
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::FreeAssAlgebra{T})() where T
    return zero(a)
end

function (a::FreeAssAlgebra{T})(b::T) where T
    iszero(b) && return zero(a)
    return FreeAssAlgElem{T}(a, T[b], [Int[]], 1)
end

function (a::FreeAssAlgebra{T})(b::Integer) where T
    iszero(b) && return zero(a)
    R = base_ring(a)
    return FreeAssAlgElem{T}(a, T[R(b)], [Int[]], 1)
end

function (a::FreeAssAlgebra{T})(b::FreeAssAlgElem{T}) where T <: RingElement
    parent(b) != a && error("Unable to coerce element")
    return b
end

function (a::FreeAssAlgebra{T})(c::Vector{T}, e::Vector{Vector{Int}}) where T
    for ei in e
        @boundscheck all(i -> (1 <= i <= nvars(a)), ei) ||
                     throw(ArgumentError("variable index out of range"))
    end
    n = length(c)
    n == length(e) ||
        error("coefficient array and exponent array should have the same length")
    z = FreeAssAlgElem{T}(a, copy(c), copy(e), n)
    return combine_like_terms!(sort_terms!(z))
end

###############################################################################
#
# Coefficients, Terms, Etc.
#
###############################################################################

function coeff(a::FreeAssAlgElem, i::Int)
    @boundscheck 1 <= i <= length(a) || throw(ArgumentError("index out of range"))
    return a.coeffs[i]
end

function term(a::FreeAssAlgElem{T}, i::Int) where T <: RingElement
    @boundscheck 1 <= i <= length(a) || throw(ArgumentError("index out of range"))
    R = parent(a)
    return FreeAssAlgElem{T}(R, [a.coeffs[i]], [a.exps[i]], 1)
end

function monomial(a::FreeAssAlgElem{T}, i::Int) where T <: RingElement
    @boundscheck 1 <= i <= length(a) || throw(ArgumentError("index out of range"))
    R = parent(a)
    return FreeAssAlgElem{T}(R, T[one(base_ring(R))], [a.exps[i]], 1)
end

@doc raw"""
    exponent_word(a::FreeAssAlgElem{T}, i::Int) where T <: RingElement

Return a vector of variable indices corresponding to the monomial of the
$i$-th term of $a$. Term numbering begins at $1$, and the variable
indices are given in the order of the variables for the ring.
"""
function exponent_word(a::FreeAssAlgElem{T}, i::Int) where T <: RingElement
    @boundscheck 1 <= i <= length(a) || throw(ArgumentError("index out of range"))
    return a.exps[i]
end

function Base.iterate(a::FreeAssAlgExponentWords, state = 0)
    state += 1
    state <= length(a.poly) || return nothing
    return exponent_word(a.poly, state), state
end

function leading_coefficient(a::FreeAssAlgElem{T}) where T
    return a.length > 0 ? coeff(a, 1) : zero(base_ring(a))
end

function leading_monomial(a::FreeAssAlgElem{T}) where T
    if length(a) < 1
        throw(ArgumentError("Zero polynomial does not have a leading monomial"))
    end
    return monomial(a, 1)
end

function leading_term(a::FreeAssAlgElem{T}) where T
    if length(a) < 1
        throw(ArgumentError("Zero polynomial does not have a leading term"))
    end
    return term(a, 1)
end

function leading_exponent_word(a::FreeAssAlgElem{T}) where T
    if length(a) < 1
        throw(ArgumentError("Zero polynomial does not have a leading exponent word"))
    end
    return exponent_word(a, 1)
end

function total_degree(a::FreeAssAlgElem{T}) where T
    # currently stored in dexlex
    return length(a) > 0 ? length(a.exps[1]) : -1
end

function Base.length(
    x::FreeAssAlgExponentWords{T},
) where {S <: RingElement, T <: FreeAssAlgElem{S}}
    return length(x.poly)
end

function Base.eltype(
    x::FreeAssAlgExponentWords{T},
) where {S <: RingElement, T <: FreeAssAlgElem{S}}
    return Vector{Int}
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

function canonical_unit(a::FreeAssAlgElem{T}) where T <: RingElement
    return canonical_unit(leading_coefficient(a))
end

###############################################################################
#
# Unsafe functions
#
###############################################################################

function fit!(a::FreeAssAlgElem{T}, n::Int) where T <: RingElement
    if length(a.coeffs) < n
        resize!(a.coeffs, n)
    end
    if length(a.exps) < n
        resize!(a.exps, n)
    end
    return nothing
end

for T in [RingElem, Integer, Rational, AbstractFloat]
    @eval begin
        function setcoeff!(a::FreeAssAlgElem{S}, i::Int, c::S) where S <: $T
            fit!(a, i)
            a.coeffs[i] = c
            if i > length(a)
                a.length = i
            end
            return a
        end
    end
end

function set_exponent_word!(
    a::FreeAssAlgElem{T},
    i::Int,
    w::Vector{Int},
) where T <: RingElement
    n = nvars(parent(a))
    @boundscheck all(x -> 1 <= x <= n, w) ||
                 throw(ArgumentError("variable index out of range"))
    fit!(a, i)
    a.exps[i] = w
    if i > length(a)
        a.length = i
    end
    return a
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::FreeAssAlgElem{T}, b::FreeAssAlgElem{T}) where T
    fl = check_parent(a, b, false)
    !fl && return false
    return a.length == b.length &&
           view(a.exps, 1:a.length) == view(b.exps, 1:b.length) &&
           view(a.coeffs, 1:a.length) == view(b.coeffs, 1:b.length)
end

function word_cmp(a::Vector{Int}, b::Vector{Int})
    if length(a) > length(b)
        return +1
    elseif length(a) < length(b)
        return -1
    else
        # deglex
        for i in 1:length(a)
            if a[i] > b[i]
                return -1
            elseif a[i] < b[i]
                return +1
            end
        end
        return 0
    end
end

function word_gt(a::Vector{Int}, b::Vector{Int})
    return word_cmp(a, b) > 0
end

function sort_terms!(z::FreeAssAlgElem{T}) where T
    n = length(z)
    if n > 1
        p = sortperm(view(z.exps, 1:n), lt = word_gt)
        z.coeffs = [z.coeffs[p[i]] for i in 1:n]
        z.exps = [z.exps[p[i]] for i in 1:n]
    end
    return z
end

function combine_like_terms!(z::FreeAssAlgElem{T}) where T
    o = 0
    i = 1
    while i <= z.length
        if o > 0 && word_cmp(z.exps[o], z.exps[i]) == 0
            z.coeffs[o] += z.coeffs[i]
        else
            o += (o < 1 || !iszero(z.coeffs[o]))
            z.exps[o] = z.exps[i]
            z.coeffs[o] = z.coeffs[i]
        end
        i += 1
    end
    o += (o < 1 || !iszero(z.coeffs[o]))
    z.length = o - 1
    return z
end


@doc """
    isless(p::FreeAssAlgElem{T}, q::FreeAssAlgElem{T}) where T

Implements the degree lexicographic ordering on terms, i.e.
first, the degrees of the largest monomials are compared, and if they
are the same, they are compared lexicographically and if they are still the same,
the coefficients are compared. 
If everything is still the same, the next largest monomial is compared
and lastly the number of monomials is compared.
Since the coefficients are also compared, this only works when the 
coefficient Ring implements isless.

The order of letters is the reverse of the order given when initialising the algebra.

# Examples
```jldoctest; setup = :(using AbstractAlgebra)
julia> R, (x, y) = free_associative_algebra(QQ, ["x", "y"]);

julia> x < y^2
true

julia> x^2 < x^2 + y
true

julia> y < x
true

julia> x^2 < 2*x^2
true
```
"""
function isless(p::FreeAssAlgElem{T}, q::FreeAssAlgElem{T}) where T
    if p == q
        return false
    end
    l = min(length(p.exps), length(q.exps))
    sort_terms!(p)
    sort_terms!(q)
    for i in 1:l
        c = word_cmp(q.exps[i], p.exps[i])
        if c > 0
            return true
        elseif c < 0
            return false
        elseif p.coeffs[i] != q.coeffs[i]
            return p.coeffs[i] < q.coeffs[i]
        end
    end
    if length(p.exps) < length(q.exps)
        return true
    else
        return false
    end
end

###############################################################################
#
#   Arithmetic
#
###############################################################################

function -(a::FreeAssAlgElem{T}) where T <: RingElement
    n = length(a)
    R = parent(a)
    zcoeffs = T[-a.coeffs[i] for i in 1:n]
    return FreeAssAlgElem{T}(R, zcoeffs, copy(a.exps), n)
end

function *(a::FreeAssAlgElem{T}, b::FreeAssAlgElem{T}) where T <: RingElement
    zcoeffs = T[]
    zexps = Vector{Int}[]
    for i in 1:a.length, j in 1:b.length
        push!(zcoeffs, a.coeffs[i] * b.coeffs[j])
        push!(zexps, vcat(a.exps[i], b.exps[j]))
    end
    z = FreeAssAlgElem{T}(parent(a), zcoeffs, zexps, length(zcoeffs))
    return combine_like_terms!(sort_terms!(z))
end

function +(a::FreeAssAlgElem{T}, b::FreeAssAlgElem{T}) where T <: RingElement
    zcoeffs = T[]
    zexps = Vector{Int}[]
    i = j = 1
    while i <= a.length && j <= b.length
        c = word_cmp(a.exps[i], b.exps[j])
        if c < 0
            push!(zcoeffs, b.coeffs[j])
            push!(zexps, b.exps[j])
            j += 1
        elseif c > 0
            push!(zcoeffs, a.coeffs[i])
            push!(zexps, a.exps[i])
            i += 1
        else
            s = a.coeffs[i] + b.coeffs[j]
            if !iszero(s)
                push!(zcoeffs, s)
                push!(zexps, a.exps[i])
            end
            i += 1
            j += 1
        end
    end
    while i <= a.length
        push!(zcoeffs, a.coeffs[i])
        push!(zexps, a.exps[i])
        i += 1
    end
    while j <= b.length
        push!(zcoeffs, b.coeffs[j])
        push!(zexps, b.exps[j])
        j += 1
    end
    return FreeAssAlgElem{T}(parent(a), zcoeffs, zexps, length(zcoeffs))
end

# a - b ignoring the first "start" terms of both
function _sub_rest(
    a::FreeAssAlgElem{T},
    b::FreeAssAlgElem{T},
    start::Int,
) where T <: RingElement
    zcoeffs = T[]
    zexps = Vector{Int}[]
    i = j = start + 1
    while i <= a.length && j <= b.length
        c = word_cmp(a.exps[i], b.exps[j])
        if c < 0
            push!(zcoeffs, -b.coeffs[j])
            push!(zexps, b.exps[j])
            j += 1
        elseif c > 0
            push!(zcoeffs, a.coeffs[i])
            push!(zexps, a.exps[i])
            i += 1
        else
            s = a.coeffs[i] - b.coeffs[j]
            if !iszero(s)
                push!(zcoeffs, s)
                push!(zexps, a.exps[i])
            end
            i += 1
            j += 1
        end
    end
    while i <= a.length
        push!(zcoeffs, a.coeffs[i])
        push!(zexps, a.exps[i])
        i += 1
    end
    while j <= b.length
        push!(zcoeffs, -b.coeffs[j])
        push!(zexps, b.exps[j])
        j += 1
    end
    return FreeAssAlgElem{T}(parent(a), zcoeffs, zexps, length(zcoeffs))
end

function -(a::FreeAssAlgElem{T}, b::FreeAssAlgElem{T}) where T <: RingElement
    return _sub_rest(a, b, 0)
end

function ^(a::FreeAssAlgElem{T}, b::Integer) where T <: RingElement
    if b == 0
        return one(parent(a))
    elseif b == 1
        return deepcopy(a)
    elseif a.length == 1
        if isempty(a.exps[1])
            e = [Int[]]
        else
            b < 0 && throw(NotInvertibleError(a))
            e = Vector{Int}[reduce(vcat, [a.exps[1] for i in 1:b])]
        end
        return FreeAssAlgElem{T}(parent(a), [a.coeffs[1]^b], e, 1)
    else
        b < 0 && throw(NotInvertibleError(a))
        return AbstractAlgebra.internal_power(a, b)
    end
end

###############################################################################
#
# Division
#
###############################################################################

# return c*w*a*wp
function mul_term(c::T, w::Vector{Int}, a::FreeAssAlgElem{T}, wp::Vector{Int}) where T
    zcoeffs =
        isone(c) ? T[a.coeffs[i] for i in 1:a.length] :
        T[c * a.coeffs[i] for i in 1:a.length]
    zexps = Vector{Int}[vcat(w, a.exps[i], wp) for i in 1:a.length]
    return FreeAssAlgElem{T}(parent(a), zcoeffs, zexps, a.length)
end

# return (true, l, r) with a = l*b*r and length(l) minimal
#     or (false, junk, junk) if a is not two-sided divisible by b
function word_divides_leftmost(a::Vector{Int}, b::Vector{Int})
    n = length(b)
    for i in 0:length(a)-n
        match = true
        for j in 1:n
            if b[j] != a[i+j]
                match = false
                break
            end
        end
        if match
            return (true, Int[a[k] for k in 1:i], Int[a[k] for k in 1+i+n:length(a)])
        end
    end
    return (false, Int[], Int[])
end

# return (true, l, r) with a = l*b*r and length(r) minimal
#     or (false, junk, junk) if a is not two-sided divisible by b
function word_divides_rightmost(a::Vector{Int}, b::Vector{Int})
    n = length(b)
    for i in length(a)-n:-1:0
        match = true
        for j in 1:n
            if b[j] != a[i+j]
                match = false
                break
            end
        end
        if match
            return (true, Int[a[k] for k in 1:i], Int[a[k] for k in 1+i+n:length(a)])
        end
    end
    return (false, Int[], Int[])
end


function AbstractAlgebra.divexact_left(
    f::FreeAssAlgElem{T},
    g::FreeAssAlgElem{T};
    check::Bool = true,
) where T
    R = parent(f)
    qcoeffs = T[]
    qexps = Vector{Int}[]
    while length(f) > 0
        ok, ml, mr = word_divides_leftmost(f.exps[1], g.exps[1])
        ok && isempty(ml) || throw(ArgumentError("Not an exact division"))
        qi = divexact(f.coeffs[1], g.coeffs[1])
        push!(qcoeffs, qi)
        push!(qexps, mr)
        f = _sub_rest(f, mul_term(qi, ml, g, mr), 1) # enforce lt cancellation
    end
    return FreeAssAlgElem{T}(R, qcoeffs, qexps, length(qcoeffs))
end

function AbstractAlgebra.divexact_right(
    f::FreeAssAlgElem{T},
    g::FreeAssAlgElem{T};
    check::Bool = true,
) where T
    R = parent(f)
    qcoeffs = T[]
    qexps = Vector{Int}[]
    while length(f) > 0
        ok, ml, mr = word_divides_rightmost(f.exps[1], g.exps[1])
        ok && isempty(mr) || throw(ArgumentError("Not an exact division"))
        qi = divexact(f.coeffs[1], g.coeffs[1])
        push!(qcoeffs, qi)
        push!(qexps, ml)
        f = _sub_rest(f, mul_term(qi, ml, g, mr), 1) # enforce lt cancellation
    end
    return FreeAssAlgElem{T}(R, qcoeffs, qexps, length(qcoeffs))
end


###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

function divexact(
    a::FreeAssAlgElem{T},
    b::Integer;
    check::Bool = true,
) where T <: RingElement

    n = length(a)
    R = parent(a)
    b = base_ring(R)(b)
    zcoeffs = T[divexact(a.coeffs[i], b, check = check) for i in 1:n]
    return combine_like_terms!(FreeAssAlgElem{T}(R, zcoeffs, copy(a.exps), n))
end


################################################################################
#
#  Change base ring
#
################################################################################

function _change_freeassalg_ring(R, Rx, cached)
    P, _ = AbstractAlgebra.free_associative_algebra(R, symbols(Rx); cached = cached)
    return P
end

function change_base_ring(
    R::Ring,
    a::FreeAssAlgElem{T};
    cached::Bool = true,
    parent::AbstractAlgebra.FreeAssAlgebra = _change_freeassalg_ring(R, parent(a), cached),
) where T <: RingElement
    base_ring(parent) != R && error("Base rings do not match.")
    return _map(R, a, parent)
end

function map_coefficients(
    f::S,
    a::FreeAssAlgElem{T};
    cached::Bool = true,
    parent::AbstractAlgebra.FreeAssAlgebra = _change_freeassalg_ring(
        parent(f(zero(base_ring(a)))),
        parent(a),
        cached,
    ),
) where {S, T <: RingElement}
    return _map(f, a, parent)
end

function _map(g::S, a::FreeAssAlgElem{T}, Rx) where {S, T <: RingElement}
    cvzip = zip(coefficients(a), exponent_words(a))
    M = MPolyBuildCtx(Rx)
    for (c, v) in cvzip
        push_term!(M, g(c), v)
    end

    return finish(M)
end


###############################################################################
#
#   free_associative_algebra constructor
#
###############################################################################

function free_associative_algebra(
    R::AbstractAlgebra.Ring,
    s::Vector{Symbol};
    cached::Bool = true,
)
    parent_obj = FreeAssAlgebra{elem_type(R)}(R, s, cached)
    return (parent_obj, gens(parent_obj))
end
