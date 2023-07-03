ngens(R::MPolyRing) = nvars(R)

AbstractAlgebra.promote_rule(::Type{T}, ::Type{S}) where {S<:NumFieldElem,T<:Integer} = S

AbstractAlgebra.promote_rule(::Type{S}, ::Type{T}) where {S<:NumFieldElem,T<:Integer} = S

elem_type(::Type{Generic.ResidueRing{T}}) where {T} = Generic.ResidueRingElem{T}

nrows(A::Matrix{T}) where {T} = size(A)[1]
ncols(A::Matrix{T}) where {T} = size(A)[2]

function polynomial_ring(R::Ring; cached::Bool=false)
    return polynomial_ring(R, "x", cached=cached)
end

function content(a::PolyRingElem{<:FieldElem})
    return one(base_ring(a))
end

function Base.rem(f::PolyRingElem, g::PolyRingElem)
    return mod(f, g)
end

function mod(f::AbstractAlgebra.PolyRingElem{T}, g::AbstractAlgebra.PolyRingElem{T}) where {T<:RingElem}
    check_parent(f, g)
    if length(g) == 0
        throw(DivideError())
    end
    if length(f) >= length(g)
        f = deepcopy(f)
        b = leading_coefficient(g)
        g = inv(b) * g
        c = base_ring(f)()
        while length(f) >= length(g)
            l = -leading_coefficient(f)
            for i = 1:length(g)
                c = mul!(c, coeff(g, i - 1), l)
                u = coeff(f, i + length(f) - length(g) - 1)
                u = addeq!(u, c)
                f = setcoeff!(f, i + length(f) - length(g) - 1, u)
            end
            set_length!(f, normalise(f, length(f) - 1))
        end
    end
    return f
end

function matrix(a::Vector{Vector{T}}) where {T}
    return matrix(permutedims(reduce(hcat, a), (2, 1)))
end

export is_constant

function is_constant(f::PolyRingElem)
    return f.length < 2
end

function identity_matrix(::Type{MatElem}, R::Ring, n::Int)
    return identity_matrix(R, n)
end

zero_matrix(::Type{Int}, r, c) = zeros(Int, r, c)

base_ring(::Vector{Int}) = Int

function AbstractAlgebra.is_symmetric(M::MatElem)
    for i in 1:nrows(M)
        for j in i:ncols(M)
            if M[i, j] != M[j, i]
                return false
            end
        end
    end
    return true
end

################################################################################
#
#  Create a matrix from rows
#
################################################################################

function matrix(K::Ring, R::Vector{<:Vector})
    if length(R) == 0
        return zero_matrix(K, 0, 0)
    else
        n = length(R)
        m = length(R[1])
        z = zero_matrix(K, n, m)
        for i in 1:n
            @assert length(R[i]) == m
            for j in 1:m
                z[i, j] = R[i][j]
            end
        end
        return z
    end
end

export neg!

sub!(z::Rational{Int}, x::Rational{Int}, y::Int) = x - y

neg!(z::Rational{Int}, x::Rational{Int}) = -x

add!(z::Rational{Int}, x::Rational{Int}, y::Int) = x + y

mul!(z::Rational{Int}, x::Rational{Int}, y::Int) = x * y

is_negative(x::Rational) = x.num < 0

function is_upper_triangular(A::Generic.Mat)
    m = nrows(A)
    n = ncols(A)
    d = 0
    for r = 1:m
        for c = 1:n
            if !iszero(A[r, c])
                if c <= d
                    return false
                end
                d = c
                break
            end
        end
    end
    return true
end

function sub(M::Generic.Mat, rows::UnitRange{Int}, cols::UnitRange{Int})
    @assert step(rows) == 1 && step(cols) == 1
    z = zero_matrix(base_ring(M), length(rows), length(cols))
    for i in rows
        for j in cols
            z[i-first(rows)+1, j-first(cols)+1] = M[i, j]
        end
    end
    return z
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

function mulmod(a::S, b::S, mod::Vector{S}) where {S<:MPolyRingElem{T}} where {T<:RingElem}
    return Base.divrem(a * b, mod)[2]
end

@inline ngens(R::AbstractAlgebra.Generic.MPolyRing) = R.num_vars

#trivia to make life easier

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

ngens(L::SimpleNumField{T}) where {T} = 1

is_unit(a::NumFieldElem) = !iszero(a)

canonical_unit(a::NumFieldElem) = a

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

function addmul!(z::T, x::T, y::T) where {T<:RingElement}
    zz = parent(z)()
    zz = mul!(zz, x, y)
    return addeq!(z, zz)
end

function gcd(a::ResElem{T}, b::ResElem{T}) where {T<:Integer}
    m = modulus(a)
    return parent(a)(gcd(gcd(a.data, m), b.data))
end

function Base.minimum(::typeof(precision), a::Vector{<:SeriesElem})
    return minimum(map(precision, a))
end

function Base.maximum(::typeof(precision), a::Vector{<:SeriesElem})
    return maximum(map(precision, a))
end

function canonical_unit(a::SeriesElem)
    iszero(a) && return one(parent(a))
    v = valuation(a)
    v == 0 && return a
    v > 0 && return shift_right(a, v)
    return shift_left(a, -v)
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

# should be Nemo/AA
# TODO: symbols vs strings
#       lift(PolyRing, Series)
#       lift(FracField, Series)
#       (to be in line with lift(ZZ, padic) and lift(QQ, padic)
#TODO: some of this would only work for Abs, not Rel, however, this should be fine here
function map_coefficients(f, a::RelPowerSeriesRingElem; parent::SeriesRing)
    c = typeof(f(coeff(a, 0)))[]
    for i = 0:pol_length(a)-1
        push!(c, f(polcoeff(a, i)))
    end
    b = parent(c, length(c), precision(a), valuation(a))
    return b
end

function lift(R::PolyRing{S}, s::SeriesElem{S}) where {S}
    t = R()
    for x = 0:pol_length(s)
        setcoeff!(t, x, polcoeff(s, x))
    end
    return shift_left(t, valuation(s))
end

function gen(R::Union{Generic.ResidueRing{T},Generic.ResidueField{T}}) where {T<:PolyRingElem}
    return R(gen(base_ring(R)))
end

function characteristic(R::Union{Generic.ResidueRing{T},Generic.ResidueField{T}}) where {T<:PolyRingElem}
    return characteristic(base_ring(base_ring(R)))
end

function size(R::Union{Generic.ResidueRing{T},Generic.ResidueField{T}}) where {T<:ResElem}
    return size(base_ring(base_ring(R)))^degree(modulus(R))
end

function size(R::Union{Generic.ResidueRing{T},Generic.ResidueField{T}}) where {T<:PolyRingElem}
    return size(base_ring(base_ring(R)))^degree(R.modulus)
end

function rand(R::Union{Generic.ResidueRing{T},Generic.ResidueField{T}}) where {T<:PolyRingElem}
    r = rand(base_ring(base_ring(R)))
    g = gen(R)
    for i = 1:degree(R.modulus)
        r = r * g + rand(base_ring(base_ring(R)))
    end
    return r
end

function gens(R::Union{Generic.ResidueRing{T},Generic.ResidueField{T}}) where {T<:PolyRingElem} ## probably needs more cases
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

function lift(a::Generic.ResidueRingElem)
    return a.data
end

function lift(a::Generic.ResidueFieldElem)
    return a.data
end

function leading_monomial(f::Generic.MPoly)
    R = parent(f)
    l = length(f)
    if l == 0
        return f
    end
    A = f.exps
    r, c = size(A)
    e = A[1:r, 1:1]
    return R([one(base_ring(R))], e)
end

function leading_coefficient(f::Generic.MPoly)
    return f.coeffs[1]
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

function show(io::IO, M::Map)
    @show_name(io, M)
    if get(io, :compact, false)
        print(io, domain(M), " --> ", codomain(M), "\n")
        return
    end
    io = Base.IOContext(io, :compact => true)
    print(io, "Map with following data\n")
    print(io, "Domain:\n")
    print(io, "=======\n")
    print(io, domain(M))
    print(io, "\nCodomain:\n")
    print(io, "=========\n")
    print(io, codomain(M))
end

function image(M::Map{D,C}, a) where {D,C}
    if isdefined(M, :header)
        if isdefined(M.header, :image)
            return M.header.image(a)::elem_type(C)
        else
            error("No image function known")
        end
    else
        return M(a)
    end
end

function preimage(M::Map{D,C}, a) where {D,C}
    if isdefined(M.header, :preimage)
        p = M.header.preimage(a)::elem_type(D)
        @assert parent(p) === domain(M)
        return p
    end
    error("No preimage function known")
end

\(f::Map, x) = preimage(f, x)

function preimage(f::AbstractAlgebra.Generic.CompositeMap, a)
    return preimage(f.map1, preimage(f.map2, a))
end

dense_matrix_type(::Type{T}) where {T} = Generic.MatSpaceElem{T}

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

transpose!(A::MatrixElem) = transpose(A)

function zero_matrix(::Type{MatElem}, R::Ring, n::Int)
    return zero_matrix(R, n)
end

function zero_matrix(::Type{MatElem}, R::Ring, n::Int, m::Int)
    return zero_matrix(R, n, m)
end

function matrix(A::Matrix{T}) where {T<:RingElem}
    r, c = size(A)
    (r < 0 || c < 0) && error("Array must be non-empty")
    m = matrix(parent(A[1, 1]), A)
    return m
end

function matrix(A::Vector{T}) where {T<:RingElem}
    return matrix(reshape(A, length(A), 1))
end

function is_zero_row(M::MatElem{T}, i::Int) where {T}
    for j in 1:ncols(M)
        if !iszero(M[i, j])
            return false
        end
    end
    return true
end

function is_zero_row(M::Matrix{T}, i::Int) where {T<:Integer}
    for j = 1:Base.size(M, 2)
        if M[i, j] != 0
            return false
        end
    end
    return true
end

function is_zero_row(M::Matrix{T}, i::Int) where {T<:RingElem}
    for j in 1:Base.size(M, 2)
        if !iszero(M[i, j])
            return false
        end
    end
    return true
end

################################################################################
#
#  Kernel function
#
################################################################################

@doc raw"""
    kernel(a::MatElem{T}; side::Symbol = :right) -> Int, MatElem{T}

It returns a tuple $(n, M)$, where $n$ is the rank of the kernel and $M$ is a basis for it. If side is $:right$ or not
specified, the right kernel is computed. If side is $:left$, the left kernel is computed.
"""
function kernel(A::MatElem; side::Symbol=:right)
    if side == :right
        return right_kernel(A)
    elseif side == :left
        return left_kernel(A)
    else
        error("Unsupported argument: :$side for side: Must be :left or :right")
    end
end

right_kernel(M::MatElem) = nullspace(M)

function left_kernel(M::MatElem)
    rk, M1 = nullspace(transpose(M))
    return rk, transpose(M1)
end

################################################################################
#
#  Kernel over different rings
#
################################################################################

@doc raw"""
kernel(a::MatrixElem{T}, R::Ring; side::Symbol = :right) -> n, MatElem{elem_type(R)}

It returns a tuple $(n, M)$, where $n$ is the rank of the kernel over $R$ and $M$ is a basis for it. If side is $:right$ or not
specified, the right kernel is computed. If side is $:left$, the left kernel is computed.
"""
function kernel(M::MatrixElem, R::Ring; side::Symbol=:right)
    MP = change_base_ring(R, M)
    return kernel(MP, side=side)
end

################################################################################
#
#  Diagonal (block) matrix creation
#
################################################################################

@doc raw"""
    diagonal_matrix(x::T...) where T <: RingElem -> MatElem{T}
    diagonal_matrix(x::Vector{T}) where T <: RingElem -> MatElem{T}
    diagonal_matrix(Q, x::Vector{T}) where T <: RingElem -> MatElem{T}

Returns a diagonal matrix whose diagonal entries are the elements of $x$.

# Examples

```jldoctest
julia> diagonal_matrix(QQ(1), QQ(2))
[1   0]
[0   2]

julia> diagonal_matrix([QQ(3), QQ(4)])
[3   0]
[0   4]

julia> diagonal_matrix(QQ, [5, 6])
[5   0]
[0   6]
```
"""
function diagonal_matrix(R::Ring, x::Vector{<:RingElement})
    x = R.(x)
    M = zero_matrix(R, length(x), length(x))
    for i = 1:length(x)
        M[i, i] = x[i]
    end
    return M
end

function diagonal_matrix(x::T, xs::T...) where {T<:RingElem}
    return diagonal_matrix(collect((x, xs...)))
end

diagonal_matrix(x::Vector{<:RingElement}) = diagonal_matrix(parent(x[1]), x)

@doc raw"""
    diagonal_matrix(x::Vector{T}) where T <: MatElem -> MatElem

Returns a block diagonal matrix whose diagonal blocks are the matrices in $x$.
"""
function diagonal_matrix(x::Vector{T}) where {T<:MatElem}
    return cat(x..., dims=(1, 2))::T
end

function diagonal_matrix(x::T, xs::T...) where {T<:MatElem}
    return cat(x, xs..., dims=(1, 2))::T
end

function diagonal_matrix(R::Ring, x::Vector{<:MatElem})
    if length(x) == 0
        return zero_matrix(R, 0, 0)
    end
    x = [change_base_ring(R, i) for i in x]
    return diagonal_matrix(x)
end

@doc raw"""
vcat(A::Vector{Generic.Mat}) -> Generic.Mat
vcat(A::Vector{ZZMatrix}) -> ZZMatrix

Forms a big matrix by vertically concatenating the matrices in $A$.
All component matrices need to have the same number of columns.
"""
function Base.vcat(A::Vector{T}) where {S<:RingElem,T<:MatElem{S}}
    if any(x -> ncols(x) != ncols(A[1]), A)
        error("Matrices must have same number of columns")
    end
    M = zero_matrix(base_ring(A[1]), sum(nrows, A), ncols(A[1]))
    s = 0
    for i = A
        for j = 1:nrows(i)
            for k = 1:ncols(i)
                M[s+j, k] = i[j, k]
            end
        end
        s += nrows(i)
    end
    return M
end

function Base.vcat(A::MatElem...)
    r = nrows(A[1])
    c = ncols(A[1])
    R = base_ring(A[1])
    for i = 2:length(A)
        @assert ncols(A[i]) == c
        @assert base_ring(A[i]) == R
        r += nrows(A[i])
    end
    X = zero_matrix(R, r, c)
    o = 1
    for i = 1:length(A)
        for j = 1:nrows(A[i])
            X[o, :] = A[i][j, :]
            o += 1
        end
    end
    return X
end

function Base.hcat(A::Vector{T}) where {S<:RingElem,T<:MatElem{S}}
    if any(x -> nrows(x) != nrows(A[1]), A)
        error("Matrices must have same number of rows")
    end
    M = zero_matrix(base_ring(A[1]), nrows(A[1]), sum(ncols, A))
    s = 0
    for i = A
        for j = 1:ncols(i)
            for k = 1:nrows(i)
                M[k, s+j] = i[k, j]
            end
        end
        s += ncols(i)
    end
    return M
end

function Base.hcat(A::MatElem...)
    r = nrows(A[1])
    c = ncols(A[1])
    R = base_ring(A[1])
    for i = 2:length(A)
        @assert nrows(A[i]) == r
        @assert base_ring(A[i]) == R
        c += ncols(A[i])
    end
    X = zero_matrix(R, r, c)
    o = 1
    for i = 1:length(A)
        for j = 1:ncols(A[i])
            X[:, o] = A[i][:, j]
            o += 1
        end
    end
    return X
end

function Base.cat(A::MatElem...; dims)
    @assert dims == (1, 2) || isa(dims, Int)

    if isa(dims, Int)
        if dims == 1
            return hcat(A...)
        elseif dims == 2
            return vcat(A...)
        else
            error("dims must be 1, 2, or (1,2)")
        end
    end

    local X
    for i = 1:length(A)
        if i == 1
            X = hcat(A[1], zero_matrix(base_ring(A[1]), nrows(A[1]), sum(Int[ncols(A[j]) for j = 2:length(A)])))
        else
            X = vcat(X, hcat(zero_matrix(base_ring(A[1]), nrows(A[i]), sum(ncols(A[j]) for j = 1:i-1)), A[i], zero_matrix(base_ring(A[1]), nrows(A[i]), sum(Int[ncols(A[j]) for j = i+1:length(A)]))))
        end
    end
    return X
end

#= seems to be in AA now
function Base.hvcat(rows::Tuple{Vararg{Int}}, A::MatElem...)
  B = hcat([A[i] for i=1:rows[1]]...)
  o = rows[1]
  for j=2:length(rows)
    C = hcat([A[i+o] for i=1:rows[j]]...)
    o += rows[j]
    B = vcat(B, C)
  end
  return B
end
=#

function is_upper_triangular(M::MatElem)
    n = nrows(M)
    for i = 2:n
        for j = 1:min(i - 1, ncols(M))
            if !iszero(M[i, j])
                return false
            end
        end
    end
    return true
end

function is_lower_triangular(M::MatElem)
    for i = 1:nrows(M)
        for j = i+1:ncols(M)
            if !iszero(M[i, j])
                return false
            end
        end
    end
    return true
end

################################################################################
#
#  Minpoly and Charpoly
#
################################################################################

function minpoly(M::MatElem)
    k = base_ring(M)
    kx, x = polynomial_ring(k, cached=false)
    return minpoly(kx, M)
end

function charpoly(M::MatElem)
    k = base_ring(M)
    kx, x = polynomial_ring(k, cached=false)
    return charpoly(kx, M)
end

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

function sub(M::MatElem{T}, r::UnitRange{<:Integer}, c::UnitRange{<:Integer}) where {T}
    z = similar(M, length(r), length(c))
    for i in 1:length(r)
        for j in 1:length(c)
            z[i, j] = M[r[i], c[j]]
        end
    end
    return z
end

sub!(z::T, x::T, y::T) where {T} = x - y

# TODO (CF):
# should be Bernstein'ed: this is slow for large valuations
# returns the maximal v s.th. z mod p^v == 0 and z div p^v
#   also useful if p is not prime....
#
# TODO: what happens to z = 0???

# function remove(z::T, p::T) where {T<:Integer}
#     z == 0 && return (0, z)
#     v = 0
#     @assert p > 1
#     while mod(z, p) == 0
#         z = Base.div(z, p)
#         v += 1
#     end
#     return (v, z)
# end

function remove(z::Rational{T}, p::T) where {T<:Integer}
    z == 0 && return (0, z)
    v, d = remove(denominator(z), p)
    w, n = remove(numerator(z), p)
    return w - v, n // d
end

# function valuation(z::T, p::T) where {T<:Integer}
#     iszero(z) && error("Not yet implemented")
#     v = 0
#     @assert p > 1
#     while mod(z, p) == 0
#         z = Base.div(z, p)
#         v += 1
#     end
#     return v
# end

function valuation(z::Rational{T}, p::T) where {T<:Integer}
    z == 0 && error("Not yet implemented")
    v = valuation(denominator(z), p)
    w = valuation(numerator(z), p)
    return w - v
end

is_negative(n::Integer) = cmp(n, 0) < 0
is_positive(n::Integer) = cmp(n, 0) > 0

function divisible(x::Integer, y::Integer)
    return iszero(rem(x, y))
end
