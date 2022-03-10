###############################################################################
#
#   LaurentMPoly.jl: Multivariate Laurent polynomials
#
###############################################################################

###############################################################################
#
#   String I/O
#
###############################################################################

function expressify(a::LaurentMPolyElem, x = symbols(parent(a)); context = nothing)
    sum = Expr(:call, :+)
    n = nvars(parent(a))
    for (c, v) in zip(coefficients(a), exponent_vectors(a))
        prod = Expr(:call, :*)
        if !isone(c)
            push!(prod.args, expressify(c, context = context))
        end
        for i in 1:n
            if isone(v[i])
                push!(prod.args, x[i])
            elseif !iszero(v[i])
                push!(prod.args, Expr(:call, :^, x[i], v[i]))
            end
        end
        push!(sum.args, prod)
    end
    return sum
end

@enable_all_show_via_expressify LaurentMPolyElem

function expressify(a::LaurentMPolyRing; context = nothing)
    return Expr(:sequence, Expr(:text, "Multivariate Laurent Polynomial Ring over "),
                           expressify(base_ring(a); context = context),
                           Expr(:text, " in "),
                           Expr(:series, symbols(a)...))
end

@enable_all_show_via_expressify LaurentMPolyRing

###############################################################################
#
#   LaurentPolynomialRing constructor
#
###############################################################################

@doc Markdown.doc"""
    LaurentPolynomialRing(R::AbstractAlgebra.Ring, s::Vector{T}; cached::Bool = true) where T <: Union{String, Char, Symbol}

Given a base ring `R` and an array of strings `s` specifying how the
generators (variables) should be printed, return a tuple `T, (x1, x2, ...)`
representing the new ring $T = R[x1, 1/x1, x2, 1/x2, ...]$ and the generators
$x1, x2, ...$ of the  ring. By default the parent object `T` will
depend only on `R` and `x1, x2, ...` and will be cached. Setting the optional
argument `cached` to `false` will prevent the parent object `T` from being
cached.
"""
function LaurentPolynomialRing(R::AbstractAlgebra.Ring, s::Vector{String}; cached::Bool = true)
   return Generic.LaurentPolynomialRing(R, [Symbol(v) for v in s], cached=cached)
end

function LaurentPolynomialRing(R::AbstractAlgebra.Ring, s::Vector{Char}; cached::Bool = true)
   return Generic.LaurentPolynomialRing(R, [Symbol(v) for v in s], cached=cached)
end

function LaurentPolynomialRing(R::AbstractAlgebra.Ring, s::Vector{Symbol}; cached::Bool = true)
   return Generic.LaurentPolynomialRing(R, s; cached=cached)
end

function LaurentPolynomialRing(R::AbstractAlgebra.Ring, n::Int, s::String; cached::Bool = false)
   return Generic.LaurentPolynomialRing(R, [Symbol(s, i) for i=1:n]; cached=cached)
end

function LaurentPolynomialRing(R::AbstractAlgebra.Ring, n::Int, s::Char; cached::Bool = false)
   return Generic.LaurentPolynomialRing(R, [Symbol(s, i) for i=1:n]; cached=cached)
end

function LaurentPolynomialRing(R::AbstractAlgebra.Ring, n::Int, s::Symbol=:x; cached::Bool = false)
   return Generic.LaurentPolynomialRing(R, [Symbol(s, i) for i=1:n], cached = cached)
end

