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

function expressify(a::LaurentMPolyRingElem, x = symbols(parent(a)); context = nothing)
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

@enable_all_show_via_expressify LaurentMPolyRingElem

function expressify(a::LaurentMPolyRing; context = nothing)
    return Expr(:sequence, Expr(:text, "Multivariate Laurent Polynomial Ring in "),
                           Expr(:series, symbols(a)...),
                           Expr(:text, " over "),
                           expressify(base_ring(a); context = context))
end

@enable_all_show_via_expressify LaurentMPolyRing

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function gens(R::LaurentMPolyRing)
    return [gen(R, i) for i in 1:nvars(R)]
end

coefficient_ring(a::LaurentMPolyRingElem) = coefficient_ring(parent(a))
base_ring(a::LaurentMPolyRingElem) = base_ring(parent(a))

###############################################################################
#
#   Derivative
#
###############################################################################

function derivative(a::LaurentMPolyRingElem{T}, x::LaurentMPolyRingElem{T}) where T <: RingElement
   check_parent(a, x)
   return derivative(a, var_index(x))
end

###############################################################################
#
#   Random elements
#
###############################################################################

function RandomExtensions.maketype(S::AbstractAlgebra.LaurentMPolyRing, _, _, _)
    return elem_type(S)
end

function RandomExtensions.make(S::AbstractAlgebra.LaurentMPolyRing,
                               term_range::UnitRange{Int},
                               exp_bound::UnitRange{Int},
                               vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, term_range, exp_bound, vs[1])
   else
      make(S, term_range, exp_bound, make(R, vs...))
   end
end

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make4{<:RingElement,
                                         <:AbstractAlgebra.LaurentMPolyRing,
                                         UnitRange{Int},
                                         UnitRange{Int}}})
   S, term_range, exp_bound, v = sp[][1:end]
   f = zero(S)
   g = gens(S)
   R = base_ring(S)
   for i = 1:rand(rng, term_range)
      term = one(S)
      for j = 1:length(g)
         term *= g[j]^rand(rng, exp_bound)
      end
      term *= rand(rng, v)
      f += term
   end
   return f
end

function rand(rng::AbstractRNG, S::AbstractAlgebra.LaurentMPolyRing,
              term_range::UnitRange{Int}, exp_bound::UnitRange{Int}, v...)
   rand(rng, make(S, term_range, exp_bound, v...))
end

function rand(S::AbstractAlgebra.LaurentMPolyRing, term_range, exp_bound, v...)
   rand(GLOBAL_RNG, S, term_range, exp_bound, v...)
end

###############################################################################
#
#   LaurentPolynomialRing constructor
#
###############################################################################

@doc raw"""
    LaurentPolynomialRing(R::AbstractAlgebra.Ring, s::Vector{T}; cached::Bool = true) where T <: VarName

Given a base ring `R` and an array of strings `s` specifying how the
generators (variables) should be printed, return a tuple `T, (x1, x2, ...)`
representing the new ring $T = R[x1, 1/x1, x2, 1/x2, ...]$ and the generators
$x1, x2, ...$ of the  ring. By default the parent object `T` will
depend only on `R` and `x1, x2, ...` and will be cached. Setting the optional
argument `cached` to `false` will prevent the parent object `T` from being
cached.
"""
function LaurentPolynomialRing(R::AbstractAlgebra.Ring, s::AbstractVector{<:VarName}; cached::Bool = true)
   return Generic.LaurentPolynomialRing(R, [Symbol(v) for v in s], cached=cached)
end

function LaurentPolynomialRing(R::AbstractAlgebra.Ring, s::Vector{Symbol}; cached::Bool = true)
   return Generic.LaurentPolynomialRing(R, s; cached=cached)
end

function LaurentPolynomialRing(R::AbstractAlgebra.Ring, n::Int, s::VarName=:x; cached::Bool = false)
   return Generic.LaurentPolynomialRing(R, [Symbol(s, i) for i=1:n]; cached=cached)
end
