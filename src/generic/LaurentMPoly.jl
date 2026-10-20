###############################################################################
#
#   LaurentMPoly.jl : Generic multivariate Laurent polynomials over rings
#
###############################################################################

# The rest of the implementation is in src/LaurentMPoly.jl. This part extends
# Generic.MPolyBuildCtx, which does not exist yet when that file is loaded.

###############################################################################
#
#   Build Context
#
###############################################################################

mutable struct LaurentMPolyBuildCtx{T, S}
    coeffs::Vector{T}
    exps::Vector{Vector{Int}}
    parent::S

   function LaurentMPolyBuildCtx{T, S}(
         coeffs::Vector{T},
         exps::Vector{Vector{Int}},
         parent::S) where {T, S}
      return new{T, S}(coeffs, exps, parent)
   end
end

function MPolyBuildCtx(R::AbstractAlgebra.LaurentMPolyRing{T}) where T
    return LaurentMPolyBuildCtx{T, typeof(R)}(T[], Vector{Int}[], R)
end

function MPolyBuildCtx(R::AbstractAlgebra.UniversalRing{<:AbstractAlgebra.LaurentMPolyRingElem, T}) where T
    return LaurentMPolyBuildCtx{T, typeof(R)}(T[], Vector{Int}[], R)
end

function push_term!(B::LaurentMPolyBuildCtx{T, S}, c::U, expv::Vector{Int}) where {S, T, U}
    length(expv) == nvars(B.parent) || error("length of exponent vector should match the number of variables")
    push!(B.coeffs, c)
    push!(B.exps, expv)
    return B
end

function finish(B::LaurentMPolyBuildCtx{T, S}) where {T, S}
    res = B.parent(B.coeffs, B.exps)
    B.coeffs = T[]
    B.exps = Vector{Int}[]
    return res
end
