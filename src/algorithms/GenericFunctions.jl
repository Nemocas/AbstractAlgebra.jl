###############################################################################
#
#   GenericFunctions.jl : Functions for Generic types
#
###############################################################################

# TODO: Move from CRT from Hecke to AbstractAlgebra?
# Currently no implementation, only example on how the arbitrary inputs `crt`
# should look like.
# @doc Markdown.doc"""
#     crt(r::AbstractVector{T}, m::AbstractVector{T}) where T
#     crt(r::T, m::T...) where T

# Return $x$ in the Euclidean domain $T$ such that $x \equiv r_i \mod m_i$
# for all $i$.
# """
function crt end

@doc Markdown.doc"""
    factor(a::T)

Return a factorization of the element $a$ as a `Fac{T}`.
"""
function factor end

@doc Markdown.doc"""
    factor_squarefree(a::T)

Return a squarefree factorization of the element $a$ as a `Fac{T}`.
"""
function factor_squarefree end

@doc Markdown.doc"""
    isirreducible(a)

Return `true` if $a$ is irreducible, else return `false`.
"""
function isirreducible end

@doc Markdown.doc"""
    issquarefree(a)

Return `true` if $a$ is squarefree, else return `false`.
"""
function issquarefree end

