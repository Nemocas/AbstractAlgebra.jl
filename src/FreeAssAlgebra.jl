###############################################################################
#
#   FreeAssAlgebra.jl : free associative algebra R<x1,...,xn>
#
###############################################################################

export exponent_words, FreeAssociativeAlgebra, leading_exponent_word

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

base_ring(a::FreeAssAlgElem{T}) where T <: RingElement = base_ring(parent(a))

coefficient_ring(a::FreeAssAlgElem{T}) where T <: RingElement = base_ring(a)

coefficient_ring(R::FreeAssAlgebra{T}) where T <: RingElement = base_ring(R)

function isdomain_type(::Type{S}) where {T <: RingElement, S <: FreeAssAlgElem{T}}
   return isdomain_type(T)
end

function isexact_type(a::Type{S}) where {T <: RingElement, S <: FreeAssAlgElem{T}}
   return isexact_type(T)
end

###############################################################################
#
#   String IO
#
###############################################################################

function _expressify_word!(prod::Expr, x, v::Vector)
   j = -1
   e = 0
   for i in v
      if j != i
         if j > 0 && !iszero(e)
            push!(prod.args, e == 1 ? x[j] : Expr(:call, :^, x[j], e))
         end
         e = 0
      end
      j = i
      e += 1
   end
   if j > 0 && !iszero(e)
      push!(prod.args, e == 1 ? x[j] : Expr(:call, :^, x[j], e))
   end
end

function expressify(a::FreeAssAlgElem, x = symbols(parent(a)); context = nothing)
   sum = Expr(:call, :+)
   for (c, v) in zip(coefficients(a), exponent_words(a))
      prod = Expr(:call, :*)
      if !isone(c)
         push!(prod.args, expressify(c, context = context))
      end
      _expressify_word!(prod, x, v)
      push!(sum.args, prod)
   end
   return sum
end

@enable_all_show_via_expressify FreeAssAlgElem

function expressify(a::FreeAssAlgebra; context = nothing)
   return Expr(:sequence, Expr(:text, "Free associative algebra over "),
                          expressify(base_ring(a)),
                          Expr(:text, " on "),
                          Expr(:series, symbols(a)...))
end

@enable_all_show_via_expressify FreeAssAlgebra

###############################################################################
#
#   Coeffs, Terms, etc.
#
###############################################################################

function coefficients(a::AbstractAlgebra.FreeAssAlgElem)
   return Generic.MPolyCoeffs(a)
end

function terms(a::AbstractAlgebra.FreeAssAlgElem)
   return Generic.MPolyTerms(a)
end

function monomials(a::AbstractAlgebra.FreeAssAlgElem)
   return Generic.MPolyMonomials(a)
end

@doc Markdown.doc"""
    exponent_words(a::AbstractAlgebra.FreeAssAlgElem{T}) where T <: RingElement

Return an iterator for the exponent words of the given polynomial. To
retrieve an array of the exponent words, use `collect(exponent_words(a))`.
"""
function exponent_words(a::AbstractAlgebra.FreeAssAlgElem{T}) where T <: RingElement
   return Generic.FreeAssAlgExponentWords(a)
end

###############################################################################
#
# Hashing
#
###############################################################################

function Base.hash(x::FreeAssAlgElem{T}, h::UInt) where T <: RingElement
   b = 0x6220ed52502c8d9f%UInt
   for (c, e) in zip(coefficients(x), exponent_words(x))
      b = xor(b, xor(hash(c, h), h))
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
      b = xor(b, xor(hash(e, h), h))
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

###############################################################################
#
#   Random elements
#
###############################################################################

RandomExtensions.maketype(S::AbstractAlgebra.FreeAssAlgebra, _, _, _) = elem_type(S)

function RandomExtensions.make(S::AbstractAlgebra.FreeAssAlgebra,
                               term_range::UnitRange{Int},
                               exp_bound::UnitRange{Int}, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, term_range, exp_bound, vs[1])
   else
      make(S, term_range, exp_bound, make(R, vs...))
   end
end

function rand(rng::AbstractRNG, sp::SamplerTrivial{<:Make4{
                 <:NCRingElement,<:AbstractAlgebra.FreeAssAlgebra,UnitRange{Int},UnitRange{Int}}})
   S, term_range, exp_bound, v = sp[][1:end]
   f = S()
   g = gens(S)
   R = base_ring(S)
   isempty(g) && return S(rand(rng, v))
   for i = 1:rand(rng, term_range)
      term = S(1)
      for j = 1:rand(rng, exp_bound)
         term *= rand(g)
      end
      term *= rand(rng, v)
      f += term
   end
   return f
end

function rand(rng::AbstractRNG, S::AbstractAlgebra.FreeAssAlgebra,
              term_range::UnitRange{Int}, exp_bound::UnitRange{Int}, v...)
   m = make(S, term_range, exp_bound, v...)
   rand(rng, m)
end

function rand(S::AbstractAlgebra.FreeAssAlgebra, term_range, exp_bound, v...)
   rand(GLOBAL_RNG, S, term_range, exp_bound, v...)
end

###############################################################################
#
#   FreeAssociativeAlgebra constructor
#
###############################################################################

function FreeAssociativeAlgebra(
   R::AbstractAlgebra.Ring,
   s::Union{AbstractVector{<:AbstractString}, AbstractVector{Symbol}, AbstractVector{Char}};
   cached::Bool = true)

   S = [Symbol(v) for v in s]
   return Generic.FreeAssociativeAlgebra(R, S, cached=cached)
end

function FreeAssociativeAlgebra(
   R::AbstractAlgebra.Ring,
   s::Vector{Symbol};
   cached::Bool = true)

   return Generic.FreeAssociativeAlgebra(R, s, cached=cached)
end

function FreeAssociativeAlgebra(
   R::AbstractAlgebra.Ring,
   n::Int,
   s::Union{AbstractString, Symbol, Char};
   cached::Bool = false)

   S = [Symbol(s, i) for i in 1:n]
   return Generic.FreeAssociativeAlgebra(R, S; cached=cached)
end

function FreeAssociativeAlgebra(
   R::AbstractAlgebra.Ring,
   n::Int,
   s::Symbol=:x;
   cached::Bool = false)

   S = [Symbol(s, i) for i in 1:n]
   return Generic.FreeAssociativeAlgebra(R, S; cached=cached)
end

