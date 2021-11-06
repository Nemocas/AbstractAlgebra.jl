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

function isdomain_type(::Type{S}) where {T <: RingElement, S <: FreeAssAlgElem{T}}
   return isdomain_type(S)
end

function isexact_type(a::Type{S}) where {T <: RingElement, S <: FreeAssAlgElem{T}}
   return isexact_type(S)
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

function coefficients(a::FreeAssAlgElem)
   return Generic.MPolyCoeffs(a)
end

function terms(a::FreeAssAlgElem)
   return Generic.MPolyTerms(a)
end

function monomials(a::FreeAssAlgElem)
   return Generic.MPolyMonomials(a)
end

function exponent_words(a::FreeAssAlgElem)
   return Generic.FreeAssAlgExponentWords(a)
end

###############################################################################
#
#   FreeAssociativeAlgebra constructor
#
###############################################################################

function FreeAssociativeAlgebra(
   R::AbstractAlgebra.Ring,
   s::Vector{T};
   cached::Bool = true) where T <: Union{String, Char}

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
   s::Union{String, Char, Symbol};
   cached::Bool = false)

   S = [Symbol(s, i) for i in 1:n]
   return Generic.FreeAssociativeAlgebra(R, n, S; cached=cached)
end

