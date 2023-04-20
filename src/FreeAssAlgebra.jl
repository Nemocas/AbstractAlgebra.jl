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

base_ring(a::FreeAssAlgElem{T}) where T <: RingElement = base_ring(parent(a))

coefficient_ring(a::FreeAssAlgElem{T}) where T <: RingElement = base_ring(a)

coefficient_ring(R::FreeAssAlgebra{T}) where T <: RingElement = base_ring(R)

function is_domain_type(::Type{S}) where {T <: RingElement, S <: FreeAssAlgElem{T}}
   return is_domain_type(T)
end

function is_exact_type(a::Type{S}) where {T <: RingElement, S <: FreeAssAlgElem{T}}
   return is_exact_type(T)
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

function show(io::IO, ::MIME"text/plain", a::FreeAssAlgebra)
  max_vars = 5 # largest number of variables to print
  n = nvars(a)
  print(io, "Free associative algebra")
  print(io, " on ", ItemQuantity(nvars(a), "indeterminate"), " ")
  if n > max_vars
    join(io, symbols(a)[1:max_vars - 1], ", ")
    println(io, "..., ", symbols(a)[n])
  else
    join(io, symbols(a), ", ")
    println(io)
  end
  io = pretty(io)
  print(io, Indent(), "over ", Lowercase(), base_ring(a))
end

function show(io::IO, a::FreeAssAlgebra)
  if get(io, :supercompact, false)
    # no nested printing
    print(io, "Free associative algebra")
  else
    # nested printing allowed, preferably supercompact
    io = pretty(io)
    print(io, "Free associative algebra on ", ItemQuantity(nvars(a), "indeterminate"))
    print(IOContext(io, :supercompact => true), " over ", Lowercase(), base_ring(a))
  end
end

###############################################################################
#
#   Basic Manipulation
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

@doc raw"""
    exponent_words(a::AbstractAlgebra.FreeAssAlgElem{T}) where T <: RingElement

Return an iterator for the exponent words of the given polynomial. To
retrieve an array of the exponent words, use `collect(exponent_words(a))`.
"""
function exponent_words(a::AbstractAlgebra.FreeAssAlgElem{T}) where T <: RingElement
   return Generic.FreeAssAlgExponentWords(a)
end

function is_unit(a::FreeAssAlgElem{T}) where T
   if is_constant(a)
      return is_unit(leading_coefficient(a))
   elseif is_domain_type(elem_type(coefficient_ring(a)))
      return false
   elseif length(a) == 1
      return false
   else
      throw(NotImplementedError(:is_unit, a))
   end
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
#   free_associative_algebra constructor
#
###############################################################################

function free_associative_algebra(
   R::AbstractAlgebra.Ring,
   s::AbstractVector{<:VarName};
   cached::Bool = true)

   S = [Symbol(v) for v in s]
   return Generic.free_associative_algebra(R, S, cached=cached)
end

function free_associative_algebra(
   R::AbstractAlgebra.Ring,
   s::Vector{Symbol};
   cached::Bool = true)

   return Generic.free_associative_algebra(R, s, cached=cached)
end

function free_associative_algebra(
   R::AbstractAlgebra.Ring,
   n::Int,
   s::VarName;
   cached::Bool = false)

   S = [Symbol(s, i) for i in 1:n]
   return Generic.free_associative_algebra(R, S; cached=cached)
end

function free_associative_algebra(
   R::AbstractAlgebra.Ring,
   n::Int,
   s::Symbol=:x;
   cached::Bool = false)

   S = [Symbol(s, i) for i in 1:n]
   return Generic.free_associative_algebra(R, S; cached=cached)
end

