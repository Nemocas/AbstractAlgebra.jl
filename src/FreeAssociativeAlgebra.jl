###############################################################################
#
#   FreeAssociativeAlgebra.jl : free associative algebra R<x1,...,xn>
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

coefficient_ring(R::FreeAssociativeAlgebra{T}) where T <: RingElement = base_ring(R)

function is_domain_type(::Type{S}) where {T <: RingElement, S <: FreeAssociativeAlgebraElem{T}}
   return is_domain_type(T)
end

function is_exact_type(a::Type{S}) where {T <: RingElement, S <: FreeAssociativeAlgebraElem{T}}
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

function expressify(a::FreeAssociativeAlgebraElem, x = symbols(parent(a)); context = nothing)
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

@enable_all_show_via_expressify FreeAssociativeAlgebraElem

function show(io::IO, mime::MIME"text/plain", a::FreeAssociativeAlgebra)
  @show_name(io, a)
  @show_special(io, mime, a)

  max_vars = 5 # largest number of variables to print
  n = nvars(a)
  print(io, "Free associative algebra")
  print(io, " on ", ItemQuantity(nvars(a), "indeterminate"), " ")
  if n > max_vars
    join(io, symbols(a)[1:max_vars - 1], ", ")
    println(io, ", ..., ", symbols(a)[n])
  else
    join(io, symbols(a), ", ")
    println(io)
  end
  io = pretty(io)
  print(io, Indent(), "over ", Lowercase(), base_ring(a))
  print(io, Dedent())
end

function show(io::IO, a::FreeAssociativeAlgebra)
  @show_name(io, a)
  @show_special(io, a)
  if is_terse(io)
    print(io, "Free associative algebra")
  else
    io = pretty(io)
    print(io, "Free associative algebra on ", ItemQuantity(nvars(a), "indeterminate"))
    print(terse(io), " over ", Lowercase(), base_ring(a))
  end
end

###############################################################################
#
#   Basic Manipulation
#
###############################################################################

function coefficients(a::FreeAssociativeAlgebraElem)
   return Generic.MPolyCoeffs(a)
end

function terms(a::FreeAssociativeAlgebraElem)
   return Generic.MPolyTerms(a)
end

function monomials(a::FreeAssociativeAlgebraElem)
   return Generic.MPolyMonomials(a)
end

@doc raw"""
    exponent_words(a::FreeAssociativeAlgebraElem{T}) where T <: RingElement

Return an iterator for the exponent words of the given polynomial. To
retrieve an array of the exponent words, use `collect(exponent_words(a))`.
"""
function exponent_words(a::FreeAssociativeAlgebraElem{T}) where T <: RingElement
   return Generic.FreeAssAlgExponentWords(a)
end

function is_unit(a::FreeAssociativeAlgebraElem{T}) where T
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

function Base.hash(x::FreeAssociativeAlgebraElem{T}, h::UInt) where T <: RingElement
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
#   Evaluation
#
###############################################################################

@doc raw"""
    evaluate(a::FreeAssociativeAlgebraElem{T}, vals::Vector{U}) where {T <: RingElement, U <: NCRingElem}

Evaluate `a` by substituting in the array of values for each of the variables.
The evaluation will succeed if multiplication is defined between elements of
the coefficient ring of `a` and elements of `vals`.

The syntax `a(vals...)` is also supported.

# Examples

```jldoctest; setup = :(using AbstractAlgebra)
julia> R, (x, y) = free_associative_algebra(ZZ, ["x", "y"]);

julia> f = x*y - y*x
x*y - y*x

julia> S = matrix_ring(ZZ, 2);

julia> m1 = S([1 2; 3 4])
[1   2]
[3   4]

julia> m2 = S([0 1; 1 0])
[0   1]
[1   0]

julia> evaluate(f, [m1, m2])
[-1   -3]
[ 3    1]

julia> m1*m2 - m2*m1 == evaluate(f, [m1, m2])
true

julia> m1*m2 - m2*m1 == f(m1, m2)
true
```
"""
function evaluate(a::FreeAssociativeAlgebraElem{T}, vals::Vector{U}) where {T <: RingElement, U <: NCRingElem}
   length(vals) != nvars(parent(a)) && error("Number of variables does not match number of values")
   R = base_ring(parent(a))
   S = parent(one(R)*one(parent(vals[1])))
   r = zero(S)
   o = one(S)
   for (c, v) in zip(coefficients(a), exponent_words(a))
      r = add!(r, c*prod((vals[i] for i in v), init = o))
   end
   return r
end

function (a::FreeAssociativeAlgebraElem{T})(val::U, vals::U...) where {T <: RingElement, U <: NCRingElem}
   return evaluate(a, [val, vals...])
end

###############################################################################
#
#   Random elements
#
###############################################################################

RandomExtensions.maketype(S::FreeAssociativeAlgebra, _, _, _) = elem_type(S)

function RandomExtensions.make(S::FreeAssociativeAlgebra,
                               term_range::AbstractUnitRange{Int},
                               exp_bound::AbstractUnitRange{Int}, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, term_range, exp_bound, vs[1])
   else
      Make(S, term_range, exp_bound, make(R, vs...))
   end
end

function rand(rng::AbstractRNG, sp::SamplerTrivial{<:Make4{
                 <:NCRingElement, <:FreeAssociativeAlgebra, <:AbstractUnitRange{Int}, <:AbstractUnitRange{Int}}})
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

function rand(rng::AbstractRNG, S::FreeAssociativeAlgebra,
              term_range::AbstractUnitRange{Int}, exp_bound::AbstractUnitRange{Int}, v...)
   m = make(S, term_range, exp_bound, v...)
   rand(rng, m)
end

function rand(S::FreeAssociativeAlgebra, term_range, exp_bound, v...)
   rand(GLOBAL_RNG, S, term_range, exp_bound, v...)
end

@varnames_interface Generic.free_associative_algebra(R::Ring, s)
