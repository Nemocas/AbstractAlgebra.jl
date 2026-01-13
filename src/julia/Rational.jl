###############################################################################
#
#   Rational.jl : Additional AbstractAlgebra functionality for Julia Rationals
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

const JuliaQQ = Rationals{BigInt}()

const qq = Rationals{Int}()

parent(a::Rational{T}) where T <: Integer = Rationals{T}()

elem_type(::Type{Rationals{T}}) where T <: Integer = Rational{T}

parent_type(::Type{Rational{T}}) where T <: Integer = Rationals{T}

base_ring_type(::Type{Rationals{T}}) where T <: Integer = Integers{T}

base_ring(a::Rationals{T}) where T <: Integer = Integers{T}()

is_exact_type(::Type{Rational{T}}) where T <: Integer = true

is_domain_type(::Type{Rational{T}}) where T <: Integer = true

###############################################################################
#
#   Basic manipulation
#
###############################################################################

zero(::Rationals{T}) where T <: Integer = Rational{T}(0)

one(::Rationals{T}) where T <: Integer = Rational{T}(1)

is_unit(a::Rational) = !is_zero(a)

is_zero_divisor(a::Rational) = is_zero(a)

canonical_unit(a::Rational) = iszero(a) ? one(a) : a

function numerator(a::Rational, canonicalise::Bool=true)
   return Base.numerator(a) # all other types ignore canonicalise
end

function denominator(a::Rational, canonicalise::Bool=true)
   return Base.denominator(a) # all other types ignore canonicalise
end

characteristic(::Rationals) = 0
is_known(::typeof(characteristic), ::Rationals) = true

###############################################################################
#
#   String I/O
#
###############################################################################

function expressify(a::Rational; context = nothing)
    n = numerator(a)
    d = denominator(a)
    if isone(d)
        return n
    else
        return Expr(:call, ://, n, d)
    end
end

function show(io::IO, R::Rationals{BigInt})
   print(io, "Rationals")
end

function show(io::IO, R::Rationals{T}) where T
   print(io, "Rationals{$T}()")
end

###############################################################################
#
#   Divides
#
###############################################################################

function divides(a::Rational{T}, b::Rational{T}) where T <: Integer
   return true, a//b
end

###############################################################################
#
#   Division with remainder
#
###############################################################################

function divrem(a::Rational{T}, b::Rational{T}) where T <: Integer
   return a//b, zero(Rational{T})
end

function div(a::Rational{T}, b::Rational{T}) where T <: Integer
   return a//b
end

###############################################################################
#
#   Exact division
#
###############################################################################

divexact(a::Rational, b::Integer; check::Bool=true) = a//b

divexact(a::Integer, b::Rational; check::Bool=true) = a//b

divexact(a::Rational, b::Rational; check::Bool=true) = a//b

function divides(a::T, b::T) where T <: Rational
   if iszero(b)
      return false, T(0)
   else
      return true, divexact(a, b; check=false)
   end
end

###############################################################################
#
#   Square root
#
###############################################################################

function sqrt(a::Rational{T}; check::Bool=true) where T <: Integer
   return sqrt(numerator(a, false); check=check)//sqrt(denominator(a, false); check=check)
end

function is_square(a::Rational{T}) where T <: Integer
   return is_square(numerator(a)) && is_square(denominator(a))
end

function is_square_with_sqrt(a::Rational{T}) where T <: Integer
   f1, s1 = is_square_with_sqrt(numerator(a))
   if !f1
      return false, zero(T)
   end
   f2, s2 = is_square_with_sqrt(denominator(a))
   if !f2
      return false, zero(T)
   end
   return true, s1//s2
end

###############################################################################
#
#   Root
#
###############################################################################

function root(a::Rational{T}, n::Int; check::Bool=true) where T <: Integer
   num = root(numerator(a, false), n; check=check)
   den = root(denominator(a, false), n; check=check)
   return num//den
end

function is_power(a::Rational{T}, n::Int) where T <: Integer
   f1, r1 = is_power(numerator(a), n)
   if !f1
      return false, zero(T)
   end
   f2, r2 = is_power(denominator(a), n)
   if !f2
      return false, zero(T)
   end
   return true, r1//r2
end

###############################################################################
#
#   Exponential
#
###############################################################################

@doc raw"""
    exp(a::Rational{T}) where T <: Integer

Return $1$ if $a = 0$, otherwise throw an exception.
"""
function exp(a::Rational{T}) where T <: Integer
   a != 0 && throw(DomainError(a, "a must be 0"))
   return Rational{T}(1)
end

@doc raw"""
    log(a::Rational{T}) where T <: Integer

Return $0$ if $a = 1$, otherwise throw an exception.
"""
function log(a::Rational{T}) where T <: Integer
   a != 1 && throw(DomainError(a, "a must be 1"))
end

###############################################################################
#
#   Random generation
#
###############################################################################

RandomExtensions.maketype(R::Rationals{T}, _) where {T} = Rational{T}

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make2{Rational{T}, Rationals{T}, <:AbstractArray{<:Integer}}}
              ) where {T}
   R, n = sp[][1:end]
   d = T(0)
   while iszero(d)
      d = T(rand(rng, n))
   end
   n = T(rand(rng, n))
   return R(n, d)
end


rand(rng::AbstractRNG, R::Rationals, n) = rand(rng, make(R, n))

rand(R::Rationals, n) = rand(Random.default_rng(), R, n)

###############################################################################
#
#   Conformance test element generation
#
###############################################################################

function ConformanceTests.generate_element(R::Rationals)
  B = base_ring(R)
  n = ConformanceTests.generate_element(B)
  d = ConformanceTests.generate_element(B)
  return is_zero(d) ? R(n) : R(n, d)
end

###############################################################################
#
#   valutaion / remove
#
###############################################################################

# TODO (CF):
# should be Bernstein'ed: this is slow for large valuations
# returns the maximal v s.th. z mod p^v == 0 and z div p^v
#   also useful if p is not prime....
#
# TODO: what happens to z = 0???
function remove(z::Rational{T}, p::T) where {T<:Integer}
    iszero(z) && return (0, z)
    v, d = remove(denominator(z), p)
    w, n = remove(numerator(z), p)
    return w - v, n // d
end

function valuation(z::Rational{T}, p::T) where {T<:Integer}
    iszero(z) && error("Not yet implemented")
    v = valuation(denominator(z), p)
    w = valuation(numerator(z), p)
    return w - v
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (R::Rationals{T})() where T <: Integer
   return Rational{T}(0)
end

function (R::Rationals{T})(b) where T <: Integer
   return Rational{T}(b)
end

function (R::Rationals{T})(b::Integer, c::Integer) where T <: Integer
   return Rational{T}(b, c)
end

###############################################################################
#
#   fraction_field constructor
#
###############################################################################

fraction_field(R::Integers{T}) where T <: Integer = Rationals{T}()
