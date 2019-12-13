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

base_ring(a::Rational{Int}) = zz

base_ring(a::Rational{BigInt}) = JuliaZZ

base_ring(a::Rationals{Int}) = zz

base_ring(a::Rationals{BigInt}) = JuliaZZ

base_ring(a::Rationals{T}) where T <: Integer = Integers{T}()

base_ring(a::Rational{T}) where T <: Integer = Integers{T}()

isexact_type(::Type{Rational{T}}) where T <: Integer = true

isdomain_type(::Type{Rational{T}}) where T <: Integer = true

###############################################################################
#
#   Basic manipulation
#
###############################################################################

zero(::Rationals{T}) where T <: Integer = Rational{T}(0)

one(::Rationals{T}) where T <: Integer = Rational{T}(1)

isunit(a::Rational) = a != 0

canonical_unit(a::Rational)  = a

function numerator(a::Rational, canonicalise::Bool=true)
   return Base.numerator(a) # all other types ignore canonicalise
end

function denominator(a::Rational, canonicalise::Bool=true)
   return Base.denominator(a) # all other types ignore canonicalise
end

characteristic(a::Rational{T}) where T <: Integer = 0

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::Rationals)
   print(io, "Rationals")
end

needs_parentheses(::Rational) = false

displayed_with_minus_in_front(a::Rational) = a < 0

show_minus_one(::Type{Rational{T}}) where T <: Integer = false

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

divexact(a::Rational, b::Integer) = a//b

divexact(a::Integer, b::Rational) = a//b

divexact(a::Rational, b::Rational) = a//b

function divides(a::T, b::T) where T <: Rational
   if b == 0
      return false, T(0)
   else
      return true, divexact(a, b)
   end
end

###############################################################################
#
#   GCD
#
###############################################################################

function gcd(p::Rational{T}, q::Rational{T}) where T <: Integer
   a = p.num*q.den
   b = p.den*q.num
   n = gcd(a, b)
   d = p.den*q.den
   if d != 1 && n != 0
      g = gcd(n, d)
      n = divexact(n, g)
      d = divexact(d, g)
   end
   if n == 0
      return Rational{T}(n, T(1))
   else
      return Rational{T}(n, d)
   end
end

function gcdx(p::Rational{T}, q::Rational{T}) where {T <: Integer}
   g = gcd(p, q)
   if !iszero(p)
      return (g, g//p, zero(q))
   elseif !iszero(q)
      return (g, zero(p), g//q)
   else
      @assert iszero(g)
      return (g, zero(p), zero(p))
   end
end


###############################################################################
#
#   Square root
#
###############################################################################

@doc Markdown.doc"""
    sqrt(a::Rational{T}) where T <: Integer
> Return the square root of $a$ if it is the square of a rational, otherwise
> throw an error.
"""
function sqrt(a::Rational{T}) where T <: Integer
   return sqrt(numerator(a))//sqrt(denominator(a))
end

###############################################################################
#
#   Exponential
#
###############################################################################

@doc Markdown.doc"""
    exp(a::Rational{T}) where T <: Integer
> Return $1$ if $a = 0$, otherwise throw an exception.
"""
function exp(a::Rational{T}) where T <: Integer
   a != 0 && throw(DomainError(a, "a must be 0"))
   return Rational{T}(1)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(a::Rational{T}) where T <: Integer
   n = a.num
   n = zero!(n)
   if a.den == 1
      return Rational{T}(n, a.den)
   else
      return Rational{T}(n, T(1))
   end
end

function mul!(a::Rational{T}, b::Rational{T}, c::Rational{T}) where T <: Integer
   n = a.num
   d = a.den
   n = mul!(n, b.num, c.num)
   d = mul!(d, b.den, c.den)
   if d != 1 && n != 0
      g = gcd(n, d)
      n = divexact(n, g)
      d = divexact(d, g)
   end
   if n == 0
      return Rational{T}(n, T(1))
   else
      return Rational{T}(n, d)
   end
end

function add!(a::Rational{T}, b::Rational{T}, c::Rational{T}) where T <: Integer
   if a === b
      return addeq!(a, c)
   elseif a == c
      return addeq!(a, b)
   else # no aliasing
      n = a.num
      d = a.den
      d = mul!(d, b.den, c.den)
      n = mul!(n, b.num, c.den)
      n = addmul!(n, b.den, c.num)
      if d != 1 && n != 0
         g = gcd(n, d)
         n = divexact(n, g)
         d = divexact(d, g)
      end
      if n == 0
         return Rational{T}(n, T(1))
      else
         return Rational{T}(n, d)
      end
   end
end

function addeq!(a::Rational{T}, b::Rational{T}) where T <: Integer
   if a === b
      if iseven(a.den)
         return Rational{T}(a.num, div(b.den, 2))
      else
         return Rational{T}(2*a.num, b.den)
      end
   else
      n = a.num
      n = mul!(n, n, b.den)
      n = addmul!(n, b.num, a.den)
      d = a.den
      d = mul!(d, d, b.den)
      if d != 1 && n != 0
         g = gcd(n, d)
         n = divexact(n, g)
         d = divexact(d, g)
      end
      if n == 0
         return Rational{T}(n, T(1))
      else
         return Rational{T}(n, d)
      end
   end
end

function addmul!(a::Rational{T}, b::Rational{T}, c::Rational{T}, d::Rational{T}) where T <: Integer
   d = mul!(d, b, c)
   a = addeq!(a, d)
   return a
end

###############################################################################
#
#   Random generation
#
###############################################################################

function rand(rng::AbstractRNG, R::Rationals{T}, n::UnitRange{Int}) where T <: Integer
   d = T(0)
   while d == 0
      d = T(rand(rng, n))
   end
   n = T(rand(rng, n))
   return R(n, d)
end

rand(R::Rationals, n) = rand(Random.GLOBAL_RNG, R, n)

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
#   FractionField constructor
#
###############################################################################

FractionField(R::Integers{T}) where T <: Integer = Rationals{T}()
