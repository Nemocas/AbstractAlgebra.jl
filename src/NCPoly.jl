###############################################################################
#
#   NCPoly.jl : polynomials over noncommutative rings
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

base_ring(R::NCPolyRing{T}) where T <: NCRingElem = R.base_ring::parent_type(T)

coefficient_ring(R::NCPolyRing) = base_ring(R)

function isexact_type(a::Type{T}) where {S <: NCRingElem, T <: NCPolyElem{S}}
   return isexact_type(S)
end

@doc Markdown.doc"""
    var(a::NCPolyRing)

Return the internal name of the generator of the polynomial ring. Note that
this is returned as a `Symbol` not a `String`.
"""
var(a::NCPolyRing) = a.S

@doc Markdown.doc"""
    symbols(a::NCPolyRing)

Return an array of the variable names for the polynomial ring. Note that
this is returned as an array of `Symbol` not `String`.
"""
symbols(a::NCPolyRing) = [a.S]

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::NCPolyElem, h::UInt)
   b = 0xd3f41ffbf953cbd8%UInt
   for i in 0:length(a) - 1
      b = xor(b, xor(hash(coeff(a, i), h), h))
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

zero(R::NCPolyRing) = R(0)

one(R::NCPolyRing) = R(1)

@doc Markdown.doc"""
    gen(R::NCPolyRing)

Return the generator of the given polynomial ring.
"""
gen(R::NCPolyRing) = R([zero(base_ring(R)), one(base_ring(R))])

isterm(a::T) where T <: NCRingElem = true

ismonomial_monomial(a::T) where T <: NCRingElem = isone(a)

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, p::NCPolyRing)
   print(io, "Univariate Polynomial Ring in ")
   print(io, string(var(p)))
   print(io, " over ")
   print(IOContext(io, :compact => true), base_ring(p))
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function +(a::NCPolyElem{T}, b::NCPolyElem{T}) where T <: NCRingElem
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)
   lenz = max(lena, lenb)
   z = parent(a)()
   fit!(z, lenz)
   i = 1
   while i <= min(lena, lenb)
      z = setcoeff!(z, i - 1, coeff(a, i - 1) + coeff(b, i - 1))
      i += 1
   end
   while i <= lena
      z = setcoeff!(z, i - 1, deepcopy(coeff(a, i - 1)))
      i += 1
   end
   while i <= lenb
      z = setcoeff!(z, i - 1, deepcopy(coeff(b, i - 1)))
      i += 1
   end
   z = set_length!(z, normalise(z, i - 1))
   return z
end

function -(a::NCPolyElem{T}, b::NCPolyElem{T}) where T <: NCRingElem
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)
   lenz = max(lena, lenb)
   z = parent(a)()
   fit!(z, lenz)
   i = 1
   while i <= min(lena, lenb)
      z = setcoeff!(z, i - 1, coeff(a, i - 1) - coeff(b, i - 1))
      i += 1
   end
   while i <= lena
      z = setcoeff!(z, i - 1, deepcopy(coeff(a, i - 1)))
      i += 1
   end
   while i <= lenb
      z = setcoeff!(z, i - 1, -coeff(b, i - 1))
      i += 1
   end
   z = set_length!(z, normalise(z, i - 1))
   return z
end

function *(a::NCPolyElem{T}, b::NCPolyElem{T}) where T <: NCRingElem
   lena = length(a)
   lenb = length(b)
   if lena == 0 || lenb == 0
      return parent(a)()
   end
   t = base_ring(a)()
   lenz = lena + lenb - 1
   d = Array{T}(undef, lenz)
   for i = 1:lena
      d[i] = coeff(a, i - 1)*coeff(b, 0)
   end
   for i = 2:lenb
      d[lena + i - 1] = a.coeffs[lena]*coeff(b, i - 1)
   end
   for i = 1:lena - 1
      for j = 2:lenb
         t = mul!(t, coeff(a, i - 1), b.coeffs[j])
         d[i + j - 1] = addeq!(d[i + j - 1], t)
      end
   end
   z = parent(a)(d)
   z = set_length!(z, normalise(z, lenz))
   return z
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(a::T, b::NCPolyElem{T}) where T <: NCRingElem
   len = length(b)
   z = parent(b)()
   fit!(z, len)
   for i = 1:len
      z = setcoeff!(z, i - 1, a*coeff(b, i - 1))
   end
   z = set_length!(z, normalise(z, len))
   return z
end

function *(a::NCPolyElem{T}, b::T) where T <: NCRingElem
   len = length(a)
   z = parent(a)()
   fit!(z, len)
   for i = 1:len
      z = setcoeff!(z, i - 1, coeff(a, i - 1)*b)
   end
   z = set_length!(z, normalise(z, len))
   return z
end

###############################################################################
#
#   Powering
#
###############################################################################

@doc Markdown.doc"""
    ^(a::NCPolyElem{T}, b::Int) where T <: NCRingElem

Return $a^b$. We require $b \geq 0$.
"""
function ^(a::NCPolyElem{T}, b::Int) where T <: NCRingElem
   b < 0 && throw(DomainError(b, "exponent must be >= 0"))
   # special case powers of x for constructing polynomials efficiently
   R = parent(a)
   if isgen(a)
      z = R()
      fit!(z, b + 1)
      z = setcoeff!(z, b, deepcopy(coeff(a, 1)))
      for i = 1:b
         z = setcoeff!(z, i - 1, deepcopy(coeff(a, 0)))
      end
      z = set_length!(z, b + 1)
      return z
   elseif length(a) == 0
      return zero(R)
   elseif length(a) == 1
      return R(coeff(a, 0)^b)
   elseif b == 0
      return one(R)
   elseif b == 1
      return deepcopy(a)
   else
      bit = ~((~UInt(0)) >> 1)
      while (UInt(bit) & b) == 0
         bit >>= 1
      end
      z = a
      bit >>= 1
      while bit != 0
         z = z*z
         if (UInt(bit) & b) != 0
            z *= a
         end
         bit >>= 1
      end
      return z
   end
end

###############################################################################
#
#   Comparisons
#
###############################################################################

@doc Markdown.doc"""
    ==(x::NCPolyElem{T}, y::NCPolyElem{T}) where T <: NCRingElem

Return `true` if $x == y$ arithmetically, otherwise return `false`. Recall
that power series to different precisions may still be arithmetically
equal to the minimum of the two precisions.
"""
function ==(x::NCPolyElem{T}, y::NCPolyElem{T}) where T <: NCRingElem
   b = check_parent(x, y, false)
   !b && return false
   if length(x) != length(y)
      return false
   else
      for i = 1:length(x)
         if coeff(x, i - 1) != coeff(y, i - 1)
            return false
         end
      end
   end
   return true
end

@doc Markdown.doc"""
    isequal(x::NCPolyElem{T}, y::NCPolyElem{T}) where T <: NCRingElem

Return `true` if $x == y$ exactly, otherwise return `false`. This function is
useful in cases where the coefficients of the polynomial are inexact, e.g.
power series. Only if the power series are precisely the same, to the same
precision, are they declared equal by this function.
"""
function isequal(x::NCPolyElem{T}, y::NCPolyElem{T}) where T <: NCRingElem
   if parent(x) != parent(y)
      return false
   end
   if length(x) != length(y)
      return false
   end
   for i = 1:length(x)
      if !isequal(coeff(x, i - 1), coeff(y, i - 1))
         return false
      end
   end
   return true
end

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

@doc Markdown.doc"""
    ==(x::NCPolyElem{T}, y::T) where T <: NCRingElem

Return `true` if $x == y$.
"""
==(x::NCPolyElem{T}, y::T) where T <: NCRingElem = ((length(x) == 0 && y == 0)
                        || (length(x) == 1 && coeff(x, 0) == y))

@doc Markdown.doc"""
    ==(x::T, y::NCPolyElem{T}) where T <: NCRingElem

Return `true` if $x = y$.
"""
==(x::T, y::NCPolyElem{T}) where T <: NCRingElem = y == x

@doc Markdown.doc"""
    ==(x::Union{Integer, Rational, AbstractFloat}, y::NCPolyElem)

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::Union{Integer, Rational, AbstractFloat}, y::NCPolyElem) = y == x

###############################################################################
#
#   Truncation
#
###############################################################################

@doc Markdown.doc"""
    mullow(a::NCPolyElem{T}, b::NCPolyElem{T}, n::Int) where T <: NCRingElem

Return $a\times b$ truncated to $n$ terms.
"""
function mullow(a::NCPolyElem{T}, b::NCPolyElem{T}, n::Int) where T <: NCRingElem
   check_parent(a, b)
   lena = length(a)
   lenb = length(b)
   if lena == 0 || lenb == 0
      return zero(parent(a))
   end
   if n < 0
      n = 0
   end
   t = base_ring(a)()
   lenz = min(lena + lenb - 1, n)
   d = Array{T}(undef, lenz)
   for i = 1:min(lena, lenz)
      d[i] = coeff(a, i - 1)*coeff(b, 0)
   end
   if lenz > lena
      for j = 2:min(lenb, lenz - lena + 1)
          d[lena + j - 1] = coeff(a, lena - 1)*coeff(b, j - 1)
      end
   end
   for i = 1:lena - 1
      if lenz > i
         for j = 2:min(lenb, lenz - i + 1)
            t = mul!(t, coeff(a, i - 1), coeff(b, j - 1))
            d[i + j - 1] = addeq!(d[i + j - 1], t)
         end
      end
   end
   z = parent(a)(d)
   z = set_length!(z, normalise(z, lenz))
   return z
end

###############################################################################
#
#   Exact division
#
###############################################################################

@doc Markdown.doc"""
    divexact_right(f::NCPolyElem{T}, g::NCPolyElem{T}) where T <: NCRingElem

Assuming $f = qg$, return $q$.
"""
function divexact_right(f::NCPolyElem{T}, g::NCPolyElem{T}) where T <: NCRingElem
   check_parent(f, g)
   iszero(g) && throw(DivideError())
   if iszero(f)
      return zero(parent(f))
   end
   lenq = length(f) - length(g) + 1
   d = Array{T}(undef, lenq)
   for i = 1:lenq
      d[i] = zero(base_ring(f))
   end
   x = gen(parent(f))
   leng = length(g)
   while length(f) >= leng
      lenf = length(f)
      q1 = d[lenf - leng + 1] = divexact_right(coeff(f, lenf - 1), coeff(g, leng - 1))
      f = f - shift_left(q1*g, lenf - leng)
      if length(f) == lenf # inexact case
         f = set_length!(f, normalise(f, lenf - 1))
      end
   end
   q = parent(f)(d)
   q = set_length!(q, lenq)
   return q
end

@doc Markdown.doc"""
    divexact_left(f::NCPolyElem{T}, g::NCPolyElem{T}) where T <: NCRingElem

Assuming $f = gq$, return $q$.
"""
function divexact_left(f::NCPolyElem{T}, g::NCPolyElem{T}) where T <: NCRingElem
   check_parent(f, g)
   iszero(g) && throw(DivideError())
   if iszero(f)
      return zero(parent(f))
   end
   lenq = length(f) - length(g) + 1
   d = Array{T}(undef, lenq)
   for i = 1:lenq
      d[i] = zero(base_ring(f))
   end
   x = gen(parent(f))
   leng = length(g)
   while length(f) >= leng
      lenf = length(f)
      q1 = d[lenf - leng + 1] = divexact_left(coeff(f, lenf - 1), coeff(g, leng - 1))
      f = f - shift_left(g*q1, lenf - leng)
      if length(f) == lenf # inexact case
         f = set_length!(f, normalise(f, lenf - 1))
      end
   end
   q = parent(f)(d)
   q = set_length!(q, lenq)
   return q
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

@doc Markdown.doc"""
    divexact_right(a::NCPolyElem{T}, b::T) where T <: NCRingElem

Assuming $a = qb$, return $q$.
"""
function divexact_right(a::NCPolyElem{T}, b::T) where T <: NCRingElem
   iszero(b) && throw(DivideError())
   z = parent(a)()
   fit!(z, length(a))
   for i = 1:length(a)
      z = setcoeff!(z, i - 1, divexact_right(coeff(a, i - 1), b))
   end
   z = set_length!(z, length(a))
   return z
end

@doc Markdown.doc"""
    divexact_left(a::NCPolyElem{T}, b::T) where T <: NCRingElem

Assuming $a = bq$, return $q$.
"""
function divexact_left(a::NCPolyElem{T}, b::T) where T <: NCRingElem
   iszero(b) && throw(DivideError())
   z = parent(a)()
   fit!(z, length(a))
   for i = 1:length(a)
      z = setcoeff!(z, i - 1, divexact_left(coeff(a, i - 1), b))
   end
   z = set_length!(z, length(a))
   return z
end

@doc Markdown.doc"""
    divexact_right(a::NCPolyElem, b::Union{Integer, Rational, AbstractFloat})

Assuming $a = qb$, return $q$.
"""
function divexact_right(a::NCPolyElem, b::Union{Integer, Rational, AbstractFloat})
   iszero(b) && throw(DivideError())
   z = parent(a)()
   fit!(z, length(a))
   for i = 1:length(a)
      z = setcoeff!(z, i - 1, divexact_right(coeff(a, i - 1), b))
   end
   z = set_length!(z, length(a))
   return z
end

@doc Markdown.doc"""
    divexact_left(a::NCPolyElem, b::Union{Integer, Rational, AbstractFloat})

Assuming $a = bq$, return $q$.
"""
divexact_left(a::NCPolyElem, b::Union{Integer, Rational, AbstractFloat}) = divexact_right(a, b)

###############################################################################
#
#   Evaluation
#
###############################################################################

@doc Markdown.doc"""
    evaluate(a::NCPolyElem, b::T) where T <: NCRingElem

Evaluate the polynomial $a$ at the value $b$ and return the result.
"""
function evaluate(a::NCPolyElem, b::T) where T <: NCRingElem
   i = length(a)
   R = base_ring(a)
   if i == 0
       return zero(R)
   end
   z = R(coeff(a, i - 1))
   while i > 1
      i -= 1
      z = R(coeff(a, i - 1)) + z*b
      parent(z) # To work around a bug in julia
   end
   return z
end

@doc Markdown.doc"""
    evaluate(a::NCPolyElem, b::Union{Integer, Rational, AbstractFloat})

Evaluate the polynomial $a$ at the value $b$ and return the result.
"""
function evaluate(a::NCPolyElem, b::Union{Integer, Rational, AbstractFloat})
   i = length(a)
   R = base_ring(a)
   if i == 0
       return zero(R)
   end
   z = R(coeff(a, i - 1))
   while i > 1
      i -= 1
      z = R(coeff(a, i - 1)) + z*b
      parent(z) # To work around a bug in julia
   end
   return z
end

# Note: composition is not associative, e.g. consider fo(goh) vs (fog)oh
# for f and g of degree 2 and h of degree 1 -- and recall coeffs don't commute

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function addmul!(z::NCPolyElem{T}, x::NCPolyElem{T}, y::NCPolyElem{T}, c::NCPolyElem{T}) where T <: NCRingElem
   c = mul!(c, x, y)
   z = addeq!(z, c)
   return z
end

###############################################################################
#
#   Random elements
#
###############################################################################

RandomExtensions.maketype(S::NCPolyRing, dr::UnitRange{Int}, _) = elem_type(S)

function RandomExtensions.make(S::NCPolyRing, deg_range::UnitRange{Int}, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, deg_range, vs[1]) # forward to default Make constructor
   else
      make(S, deg_range, make(R, vs...))
   end
end

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make3{<:NCPolyElem,
                                         <:NCPolyRing,
                                         UnitRange{Int}}})
   S, deg_range, v = sp[][1:end]
   R = base_ring(S)
   f = S()
   x = gen(S)
   for i = 0:rand(rng, deg_range)
      f += rand(rng, v)*x^i
   end
   return f
end

rand(rng::AbstractRNG, S::NCPolyRing, deg_range::UnitRange{Int}, v...) =
   rand(rng, make(S, deg_range, v...))

rand(S::NCPolyRing, deg_range, v...) = rand(Random.GLOBAL_RNG, S, deg_range, v...)


###############################################################################
#
#   PolynomialRing constructor
#
###############################################################################

@doc Markdown.doc"""
    PolynomialRing(R::NCRing, s::Union{AbstractString, Char, Symbol}; cached::Bool = true)

Given a base ring `R` and string `s` specifying how the generator (variable)
should be printed, return a tuple `S, x` representing the new polynomial
ring $S = R[x]$ and the generator $x$ of the ring. By default the parent
object `S` will depend only on `R` and `x` and will be cached. Setting the
optional argument `cached` to `false` will prevent the parent object `S` from
being cached.
"""
PolynomialRing(R::NCRing, s::Union{AbstractString, Char, Symbol}; cached::Bool = true)

function PolynomialRing(R::NCRing, s::Symbol; cached::Bool = true)
   return Generic.PolynomialRing(R, s, cached=cached)
end

function PolynomialRing(R::NCRing, s::AbstractString; cached::Bool = true)
   return Generic.PolynomialRing(R, Symbol(s), cached=cached)
end

function PolynomialRing(R::NCRing, s::Char; cached::Bool = true)
   return Generic.PolynomialRing(R, Symbol(s); cached=cached)
end
