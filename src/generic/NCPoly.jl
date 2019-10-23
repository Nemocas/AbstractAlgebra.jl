###############################################################################
#
#   NCPoly.jl : Generic polynomials over noncommutative rings
#
###############################################################################

export NCPolyRing, NCPoly

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{NCPoly{T}}) where T <: NCRingElem = NCPolyRing{T}

elem_type(::Type{NCPolyRing{T}}) where T <: NCRingElem = NCPoly{T}

@doc Markdown.doc"""
    base_ring(R::AbstractAlgebra.NCPolyRing{T}) where T <: NCRingElem
> Return the base ring of the given polynomial ring.
"""
base_ring(R::AbstractAlgebra.NCPolyRing{T}) where T <: NCRingElem = R.base_ring::parent_type(T)

function isexact_type(a::Type{T}) where {S <: NCRingElem, T <: AbstractAlgebra.NCPolyElem{S}}
   return isexact_type(S)
end

@doc Markdown.doc"""
    var(a::AbstractAlgebra.NCPolyRing)
> Return the internal name of the generator of the polynomial ring. Note that
> this is returned as a `Symbol` not a `String`.
"""
var(a::AbstractAlgebra.NCPolyRing) = a.S

@doc Markdown.doc"""
    symbols(a::AbstractAlgebra.NCPolyRing)
> Return an array of the variable names for the polynomial ring. Note that
> this is returned as an array of `Symbol` not `String`.
"""
symbols(a::AbstractAlgebra.NCPolyRing) = [a.S]

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::AbstractAlgebra.NCPolyElem, h::UInt)
   b = 0xd3f41ffbf953cbd8%UInt
   for i in 0:length(a) - 1
      b = xor(b, xor(hash(coeff(a, i), h), h))
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

function setcoeff!(c::NCPoly{T}, n::Int, a::T) where {T <: NCRingElem}
   if !iszero(a) || n + 1 <= length(c)
      fit!(c, n + 1)
      c.coeffs[n + 1] = a
      c.length = max(length(c), n + 1)
      # don't normalise
   end
   return c
end

function normalise(a::NCPoly, n::Int)
   while n > 0 && iszero(a.coeffs[n])
      n -= 1
   end
   return n
end

coeff(a::NCPoly, n::Int) = n >= length(a) ? base_ring(a)(0) : a.coeffs[n + 1]

@doc Markdown.doc"""
    zero(R::AbstractAlgebra.NCPolyRing)
> Return the zero polynomial in the given polynomial ring.
"""
zero(R::AbstractAlgebra.NCPolyRing) = R(0)

@doc Markdown.doc"""
    one(R::AbstractAlgebra.NCPolyRing)
> Return the constant polynomial $1$ in the given polynomial ring.
"""
one(R::AbstractAlgebra.NCPolyRing) = R(1)

@doc Markdown.doc"""
    gen(R::AbstractAlgebra.NCPolyRing)
> Return the generator of the given polynomial ring.
"""
gen(R::AbstractAlgebra.NCPolyRing) = R([zero(base_ring(R)), one(base_ring(R))])

isterm(a::T) where {T <: NCRingElem} = true

ismonomial(a::T) where {T <: NCRingElem} = isone(a)

function deepcopy_internal(a::NCPoly{T}, dict::IdDict) where {T <: NCRingElem}
   coeffs = Array{T}(undef, length(a))
   for i = 1:length(a)
      coeffs[i] = deepcopy(a.coeffs[i])
   end
   return parent(a)(coeffs)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, p::AbstractAlgebra.NCPolyRing)
   print(io, "Univariate Polynomial Ring in ")
   print(io, string(var(p)))
   print(io, " over ")
   print(IOContext(io, :compact => true), base_ring(p))
end

show_minus_one(::Type{NCPoly{T}}) where {T <: NCRingElem} = show_minus_one(T)

###############################################################################
#
#   Binary operations
#
###############################################################################

@doc Markdown.doc"""
    +(a::AbstractAlgebra.NCPolyElem{T}, b::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
> Return $a + b$.
"""
function +(a::AbstractAlgebra.NCPolyElem{T}, b::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
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
   set_length!(z, normalise(z, i - 1))
   return z
end

@doc Markdown.doc"""
    -(a::AbstractAlgebra.NCPolyElem{T}, b::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
> Return $a - b$.
"""
function -(a::AbstractAlgebra.NCPolyElem{T}, b::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
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
   set_length!(z, normalise(z, i - 1))
   return z
end

function *(a::AbstractAlgebra.NCPolyElem{T}, b::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
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
   set_length!(z, normalise(z, lenz))
   return z
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(a::T, b::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
   len = length(b)
   z = parent(b)()
   fit!(z, len)
   for i = 1:len
      z = setcoeff!(z, i - 1, a*coeff(b, i - 1))
   end
   set_length!(z, normalise(z, len))
   return z
end

function *(a::AbstractAlgebra.NCPolyElem{T}, b::T) where {T <: NCRingElem}
   len = length(a)
   z = parent(a)()
   fit!(z, len)
   for i = 1:len
      z = setcoeff!(z, i - 1, coeff(a, i - 1)*b)
   end
   set_length!(z, normalise(z, len))
   return z
end

###############################################################################
#
#   Powering
#
###############################################################################

@doc Markdown.doc"""
    ^(a::AbstractAlgebra.NCPolyElem{T}, b::Int) where {T <: NCRingElem}
> Return $a^b$. We require $b \geq 0$.
"""
function ^(a::AbstractAlgebra.NCPolyElem{T}, b::Int) where {T <: NCRingElem}
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
      set_length!(z, b + 1)
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
    ==(x::AbstractAlgebra.NCPolyElem{T}, y::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
> Return `true` if $x == y$ arithmetically, otherwise return `false`. Recall
> that power series to different precisions may still be arithmetically
> equal to the minimum of the two precisions.
"""
function ==(x::AbstractAlgebra.NCPolyElem{T}, y::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
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
    isequal(x::AbstractAlgebra.NCPolyElem{T}, y::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
> Return `true` if $x == y$ exactly, otherwise return `false`. This function is
> useful in cases where the coefficients of the polynomial are inexact, e.g.
> power series. Only if the power series are precisely the same, to the same
> precision, are they declared equal by this function.
"""
function isequal(x::AbstractAlgebra.NCPolyElem{T}, y::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
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
    ==(x::AbstractAlgebra.NCPolyElem{T}, y::T) where {T <: NCRingElem}
> Return `true` if $x == y$.
"""
==(x::AbstractAlgebra.NCPolyElem{T}, y::T) where T <: NCRingElem = ((length(x) == 0 && y == 0)
                        || (length(x) == 1 && coeff(x, 0) == y))

@doc Markdown.doc"""
    ==(x::T, y::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
> Return `true` if $x = y$.
"""
==(x::T, y::AbstractAlgebra.NCPolyElem{T}) where T <: NCRingElem = y == x

@doc Markdown.doc"""
    ==(x::Union{Integer, Rational, AbstractFloat}, y::AbstractAlgebra.NCPolyElem)
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::Union{Integer, Rational, AbstractFloat}, y::AbstractAlgebra.NCPolyElem) = y == x

###############################################################################
#
#   Truncation
#
###############################################################################

@doc Markdown.doc"""
    mullow(a::AbstractAlgebra.NCPolyElem{T}, b::AbstractAlgebra.NCPolyElem{T}, n::Int) where {T <: NCRingElem}
> Return $a\times b$ truncated to $n$ terms.
"""
function mullow(a::AbstractAlgebra.NCPolyElem{T}, b::AbstractAlgebra.NCPolyElem{T}, n::Int) where {T <: NCRingElem}
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
   set_length!(z, normalise(z, lenz))
   return z
end

###############################################################################
#
#   Exact division
#
###############################################################################

@doc Markdown.doc"""
    divexact_right(f::AbstractAlgebra.NCPolyElem{T}, g::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
> Assuming $f = qg$, return $q$.
"""
function divexact_right(f::AbstractAlgebra.NCPolyElem{T}, g::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
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
         set_length!(f, normalise(f, lenf - 1))
      end
   end
   q = parent(f)(d)
   set_length!(q, lenq)
   return q
end

@doc Markdown.doc"""
    divexact_left(f::AbstractAlgebra.NCPolyElem{T}, g::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
> Assuming $f = gq$, return $q$.
"""
function divexact_left(f::AbstractAlgebra.NCPolyElem{T}, g::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
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
         set_length!(f, normalise(f, lenf - 1))
      end
   end
   q = parent(f)(d)
   set_length!(q, lenq)
   return q
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

@doc Markdown.doc"""
    divexact_right(a::AbstractAlgebra.NCPolyElem{T}, b::T) where {T <: NCRingElem}
> Assuming $a = qb$, return $q$.
"""
function divexact_right(a::AbstractAlgebra.NCPolyElem{T}, b::T) where {T <: NCRingElem}
   iszero(b) && throw(DivideError())
   z = parent(a)()
   fit!(z, length(a))
   for i = 1:length(a)
      z = setcoeff!(z, i - 1, divexact_right(coeff(a, i - 1), b))
   end
   set_length!(z, length(a))
   return z
end

@doc Markdown.doc"""
    divexact_left(a::AbstractAlgebra.NCPolyElem{T}, b::T) where {T <: NCRingElem}
> Assuming $a = bq$, return $q$.
"""
function divexact_left(a::AbstractAlgebra.NCPolyElem{T}, b::T) where {T <: NCRingElem}
   iszero(b) && throw(DivideError())
   z = parent(a)()
   fit!(z, length(a))
   for i = 1:length(a)
      z = setcoeff!(z, i - 1, divexact_left(coeff(a, i - 1), b))
   end
   set_length!(z, length(a))
   return z
end

@doc Markdown.doc"""
    divexact_right(a::AbstractAlgebra.NCPolyElem, b::Union{Integer, Rational, AbstractFloat})
> Assuming $a = qb$, return $q$.
"""
function divexact_right(a::AbstractAlgebra.NCPolyElem, b::Union{Integer, Rational, AbstractFloat})
   iszero(b) && throw(DivideError())
   z = parent(a)()
   fit!(z, length(a))
   for i = 1:length(a)
      z = setcoeff!(z, i - 1, divexact_right(coeff(a, i - 1), b))
   end
   set_length!(z, length(a))
   return z
end

@doc Markdown.doc"""
    divexact_left(a::AbstractAlgebra.NCPolyElem, b::Union{Integer, Rational, AbstractFloat})
> Assuming $a = bq$, return $q$.
"""
divexact_left(a::AbstractAlgebra.NCPolyElem, b::Union{Integer, Rational, AbstractFloat}) = divexact_right(a, b)

###############################################################################
#
#   Evaluation
#
###############################################################################

@doc Markdown.doc"""
    evaluate(a::AbstractAlgebra.NCPolyElem, b::T) where {T <: NCRingElem}
> Evaluate the polynomial $a$ at the value $b$ and return the result.
"""
function evaluate(a::AbstractAlgebra.NCPolyElem, b::T) where {T <: NCRingElem}
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
    evaluate(a::AbstractAlgebra.NCPolyElem, b::Union{Integer, Rational, AbstractFloat})
> Evaluate the polynomial $a$ at the value $b$ and return the result.
"""
function evaluate(a::AbstractAlgebra.NCPolyElem, b::Union{Integer, Rational, AbstractFloat})
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

function set_length!(c::NCPoly{T}, n::Int) where T <: NCRingElem
   if n < c.length
      for i = n + 1:c.length
         c.coeffs[i] = zero!(c.coeffs[i])
      end
   end
   c.length = n
end

function fit!(c::NCPoly{T}, n::Int) where {T <: NCRingElem}
   if length(c.coeffs) < n
      t = c.coeffs
      c.coeffs = Array{T}(undef, n)
      for i = 1:length(c)
         c.coeffs[i] = t[i]
      end
      for i = length(c) + 1:n
         c.coeffs[i] = zero(base_ring(c))
      end
   end
   return nothing
end

function zero!(c::NCPoly{T}) where {T <: NCRingElem}
   set_length!(c, 0)
   return c
end

function mul!(c::AbstractAlgebra.NCPolyElem{T}, a::AbstractAlgebra.NCPolyElem{T}, b::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
   lena = length(a)
   lenb = length(b)

   if lena == 0 || lenb == 0
      set_length!(c, 0)
   else
      if a === c
         a = deepcopy(a)
      end
      if b === c
         b = deepcopy(b)
      end

      t = base_ring(a)()

      lenc = lena + lenb - 1
      fit!(c, lenc)

      for i = 1:lena
         c.coeffs[i] = mul!(c.coeffs[i], coeff(a, i - 1), coeff(b, 0))
      end

      for i = 2:lenb
         c.coeffs[lena + i - 1] = mul!(c.coeffs[lena + i - 1], coeff(a, lena - 1), coeff(b, i - 1))
      end

      for i = 1:lena - 1
         for j = 2:lenb
            t = mul!(t, coeff(a, i - 1), coeff(b, j - 1))
            c.coeffs[i + j - 1] = addeq!(c.coeffs[i + j - 1], t)
         end
      end

      set_length!(c, normalise(c, lenc))
   end
   return c
end

function addeq!(c::AbstractAlgebra.NCPolyElem{T}, a::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
   lenc = length(c)
   lena = length(a)
   len = max(lenc, lena)
   fit!(c, len)
   for i = 1:lena
      c.coeffs[i] = addeq!(c.coeffs[i], coeff(a, i - 1))
   end
   set_length!(c, normalise(c, len))
   return c
end

function add!(c::AbstractAlgebra.NCPolyElem{T}, a::AbstractAlgebra.NCPolyElem{T}, b::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
   lena = length(a)
   lenb = length(b)
   len = max(lena, lenb)
   fit!(c, len)
   i = 1
   while i <= min(lena, lenb)
      c.coeffs[i] = add!(c.coeffs[i], coeff(a, i - 1), coeff(b, i - 1))
      i += 1
   end
   while i <= lena
      # mutating operators must ensure they don't introduce new aliasing
      c = setcoeff!(c, i - 1, deepcopy(coeff(a, i - 1)))
      i += 1
   end
   while i <= lenb
      # mutating operators must ensure they don't introduce new aliasing
      c = setcoeff!(c, i - 1, deepcopy(coeff(b, i - 1)))
      i += 1
   end
   set_length!(c, normalise(c, len))
   return c
end

function addmul!(z::AbstractAlgebra.NCPolyElem{T}, x::AbstractAlgebra.NCPolyElem{T}, y::AbstractAlgebra.NCPolyElem{T}, c::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
   c = mul!(c, x, y)
   z = addeq!(z, c)
   return z
end

###############################################################################
#
#   Random elements
#
###############################################################################

function rand(rng::AbstractRNG, S::AbstractAlgebra.NCPolyRing, deg_range::UnitRange{Int}, v...)
   R = base_ring(S)
   f = S()
   x = gen(S)
   for i = 0:rand(rng, deg_range)
      f += rand(rng, R, v...)*x^i
   end
   return f
end

function rand(S::AbstractAlgebra.NCPolyRing, deg_range, v...)
   rand(Random.GLOBAL_RNG, S, deg_range, v...)
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{NCPoly{T}}, ::Type{NCPoly{T}}) where T <: NCRingElem = NCPoly{T}

function promote_rule(::Type{NCPoly{T}}, ::Type{U}) where {T <: NCRingElem, U <: NCRingElem}
   promote_rule(T, U) == T ? NCPoly{T} : Union{}
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::NCPolyRing{T})(b::NCRingElem) where {T <: NCRingElem}
   return a(base_ring(a)(b))
end

function (a::NCPolyRing{T})() where {T <: NCRingElem}
   z = NCPoly{T}()
   z.parent = a
   return z
end

function (a::NCPolyRing{T})(b::Union{Integer, Rational, AbstractFloat}) where {T <: NCRingElem}
   z = NCPoly{T}(base_ring(a)(b))
   z.parent = a
   return z
end

function (a::NCPolyRing{T})(b::T) where {T <: NCRingElem}
   parent(b) != base_ring(a) && error("Unable to coerce to polynomial")
   z = NCPoly{T}(b)
   z.parent = a
   return z
end

function (a::NCPolyRing{T})(b::AbstractAlgebra.NCPolyElem{T}) where {T <: NCRingElem}
   parent(b) != a && error("Unable to coerce polynomial")
   return b
end

function (a::NCPolyRing{T})(b::Array{T, 1}) where T <: NCRingElem
   R = base_ring(a)
   for i = 1:length(b)
      b[i] = R(b[i])
   end
   z = NCPoly{T}(b)
   z.parent = a
   return z
end

function (a::NCPolyRing{T})(b::Array{S, 1}) where {S <: RingElement, T <: NCRingElem}
   R = base_ring(a)
   len = length(b)
   entries = Array{T}(undef, len)
   for i = 1:length(b)
      entries[i] = R(b[i])
   end
   z = NCPoly{T}(entries)
   z.parent = a
   return z
end

function (a::NCPolyRing{T})(b::Array{S, 1}) where {S <: NCRingElem, T <: NCRingElem}
   R = base_ring(a)
   len = length(b)
   entries = Array{T}(undef, len)
   for i = 1:length(b)
      entries[i] = R(b[i])
   end
   z = NCPoly{T}(entries)
   z.parent = a
   return z
end

# Functions to remove ambiguities on julia 0.7
function (a::NCPolyRing{T})(b::T) where {T <: Rational}
   parent(b) != base_ring(a) && error("Unable to coerce to polynomial")
   z = NCPoly{T}(b)
   z.parent = a
   return z
end

function (a::NCPolyRing{T})(b::T) where {T <: AbstractFloat}
   parent(b) != base_ring(a) && error("Unable to coerce to polynomial")
   z = NCPoly{T}(b)
   z.parent = a
   return z
end

function (a::NCPolyRing{T})(b::T) where {T <: Integer}
   parent(b) != base_ring(a) && error("Unable to coerce to polynomial")
   z = NCPoly{T}(b)
   z.parent = a
   return z
end

###############################################################################
#
#   PolynomialRing constructor
#
###############################################################################

@doc Markdown.doc"""
    PolynomialRing(R::AbstractAlgebra.NCRing, s::AbstractString; cached::Bool = true)
> Given a base ring `R` and string `s` specifying how the generator (variable)
> should be printed, return a tuple `S, x` representing the new polynomial
> ring $S = R[x]$ and the generator $x$ of the ring. By default the parent
> object `S` will depend only on `R` and `x` and will be cached. Setting the
> optional argument `cached` to `false` will prevent the parent object `S` from
> being cached.
"""
function PolynomialRing(R::AbstractAlgebra.NCRing, s::AbstractString; cached::Bool = true)
   S = Symbol(s)
   T = elem_type(R)
   parent_obj = NCPolyRing{T}(R, S, cached)

   return parent_obj, parent_obj([R(0), R(1)])
end
