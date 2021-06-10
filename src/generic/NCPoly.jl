###############################################################################
#
#   NCPoly.jl : Generic polynomials over noncommutative rings
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{NCPoly{T}}) where T <: NCRingElem = NCPolyRing{T}

elem_type(::Type{NCPolyRing{T}}) where T <: NCRingElem = NCPoly{T}

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function setcoeff!(c::NCPoly{T}, n::Int, a::T) where T <: NCRingElem
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

coeff(a::NCPoly, n::Int) = n >= length(a) ? coefficient_ring(a)(0) : a.coeffs[n + 1]

function deepcopy_internal(a::NCPoly{T}, dict::IdDict) where T <: NCRingElem
   coeffs = Array{T}(undef, length(a))
   for i = 1:length(a)
      coeffs[i] = deepcopy(a.coeffs[i])
   end
   return parent(a)(coeffs)
end

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
   return c
end

function fit!(c::NCPoly{T}, n::Int) where T <: NCRingElem
   if length(c.coeffs) < n
      resize!(c.coeffs, n)
      for i = length(c) + 1:n
         c.coeffs[i] = zero(coefficient_ring(c))
      end
   end
   return nothing
end

function zero!(c::NCPoly{T}) where T <: NCRingElem
   c = set_length!(c, 0)
   return c
end

function mul!(c::NCPoly{T}, a::NCPoly{T}, b::NCPoly{T}) where T <: NCRingElem
   lena = length(a)
   lenb = length(b)

   if lena == 0 || lenb == 0
      c = set_length!(c, 0)
   else
      if a === c
         a = deepcopy(a)
      end
      if b === c
         b = deepcopy(b)
      end

      t = coefficient_ring(a)()

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

      c = set_length!(c, normalise(c, lenc))
   end
   return c
end

function addeq!(c::NCPoly{T}, a::NCPoly{T}) where T <: NCRingElem
   lenc = length(c)
   lena = length(a)
   len = max(lenc, lena)
   fit!(c, len)
   for i = 1:lena
      c.coeffs[i] = addeq!(c.coeffs[i], coeff(a, i - 1))
   end
   c = set_length!(c, normalise(c, len))
   return c
end

function add!(c::NCPoly{T}, a::NCPoly{T}, b::NCPoly{T}) where T <: NCRingElem
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
   c = set_length!(c, normalise(c, len))
   return c
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

function (a::NCPolyRing{T})(b::NCRingElem) where T <: NCRingElem
   return a(coefficient_ring(a)(b))
end

function (a::NCPolyRing{T})() where T <: NCRingElem
   z = NCPoly{T}()
   z.parent = a
   return z
end

function (a::NCPolyRing{T})(b::Union{Integer, Rational, AbstractFloat}) where T <: NCRingElem
   z = NCPoly{T}(coefficient_ring(a)(b))
   z.parent = a
   return z
end

function (a::NCPolyRing{T})(b::T) where T <: NCRingElem
   parent(b) != coefficient_ring(a) && error("Unable to coerce to polynomial")
   z = NCPoly{T}(b)
   z.parent = a
   return z
end

function (a::NCPolyRing{T})(b::NCPolyElem{T}) where T <: NCRingElem
   parent(b) != a && error("Unable to coerce polynomial")
   return b
end

function (a::NCPolyRing{T})(b::Array{T, 1}) where T <: NCRingElem
   R = coefficient_ring(a)
   for i = 1:length(b)
      b[i] = R(b[i])
   end
   z = NCPoly{T}(b)
   z.parent = a
   return z
end

function (a::NCPolyRing{T})(b::Array{S, 1}) where {S <: RingElement, T <: NCRingElem}
   R = coefficient_ring(a)
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
   R = coefficient_ring(a)
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
   parent(b) != coefficient_ring(a) && error("Unable to coerce to polynomial")
   z = NCPoly{T}(b)
   z.parent = a
   return z
end

function (a::NCPolyRing{T})(b::T) where {T <: AbstractFloat}
   parent(b) != coefficient_ring(a) && error("Unable to coerce to polynomial")
   z = NCPoly{T}(b)
   z.parent = a
   return z
end

function (a::NCPolyRing{T})(b::T) where {T <: Integer}
   parent(b) != coefficient_ring(a) && error("Unable to coerce to polynomial")
   z = NCPoly{T}(b)
   z.parent = a
   return z
end

###############################################################################
#
#   PolynomialRing constructor
#
###############################################################################

function PolynomialRing(R::AbstractAlgebra.NCRing, s::Symbol; cached::Bool = true)
   T = elem_type(R)
   parent_obj = NCPolyRing{T}(R, s, cached)

   return parent_obj, parent_obj([R(0), R(1)])
end
