###############################################################################
#
#   Poly.jl : Generic polynomials over rings
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{Poly{T}}) where T <: RingElement = PolyRing{T}

elem_type(::Type{PolyRing{T}}) where T <: RingElement = Poly{T}

@doc Markdown.doc"""
    dense_poly_type(::Type{T}) where T <: RingElement

Return the type of a polynomial whose coefficients have the given type.
"""
dense_poly_type(::Type{T}) where T <: RingElement = Generic.Poly{T}

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function setcoeff!(c::Poly{T}, n::Int, a::T) where T <: RingElement
   if !iszero(a) || n + 1 <= length(c)
      fit!(c, n + 1)
      c.coeffs[n + 1] = a
      c.length = max(length(c), n + 1)
      # don't normalise
   end
   return c
end

function normalise(a::Poly, n::Int)
   while n > 0 && iszero(a.coeffs[n])
      n -= 1
   end
   return n
end

coeff(a::Poly, n::Int) = n >= length(a) ? base_ring(a)(0) : a.coeffs[n + 1]

function deepcopy_internal(a::Poly{T}, dict::IdDict) where T <: RingElement
   coeffs = Array{T}(undef, length(a))
   for i = 1:length(a)
      coeffs[i] = deepcopy(a.coeffs[i])
   end
   return parent(a)(coeffs)
end

###############################################################################
#
#   Karatsuba multiplication
#
###############################################################################

@doc Markdown.doc"""
    mul_karatsuba(a::Poly{T}, b::Poly{T}, cutoff::Int) where T <: RingElement

Return $a \times b$ using the Karatsuba algorithm recursively for problems of
size roughly greater than `cutoff`.
"""
function mul_karatsuba(a::Poly{T}, b::Poly{T}, cutoff::Int) where T <: RingElement
   alen = length(a)
   blen = length(b)
   (alen < 1 || blen < 1) && return zero(parent(a))
   zlen = alen + blen - 1
   zcoeffs = Vector{T}(undef, zlen)
   AbstractAlgebra.DensePoly.mullow_fast!(zcoeffs, zlen,
                          a.coeffs, alen, b.coeffs, blen, base_ring(a), cutoff)
   z = parent(a)(zcoeffs)
   z = set_length!(z, normalise(z, zlen))
   return z
end

@doc Markdown.doc"""
    mullow_karatsuba(a::Poly{T}, b::Poly{T}, n::Int, cutoff::Int) where T <: RingElement

Return $a \times b$ truncated to $n$ terms using the Karatsuba algorithm
recursively for problems of size roughly greater than `cutoff`.
"""
function mullow_karatsuba(a::Poly{T}, b::Poly{T}, n::Int, cutoff::Int) where T <: RingElement
   alen = length(a)
   blen = length(b)
   (n < 1 || alen < 1 || blen < 1) && return zero(parent(a))
   zlen = min(alen + blen - 1, n)
   zcoeffs = Vector{T}(undef, zlen)
   AbstractAlgebra.DensePoly.mullow_fast!(zcoeffs, zlen,
                          a.coeffs, alen, b.coeffs, blen, base_ring(a), cutoff)
   z = parent(a)(zcoeffs)
   z = set_length!(z, normalise(z, zlen))
   return z
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function set_length!(c::Poly{T}, n::Int) where T <: RingElement
   if n < c.length
      for i = n + 1:c.length
         c.coeffs[i] = zero!(c.coeffs[i])
      end
   end
   c.length = n
   return c
end

function fit!(c::Poly{T}, n::Int) where T <: RingElement
   if length(c.coeffs) < n
      resize!(c.coeffs, n)
      for i = length(c) + 1:n
         c.coeffs[i] = zero(base_ring(c))
      end
   end
   return nothing
end

function zero!(c::Poly{T}) where T <: RingElement
   c = set_length!(c, 0)
   return c
end

function mul!(c::Poly{T}, a::Poly{T}, b::Poly{T}) where T <: RingElement
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

      c = set_length!(c, normalise(c, lenc))
   end
   return c
end

function addeq!(c::Poly{T}, a::Poly{T}) where T <: RingElement
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

function add!(c::Poly{T}, a::Poly{T}, b::Poly{T}) where T <: RingElement
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

promote_rule(::Type{Poly{T}}, ::Type{Poly{T}}) where T <: RingElement = Poly{T}

function promote_rule(::Type{Poly{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? Poly{T} : Union{}
end



###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (a::PolyRing{T})(b::RingElement) where T <: RingElement
   return a(base_ring(a)(b))
end

function (a::PolyRing{T})() where T <: RingElement
   z = Poly{T}()
   z.parent = a
   return z
end

function (a::PolyRing{T})(b::Union{Integer, Rational, AbstractFloat}) where T <: RingElement
   z = Poly{T}(base_ring(a)(b))
   z.parent = a
   return z
end

function (a::PolyRing{T})(b::T) where T <: RingElement
   parent(b) != base_ring(a) && error("Unable to coerce to polynomial")
   z = Poly{T}(b)
   z.parent = a
   return z
end

function (a::PolyRing{T})(b::AbstractAlgebra.PolyElem{T}) where T <: RingElement
   parent(b) != a && error("Unable to coerce polynomial")
   return b
end

function (a::PolyRing{T})(b::Array{T, 1}) where T <: RingElement
   R = base_ring(a)
   for i = 1:length(b)
      b[i] = R(b[i])
   end
   z = Poly{T}(b)
   z.parent = a
   return z
end

function (a::PolyRing{T})(b::Array{S, 1}) where {S <: RingElement, T <: RingElement}
   R = base_ring(a)
   len = length(b)
   entries = Array{T}(undef, len)
   for i = 1:length(b)
      entries[i] = R(b[i])
   end
   z = Poly{T}(entries)
   z.parent = a
   return z
end

# Functions to remove ambiguities on julia 0.7
function (a::PolyRing{T})(b::T) where {T <: Rational}
   parent(b) != base_ring(a) && error("Unable to coerce to polynomial")
   z = Poly{T}(b)
   z.parent = a
   return z
end

function (a::PolyRing{T})(b::T) where {T <: AbstractFloat}
   parent(b) != base_ring(a) && error("Unable to coerce to polynomial")
   z = Poly{T}(b)
   z.parent = a
   return z
end

function (a::PolyRing{T})(b::T) where {T <: Integer}
   parent(b) != base_ring(a) && error("Unable to coerce to polynomial")
   z = Poly{T}(b)
   z.parent = a
   return z
end

###############################################################################
#
#   PolynomialRing constructor
#
###############################################################################

function PolynomialRing(R::AbstractAlgebra.Ring, s::Symbol; cached::Bool = true)
   T = elem_type(R)
   parent_obj = PolyRing{T}(R, s, cached)

   return parent_obj, parent_obj([R(0), R(1)])
end
