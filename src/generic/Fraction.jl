###############################################################################
#
#   Fraction.jl : generic fraction fields
#
###############################################################################

export FractionField

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{Frac{T}}) where T <: RingElem = FracField{T}

elem_type(::Type{FracField{T}}) where {T <: RingElem} = Frac{T}

@doc Markdown.doc"""
    base_ring(a::AbstractAlgebra.FracField{T}) where T <: RingElem
> Return the base ring $R$ of the given fraction field.
"""
base_ring(a::AbstractAlgebra.FracField{T}) where T <: RingElem = a.base_ring::parent_type(T)

@doc Markdown.doc"""
    base_ring(a::AbstractAlgebra.FracElem)
> Return the base ring $R$ of the fraction field that the supplied
> element $a$ belongs to.
"""
base_ring(a::AbstractAlgebra.FracElem) = base_ring(parent(a))

@doc Markdown.doc"""
    parent(a::AbstractAlgebra.FracElem)
> Return the parent object of the given fraction element.
"""
parent(a::AbstractAlgebra.FracElem) = a.parent

function isdomain_type(::Type{T}) where {S <: RingElement, T <: AbstractAlgebra.FracElem{S}}
   return isdomain_type(S)
end

function isexact_type(a::Type{T}) where {S <: RingElement, T <: AbstractAlgebra.FracElem{S}}
   return isexact_type(S)
end

@doc Markdown.doc"""
    characteristic(R::AbstractAlgebra.FracField{T}) where T <: RingElem
> Return the characteristic of the given field.
"""
function characteristic(R::AbstractAlgebra.FracField{T}) where T <: RingElem
   return characteristic(base_ring(R))
end

function check_parent(a::AbstractAlgebra.FracElem, b::AbstractAlgebra.FracElem, throw::Bool = true)
   fl = parent(a) != parent(b)
   fl && throw && error("Incompatible rings in fraction field operation")
   return !fl
end

###############################################################################
#
#   Constructors
#
###############################################################################

function //(x::T, y::T) where {T <: RingElem}
   iszero(y) && throw(DivideError())
   g = gcd(x, y)
   z = Frac{T}(divexact(x, g), divexact(y, g))
   try
      z.parent = FracDict[R]
   catch
      z.parent = FractionField(parent(x))
   end
   return z
end

//(x::T, y::AbstractAlgebra.FracElem{T}) where {T <: RingElem} = parent(y)(x)//y

//(x::AbstractAlgebra.FracElem{T}, y::T) where {T <: RingElem} = x//parent(x)(y)

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::AbstractAlgebra.FracElem, h::UInt)
   b = 0x8a30b0d963237dd5%UInt
   # We canonicalise before hashing
   return xor(b, hash(AbstractAlgebra.numerator(a, true), h), hash(AbstractAlgebra.denominator(a, true), h), h)
end

function Base.numerator(a::Frac, canonicalise::Bool=true)
   if canonicalise
      u = canonical_unit(a.den)
      return divexact(a.num, u)
   else
      return a.num
   end
end

function Base.denominator(a::Frac, canonicalise::Bool=true)
   if canonicalise
      u = canonical_unit(a.den)
      return divexact(a.den, u)
   else
      return a.den
   end
end

# Fall back method for all other fraction types in system
function Base.numerator(a::AbstractAlgebra.FracElem, canonicalise::Bool=true)
   return Base.numerator(a) # all other types ignore canonicalise
end

# Fall back method for all other fraction types in system
function Base.denominator(a::AbstractAlgebra.FracElem, canonicalise::Bool=true)
   return Base.denominator(a) # all other types ignore canonicalise
end

@doc Markdown.doc"""
    zero(R::AbstractAlgebra.FracField)
> Return $0/1$ in the given fraction field.
"""
zero(R::AbstractAlgebra.FracField) = R(0)

@doc Markdown.doc"""
    one(R::AbstractAlgebra.FracField)
> Return $1/1$ in the given fraction field.
"""
one(R::AbstractAlgebra.FracField) = R(1)

@doc Markdown.doc"""
    iszero(a::AbstractAlgebra.FracElem)
> Return `true` if the supplied element $a$ is zero in the fraction field it
> belongs to, otherwise return `false`.
"""
iszero(a::AbstractAlgebra.FracElem) = iszero(AbstractAlgebra.numerator(a, false))

@doc Markdown.doc"""
    isone(a::AbstractAlgebra.FracElem)
> Return `true` if the supplied element $a$ is one in the fraction field it
> belongs to, otherwise return `false`.
"""
isone(a::AbstractAlgebra.FracElem) = AbstractAlgebra.numerator(a, false) == AbstractAlgebra.denominator(a, false)

@doc Markdown.doc"""
    isunit(a::AbstractAlgebra.FracElem)
> Return `true` if the supplied element $a$ is invertible in the fraction field
> it belongs to, i.e. the AbstractAlgebra.numerator is nonzero, otherwise return `false`.
"""
isunit(a::AbstractAlgebra.FracElem) = !iszero(AbstractAlgebra.numerator(a, false))

function deepcopy_internal(a::Frac{T}, dict::IdDict) where {T <: RingElem}
   v = Frac{T}(deepcopy(AbstractAlgebra.numerator(a, false)), deepcopy(AbstractAlgebra.denominator(a, false)))
   v.parent = parent(a)
   return v
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(a::AbstractAlgebra.FracElem) = a

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function show(io::IO, x::AbstractAlgebra.FracElem)
   # Canonicalise for display
   n = AbstractAlgebra.numerator(x, true)
   d = AbstractAlgebra.denominator(x, true)
   if !isone(d) && needs_parentheses(n)
      print(io, "(")
   end
   print(IOContext(io, :compact => true), n)
   if !isone(d)
      if needs_parentheses(n)
         print(io, ")")
      end
      print(io, "//")
      print(io, "(") # always print parentheses for denoninators e.g. x//(x*y*z)
      print(IOContext(io, :compact => true), d)
      print(io, ")")
   end
end

function show(io::IO, a::AbstractAlgebra.FracField)
   print(IOContext(io, :compact => true), "Fraction field of ", base_ring(a))
end

# Parentheses are only needed for fractions if we didn't print them already
needs_parentheses(x::AbstractAlgebra.FracElem) = isone(AbstractAlgebra.denominator(x, true)) &&
                                     needs_parentheses(AbstractAlgebra.numerator(x, true))

function displayed_with_minus_in_front(x::AbstractAlgebra.FracElem)
   n = AbstractAlgebra.numerator(x, true)
   return !needs_parentheses(n) && displayed_with_minus_in_front(n)
end

function show_minus_one(::Type{AbstractAlgebra.FracElem{T}}) where {T <: RingElem}
   return show_minus_one(T)
end

###############################################################################
#
#   Unary operators
#
###############################################################################

@doc Markdown.doc"""
    -(a::AbstractAlgebra.FracElem)
> Return $-a$.
"""
function -(a::AbstractAlgebra.FracElem)
   return parent(a)(-AbstractAlgebra.numerator(a, false), deepcopy(AbstractAlgebra.denominator(a, false)))
end

###############################################################################
#
#   Binary operators
#
###############################################################################

@doc Markdown.doc"""
    +(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
> Return $a + b$.
"""
function +(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   d1 = AbstractAlgebra.denominator(a, false)
   d2 = AbstractAlgebra.denominator(b, false)
   n1 = AbstractAlgebra.numerator(a, false)
   n2 = AbstractAlgebra.numerator(b, false)
   gd = gcd(d1, d2)
   if d1 == d2
      rnum = n1 + n2
      if isone(d1)
         rden = deepcopy(d1)
      else
         gd = gcd(rnum, d1)
         if isone(gd)
            rden = deepcopy(d1)
         else
            rnum = divexact(rnum, gd)
            rden = divexact(d1, gd)
         end
      end
   elseif isone(d1)
      rnum = n1*d2 + n2
      rden = deepcopy(d2)
   elseif isone(d2)
      rnum = n1 + n2*d1
      rden = deepcopy(d1)
   else
      if isone(gd)
         rnum = n1*d2 + n2*d1
         rden = d1*d2
      else
         q1 = divexact(d1, gd)
         q2 = divexact(d2, gd)
         rnum = q1*n2 + q2*n1
         t = gcd(rnum, gd)
         if isone(t)
            rden = q2*d1
         else
            rnum = divexact(rnum, t)
            gd = divexact(d1, t)
            rden = gd*q2
         end
      end
   end
   return parent(a)(rnum, rden)
end

@doc Markdown.doc"""
    -(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
> Return $a - b$.
"""
function -(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   d1 = AbstractAlgebra.denominator(a, false)
   d2 = AbstractAlgebra.denominator(b, false)
   n1 = AbstractAlgebra.numerator(a, false)
   n2 = AbstractAlgebra.numerator(b, false)
   if d1 == d2
      rnum = n1 - n2
      if isone(d1)
         rden = deepcopy(d1)
      else
         gd = gcd(rnum, d1)
         if isone(gd)
            rden = deepcopy(d1)
         else
            rnum = divexact(rnum, gd)
            rden = divexact(d1, gd)
         end
      end
   elseif isone(d1)
      rnum = n1*d2 - n2
      rden = deepcopy(d2)
   elseif isone(d2)
      rnum = n1 - n2*d1
      rden = deepcopy(d1)
   else
      gd = gcd(d1, d2)
      if isone(gd)
         rnum = n1*d2 - n2*d1
         rden = d1*d2
      else
         q1 = divexact(d1, gd)
         q2 = divexact(d2, gd)
         rnum = q2*n1 - q1*n2
         t = gcd(rnum, gd)
         if isone(t)
            rden = q2*d1
         else
            rnum = divexact(rnum, t)
            gd = divexact(d1, t)
            rden = gd*q2
         end
      end
   end
   return parent(a)(rnum, rden)
end

@doc Markdown.doc"""
    *(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
> Return $a\times b$.
"""
function *(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   n1 = AbstractAlgebra.numerator(a, false)
   d2 = AbstractAlgebra.denominator(b, false)
   n2 = AbstractAlgebra.numerator(b, false)
   d1 = AbstractAlgebra.denominator(a, false)
   if d1 == d2
      n = n1*n2
      d = d1*d2
   elseif isone(d1)
      gd = gcd(n1, d2)
      if isone(gd)
         n = n1*n2
         d = deepcopy(d2)
      else
         n = divexact(n1, gd)*n2
         d = divexact(d2, gd)
      end
   elseif isone(d2)
      gd = gcd(n2, d1)
      if isone(gd)
         n = n2*n1
         d = deepcopy(d1)
      else
         n = divexact(n2, gd)*n1
         d = divexact(d1, gd)
      end
   else
      g1 = gcd(n1, d2)
      g2 = gcd(n2, d1)
      if !isone(g1)
         n1 = divexact(n1, g1)
         d2 = divexact(d2, g1)
      end
      if !isone(g2)
         n2 = divexact(n2, g2)
         d1 = divexact(d1, g2)
      end
      n = n1*n2
      d = d1*d2
   end
   return parent(a)(n, d)
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

@doc Markdown.doc"""
    *(a::AbstractAlgebra.FracElem, b::Union{Integer, Rational, AbstractFloat})
> Return $a\times b$.
"""
function *(a::AbstractAlgebra.FracElem, b::Union{Integer, Rational, AbstractFloat})
   c = base_ring(a)(b)
   g = gcd(AbstractAlgebra.denominator(a, false), c)
   n = AbstractAlgebra.numerator(a, false)*divexact(c, g)
   d = divexact(AbstractAlgebra.denominator(a, false), g)
   return parent(a)(n, d)
end

@doc Markdown.doc"""
    *(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.FracElem)
> Return $a\times b$.
"""
function *(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.FracElem)
   c = base_ring(b)(a)
   g = gcd(AbstractAlgebra.denominator(b, false), c)
   n = AbstractAlgebra.numerator(b, false)*divexact(c, g)
   d = divexact(AbstractAlgebra.denominator(b, false), g)
   return parent(b)(n, d)
end

@doc Markdown.doc"""
    *(a::AbstractAlgebra.FracElem{T}, b::T) where {T <: RingElem}
> Return $a\times b$.
"""
function *(a::AbstractAlgebra.FracElem{T}, b::T) where {T <: RingElem}
   g = gcd(AbstractAlgebra.denominator(a, false), b)
   n = AbstractAlgebra.numerator(a, false)*divexact(b, g)
   d = divexact(AbstractAlgebra.denominator(a, false), g)
   return parent(a)(n, d)
end

@doc Markdown.doc"""
    *(a::T, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
> Return $a\times b$.
"""
function *(a::T, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   g = gcd(AbstractAlgebra.denominator(b, false), a)
   n = AbstractAlgebra.numerator(b, false)*divexact(a, g)
   d = divexact(AbstractAlgebra.denominator(b, false), g)
   return parent(b)(n, d)
end

@doc Markdown.doc"""
    +(a::AbstractAlgebra.FracElem, b::Union{Integer, Rational, AbstractFloat})
> Return $a + b$.
"""
function +(a::AbstractAlgebra.FracElem, b::Union{Integer, Rational, AbstractFloat})
   n = AbstractAlgebra.numerator(a, false) + AbstractAlgebra.denominator(a, false)*b
   d = AbstractAlgebra.denominator(a, false)
   return parent(a)(n, deepcopy(d))
end

@doc Markdown.doc"""
    -(a::AbstractAlgebra.FracElem, b::Union{Integer, Rational, AbstractFloat})
> Return $a - b$.
"""
function -(a::AbstractAlgebra.FracElem, b::Union{Integer, Rational, AbstractFloat})
   n = AbstractAlgebra.numerator(a, false) - AbstractAlgebra.denominator(a, false)*b
   d = AbstractAlgebra.denominator(a, false)
   return parent(a)(n, deepcopy(d))
end

@doc Markdown.doc"""
    +(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.FracElem)
> Return $a + b$.
"""
+(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.FracElem) = b + a

@doc Markdown.doc"""
    -(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.FracElem)
> Return $a - b$.
"""
function -(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.FracElem)
   n = a*AbstractAlgebra.denominator(b, false) - AbstractAlgebra.numerator(b, false)
   d = AbstractAlgebra.denominator(b, false)
   return parent(b)(n, deepcopy(d))
end

@doc Markdown.doc"""
    +(a::AbstractAlgebra.FracElem{T}, b::T) where {T <: RingElem}
> Return $a + b$.
"""
function +(a::AbstractAlgebra.FracElem{T}, b::T) where {T <: RingElem}
   n = AbstractAlgebra.numerator(a, false) + AbstractAlgebra.denominator(a, false)*b
   d = AbstractAlgebra.denominator(a, false)
   return parent(a)(n, deepcopy(d))
end

@doc Markdown.doc"""
    -(a::AbstractAlgebra.FracElem{T}, b::T) where {T <: RingElem}
> Return $a - b$.
"""
function -(a::AbstractAlgebra.FracElem{T}, b::T) where {T <: RingElem}
   n = AbstractAlgebra.numerator(a, false) - AbstractAlgebra.denominator(a, false)*b
   d = AbstractAlgebra.denominator(a, false)
   return parent(a)(n, deepcopy(d))
end

@doc Markdown.doc"""
    +(a::T, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
> Return $a + b$.
"""
+(a::T, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem} = b + a

@doc Markdown.doc"""
    -(a::T, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
> Return $a - b$.
"""
function -(a::T, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   n = a*AbstractAlgebra.denominator(b, false) - AbstractAlgebra.numerator(b, false)
   d = AbstractAlgebra.denominator(b, false)
   return parent(b)(n, deepcopy(d))
end

###############################################################################
#
#   Comparisons
#
###############################################################################

@doc Markdown.doc"""
    ==(x::AbstractAlgebra.FracElem{T}, y::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
> Return `true` if $x == y$ arithmetically, otherwise return `false`. Recall
> that power series to different precisions may still be arithmetically
> equal to the minimum of the two precisions.
"""
function ==(x::AbstractAlgebra.FracElem{T}, y::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   b  = check_parent(x, y, false)
   !b && return false

   return (AbstractAlgebra.denominator(x, false) == AbstractAlgebra.denominator(y, false) &&
           AbstractAlgebra.numerator(x, false) == AbstractAlgebra.numerator(y, false)) ||
          (AbstractAlgebra.denominator(x, true) == AbstractAlgebra.denominator(y, true) &&
           AbstractAlgebra.numerator(x, true) == AbstractAlgebra.numerator(y, true)) ||
          (AbstractAlgebra.numerator(x, false)*AbstractAlgebra.denominator(y, false) ==
           AbstractAlgebra.denominator(x, false)*AbstractAlgebra.numerator(y, false))
end

@doc Markdown.doc"""
    isequal(x::AbstractAlgebra.FracElem{T}, y::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
> Return `true` if $x == y$ exactly, otherwise return `false`. This function is
> useful in cases where the AbstractAlgebra.numerators and AbstractAlgebra.denominators of the fractions are
> inexact, e.g. power series. Only if the power series are precisely the same,
> to the same precision, are they declared equal by this function.
"""
function isequal(x::AbstractAlgebra.FracElem{T}, y::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   if parent(x) != parent(y)
      return false
   end
   return isequal(AbstractAlgebra.numerator(x, false)*AbstractAlgebra.denominator(y, false),
                  AbstractAlgebra.denominator(x, false)*AbstractAlgebra.numerator(y, false))
end

###############################################################################
#
#   Ad hoc comparisons
#
###############################################################################

@doc Markdown.doc"""
    ==(x::AbstractAlgebra.FracElem, y::Union{Integer, Rational, AbstractFloat})
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
function ==(x::AbstractAlgebra.FracElem, y::Union{Integer, Rational, AbstractFloat})
   return (isone(AbstractAlgebra.denominator(x, false)) && AbstractAlgebra.numerator(x, false) == y) ||
          (isone(AbstractAlgebra.denominator(x, true)) && AbstractAlgebra.numerator(x, true) == y) ||
          (AbstractAlgebra.numerator(x, false) == AbstractAlgebra.denominator(x, false)*y)
end

@doc Markdown.doc"""
    ==(x::Union{Integer, Rational, AbstractFloat}, y::AbstractAlgebra.FracElem)
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::Union{Integer, Rational, AbstractFloat}, y::AbstractAlgebra.FracElem) = y == x

@doc Markdown.doc"""
    ==(x::AbstractAlgebra.FracElem{T}, y::T) where {T <: RingElem}
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
function ==(x::AbstractAlgebra.FracElem{T}, y::T) where {T <: RingElem}
   return (isone(AbstractAlgebra.denominator(x, false)) && AbstractAlgebra.numerator(x, false) == y) ||
          (isone(AbstractAlgebra.denominator(x, true)) && AbstractAlgebra.numerator(x, true) == y) ||
          (AbstractAlgebra.numerator(x, false) == AbstractAlgebra.denominator(x, false)*y)
end

@doc Markdown.doc"""
    ==(x::T, y::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::T, y::AbstractAlgebra.FracElem{T}) where {T <: RingElem} = y == x

###############################################################################
#
#   Inversion
#
###############################################################################

@doc Markdown.doc"""
    inv(a::AbstractAlgebra.FracElem)
> Return the inverse of the fraction $a$.
"""
function inv(a::AbstractAlgebra.FracElem)
   iszero(AbstractAlgebra.numerator(a, false)) && throw(DivideError())
   return parent(a)(deepcopy(AbstractAlgebra.denominator(a, false)),
                    deepcopy(AbstractAlgebra.numerator(a, false)))
end

###############################################################################
#
#   Exact division
#
###############################################################################

@doc Markdown.doc"""
    divexact(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
> Return $a/b$.
"""
function divexact(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   n1 = AbstractAlgebra.numerator(a, false)
   d2 = AbstractAlgebra.denominator(b, false)
   n2 = AbstractAlgebra.numerator(b, false)
   d1 = AbstractAlgebra.denominator(a, false)
   if d1 == n2
      n = n1*d2
      d = d1*n2
   elseif isone(d1)
      gd = gcd(n1, n2)
      if isone(gd)
         n = n1*d2
         d = deepcopy(n2)
      else
         n = divexact(n1, gd)*d2
         d = divexact(n2, gd)
      end
   elseif isone(n2)
      gd = gcd(d2, d1)
      if isone(gd)
         n = d2*n1
         d = deepcopy(d1)
      else
         n = divexact(d2, gd)*n1
         d = divexact(d1, gd)
      end
   else
      g1 = gcd(n1, n2)
      g2 = gcd(d2, d1)
      if !isone(g1)
         n1 = divexact(n1, g1)
         n2 = divexact(n2, g1)
      end
      if !isone(g2)
         d2 = divexact(d2, g2)
         d1 = divexact(d1, g2)
      end
      n = n1*d2
      d = d1*n2
   end
   return parent(a)(n, d)
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

@doc Markdown.doc"""
    divexact(a::AbstractAlgebra.FracElem, b::Union{Integer, Rational, AbstractFloat})
> Return $a/b$.
"""
function divexact(a::AbstractAlgebra.FracElem, b::Union{Integer, Rational, AbstractFloat})
   b == 0 && throw(DivideError())
   c = base_ring(a)(b)
   g = gcd(AbstractAlgebra.numerator(a, false), c)
   n = divexact(AbstractAlgebra.numerator(a, false), g)
   d = AbstractAlgebra.denominator(a, false)*divexact(c, g)
   return parent(a)(n, d)
end

@doc Markdown.doc"""
    divexact(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.FracElem)
> Return $a/b$.
"""
function divexact(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.FracElem)
   iszero(b) && throw(DivideError())
   c = base_ring(b)(a)
   g = gcd(AbstractAlgebra.numerator(b, false), c)
   n = AbstractAlgebra.denominator(b, false)*divexact(c, g)
   d = divexact(AbstractAlgebra.numerator(b, false), g)
   return parent(b)(n, d)
end

@doc Markdown.doc"""
    divexact(a::AbstractAlgebra.FracElem{T}, b::T) where {T <: RingElem}
> Return $a/b$.
"""
function divexact(a::AbstractAlgebra.FracElem{T}, b::T) where {T <: RingElem}
   iszero(b) && throw(DivideError())
   g = gcd(AbstractAlgebra.numerator(a, false), b)
   n = divexact(AbstractAlgebra.numerator(a, false), g)
   d = AbstractAlgebra.denominator(a, false)*divexact(b, g)
   return parent(a)(n, d)
end

@doc Markdown.doc"""
    divexact(a::T, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
> Return $a/b$.
"""
function divexact(a::T, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   iszero(b) && throw(DivideError())
   g = gcd(AbstractAlgebra.numerator(b, false), a)
   n = AbstractAlgebra.denominator(b, false)*divexact(a, g)
   d = divexact(AbstractAlgebra.numerator(b, false), g)
   return parent(b)(n, d)
end

function divides(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   if iszero(a)
     return true, parent(a)()
   end
   if iszero(b)
     return false, parent(a)()
   end
   return true, divexact(a, b)
end

###############################################################################
#
#   Powering
#
###############################################################################

@doc Markdown.doc"""
    ^(a::AbstractAlgebra.FracElem, b::Int)
> Return $a^b$.
"""
function ^(a::AbstractAlgebra.FracElem{T}, b::Int) where {T <: RingElem}
   if b < 0
      a = inv(a)
      b = -b
   end
   return parent(a)(AbstractAlgebra.numerator(a)^b, AbstractAlgebra.denominator(a)^b)
end

###############################################################################
#
#   GCD
#
###############################################################################

@doc Markdown.doc"""
    gcd(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
> Return a greatest common divisor of $a$ and $b$ if one exists. N.B: we define
> the GCD of $a/b$ and $c/d$ to be gcd$(ad, bc)/bd$, reduced to lowest terms.
> This requires the existence of a greatest common divisor function for the
> base ring.
"""
function gcd(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   gbd = gcd(AbstractAlgebra.denominator(a, false), AbstractAlgebra.denominator(b, false))
   n = gcd(AbstractAlgebra.numerator(a, false), AbstractAlgebra.numerator(b, false))
   d = divexact(AbstractAlgebra.denominator(a, false), gbd)*AbstractAlgebra.denominator(b, false)
   n = divexact(n, canonical_unit(n))
   d = divexact(d, canonical_unit(d))
   return parent(a)(n, d)
end

################################################################################
#
#   Remove and valuation
#
################################################################################

@doc Markdown.doc"""
    remove(z::AbstractAlgebra.FracElem{T}, p::T) where {T <: RingElem}
> Return the tuple $n, x$ such that $z = p^nx$ where $x$ has valuation $0$ at
> $p$.
"""
function remove(z::AbstractAlgebra.FracElem{T}, p::T) where {T <: RingElem}
   iszero(z) && error("Not yet implemented")
   v, d = remove(AbstractAlgebra.denominator(z, false), p)
   w, n = remove(AbstractAlgebra.numerator(z, false), p)
   return w-v, parent(z)(deepcopy(n), deepcopy(d))
end

@doc Markdown.doc"""
    valuation(z::AbstractAlgebra.FracElem{T}, p::T) where {T <: RingElem}
> Return the valuation of $z$ at $p$.
"""
function valuation(z::AbstractAlgebra.FracElem{T}, p::T) where {T <: RingElem}
   v, _ = remove(z, p)
   return v
end

###############################################################################
#
#   Unsafe operators and functions
#
###############################################################################

function zero!(c::AbstractAlgebra.FracElem)
   c.num = zero!(c.num)
   if !isone(c.den)
      c.den = one(base_ring(c))
   end
   return c
end

function mul!(c::AbstractAlgebra.FracElem{T}, a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   n1 = AbstractAlgebra.numerator(a, false)
   d2 = AbstractAlgebra.denominator(b, false)
   n2 = AbstractAlgebra.numerator(b, false)
   d1 = AbstractAlgebra.denominator(a, false)
   if d1 == d2
      c.num = n1*n2
      c.den = d1*d2
   elseif isone(d1)
      gd = gcd(n1, d2)
      if isone(gd)
         c.num = n1*n2
         c.den = deepcopy(d2)
      else
         c.num = divexact(n1, gd)*n2
         c.den = divexact(d2, gd)
      end
   elseif isone(d2)
      gd = gcd(n2, d1)
      if isone(gd)
         c.num = n2*n1
         c.den = deepcopy(d1)
      else
         c.num = divexact(n2, gd)*n1
         c.den = divexact(d1, gd)
      end
   else
      g1 = gcd(n1, d2)
      g2 = gcd(n2, d1)
      if !isone(g1)
         n1 = divexact(n1, g1)
         d2 = divexact(d2, g1)
      end
      if !isone(g2)
         n2 = divexact(n2, g2)
         d1 = divexact(d1, g2)
      end
      c.num = n1*n2
      c.den = d1*d2
   end
   return c
end

function addeq!(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   d1 = AbstractAlgebra.denominator(a, false)
   d2 = AbstractAlgebra.denominator(b, false)
   n1 = AbstractAlgebra.numerator(a, false)
   n2 = AbstractAlgebra.numerator(b, false)
   gd = gcd(d1, d2)
   if d1 == d2
      a.num = addeq!(a.num, b.num)
      if !isone(d1)
         gd = gcd(a.num, d1)
         if !isone(gd)
            a.num = divexact(a.num, gd)
            a.den = divexact(d1, gd)
         end
      end
   elseif isone(d1)
      if n1 !== n2
         a.num = mul!(a.num, a.num, d2)
         a.num = addeq!(a.num, n2)
      else
         a.num = n1*d2 + n2
      end
      a.den = deepcopy(d2)
   elseif isone(d2)
      a.num = addeq!(a.num, n2*d1)
      a.den = deepcopy(d1)
   else
      if isone(gd)
         if n1 !== n2
            a.num = mul!(a.num, a.num, d2)
            a.num = addeq!(a.num, n2*d1)
         else
            a.num = n1*d2 + n2*d1
         end
         a.den = d1*d2
      else
         q1 = divexact(d1, gd)
         q2 = divexact(d2, gd)
         a.num = q1*n2 + q2*n1
         t = gcd(a.num, gd)
         if isone(t)
            a.den = mul!(a.den, a.den, q2)
         else
            gd = divexact(d1, t)
            a.num = divexact(a.num, t)
            a.den = gd*q2
         end
      end
   end
   return a
end

function add!(c::AbstractAlgebra.FracElem{T}, a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   d1 = AbstractAlgebra.denominator(a, false)
   d2 = AbstractAlgebra.denominator(b, false)
   n1 = AbstractAlgebra.numerator(a, false)
   n2 = AbstractAlgebra.numerator(b, false)
   gd = gcd(d1, d2)
   if d1 == d2
      c.num = n1 + n2
      if isone(d1)
         c.den = deepcopy(d1)
      else
         gd = gcd(c.num, d1)
         if isone(gd)
            c.den = deepcopy(d1)
         else
            c.num = divexact(c.num, gd)
            c.den = divexact(d1, gd)
         end
      end
   elseif isone(d1)
      c.num = n1*d2 + n2
      c.den = deepcopy(d2)
   elseif isone(d2)
      c.num = n1 + n2*d1
      c.den = deepcopy(d1)
   else
      if isone(gd)
         c.num = n1*d2 + n2*d1
         c.den = d1*d2
      else
         q1 = divexact(d1, gd)
         q2 = divexact(d2, gd)
         c.num = q1*n2 + q2*n1
         t = gcd(c.num, gd)
         if isone(t)
            c.den = q2*d1
         else
            gd = divexact(d1, t)
            c.num = divexact(c.num, t)
            c.den = gd*q2
         end
      end
   end
   return c
end

###############################################################################
#
#   Random functions
#
###############################################################################

function rand(S::AbstractAlgebra.FracField{T}, v...) where {T <: RingElem}
   R = base_ring(S)
   n = rand(R, v...)
   d = R()
   while iszero(d)
      d = rand(R, v...)
   end
   return S(n, d)
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{Frac{T}}, ::Type{Frac{T}}) where T <: RingElement = Frac{T}

function promote_rule(::Type{Frac{T}}, ::Type{U}) where {T <: RingElem, U <: RingElem}
   promote_rule(T, U) == T ? Frac{T} : Union{}
end

###############################################################################
#
#   Parent object call overloading
#
###############################################################################

function (a::FracField{T})(b::RingElement) where {T <: RingElement}
   return a(base_ring(a)(b))
end

function (a::FracField{T})() where {T <: RingElement}
   z = Frac{T}(zero(base_ring(a)), one(base_ring(a)))
   z.parent = a
   return z
end

function (a::FracField{T})(b::T) where {T <: RingElement}
   parent(b) != base_ring(a) && error("Could not coerce to fraction")
   z = Frac{T}(b, one(base_ring(a)))
   z.parent = a
   return z
end

function (a::FracField{T})(b::T, c::T) where {T <: RingElement}
   parent(b) != base_ring(a) && error("Could not coerce to fraction")
   parent(c) != base_ring(a) && error("Could not coerce to fraction")
   z = Frac{T}(b, c)
   z.parent = a
   return z
end

function (a::FracField{T})(b::T, c::Union{Integer, Rational, AbstractFloat}) where {T <: RingElement}
   parent(b) != base_ring(a) && error("Could not coerce to fraction")
   z = Frac{T}(b, base_ring(a)(c))
   z.parent = a
   return z
end

function (a::FracField{T})(b::Union{Integer, Rational, AbstractFloat}, c::T) where {T <: RingElement}
   parent(c) != base_ring(a) && error("Could not coerce to fraction")
   z = Frac{T}(base_ring(a)(b), c)
   z.parent = a
   return z
end

function (a::FracField{T})(b::Union{Integer, Rational, AbstractFloat}) where {T <: RingElement}
   z = Frac{T}(base_ring(a)(b), one(base_ring(a)))
   z.parent = a
   return z
end

function (a::FracField{T})(b::Integer, c::Integer) where {T <: RingElement}
   z = Frac{T}(base_ring(a)(b), base_ring(a)(c))
   z.parent = a
   return z
end

function (a::FracField{T})(b::Frac{T}) where {T <: RingElement}
   a != parent(b) && error("Could not coerce to fraction")
   return b
end

###############################################################################
#
#   FractionField constructor
#
###############################################################################

@doc Markdown.doc"""
    FractionField(R::AbstractAlgebra.Ring; cached=true)
> Return the parent object of the fraction field over the given base ring $R$.
> If `cached == true` (the default), the returned parent object is cached so
> that it will always be returned by a call to the constructor when the same
> base ring $R$ is supplied.
"""
function FractionField(R::AbstractAlgebra.Ring; cached=true)
   R2 = R
   T = elem_type(R)

   return FracField{T}(R, cached)
end
