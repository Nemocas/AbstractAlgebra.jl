###############################################################################
#
#   Fraction.jl : generic fraction fields
#
###############################################################################

export FractionField, num, den

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{Frac{T}}) where T <: RingElem = FracField{T}

elem_type(::Type{FracField{T}}) where {T <: RingElem} = Frac{T}

@doc Markdown.doc"""
    base_ring{T}(S::AbstractAlgebra.FracField{T})
> Return the base ring $R$ of the given fraction field.
"""
base_ring(a::AbstractAlgebra.FracField{T}) where T <: RingElem = a.base_ring::parent_type(T)

@doc Markdown.doc"""
    base_ring{T}(r::AbstractAlgebra.FracElem)
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
    characteristic{T <: RingElem}(R::AbstractAlgebra.FracField{T})
> Return the characteristic of the given field.
"""
function characteristic(R::AbstractAlgebra.FracField{T}) where T <: RingElem
   return characteristic(base_ring(R))
end

function check_parent(a::AbstractAlgebra.FracElem, b::AbstractAlgebra.FracElem)
   parent(a) != parent(b) && error("Incompatible rings in fraction field operation")
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
   return xor(b, hash(numerator(a), h), hash(denominator(a), h), h)
end

function numerator(a::AbstractAlgebra.FracElem)
   u = canonical_unit(a.den)
   return divexact(a.num, u)
end

function denominator(a::AbstractAlgebra.FracElem)
   u = canonical_unit(a.den)
   return divexact(a.den, u)
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
iszero(a::AbstractAlgebra.FracElem) = iszero(numerator(a))

@doc Markdown.doc"""
    isone(a::AbstractAlgebra.FracElem)
> Return `true` if the supplied element $a$ is one in the fraction field it
> belongs to, otherwise return `false`.
"""
isone(a::AbstractAlgebra.FracElem) = numerator(a) == denominator(a)

@doc Markdown.doc"""
    isunit(a::AbstractAlgebra.FracElem)
> Return `true` if the supplied element $a$ is invertible in the fraction field
> it belongs to, i.e. the numerator is nonzero, otherwise return `false`.
"""
isunit(a::AbstractAlgebra.FracElem) = !iszero(numerator(a))

function deepcopy_internal(a::Frac{T}, dict::IdDict) where {T <: RingElem}
   v = Frac{T}(deepcopy(numerator(a)), deepcopy(denominator(a)))
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
   u = canonical_unit(denominator(x))
   n = divexact(numerator(x), u)
   d = divexact(denominator(x), u);
   if d != 1 && needs_parentheses(n)
      print(io, "(")
   end
   print(io, n)
   if d != 1
      if needs_parentheses(n)
         print(io, ")")
      end
      print(io, "//")
      if needs_parentheses(d)
         print(io, "(")
      end
      print(io, d)
      if needs_parentheses(d)
         print(io, ")")
      end
   end
end

function show(io::IO, a::AbstractAlgebra.FracField)
   print(io, "Fraction field of ", base_ring(a))
end

needs_parentheses(x::AbstractAlgebra.FracElem) = isone(denominator(x)) && needs_parentheses(numerator(x))

displayed_with_minus_in_front(x::AbstractAlgebra.FracElem) = !needs_parentheses(numerator(x)) && displayed_with_minus_in_front(numerator(x))

show_minus_one(::Type{AbstractAlgebra.FracElem{T}}) where {T <: RingElem} = show_minus_one(T)

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
   return parent(a)(-numerator(a), denominator(a))
end

###############################################################################
#
#   Binary operators
#
###############################################################################

@doc Markdown.doc"""
    +{T <: RingElem}(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T})
> Return $a + b$.
"""
function +(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   d1 = denominator(a)
   d2 = denominator(b)
   n1 = numerator(a)
   n2 = numerator(b)
   gd = gcd(d1, d2)
   if d1 == d2
      rnum = n1 + n2
      if isone(d1)
         rden = d1
      else
         gd = gcd(rnum, d1)
         if isone(gd)
            rden = d1
         else
            rnum = divexact(rnum, gd)
            rden = divexact(d1, gd)
         end
      end
   elseif isone(d1)
      rnum = n1*d2 + n2
      rden = d2
   elseif isone(d2)
      rnum = n1 + n2*d1
      rden = d1
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
    -{T <: RingElem}(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T})
> Return $a - b$.
"""
function -(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   d1 = denominator(a)
   d2 = denominator(b)
   n1 = numerator(a)
   n2 = numerator(b)
   if d1 == d2
      rnum = n1 - n2
      if isone(d1)
         rden = d1
      else
         gd = gcd(rnum, d1)
         if isone(gd)
            rden = d1
         else
            rnum = divexact(rnum, gd)
            rden = divexact(d1, gd)
         end
      end
   elseif isone(d1)
      rnum = n1*d2 - n2
      rden = d2
   elseif isone(d2)
      rnum = n1 - n2*d1
      rden = d1
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
    *{T <: RingElem}(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T})
> Return $a\times b$.
"""
function *(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   n1 = numerator(a)
   d2 = denominator(b)
   n2 = numerator(b)
   d1 = denominator(a)
   if d1 == d2
      n = n1*n2
      d = d1*d2
   elseif isone(d1)
      gd = gcd(n1, d2)
      if isone(gd)
         n = n1*n2
         d = d2
      else
         n = divexact(n1, gd)*n2
         d = divexact(d2, gd)
      end
   elseif isone(d2)
      gd = gcd(n2, d1)
      if isone(gd)
         n = n2*n1
         d = d1
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
   g = gcd(denominator(a), c)
   n = numerator(a)*divexact(c, g)
   d = divexact(denominator(a), g)
   return parent(a)(n, d)
end

@doc Markdown.doc"""
    *(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.FracElem)
> Return $a\times b$.
"""
function *(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.FracElem)
   c = base_ring(b)(a)
   g = gcd(denominator(b), c)
   n = numerator(b)*divexact(c, g)
   d = divexact(denominator(b), g)
   return parent(b)(n, d)
end

@doc Markdown.doc"""
    *{T <: RingElem}(a::AbstractAlgebra.FracElem{T}, b::T)
> Return $a\times b$.
"""
function *(a::AbstractAlgebra.FracElem{T}, b::T) where {T <: RingElem}
   g = gcd(denominator(a), b)
   n = numerator(a)*divexact(b, g)
   d = divexact(denominator(a), g)
   return parent(a)(n, d)
end

@doc Markdown.doc"""
    *{T <: RingElem}(a::T, b::AbstractAlgebra.FracElem{T})
> Return $a\times b$.
"""
function *(a::T, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   g = gcd(denominator(b), a)
   n = numerator(b)*divexact(a, g)
   d = divexact(denominator(b), g)
   return parent(b)(n, d)
end

@doc Markdown.doc"""
    +(a::AbstractAlgebra.FracElem, b::Union{Integer, Rational, AbstractFloat})
> Return $a + b$.
"""
function +(a::AbstractAlgebra.FracElem, b::Union{Integer, Rational, AbstractFloat})
   n = numerator(a) + denominator(a)*b
   d = denominator(a)
   return parent(a)(n, d)
end

@doc Markdown.doc"""
    -(a::AbstractAlgebra.FracElem, b::Union{Integer, Rational, AbstractFloat})
> Return $a - b$.
"""
function -(a::AbstractAlgebra.FracElem, b::Union{Integer, Rational, AbstractFloat})
   n = numerator(a) - denominator(a)*b
   d = denominator(a)
   return parent(a)(n, d)
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
   n = a*denominator(b) - numerator(b)
   d = denominator(b)
   return parent(b)(n, d)
end

@doc Markdown.doc"""
    +{T <: RingElem}(a::AbstractAlgebra.FracElem{T}, b::T)
> Return $a + b$.
"""
function +(a::AbstractAlgebra.FracElem{T}, b::T) where {T <: RingElem}
   n = numerator(a) + denominator(a)*b
   d = denominator(a)
   return parent(a)(n, d)
end

@doc Markdown.doc"""
    -{T <: RingElem}(a::AbstractAlgebra.FracElem{T}, b::T)
> Return $a - b$.
"""
function -(a::AbstractAlgebra.FracElem{T}, b::T) where {T <: RingElem}
   n = numerator(a) - denominator(a)*b
   d = denominator(a)
   return parent(a)(n, d)
end

@doc Markdown.doc"""
    +{T <: RingElem}(a::T, b::AbstractAlgebra.FracElem{T})
> Return $a + b$.
"""
+(a::T, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem} = b + a

@doc Markdown.doc"""
    -{T <: RingElem}(a::T, b::AbstractAlgebra.FracElem{T})
> Return $a - b$.
"""
function -(a::T, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   n = a*denominator(b) - numerator(b)
   d = denominator(b)
   return parent(b)(n, d)
end

###############################################################################
#
#   Comparisons
#
###############################################################################

@doc Markdown.doc"""
    =={T <: RingElem}(x::AbstractAlgebra.FracElem{T}, y::AbstractAlgebra.FracElem{T})
> Return `true` if $x == y$ arithmetically, otherwise return `false`. Recall
> that power series to different precisions may still be arithmetically
> equal to the minimum of the two precisions.
"""
function ==(x::AbstractAlgebra.FracElem{T}, y::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   check_parent(x, y)
   return (denominator(x) == denominator(y) && numerator(x) == numerator(y)) || (numerator(x)*denominator(y) == denominator(x)*numerator(y))
end

@doc Markdown.doc"""
    isequal{T <: RingElem}(x::AbstractAlgebra.FracElem{T}, y::AbstractAlgebra.FracElem{T})
> Return `true` if $x == y$ exactly, otherwise return `false`. This function is
> useful in cases where the numerators and denominators of the fractions are
> inexact, e.g. power series. Only if the power series are precisely the same,
> to the same precision, are they declared equal by this function.
"""
function isequal(x::AbstractAlgebra.FracElem{T}, y::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   if parent(x) != parent(y)
      return false
   end
   return isequal(numerator(x)*denominator(y), denominator(x)*numerator(y))
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
   return (isone(denominator(x)) && numerator(x) == y) || (numerator(x) == denominator(x)*y)
end

@doc Markdown.doc"""
    ==(x::Union{Integer, Rational, AbstractFloat}, y::AbstractAlgebra.FracElem)
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::Union{Integer, Rational, AbstractFloat}, y::AbstractAlgebra.FracElem) = y == x

@doc Markdown.doc"""
    =={T <: RingElem}(x::AbstractAlgebra.FracElem{T}, y::T)
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
function ==(x::AbstractAlgebra.FracElem{T}, y::T) where {T <: RingElem}
   return (isone(denominator(x)) && numerator(x) == y) || (numerator(x) == denominator(x)*y)
end

@doc Markdown.doc"""
    =={T <: RingElem}(x::T, y::AbstractAlgebra.FracElem{T})
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
   iszero(numerator(a)) && throw(DivideError())
   return parent(a)(denominator(a), numerator(a))
end

###############################################################################
#
#   Exact division
#
###############################################################################

@doc Markdown.doc"""
    divexact{T <: RingElem}(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T})
> Return $a/b$.
"""
function divexact(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   n1 = numerator(a)
   d2 = denominator(b)
   n2 = numerator(b)
   d1 = denominator(a)
   if d1 == n2
      n = n1*d2
      d = d1*n2
   elseif isone(d1)
      gd = gcd(n1, n2)
      if isone(gd)
         n = n1*d2
         d = n2
      else
         n = divexact(n1, gd)*d2
         d = divexact(n2, gd)
      end
   elseif isone(n2)
      gd = gcd(d2, d1)
      if isone(gd)
         n = d2*n1
         d = d1
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
   g = gcd(numerator(a), c)
   n = divexact(numerator(a), g)
   d = denominator(a)*divexact(c, g)
   return parent(a)(n, d)
end

@doc Markdown.doc"""
    divexact(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.FracElem)
> Return $a/b$.
"""
function divexact(a::Union{Integer, Rational, AbstractFloat}, b::AbstractAlgebra.FracElem)
   iszero(b) && throw(DivideError())
   c = base_ring(b)(a)
   g = gcd(numerator(b), c)
   n = denominator(b)*divexact(c, g)
   d = divexact(numerator(b), g)
   return parent(b)(n, d)
end

@doc Markdown.doc"""
    divexact{T <: RingElem}(a::AbstractAlgebra.FracElem{T}, b::T)
> Return $a/b$.
"""
function divexact(a::AbstractAlgebra.FracElem{T}, b::T) where {T <: RingElem}
   iszero(b) && throw(DivideError())
   g = gcd(numerator(a), b)
   n = divexact(numerator(a), g)
   d = denominator(a)*divexact(b, g)
   return parent(a)(n, d)
end

@doc Markdown.doc"""
    divexact{T <: RingElem}(a::T, b::AbstractAlgebra.FracElem{T})
> Return $a/b$.
"""
function divexact(a::T, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   iszero(b) && throw(DivideError())
   g = gcd(numerator(b), a)
   n = denominator(b)*divexact(a, g)
   d = divexact(numerator(b), g)
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
   return parent(a)(numerator(a)^b, denominator(a)^b)
end

###############################################################################
#
#   GCD
#
###############################################################################

@doc Markdown.doc"""
    gcd{T <: RingElem}(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T})
> Return a greatest common divisor of $a$ and $b$ if one exists. N.B: we define
> the GCD of $a/b$ and $c/d$ to be gcd$(ad, bc)/bd$, reduced to lowest terms.
> This requires the existence of a greatest common divisor function for the
> base ring.
"""
function gcd(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   n = gcd(numerator(a)*denominator(b), denominator(a)*numerator(b))
   d = denominator(a)*denominator(b)
   g = gcd(n, d)
   return parent(a)(divexact(n, g), divexact(d, g))
end

################################################################################
#
#   Remove and valuation
#
################################################################################

@doc Markdown.doc"""
    remove{T <: RingElem}(z::AbstractAlgebra.FracElem{T}, p::T)
> Return the tuple $n, x$ such that $z = p^nx$ where $x$ has valuation $0$ at
> $p$.
"""
function remove(z::AbstractAlgebra.FracElem{T}, p::T) where {T <: RingElem}
   iszero(z) && error("Not yet implemented")
   v, d = remove(denominator(z), p)
   w, n = remove(numerator(z), p)
   return w-v, n//d
end

@doc Markdown.doc"""
    valuation{T <: RingElem}(z::AbstractAlgebra.FracElem{T}, p::T)
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
   g1 = gcd(numerator(a), denominator(b))
   g2 = gcd(numerator(b), denominator(a))
   c.num = divexact(numerator(a), g1)*divexact(numerator(b), g2)
   c.den = divexact(denominator(a), g2)*divexact(denominator(b), g1)
   return c
end

function addeq!(c::AbstractAlgebra.FracElem{T}, a::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   n = c.num*denominator(a) + numerator(a)*c.den
   c.den = mul!(c.den, c.den, denominator(a))
   g = gcd(n, c.den)
   c.num = divexact(n, g)
   c.den = divexact(c.den, g)
   return c
end

function add!(c::AbstractAlgebra.FracElem{T}, a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   n = c.num*denominator(a) + numerator(a)*c.den
   d = c.den*denominator(a)
   g = gcd(n, d)
   c.num = divexact(n, g)
   c.den = divexact(d, g)
   return c
end

function addeq!(c::AbstractAlgebra.FracElem{T}, a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}
   n = numerator(b)*denominator(a) + numerator(a)*denominator(b)
   c.den = mul!(c.den, denominator(b), denominator(a))
   g = gcd(n, d)
   c.num = divexact(n, g)
   c.den = divexact(c.den, g)
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
   while d == 0
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
