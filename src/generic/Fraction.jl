###############################################################################
#
#   Fraction.jl : generic fraction fields
#
###############################################################################

export FractionField, GenFrac, GenFracField, num, den

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{GenFrac{T}}) where T <: RingElem = GenFracField{T}

elem_type(::Type{GenFracField{T}}) where {T <: RingElem} = GenFrac{T}

doc"""
    base_ring{T}(S::FracField{T})
> Return the base ring $R$ of the given fraction field.
"""
base_ring(a::FracField{T}) where T <: RingElem = a.base_ring::parent_type(T)

doc"""
    base_ring{T}(r::FracElem)
> Return the base ring $R$ of the fraction field that the supplied
> element $a$ belongs to.
"""
base_ring(a::FracElem) = base_ring(parent(a))

doc"""
    parent(a::FracElem)
> Return the parent object of the given fraction element.
"""
parent(a::FracElem) = a.parent

function check_parent(a::FracElem, b::FracElem)
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
   z = GenFrac{T}(divexact(x, g), divexact(y, g))
   try
      z.parent = GenFracDict[R]
   catch
      z.parent = FractionField(parent(x))
   end
   return z
end

//(x::T, y::Integer) where {T <: RingElem} = x//parent(x)(y)
                                          
//(x::Integer, y::T) where {T <: RingElem} = parent(y)(x)//y

//(x::T, y::FracElem{T}) where {T <: RingElem} = parent(y)(x)//y

//(x::FracElem{T}, y::T) where {T <: RingElem} = x//parent(x)(y)

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::FracElem, h::UInt)
   b = 0x8a30b0d963237dd5%UInt
   return xor(b, hash(num(a), h), hash(den(a), h), h)
end

function num(a::FracElem)
   u = canonical_unit(a.den)
   return divexact(a.num, u)
end

function den(a::FracElem)
   u = canonical_unit(a.den)
   return divexact(a.den, u)
end

doc"""
    zero(R::FracField)
> Return $0/1$ in the given fraction field.
"""
zero(R::FracField) = R(0)

doc"""
    one(R::FracField)
> Return $1/1$ in the given fraction field.
"""
one(R::FracField) = R(1)

doc"""
    iszero(a::FracElem)
> Return `true` if the supplied element $a$ is zero in the fraction field it
> belongs to, otherwise return `false`.
"""
iszero(a::FracElem) = iszero(num(a))

doc"""
    isone(a::FracElem)
> Return `true` if the supplied element $a$ is one in the fraction field it
> belongs to, otherwise return `false`.
"""
isone(a::FracElem) = num(a) == den(a)

doc"""
    isunit(a::FracElem)
> Return `true` if the supplied element $a$ is invertible in the fraction field
> it belongs to, i.e. the numerator is nonzero, otherwise return `false`.
"""
isunit(a::FracElem) = !iszero(num(a))

function deepcopy_internal(a::GenFrac{T}, dict::ObjectIdDict) where {T <: RingElem}
   v = GenFrac{T}(deepcopy(num(a)), deepcopy(den(a)))
   v.parent = parent(a)
   return v
end 

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(a::FracElem) = a

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function show(io::IO, x::FracElem)
   u = canonical_unit(den(x))
   n = divexact(num(x), u)
   d = divexact(den(x), u);
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

function show(io::IO, a::FracField)
   print(io, "Fraction field of ", base_ring(a))
end

needs_parentheses(x::FracElem) = isone(den(x)) && needs_parentheses(num(x))

isnegative(x::FracElem) = !needs_parentheses(num(x)) && isnegative(num(x))

show_minus_one(::Type{FracElem{T}}) where {T <: RingElem} = show_minus_one(T)

###############################################################################
#
#   Unary operators
#
###############################################################################

doc"""
    -(a::FracElem)
> Return $-a$.
"""
function -(a::FracElem)
   return parent(a)(-num(a), den(a))
end

###############################################################################
#
#   Binary operators
#
###############################################################################

doc"""
    +{T <: RingElem}(a::FracElem{T}, b::FracElem{T})
> Return $a + b$.
"""
function +(a::FracElem{T}, b::FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   n = num(a)*den(b) + num(b)*den(a)
   d = den(a)*den(b)
   g = gcd(n, d)
   return parent(a)(divexact(n, g), divexact(d, g))
end

doc"""
    -{T <: RingElem}(a::FracElem{T}, b::FracElem{T})
> Return $a - b$.
"""
function -(a::FracElem{T}, b::FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   n = num(a)*den(b) - num(b)*den(a)
   d = den(a)*den(b)
   g = gcd(n, d)
   return parent(a)(divexact(n, g), divexact(d, g))
end

doc"""
    *{T <: RingElem}(a::FracElem{T}, b::FracElem{T})
> Return $a\times b$.
"""
function *(a::FracElem{T}, b::FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   g1 = gcd(num(a), den(b))
   g2 = gcd(num(b), den(a))
   n = divexact(num(a), g1)*divexact(num(b), g2)
   d = divexact(den(a), g2)*divexact(den(b), g1)
   return parent(a)(n, d)
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

doc"""
    *(a::FracElem, b::Integer)
> Return $a\times b$.
"""
function *(a::FracElem, b::Integer)
   c = base_ring(a)(b)
   g = gcd(den(a), c)
   n = num(a)*divexact(c, g)
   d = divexact(den(a), g)
   return parent(a)(n, d)
end

doc"""
    *(a::Integer, b::FracElem)
> Return $a\times b$.
"""
function *(a::Integer, b::FracElem)
   c = base_ring(b)(a)
   g = gcd(den(b), c)
   n = num(b)*divexact(c, g)
   d = divexact(den(b), g)
   return parent(b)(n, d)
end

doc"""
    *(a::FracElem, b::fmpz)
> Return $a\times b$.
"""
function *(a::FracElem, b::fmpz)
   c = base_ring(a)(b)
   g = gcd(den(a), c)
   n = num(a)*divexact(c, g)
   d = divexact(den(a), g)
   return parent(a)(n, d)
end

doc"""
    *(a::fmpz, b::FracElem)
> Return $a\times b$.
"""
function *(a::fmpz, b::FracElem)
   c = base_ring(b)(a)
   g = gcd(den(b), c)
   n = num(b)*divexact(c, g)
   d = divexact(den(b), g)
   return parent(b)(n, d)
end

doc"""
    *{T <: RingElem}(a::FracElem{T}, b::T)
> Return $a\times b$.
"""
function *(a::FracElem{T}, b::T) where {T <: RingElem}
   g = gcd(den(a), b)
   n = num(a)*divexact(b, g)
   d = divexact(den(a), g)
   return parent(a)(n, d)
end

doc"""
    *{T <: RingElem}(a::T, b::FracElem{T})
> Return $a\times b$.
"""
function *(a::T, b::FracElem{T}) where {T <: RingElem}
   g = gcd(den(b), a)
   n = num(b)*divexact(a, g)
   d = divexact(den(b), g)
   return parent(b)(n, d)
end

doc"""
    +(a::FracElem, b::Integer)
> Return $a + b$.
"""
function +(a::FracElem, b::Integer)
   n = num(a) + den(a)*b
   d = den(a)
   g = gcd(n, d)
   return parent(a)(divexact(n, g), divexact(d, g))
end

doc"""
    +(a::FracElem, b::fmpz)
> Return $a + b$.
"""
function +(a::FracElem, b::fmpz)
   n = num(a) + den(a)*b
   d = den(a)
   g = gcd(n, d)
   return parent(a)(divexact(n, g), divexact(d, g))
end

doc"""
    -(a::FracElem, b::Integer)
> Return $a - b$.
"""
function -(a::FracElem, b::Integer)
   n = num(a) - den(a)*b
   d = den(a)
   g = gcd(n, d)
   return parent(a)(divexact(n, g), divexact(d, g))
end

doc"""
    -(a::FracElem, b::fmpz)
> Return $a - b$.
"""
function -(a::FracElem, b::fmpz)
   n = num(a) - den(a)*b
   d = den(a)
   g = gcd(n, d)
   return parent(a)(divexact(n, g), divexact(d, g))
end

doc"""
    +(a::Integer, b::FracElem)
> Return $a + b$.
"""
+(a::Integer, b::FracElem) = b + a

doc"""
    +(a::fmpz, b::FracElem)
> Return $a + b$.
"""
+(a::fmpz, b::FracElem) = b + a

doc"""
    -(a::Integer, b::FracElem)
> Return $a - b$.
"""
function -(a::Integer, b::FracElem)
   n = a*den(b) - num(b)
   d = den(b)
   g = gcd(n, d)
   return parent(b)(divexact(n, g), divexact(d, g))
end

doc"""
    -(a::fmpz, b::FracElem)
> Return $a - b$.
"""
function -(a::fmpz, b::FracElem)
   n = a*den(b) - num(b)
   d = den(b)
   g = gcd(n, d)
   return parent(b)(divexact(n, g), divexact(d, g))
end

doc"""
    +{T <: RingElem}(a::FracElem{T}, b::T)
> Return $a + b$.
"""
function +(a::FracElem{T}, b::T) where {T <: RingElem}
   n = num(a) + den(a)*b
   d = den(a)
   g = gcd(n, d)
   return parent(a)(divexact(n, g), divexact(d, g))
end

doc"""
    -{T <: RingElem}(a::FracElem{T}, b::T)
> Return $a - b$.
"""
function -(a::FracElem{T}, b::T) where {T <: RingElem}
   n = num(a) - den(a)*b
   d = den(a)
   g = gcd(n, d)
   return parent(a)(divexact(n, g), divexact(d, g))
end

doc"""
    +{T <: RingElem}(a::T, b::FracElem{T})
> Return $a + b$.
"""
+(a::T, b::FracElem{T}) where {T <: RingElem} = b + a

doc"""
    -{T <: RingElem}(a::T, b::FracElem{T})
> Return $a - b$.
"""
function -(a::T, b::FracElem{T}) where {T <: RingElem}
   n = a*den(b) - num(b)
   d = den(b)
   g = gcd(n, d)
   return parent(b)(divexact(n, g), divexact(d, g))
end

###############################################################################
#
#   Comparisons
#
###############################################################################

doc"""
    =={T <: RingElem}(x::FracElem{T}, y::FracElem{T})
> Return `true` if $x == y$ arithmetically, otherwise return `false`. Recall
> that power series to different precisions may still be arithmetically
> equal to the minimum of the two precisions.
"""
function ==(x::FracElem{T}, y::FracElem{T}) where {T <: RingElem}
   check_parent(x, y)
   return (den(x) == den(y) && num(x) == num(y)) || (num(x)*den(y) == den(x)*num(y))
end

doc"""
    isequal{T <: RingElem}(x::FracElem{T}, y::FracElem{T})
> Return `true` if $x == y$ exactly, otherwise return `false`. This function is
> useful in cases where the numerators and denominators of the fractions are
> inexact, e.g. power series. Only if the power series are precisely the same,
> to the same precision, are they declared equal by this function.
"""
function isequal(x::FracElem{T}, y::FracElem{T}) where {T <: RingElem}
   if parent(x) != parent(y)
      return false
   end
   return isequal(num(x)*den(y), den(x)*num(y))
end

###############################################################################
#
#   Ad hoc comparisons
#
###############################################################################

doc"""
    ==(x::FracElem, y::Integer)
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
function ==(x::FracElem, y::Integer)
   return (isone(den(x)) && num(x) == y) || (num(x) == den(x)*y)
end

doc"""
    ==(x::Integer, y::FracElem)
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::Integer, y::FracElem) = y == x

doc"""
    ==(x::FracElem, y::fmpz)
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
function ==(x::FracElem, y::fmpz)
   return (isone(den(x)) && num(x) == y) || (num(x) == den(x)*y)
end

doc"""
    ==(x::fmpz, y::FracElem)
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::fmpz, y::FracElem) = y == x

doc"""
    =={T <: RingElem}(x::FracElem{T}, y::T)
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
function ==(x::FracElem{T}, y::T) where {T <: RingElem}
   return (isone(den(x)) && num(x) == y) || (num(x) == den(x)*y)
end

doc"""
    =={T <: RingElem}(x::T, y::FracElem{T})
> Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::T, y::FracElem{T}) where {T <: RingElem} = y == x

###############################################################################
#
#   Inversion
#
###############################################################################

doc"""
    inv(a::FracElem)
> Return the inverse of the fraction $a$.
"""
function inv(a::FracElem)
   iszero(num(a)) && throw(DivideError())
   return parent(a)(den(a), num(a))
end

###############################################################################
#
#   Exact division
#
###############################################################################

doc"""
    divexact{T <: RingElem}(a::FracElem{T}, b::FracElem{T})
> Return $a/b$.
"""
function divexact(a::FracElem{T}, b::FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   g1 = gcd(num(a), num(b))
   g2 = gcd(den(b), den(a))
   n = divexact(num(a), g1)*divexact(den(b), g2)
   d = divexact(den(a), g2)*divexact(num(b), g1)
   return parent(a)(n, d)
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

doc"""
    divexact(a::FracElem, b::Integer)
> Return $a/b$.
"""
function divexact(a::FracElem, b::Integer)
   b == 0 && throw(DivideError())
   c = base_ring(a)(b)
   g = gcd(num(a), c)
   n = divexact(num(a), g)
   d = den(a)*divexact(c, g)
   return parent(a)(n, d)
end

doc"""
    divexact(a::Integer, b::FracElem)
> Return $a/b$.
"""
function divexact(a::Integer, b::FracElem)
   iszero(b) && throw(DivideError())
   c = base_ring(b)(a)
   g = gcd(num(b), c)
   n = den(b)*divexact(c, g)
   d = divexact(num(b), g)
   return parent(b)(n, d)
end

doc"""
    divexact(a::FracElem, b::fmpz)
> Return $a/b$.
"""
function divexact(a::FracElem, b::fmpz)
   iszero(b) && throw(DivideError())
   c = base_ring(a)(b)
   g = gcd(num(a), c)
   n = divexact(num(a), g)
   d = den(a)*divexact(c, g)
   return parent(a)(n, d)
end

doc"""
    divexact(a::fmpz, b::FracElem)
> Return $a/b$.
"""
function divexact(a::fmpz, b::FracElem)
   iszero(b) && throw(DivideError())
   c = base_ring(b)(a)
   g = gcd(num(b), c)
   n = den(b)*divexact(c, g)
   d = divexact(num(b), g)
   return parent(b)(n, d)
end

doc"""
    divexact{T <: RingElem}(a::FracElem{T}, b::T)
> Return $a/b$.
"""
function divexact(a::FracElem{T}, b::T) where {T <: RingElem}
   iszero(b) && throw(DivideError())
   g = gcd(num(a), b)
   n = divexact(num(a), g)
   d = den(a)*divexact(b, g)
   return parent(a)(n, d)
end

doc"""
    divexact{T <: RingElem}(a::T, b::FracElem{T})
> Return $a/b$.
"""
function divexact(a::T, b::FracElem{T}) where {T <: RingElem}
   iszero(b) && throw(DivideError())
   g = gcd(num(b), a)
   n = den(b)*divexact(a, g)
   d = divexact(num(b), g)
   return parent(b)(n, d)
end

function divides(a::FracElem{T}, b::FracElem{T}) where {T <: RingElem}
   iszero(b) && error("Division by zero in divides")
   return true, divexact(a, b)
end

###############################################################################
#
#   Powering
#
###############################################################################

doc"""
    ^(a::FracElem, b::Int)
> Return $a^b$.
"""
function ^(a::FracElem{T}, b::Int) where {T <: RingElem}
   if b < 0
      a = inv(a)
      b = -b
   end
   return parent(a)(num(a)^b, den(a)^b)
end

###############################################################################
#
#   GCD
#
###############################################################################

doc"""
    gcd{T <: RingElem}(a::FracElem{T}, b::FracElem{T})
> Return a greatest common divisor of $a$ and $b$ if one exists. N.B: we define
> the GCD of $a/b$ and $c/d$ to be gcd$(ad, bc)/bd$, reduced to lowest terms.
> This requires the existence of a greatest common divisor function for the
> base ring.
"""
function gcd(a::FracElem{T}, b::FracElem{T}) where {T <: RingElem}
   check_parent(a, b)
   n = gcd(num(a)*den(b), den(a)*num(b))
   d = den(a)*den(b)
   g = gcd(n, d)
   return parent(a)(divexact(n, g), divexact(d, g))
end

################################################################################
#
#   Remove and valuation
#
################################################################################

doc"""
    remove{T <: RingElem}(z::FracElem{T}, p::T)
> Return the tuple $n, x$ such that $z = p^nx$ where $x$ has valuation $0$ at
> $p$.
"""
function remove(z::FracElem{T}, p::T) where {T <: RingElem}
   iszero(z) && error("Not yet implemented")
   v, d = remove(den(z), p)
   w, n = remove(num(z), p)
   return w-v, n//d
end 

doc"""
    valuation{T <: RingElem}(z::FracElem{T}, p::T)
> Return the valuation of $z$ at $p$.
"""
function valuation(z::FracElem{T}, p::T) where {T <: RingElem}
   v, _ = remove(z, p)
   return v
end
  
###############################################################################
#
#   Unsafe operators and functions
#
###############################################################################

function zero!(c::FracElem)
   c.num = zero!(c.num)
   if !isone(c.den)
      c.den = one(parent(c))
   end
   return c
end

function mul!(c::FracElem{T}, a::FracElem{T}, b::FracElem{T}) where {T <: RingElem}
   g1 = gcd(num(a), den(b))
   g2 = gcd(num(b), den(a))
   c.num = divexact(num(a), g1)*divexact(num(b), g2)
   c.den = divexact(den(a), g2)*divexact(den(b), g1)
   return c
end

function addeq!(c::FracElem{T}, a::FracElem{T}) where {T <: RingElem}
   n = c.num*den(a) + num(a)*c.den
   c.den = mul!(c.den, c.den, den(a))
   g = gcd(n, c.den)
   c.num = divexact(n, g)
   c.den = divexact(c.den, g)
   return c
end

function add!(c::FracElem{T}, a::FracElem{T}, b::FracElem{T}) where {T <: RingElem}
   n = c.num*den(a) + num(a)*c.den
   d = c.den*den(a)
   g = gcd(n, d)
   c.num = divexact(n, g)
   c.den = divexact(d, g)
   return c
end

function addeq!(c::FracElem{T}, a::FracElem{T}, b::FracElem{T}) where {T <: RingElem}
   n = num(b)*den(a) + num(a)*den(b)
   c.den = mul!(c.den, den(b), den(a))
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
   
function rand(S::FracField{T}, v...) where {T <: RingElem}
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

promote_rule(::Type{GenFrac{T}}, ::Type{T}) where {T <: RingElem} = GenFrac{T}

promote_rule(::Type{GenFrac{T}}, ::Type{U}) where {T <: RingElem, U <: Integer} = GenFrac{T}

function promote_rule1(::Type{GenFrac{T}}, ::Type{GenFrac{U}}) where {T <: RingElem, U <: RingElem}
   promote_rule(T, GenFrac{U}) == T ? GenFrac{T} : Union{}
end

function promote_rule(::Type{GenFrac{T}}, ::Type{U}) where {T <: RingElem, U <: RingElem}
   promote_rule(T, U) == T ? GenFrac{T} : promote_rule1(U, GenFrac{T})
end

###############################################################################
#
#   Parent object call overloading
#
###############################################################################

function (a::GenFracField{T})(b::RingElem) where {T <: RingElem}
   return a(base_ring(a)(b))
end

function (a::GenFracField{T})() where {T <: RingElem}
   z = GenFrac{T}(zero(base_ring(a)), one(base_ring(a)))
   z.parent = a
   return z
end

function (a::GenFracField{T})(b::fmpz) where {T <: RingElem}
   z = GenFrac{T}(base_ring(a)(b), one(base_ring(a)))
   z.parent = a
   return z
end

function (a::GenFracField{T})(b::T) where {T <: RingElem}
   parent(b) != base_ring(a) && error("Could not coerce to fraction")
   z = GenFrac{T}(b, one(base_ring(a)))
   z.parent = a
   return z
end

function (a::GenFracField{T})(b::T, c::T) where {T <: RingElem}
   parent(b) != base_ring(a) && error("Could not coerce to fraction")
   parent(c) != base_ring(a) && error("Could not coerce to fraction")
   z = GenFrac{T}(b, c)
   z.parent = a
   return z
end

function (a::GenFracField{T})(b::T, c::Integer) where {T <: RingElem}
   parent(b) != base_ring(a) && error("Could not coerce to fraction")
   z = GenFrac{T}(b, base_ring(a)(c))
   z.parent = a
   return z
end

function (a::GenFracField{T})(b::Integer, c::T) where {T <: RingElem}
   parent(c) != base_ring(a) && error("Could not coerce to fraction")
   z = GenFrac{T}(base_ring(a)(b), c)
   z.parent = a
   return z
end

function (a::GenFracField{T})(b::Integer) where {T <: RingElem}
   z = GenFrac{T}(base_ring(a)(b), one(base_ring(a)))
   z.parent = a
   return z
end

function (a::GenFracField{T})(b::Integer, c::Integer) where {T <: RingElem}
   z = GenFrac{T}(base_ring(a)(b), base_ring(a)(c))
   z.parent = a
   return z
end

function (a::GenFracField{T})(b::GenFrac{T}) where {T <: RingElem}
   a != parent(b) && error("Could not coerce to fraction")
   return b
end

###############################################################################
#
#   FractionField constructor
#
###############################################################################

doc"""
    FractionField(R::Ring; cached=true)
> Return the parent object of the fraction field over the given base ring $R$.
> If `cached == true` (the default), the returned parent object is cached so
> that it will always be returned by a call to the constructor when the same
> base ring $R$ is supplied.
"""
function FractionField(R::Ring; cached=true)
   R2 = R
   T = elem_type(R)
   
   return GenFracField{T}(R, cached)
end

