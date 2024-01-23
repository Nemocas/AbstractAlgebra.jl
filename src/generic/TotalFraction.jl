###############################################################################
#
#   Fraction.jl : fraction fields
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

parent_type(::Type{TotFrac{T}}) where T <: RingElem = TotFracRing{T}

elem_type(::Type{TotFracRing{T}}) where {T <: RingElem} = TotFrac{T}

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

base_ring_type(::Type{TotFracRing{T}}) where T <: RingElem = parent_type(T)

base_ring(a::TotFracRing{T}) where T <: RingElem = a.base_ring::parent_type(T)

parent(a::TotFrac) = a.parent

function is_domain_type(::Type{T}) where {S <: RingElement, T <: TotFrac{S}}
   return is_domain_type(S)
end

function is_exact_type(a::Type{T}) where {S <: RingElement, T <: TotFrac{S}}
   return is_exact_type(S)
end

function characteristic(R::TotFracRing{T}) where T <: RingElem
   return characteristic(base_ring(R))
end

function check_parent(a::TotFrac, b::TotFrac, throw::Bool = true)
   fl = parent(a) != parent(b)
   fl && throw && error("Incompatible rings in total ring of fractions operation")
   return !fl
end

###############################################################################
#
#   Constructors
#
###############################################################################

function //(x::T, y::TotFrac{T}) where {T <: RingElem}
   parent(y)(x*denominator(y, false), numerator(y, false))
end

function //(x::TotFrac{T}, y::T) where {T <: RingElem}
   parent(x)(numerator(x, false), denominator(x, false)*y)
end

function //(x::Integer, y::TotFrac{T}) where {T <: RingElem}
   parent(y)(x*denominator(y, false), numerator(y, false))
end

function //(x::TotFrac{T}, y::Integer) where {T <: RingElem}
   parent(x)(numerator(x, false), denominator(x, false)*y)
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.numerator(a::TotFrac, canonicalise::Bool=true)
   return a.num
end

function Base.denominator(a::TotFrac, canonicalise::Bool=true)
   return a.den
end

function deepcopy_internal(a::TotFrac{T}, dict::IdDict) where {T <: RingElem}
   v = TotFrac{T}(deepcopy_internal(numerator(a, false), dict),
                  deepcopy_internal(denominator(a, false), dict))
   v.parent = parent(a)
   return v
end

function Base.hash(a::TotFrac, h::UInt)
   # Elements cannot be canonicalised
   0xcdd13b94adbec4bc%UInt
end

zero(R::TotFracRing) = R(0)

one(R::TotFracRing) = R(1)

iszero(a::TotFrac) = iszero(numerator(a, false))

isone(a::TotFrac) = numerator(a, false) == denominator(a, false)

is_unit(a::TotFrac) = !is_zero_divisor(numerator(a, false))

is_zero_divisor(a::TotFrac) = is_zero_divisor(numerator(a, false))

function is_zero_divisor_with_annihilator(a::TotFrac)
   f, b = is_zero_divisor_with_annihilator(numerator(a, false))
   R = parent(a)
   return f, R(b, one(base_ring(R)), false)
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(a::TotFrac) = one(base_ring(a))

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function expressify(a::TotFrac; context = nothing)
    n = numerator(a, true)
    d = denominator(a, true)
    if isone(d)
        return expressify(n; context)
    else
        return Expr(:call, ://, expressify(n; context), expressify(d; context))
    end
end

@enable_all_show_via_expressify TotFrac

function show(io::IO, a::TotFracRing)
   print(IOContext(io, :compact => true), "Total ring of fractions of ", base_ring(a))
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::TotFrac)
   R = parent(a)
   return R(-numerator(a, false), deepcopy(denominator(a, false)), false)
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::TotFrac{T}, b::TotFrac{T}) where {T <: RingElem}
   check_parent(a, b)
   d1 = denominator(a, false)
   d2 = denominator(b, false)
   n1 = numerator(a, false)
   n2 = numerator(b, false)
   if d1 == d2
      rnum = n1 + n2
      rden = deepcopy(d1)
   elseif isone(d1)
      rnum = n1*d2 + n2
      rden = deepcopy(d2)
   elseif isone(d2)
      rnum = n1 + n2*d1
      rden = deepcopy(d1)
   else
      rnum = n1*d2 + n2*d1
      rden = d1*d2
   end
   return parent(a)(rnum, rden, false)
end

function -(a::TotFrac{T}, b::TotFrac{T}) where {T <: RingElem}
   check_parent(a, b)
   d1 = denominator(a, false)
   d2 = denominator(b, false)
   n1 = numerator(a, false)
   n2 = numerator(b, false)
   if d1 == d2
      rnum = n1 - n2
      rden = deepcopy(d1)
   elseif isone(d1)
      rnum = n1*d2 - n2
      rden = deepcopy(d2)
   elseif isone(d2)
      rnum = n1 - n2*d1
      rden = deepcopy(d1)
   else
      rnum = n1*d2 - n2*d1
      rden = d1*d2
   end
   return parent(a)(rnum, rden, false)
end

function *(a::TotFrac{T}, b::TotFrac{T}) where {T <: RingElem}
   check_parent(a, b)
   n1 = numerator(a, false)
   d2 = denominator(b, false)
   n2 = numerator(b, false)
   d1 = denominator(a, false)
   if n1 == d2
      n = deepcopy(n2)
      d = deepcopy(d1)
   elseif n2 == d1
      n = deepcopy(n1)
      d = deepcopy(d2)
   else
      n = n1*n2
      d = d1*d2
   end
   return parent(a)(n, d, false)
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(a::TotFrac, b::Union{Integer, Rational})
   n = numerator(a, false)
   d = denominator(a, false)
   if b == d
      return parent(a)(deepcopy(n))
   end
   return parent(a)(n*b, deepcopy(d), false)
end

function *(a::Union{Integer, Rational}, b::TotFrac)
   n = numerator(b, false)
   d = denominator(b, false)
   if a == d
      return parent(b)(deepcopy(n))
   end
   return parent(b)(a*n, deepcopy(d), false)
end

function *(a::TotFrac{T}, b::T) where {T <: RingElem}
   n = numerator(a, false)
   d = denominator(a, false)
   if b == d
      return parent(a)(deepcopy(n))
   end
   return parent(a)(n*b, deepcopy(d), false)
end

function *(a::T, b::TotFrac{T}) where {T <: RingElem}
   n = numerator(b, false)
   d = denominator(b, false)
   if a == d
      return parent(a)(deepcopy(n))
   end
   return parent(b)(a*n, deepcopy(d), false)
end

function +(a::TotFrac, b::Integer)
   n = numerator(a, false) + denominator(a, false)*b
   d = denominator(a, false)
   return parent(a)(n, deepcopy(d))
end

function -(a::TotFrac, b::Integer, Rational)
   n = numerator(a, false) - denominator(a, false)*b
   d = denominator(a, false)
   return parent(a)(n, deepcopy(d))
end

+(a::Integer, b::TotFrac) = b + a

function -(a::Integer, b::TotFrac)
   n = a*denominator(b, false) - numerator(b, false)
   d = denominator(b, false)
   return parent(b)(n, deepcopy(d))
end

function +(a::TotFrac{T}, b::T) where {T <: RingElem}
   n = numerator(a, false) + denominator(a, false)*b
   d = denominator(a, false)
   return parent(a)(n, deepcopy(d), false)
end

function -(a::TotFrac{T}, b::T) where {T <: RingElem}
   n = numerator(a, false) - denominator(a, false)*b
   d = denominator(a, false)
   return parent(a)(n, deepcopy(d), false)
end

+(a::T, b::TotFrac{T}) where {T <: RingElem} = b + a

function -(a::T, b::TotFrac{T}) where {T <: RingElem}
   n = a*denominator(b, false) - numerator(b, false)
   d = denominator(b, false)
   return parent(b)(n, deepcopy(d))
end

###############################################################################
#
#   Comparisons
#
###############################################################################

function ==(x::TotFrac{T}, y::TotFrac{T}) where {T <: RingElem}
   b  = check_parent(x, y, false)
   !b && return false

   return (denominator(x, false) == denominator(y, false) &&
           numerator(x, false) == numerator(y, false)) ||
          (numerator(x, false)*denominator(y, false) ==
           denominator(x, false)*numerator(y, false))
end

function isequal(x::TotFrac{T}, y::TotFrac{T}) where {T <: RingElem}
   return x == y
end

###############################################################################
#
#   Ad hoc comparisons
#
###############################################################################

function ==(x::TotFrac, y::Union{Integer, Rational})
   return (isone(denominator(x, false)) && numerator(x, false) == y) ||
          (numerator(x, false) == denominator(x, false)*y)
end

==(x::Union{Integer, Rational}, y::TotFrac) = y == x

function ==(x::TotFrac{T}, y::T) where {T <: RingElem}
   return (isone(denominator(x, false)) && numerator(x, false) == y) ||
          (numerator(x, false) == denominator(x, false)*y)
end

==(x::T, y::TotFrac{T}) where {T <: RingElem} = y == x

###############################################################################
#
#   Inversion
#
###############################################################################

@doc raw"""
    Base.inv(a::TotFrac)

Return the inverse of the fraction $a$ if it exists, otherwise raise an
exception.
"""
function Base.inv(a::TotFrac)
   is_zero_divisor(numerator(a, false)) && throw(NotInvertibleError(a))
   return parent(a)(deepcopy(denominator(a, false)),
                    deepcopy(numerator(a, false)))
end

###############################################################################
#
#   Exact division
#
###############################################################################

# Division is not possible generically due to the possibility of non-units
# that are also non-zero divisors

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

# Division is not possible generically due to the possibility of non-units
# that are also non-zero divisors

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::TotFrac{T}, b::Int) where {T <: RingElem}
   if b < 0
      a = inv(a)
      b = -b
   end
   return parent(a)(numerator(a)^b, denominator(a)^b)
end

###############################################################################
#
#   Unsafe operators and functions
#
###############################################################################

function zero!(c::TotFrac)
   c.num = zero!(c.num)
   if !isone(c.den)
      c.den = one(base_ring(c))
   end
   return c
end

function mul!(c::TotFrac{T}, a::TotFrac{T}, b::TotFrac{T}) where {T <: RingElem}
   R = parent(a)
   n1 = numerator(a, false)
   d2 = denominator(b, false)
   n2 = numerator(b, false)
   d1 = denominator(a, false)
   if n1 == d2
      n = deepcopy(n2)
      d = deepcopy(d1)
   elseif n2 == d1
      n = deepcopy(n1)
      d = deepcopy(d2)
   else
      n = n1*n2
      d = d1*d2
   end
   if n == d
      c.num = one(R)
      c.den = one(R)
   else
      c.num = n
      c.den = d
   end
   return c
end

function addeq!(a::TotFrac{T}, b::TotFrac{T}) where {T <: RingElem}
   R = base_ring(b)
   d1 = denominator(a, false)
   d2 = denominator(b, false)
   n1 = numerator(a, false)
   n2 = numerator(b, false)
   if d1 == d2
      a.num = addeq!(a.num, b.num)
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
      a.num = d1*n2 + d2*n1
      a.den = mul!(a.den, a.den, d2)
   end
   if a.num == a.den
      a.num = one(R)
      a.den = one(R)
   end
   return a
end

function add!(c::TotFrac{T}, a::TotFrac{T}, b::TotFrac{T}) where {T <: RingElem}
   R = base_ring(a)
   d1 = denominator(a, false)
   d2 = denominator(b, false)
   n1 = numerator(a, false)
   n2 = numerator(b, false)
   if d1 == d2
      c.num = add!(c.num, n1, n2)
      c.den = deepcopy(d1)
   elseif isone(d1)
      c.num = add!(c.num, n1*d2, n2)
      c.den = deepcopy(d2)
   elseif isone(d2)
      c.num = add!(c.num, n1, n2*d1)
      c.den = deepcopy(d1)
   else
      c.num = add!(c.num, n1*d2, n2*d1)
      c.den = mul!(c.den, d1, d2)
   end
   if c.num == c.den
      c.num = one(R)
      c.den = one(R)
   end
   return c
end

###############################################################################
#
#   Random functions
#
###############################################################################

RandomExtensions.maketype(R::TotFracRing, _) = elem_type(R)

function RandomExtensions.make(S::TotFracRing, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, vs[1]) # forward to default Make constructor
   else
      Make(S, make(R, vs...))
   end
end

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make2{<:RingElement, <:TotFracRing}})
   S, v = sp[][1:end]
   R = base_ring(S)
   n = rand(rng, v)
   d = R()
   while is_zero_divisor(d)
      d = rand(rng, v)
   end
   return S(n, d)
end

rand(rng::AbstractRNG, S::TotFracRing, v...) =
   rand(rng, make(S, v...))

rand(S::TotFracRing, v...) = rand(GLOBAL_RNG, S, v...)

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{TotFrac{T}}, ::Type{TotFrac{T}}) where T <: RingElement = TotFrac{T}
promote_rule(::Type{TotFrac{T}}, ::Type{TotFrac{T}}) where T <: RingElem = TotFrac{T}

function promote_rule(::Type{TotFrac{T}}, ::Type{U}) where {T <: RingElem, U <: RingElem}
   promote_rule(T, U) == T ? TotFrac{T} : Union{}
end

###############################################################################
#
#   Parent object call overloading
#
###############################################################################

function (a::TotFracRing{T})(b::RingElement) where {T <: RingElement}
   return a(base_ring(a)(b))
end

function (a::TotFracRing{T})() where {T <: RingElement}
   R = base_ring(a)
   z = TotFrac{T}(zero(R), one(R))
   z.parent = a
   return z
end

function (a::TotFracRing{T})(b::T) where {T <: RingElement}
   R = base_ring(a)
   parent(b) != R && error("Could not coerce to fraction")
   z = TotFrac{T}(b, one(R))
   z.parent = a
   return z
end

function (a::TotFracRing{T})(b::T, c::T, check=true) where {T <: RingElement}
   R = base_ring(a)
   parent(b) != R && error("Could not coerce to fraction")
   parent(c) != R && error("Could not coerce to fraction")
   check && is_zero_divisor(c) && error("Denominator is a zero divisor")
   if b == c
      return one(a)
   end
   z = TotFrac{T}(b, c)
   z.parent = a
   return z
end

function (a::TotFracRing{T})(b::T, c::Union{Integer, Rational}, check=true) where {T <: RingElement}
   a(b, base_ring(a)(c), check)
end

function (a::TotFracRing{T})(b::Union{Integer, Rational}, c::T, check=true) where {T <: RingElement}
   a(base_ring(a)(b), c, check)
end

function (a::TotFracRing{T})(b::Union{Integer, Rational}) where {T <: RingElement}
   a(base_ring(a)(b))
end

function (a::TotFracRing{T})(b::Integer, c::Integer, check=true) where {T <: RingElement}
   a(base_ring(a)(b), base_ring(a)(c), check)
end

function (a::TotFracRing{T})(b::TotFrac{T}) where {T <: RingElement}
   a != parent(b) && error("Could not coerce to fraction")
   return b
end

###############################################################################
#
#   total_ring_of_fractions constructor
#
###############################################################################

@doc raw"""
    total_ring_of_fractions(R::Ring; cached::Bool=true)

Return the parent object of the total ring of fractions over the given base
ring $R$, i.e. the localisation of `R` at the complement of the set of zero
divisors.

If `cached == true` (the default), the returned parent object is cached so
that it will always be returned by a call to the constructor when the same
base ring $R$ is supplied.
"""
function total_ring_of_fractions(R::Ring; cached::Bool=true)
   T = elem_type(R)

   return TotFracRing{T}(R, cached)
end
