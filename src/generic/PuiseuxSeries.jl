###############################################################################
#
#   PuiseuxSeries.jl : Puiseux series over rings and fields
#
###############################################################################

export laurent_ring, rescale!

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

doc"""
    laurent_ring{T <: RingElement}(R::PuiseuxSeriesRing{T})
> Return the `LaurentSeriesRing` underlying the given `PuiseuxSeriesRing`.
""" 
laurent_ring(R::PuiseuxSeriesRing{T}) where T <: RingElement = R.laurent_ring::LaurentSeriesRing{T}

doc"""
    laurent_ring{T <: FieldElement}(R::PuiseuxSeriesField{T})
> Return the `LaurentSeriesField` underlying the given `PuiseuxSeriesField`.
"""
laurent_ring(R::PuiseuxSeriesField{T}) where T <: FieldElement = R.laurent_ring::LaurentSeriesField{T}

doc"""
    O{T <: RingElement}(a::PuiseuxSeriesElem{T})
> Returns $0 + O(x^\mbox{val}(a))$. Usually this function is called with $x^n$
> as parameter for some rational $n$. Then the function returns the Puiseux series
> $0 + O(x^n)$, which can be used to set the precision of a Puiseux series when
> constructing it.
"""
function O(a::PuiseuxSeriesElem{T}) where T <: RingElement
   val = valuation(a)
   par = parent(a)
   laur = laurent_ring(par)(Array{T}(0), 0, numerator(val), numerator(val))
   return parent(a)(laur, 1//denominator(val))
end

parent_type(::Type{T}) where {S <: RingElement, T <: PuiseuxSeriesRingElem{S}} = PuiseuxSeriesRing{S}

parent_type(::Type{T}) where {S <: FieldElement, T <: PuiseuxSeriesFieldElem{S}} = PuiseuxSeriesField{S}

doc"""
    parent(a::PuiseuxSeriesElem)
> Return the parent of the given Puiseux series.
"""
parent(a::PuiseuxSeriesElem) = a.parent

elem_type(::Type{T}) where {S <: RingElement, T <: PuiseuxSeriesRing{S}} = PuiseuxSeriesRingElem{S}

elem_type(::Type{T}) where {S <: FieldElement, T <: PuiseuxSeriesField{S}} = PuiseuxSeriesFieldElem{S}

doc"""
    base_ring(R::PuiseuxSeriesRing)
> Return the base (coefficient) ring of the given Puiseux series ring.
"""
base_ring(R::PuiseuxSeriesRing{T}) where T <: RingElement = base_ring(laurent_ring(R))

doc"""
    base_ring(R::PuiseuxSeriesField)
> Return the base (coefficient) ring of the given Puiseux series field.
"""
base_ring(R::PuiseuxSeriesField{T}) where T <: FieldElement = base_ring(laurent_ring(R))

doc"""
    base_ring(a::PuiseuxSeriesElem)
> Return the base (coefficient) ring of the Puiseux series ring of the given Puiseux
> series.
"""
base_ring(a::PuiseuxSeriesElem) = base_ring(parent(a))

function isdomain_type(::Type{T}) where {S <: RingElement, T <: PuiseuxSeriesElem{S}}
   return isdomain_type(S)
end

isexact_type(a::Type{T}) where T <: PuiseuxSeriesElem = false

function check_parent(a::PuiseuxSeriesElem, b::PuiseuxSeriesElem)
   parent(a) != parent(b) &&
             error("Incompatible Puiseux series rings in Puiseux series operation")
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::PuiseuxSeriesElem, h::UInt)
   b = 0xec4c3951832c37f0%UInt
   b = xor(b, hash(a.data, h))
   b = xor(b, hash(a.scale, h))
   return b
end

doc"""
    precision(a::PuiseuxSeriesElem)
> Return the precision of the given Puiseux series in absolute terms.
"""
precision(a::PuiseuxSeriesElem) = precision(a.data)*a.scale

doc"""
    valuation(a::PuiseuxSeriesElem)
> Return the valuation of the given Puiseux series, i.e. the exponent of the first
> nonzero term (or the precision if it is arithmetically zero).
"""
valuation(a::PuiseuxSeriesElem) = valuation(a.data)*a.scale

doc"""
    zero(R::PuiseuxSeriesRing)
> Return $0 + O(x^n)$ where $n$ is the maximum precision of the Puiseux series
> ring $R$.
"""
zero(R::PuiseuxSeriesRing) = R(0)

doc"""
    zero(R::PuiseuxSeriesField)
> Return $0 + O(x^n)$ where $n$ is the maximum precision of the Puiseux series
> ring $R$.
"""
zero(R::PuiseuxSeriesField) = R(0)

doc"""
    one(R::PuiseuxSeriesRing)
> Return $1 + O(x^n)$ where $n$ is the maximum precision of the Puiseux series
> ring $R$.
"""
one(R::PuiseuxSeriesField) = R(1)

doc"""
    one(R::PuiseuxSeriesField)
> Return $1 + O(x^n)$ where $n$ is the maximum precision of the Puiseux series
> ring $R$.
"""
one(R::PuiseuxSeriesRing) = R(1)

doc"""
    gen(R::PuiseuxSeriesRing)
> Return the generator of the Puiseux series ring, i.e. $x + O(x^{n + 1})$ where
> $n$ is the maximum precision of the Puiseux series ring $R$.
"""
function gen(R::PuiseuxSeriesRing)
   S = laurent_ring(R)
   return R(gen(S), 1//1)
end

doc"""
    gen(R::PuiseuxSeriesField)
> Return the generator of the Puiseux series ring, i.e. $x + O(x^{n + 1})$ where
> $n$ is the maximum precision of the Puiseux series ring $R$.
"""
function gen(R::PuiseuxSeriesField)
   S = laurent_ring(R)
   return R(gen(S), 1//1)
end

doc"""
    iszero(a::PuiseuxSeriesElem)
> Return `true` if the given Puiseux series is arithmetically equal to zero to
> its current precision, otherwise return `false`.
"""
iszero(a::PuiseuxSeriesElem) = iszero(a.data)

doc"""
    isone(a::PuiseuxSeriesElem)
> Return `true` if the given Puiseux series is arithmetically equal to one to
> its current precision, otherwise return `false`.
"""
function isone(a::PuiseuxSeriesElem)
   return isone(a.data)
end

doc"""
    isgen(a::PuiseuxSeriesElem)
> Return `true` if the given Puiseux series is arithmetically equal to the
> generator of its Puiseux series ring to its current precision, otherwise return
> `false`.
"""
function isgen(a::PuiseuxSeriesElem)
   return valuation(a) == 1 && pol_length(a.data) == 1 && isone(polcoeff(a.data, 0))
end

doc"""
    isunit(a::PuiseuxSeriesElem)
> Return `true` if the given Puiseux series is arithmetically equal to a unit,
> i.e. is invertible, otherwise return `false`.
"""
isunit(a::PuiseuxSeriesElem) = valuation(a) == 0 && isunit(polcoeff(a.data, 0))

doc"""
    modulus{T <: ResElem}(a::PuiseuxSeriesElem{T})
> Return the modulus of the coefficients of the given Puiseux series.
"""
modulus(a::PuiseuxSeriesElem{T}) where {T <: ResElem} = modulus(base_ring(a))

doc"""
    rescale(a::PuiseuxSeriesElem)
> Rescale the polynomial underlying $a$ so that the greatest common divisor of the
> exponents of the underlying Laurent series is $1$. This function is mainly used
> internally, as the output of Puiseux series functions is assumed to be rescaled.
"""
function rescale!(a::PuiseuxSeriesElem)
   if !iszero(a)
      d = exp_gcd(a.data)
      if d == 1
         return a
      end
      a.data = deflate(a.data, d)
      a.scale *= d
      return a
   end
end

function deepcopy_internal(a::PuiseuxSeriesElem{T}, dict::ObjectIdDict) where {T <: RingElement}
    return parent(a)(deepcopy(a.data), a.scale)
end

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function show(io::IO, x::PuiseuxSeriesElem)
   print(io, x.data, ", ", x.scale)
end

function show(io::IO, a::PuiseuxSeriesRing)
   print(io, "Puiseux series ring in ", var(laurent_ring(a)), " over ")
   show(io, base_ring(a))
end

function show(io::IO, a::PuiseuxSeriesField)
   print(io, "Puiseux series field in ", var(laurent_ring(a)), " over ")
   show(io, base_ring(a))
end

needs_parentheses(x::PuiseuxSeriesElem) = pol_length(x.data) > 1

isnegative(x::PuiseuxSeriesElem) = pol_length(x) <= 1 && isnegative(polcoeff(x.data, 0))

show_minus_one(::Type{PuiseuxSeriesElem{T}}) where {T <: RingElement} = show_minus_one(T)

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::PuiseuxSeriesElem)
   R = parent(a)
   return R(-a.data, a.scale)
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::PuiseuxSeriesElem{T}, b::PuiseuxSeriesElem{T}) where T <: RingElement
    zscale = gcd(a.scale, b.scale)
    ainf = numerator(a.scale//zscale)
    binf = numerator(b.scale//zscale)
    z = parent(a)(inflate(a.data, ainf) + inflate(b.data, binf), zscale)
    z = rescale!(z)
    return z
end

function -(a::PuiseuxSeriesElem{T}, b::PuiseuxSeriesElem{T}) where T <: RingElement
    zscale = gcd(a.scale, b.scale)
    ainf = numerator(a.scale//zscale)
    binf = numerator(b.scale//zscale)
    z = parent(a)(inflate(a.data, ainf) - inflate(b.data, binf), zscale)
    z = rescale!(z)
    return z
end

function *(a::PuiseuxSeriesElem{T}, b::PuiseuxSeriesElem{T}) where T <: RingElement
    zscale = gcd(a.scale, b.scale)
    ainf = numerator(a.scale//zscale)
    binf = numerator(b.scale//zscale)
    z = parent(a)(inflate(a.data, ainf)*inflate(b.data, binf), zscale)
    z = rescale!(z)
    return z
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::PuiseuxSeriesElem{T}, b::PuiseuxSeriesElem{T}) where T <: RingElement
    zscale = gcd(a.scale, b.scale)
    ainf = numerator(a.scale//zscale)
    binf = numerator(b.scale//zscale)
    z = parent(a)(divexact(inflate(a.data, ainf), inflate(b.data, binf)), zscale)
    z = rescale!(z)
    return z
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(a::PuiseuxSeriesElem{T}, b::Int) where T <: RingElement
   # special case powers of x for constructing power series efficiently
   if pol_length(a.data) == 0
      return parent(a)(a.data^0, a.scale)
   elseif b == 0
      # in fact, the result would be exact 1 if we had exact series
      return one(parent(a))
   elseif pol_length(a.data) == 1
      return parent(a)(a.data^b, a.scale)
   elseif b == 1
      return deepcopy(a)
   elseif b == -1
      return inv(a)
   end

   if b < 0
      a = inv(a)
      b = -b
   end

   z = parent(a)(a.data^b, a.scale)
   z = rescale!(z)
   return z
end

function ^(a::PuiseuxSeriesElem{T}, b::Rational{Int}) where T <: RingElement
   (pol_length(a.data) != 1 || polcoeff(a.data, 0) != 1) && error("Rational power not implemented")
   z = parent(a)(a.data^numerator(b), a.scale//denominator(b))
   z = rescale!(z)
   return z
end
   
###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{PuiseuxSeriesElem{T}}, ::Type{PuiseuxSeriesElem{T}}) where T <: RingElement = PuiseuxSeriesElem{T}

function promote_rule(::Type{PuiseuxSeriesElem{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? PuiseuxSeriesElem{T} : Union{}
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (R::PuiseuxSeriesRing{T})(b::RingElement) where {T <: RingElement}
   return R(base_ring(R)(b))
end

function (R::PuiseuxSeriesField{T})(b::RingElement) where {T <: FieldElement}
   return R(base_ring(R)(b))
end

function (R::PuiseuxSeriesRing{T})() where {T <: RingElement}
   z = PuiseuxSeriesRingElem{T}(laurent_ring(R)(), 1//1)
   z.parent = R
   return z
end

function (R::PuiseuxSeriesField{T})() where {T <: FieldElement}
   z = PuiseuxSeriesFieldElem{T}(laurent_ring(R)(), 1//1)
   z.parent = R
   return z
end

function (R::PuiseuxSeriesRing{T})(b::LaurentSeriesRingElem{T}, scale::Rational{Int}) where T <: RingElement
   z = PuiseuxSeriesRingElem{T}(b, scale)
   z.parent = R
   return z
end

function (R::PuiseuxSeriesField{T})(b::LaurentSeriesFieldElem{T}, scale::Rational{Int}) where T <: FieldElement
   z = PuiseuxSeriesFieldElem{T}(b, scale)
   z.parent = R
   return z
end

function (R::PuiseuxSeriesRing{T})(b::Union{Integer, Rational, AbstractFloat}) where T <: RingElement
   z = PuiseuxSeriesRingElem{T}(laurent_ring(R)(b), 1//1)
   z.parent = R
   return z
end

function (R::PuiseuxSeriesField{T})(b::Union{Rational, AbstractFloat}) where T <: FieldElement
   z = PuiseuxSeriesFieldElem{T}(laurent_ring(R)(b), 1//1)
   z.parent = R
   return z
end

function (R::PuiseuxSeriesRing{T})(b::T) where T <: RingElem
   parent(b) != base_ring(R) && error("Unable to coerce to Puiseux series")
   z = PuiseuxSeriesRingElem{T}(laurent_ring(R)(b), 1//1)
   z.parent = R
   return z
end

function (R::PuiseuxSeriesField{T})(b::T) where T <: FieldElem
   parent(b) != base_ring(R) && error("Unable to coerce to Puiseux series")
   z = PuiseuxSeriesFieldElem{T}(laurent_ring(R)(b), 1//1)
   z.parent = R
   return z
end

function (R::PuiseuxSeriesRing{T})(b::PuiseuxSeriesElem{T}) where T <: RingElement
   parent(b) != R && error("Unable to coerce Puiseux series")
   return b
end

function (R::PuiseuxSeriesField{T})(b::PuiseuxSeriesElem{T}) where T <: FieldElement
   parent(b) != R && error("Unable to coerce Puiseux series")
   return b
end

###############################################################################
#
#   PuiseuxSeriesRing constructor
#
###############################################################################

doc"""
   PuiseuxSeriesRing(R::AbstractAlgebra.Ring, prec::Int, s::AbstractString; cached=true)
> Return a tuple $(S, x)$ consisting of the parent object `S` of a Puiseux series
> ring over the given base ring and a generator `x` for the Puiseux series ring.
> The maximum precision of the series in the ring is set to `prec`. This is taken as a
> maximum relative precision of the underlying Laurent series that are used to implement
> the Puiseux series in the ring. The supplied string `s` specifies the way the
> generator of the Puiseux series ring will be printed. By default, the parent
> object `S` will be cached so that supplying the same base ring, string and
> precision in future will return the same parent object and generator. If
> caching of the parent object is not required, `cached` can be set to `false`.
"""
function PuiseuxSeriesRing(R::AbstractAlgebra.Ring, prec::Int, s::AbstractString; cached=true)
   S, x = LaurentSeriesRing(R, prec, s; cached=cached)
   T = elem_type(R)

   parent_obj = PuiseuxSeriesRing{T}(S, cached)

   return parent_obj, gen(parent_obj)
end

function PuiseuxSeriesRing(R::AbstractAlgebra.Field, prec::Int, s::AbstractString; cached= true)
   S, x = LaurentSeriesField(R, prec, s; cached=cached)
   T = elem_type(R)

   parent_obj = PuiseuxSeriesField{T}(S, cached)

   return parent_obj, gen(parent_obj)
end

function PuiseuxSeriesField(R::AbstractAlgebra.Field, prec::Int, s::AbstractString; cached = true)
   S, x = LaurentSeriesField(R, prec, s; cached=cached)
   T = elem_type(R)

   parent_obj = PuiseuxSeriesField{T}(S, cached)

   return parent_obj, gen(parent_obj)
end

