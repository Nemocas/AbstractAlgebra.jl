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

@doc Markdown.doc"""
    laurent_ring(R::PuiseuxSeriesRing{T}) where T <: RingElement
> Return the `LaurentSeriesRing` underlying the given `PuiseuxSeriesRing`.
"""
laurent_ring(R::PuiseuxSeriesRing{T}) where T <: RingElement = R.laurent_ring::LaurentSeriesRing{T}

@doc Markdown.doc"""
    laurent_ring(R::PuiseuxSeriesField{T}) where T <: FieldElement
> Return the `LaurentSeriesField` underlying the given `PuiseuxSeriesField`.
"""
laurent_ring(R::PuiseuxSeriesField{T}) where T <: FieldElement = R.laurent_ring::LaurentSeriesField{T}

@doc Markdown.doc"""
    O(a::Generic.PuiseuxSeriesElem{T}) where T <: RingElement
> Return $0 + O(x^\mbox{val}(a))$. Usually this function is called with $x^n$
> as parameter for some rational $n$. Then the function returns the Puiseux series
> $0 + O(x^n)$, which can be used to set the precision of a Puiseux series when
> constructing it.
"""
function O(a::PuiseuxSeriesElem{T}) where T <: RingElement
   val = valuation(a)
   par = parent(a)
   laur = laurent_ring(par)(Array{T}(undef, 0), 0, numerator(val), numerator(val), 1)
   return parent(a)(laur, denominator(val))
end

parent_type(::Type{T}) where {S <: RingElement, T <: PuiseuxSeriesRingElem{S}} = PuiseuxSeriesRing{S}

parent_type(::Type{T}) where {S <: FieldElement, T <: PuiseuxSeriesFieldElem{S}} = PuiseuxSeriesField{S}

@doc Markdown.doc"""
    parent(a::Generic.PuiseuxSeriesElem)
> Return the parent of the given Puiseux series.
"""
parent(a::PuiseuxSeriesElem) = a.parent

elem_type(::Type{T}) where {S <: RingElement, T <: PuiseuxSeriesRing{S}} = PuiseuxSeriesRingElem{S}

elem_type(::Type{T}) where {S <: FieldElement, T <: PuiseuxSeriesField{S}} = PuiseuxSeriesFieldElem{S}

@doc Markdown.doc"""
    base_ring(R::PuiseuxSeriesRing{T}) where T <: RingElement
> Return the base (coefficient) ring of the given Puiseux series ring.
"""
base_ring(R::PuiseuxSeriesRing{T}) where T <: RingElement = base_ring(laurent_ring(R))

@doc Markdown.doc"""
    base_ring(R::PuiseuxSeriesField{T}) where T <: FieldElement
> Return the base (coefficient) ring of the given Puiseux series field.
"""
base_ring(R::PuiseuxSeriesField{T}) where T <: FieldElement = base_ring(laurent_ring(R))

@doc Markdown.doc"""
    base_ring(a::Generic.PuiseuxSeriesElem)
> Return the base (coefficient) ring of the Puiseux series ring of the given Puiseux
> series.
"""
base_ring(a::PuiseuxSeriesElem) = base_ring(parent(a))

@doc Markdown.doc"""
    max_precision(R::PuiseuxSeriesRing{T}) where T <: RingElement
> Return the maximum precision of the underlying Laurent series ring.
"""
max_precision(R::PuiseuxSeriesRing{T}) where T <: RingElement = max_precision(laurent_ring(R))

@doc Markdown.doc"""
    max_precision(R::PuiseuxSeriesField{T}) where T <: FieldElement
> Return the maximum precision of the underlying Laurent series field.
"""
max_precision(R::PuiseuxSeriesField{T}) where T <: FieldElement = max_precision(laurent_ring(R))

@doc Markdown.doc"""
    var(R::PuiseuxSeriesRing{T}) where T <: RingElement
> Return a symbol representing the variable of the given Puiseux series ring.
"""
var(R::PuiseuxSeriesRing{T}) where T <: RingElement = var(laurent_ring(R))

@doc Markdown.doc"""
    var(R::PuiseuxSeriesField{T}) where T <: FieldElement
> Return a symbol representing the variable of the given Puiseux series field.
"""
var(R::PuiseuxSeriesField{T}) where T <: FieldElement = var(laurent_ring(R))

function isdomain_type(::Type{T}) where {S <: RingElement, T <: PuiseuxSeriesElem{S}}
   return isdomain_type(S)
end

isexact_type(a::Type{T}) where T <: PuiseuxSeriesElem = false

function check_parent(a::PuiseuxSeriesElem, b::PuiseuxSeriesElem, throw::Bool = true)
   fl = parent(a) != parent(b)
   fl && throw && error("Incompatible Puiseux series rings in Puiseux series operation")
   return !fl
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

@doc Markdown.doc"""
    precision(a::Generic.PuiseuxSeriesElem)
> Return the precision of the given Puiseux series in absolute terms.
"""
precision(a::PuiseuxSeriesElem) = precision(a.data)//a.scale

@doc Markdown.doc"""
    valuation(a::Generic.PuiseuxSeriesElem)
> Return the valuation of the given Puiseux series, i.e. the exponent of the first
> nonzero term (or the precision if it is arithmetically zero).
"""
valuation(a::PuiseuxSeriesElem) = valuation(a.data)//a.scale

scale(a::PuiseuxSeriesElem) = a.scale

@doc Markdown.doc"""
    coeff(a::Generic.PuiseuxSeriesElem, n::Int)
> Return the coefficient of the term of exponent $n$ of the given Puiseux series.
"""
function coeff(a::PuiseuxSeriesElem, n::Int)
   s = scale(a)
   return coeff(a.data, n*s)
end

@doc Markdown.doc"""
    coeff(a::Generic.PuiseuxSeriesElem, r::Rational{Int})
> Return the coefficient of the term of exponent $r$ of the given Puiseux series.
"""
function coeff(a::PuiseuxSeriesElem, r::Rational{Int})
   s = scale(a)
   n = numerator(r)
   d = denominator(r)
   if mod(s, d) != 0
      return base_ring(a)()
   end
   return coeff(a.data, n*div(s, d))
end

@doc Markdown.doc"""
    zero(R::PuiseuxSeriesRing)
> Return $0 + O(x^n)$ where $n$ is the maximum precision of the Puiseux series
> ring $R$.
"""
zero(R::PuiseuxSeriesRing) = R(0)

@doc Markdown.doc"""
    zero(R::PuiseuxSeriesField)
> Return $0 + O(x^n)$ where $n$ is the maximum precision of the Puiseux series
> ring $R$.
"""
zero(R::PuiseuxSeriesField) = R(0)

@doc Markdown.doc"""
    one(R::PuiseuxSeriesField)
> Return $1 + O(x^n)$ where $n$ is the maximum precision of the Puiseux series
> ring $R$.
"""
one(R::PuiseuxSeriesField) = R(1)

@doc Markdown.doc"""
    one(R::PuiseuxSeriesRing)
> Return $1 + O(x^n)$ where $n$ is the maximum precision of the Puiseux series
> ring $R$.
"""
one(R::PuiseuxSeriesRing) = R(1)

@doc Markdown.doc"""
    gen(R::PuiseuxSeriesRing)
> Return the generator of the Puiseux series ring, i.e. $x + O(x^{n + 1})$ where
> $n$ is the maximum precision of the Puiseux series ring $R$.
"""
function gen(R::PuiseuxSeriesRing)
   S = laurent_ring(R)
   return R(gen(S), 1)
end

@doc Markdown.doc"""
    gen(R::PuiseuxSeriesField)
> Return the generator of the Puiseux series ring, i.e. $x + O(x^{n + 1})$ where
> $n$ is the maximum precision of the Puiseux series ring $R$.
"""
function gen(R::PuiseuxSeriesField)
   S = laurent_ring(R)
   return R(gen(S), 1)
end

@doc Markdown.doc"""
    iszero(a::Generic.PuiseuxSeriesElem)
> Return `true` if the given Puiseux series is arithmetically equal to zero to
> its current precision, otherwise return `false`.
"""
iszero(a::PuiseuxSeriesElem) = iszero(a.data)

@doc Markdown.doc"""
    isone(a::Generic.PuiseuxSeriesElem)
> Return `true` if the given Puiseux series is arithmetically equal to one to
> its current precision, otherwise return `false`.
"""
function isone(a::PuiseuxSeriesElem)
   return isone(a.data)
end

@doc Markdown.doc"""
    isgen(a::Generic.PuiseuxSeriesElem)
> Return `true` if the given Puiseux series is arithmetically equal to the
> generator of its Puiseux series ring to its current precision, otherwise return
> `false`.
"""
function isgen(a::PuiseuxSeriesElem)
   return valuation(a) == 1 && pol_length(a.data) == 1 && isone(polcoeff(a.data, 0))
end

@doc Markdown.doc"""
    isunit(a::Generic.PuiseuxSeriesElem)
> Return `true` if the given Puiseux series is arithmetically equal to a unit,
> i.e. is invertible, otherwise return `false`.
"""
isunit(a::PuiseuxSeriesElem) = valuation(a) == 0 && isunit(polcoeff(a.data, 0))

@doc Markdown.doc"""
    modulus(a::Generic.PuiseuxSeriesElem{T}) where {T <: ResElem}
> Return the modulus of the coefficients of the given Puiseux series.
"""
modulus(a::PuiseuxSeriesElem{T}) where {T <: ResElem} = modulus(base_ring(a))

@doc Markdown.doc"""
    rescale!(a::Generic.PuiseuxSeriesElem)
> Rescale so that the scale of the given Puiseux series and the scale of the underlying
> Laurent series are coprime. This function is used internally, as all user facing
> functions are assumed to rescale their output.
"""
function rescale!(a::PuiseuxSeriesElem)
   if !iszero(a)
      d = gcd(a.scale, gcd(scale(a.data), gcd(valuation(a.data), precision(a.data))))
      if d != 1
         set_scale!(a.data, div(scale(a.data), d))
         set_prec!(a.data, AbstractAlgebra.div(precision(a.data), d))
         set_val!(a.data, AbstractAlgebra.div(valuation(a.data), d))
         a.scale = div(a.scale, d)
      end
   else
      d = gcd(precision(a.data), a.scale)
      if d != 1
         set_prec!(a.data, AbstractAlgebra.div(precision(a.data), d))
         set_val!(a.data, AbstractAlgebra.div(valuation(a.data), d))
         a.scale = div(a.scale, d)
      end
   end
   return a
end

function deepcopy_internal(a::PuiseuxSeriesElem{T}, dict::IdDict) where {T <: RingElement}
    return parent(a)(deepcopy(a.data), a.scale)
end

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function show(io::IO, x::PuiseuxSeriesElem)
   len = pol_length(x.data)
   if len == 0
      print(IOContext(io, :compact => true), zero(base_ring(x)))
   else
      coeff_printed = false
      sc = scale(x.data)
      den = x.scale
      for i = 0:len - 1
         c = polcoeff(x.data, i)
         bracket = needs_parentheses(c)
         if !iszero(c)
            if coeff_printed && !displayed_with_minus_in_front(c)
               print(io, "+")
            end
            if i*sc + valuation(x.data) != 0
               if !isone(c) && (c != -1 || show_minus_one(elem_type(base_ring(x))))
                  if bracket
                     print(io, "(")
                  end
                  print(IOContext(io, :compact => true), c)
                  if bracket
                     print(io, ")")
                  end
                  if i*sc + valuation(x.data) != 0
                     print(io, "*")
                  end
               end
               if c == -1 && !show_minus_one(elem_type(base_ring(x)))
                  print(io, "-")
               end
               print(io, string(var(parent(x.data))))
               if (i*sc + valuation(x.data))//den != 1
                  print(io, "^")
                  q = (valuation(x.data) + i*sc)//den
                  print(IOContext(io, :compact => true), denominator(q) == 1 ? numerator(q) : q)
               end
            else
               print(IOContext(io, :compact => true), c)
            end
            coeff_printed = true
         end
      end
   end
   q = precision(x.data)//x.scale
   print(IOContext(io, :compact => true), "+O(", string(var(parent(x.data))), "^", denominator(q) == 1 ? numerator(q) :
 q, ")")
end

function show(io::IO, a::PuiseuxSeriesRing)
   print(io, "Puiseux series ring in ", var(laurent_ring(a)), " over ")
   print(IOContext(io, :compact => true), base_ring(a))
end

function show(io::IO, a::PuiseuxSeriesField)
   print(io, "Puiseux series field in ", var(laurent_ring(a)), " over ")
   print(IOContext(io, :compact => true), base_ring(a))
end

needs_parentheses(x::PuiseuxSeriesElem) = pol_length(x.data) > 1

displayed_with_minus_in_front(x::PuiseuxSeriesElem) = pol_length(x.data) <= 1 && displayed_with_minus_in_front(polcoeff(x.data, 0))

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
    s = gcd(a.scale, b.scale)
    zscale = div(a.scale*b.scale, s)
    ainf = div(a.scale, s)
    binf = div(b.scale, s)
    z = parent(a)(inflate(a.data, binf) + inflate(b.data, ainf), zscale)
    z = rescale!(z)
    return z
end

function -(a::PuiseuxSeriesElem{T}, b::PuiseuxSeriesElem{T}) where T <: RingElement
    s = gcd(a.scale, b.scale)
    zscale = div(a.scale*b.scale, s)
    ainf = div(a.scale, s)
    binf = div(b.scale, s)
    z = parent(a)(inflate(a.data, binf) - inflate(b.data, ainf), zscale)
    z = rescale!(z)
    return z
end

function *(a::PuiseuxSeriesElem{T}, b::PuiseuxSeriesElem{T}) where T <: RingElement
    s = gcd(a.scale, b.scale)
    zscale = div(a.scale*b.scale, s)
    ainf = div(a.scale, s)
    binf = div(b.scale, s)
    z = parent(a)(inflate(a.data, binf)*inflate(b.data, ainf), zscale)
    z = rescale!(z)
    return z
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(a::PuiseuxSeriesElem{T}, b::T) where T <: RingElem
   z = parent(a)(a.data*b, a.scale)
   z = rescale!(z)
   return z
end

function *(a::PuiseuxSeriesElem, b::Union{Integer, Rational, AbstractFloat})
   z = parent(a)(a.data*b, a.scale)
   z = rescale!(z)
   return z
end

*(a::T, b::PuiseuxSeriesElem{T}) where T <: RingElem = b*a

*(a::Union{Integer, Rational, AbstractFloat}, b::PuiseuxSeriesElem) = b*a

###############################################################################
#
#   Approximation
#
###############################################################################

function isapprox(a::PuiseuxSeriesElem{T}, b::PuiseuxSeriesElem{T}) where T <: RingElement
    s = gcd(a.scale, b.scale)
    zscale = div(a.scale*b.scale, s)
    ainf = div(a.scale, s)
    binf = div(b.scale, s)
    return isapprox(inflate(a.data, binf), inflate(b.data, ainf))
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::PuiseuxSeriesElem{T}, b::PuiseuxSeriesElem{T}) where T <: RingElement
    s = gcd(a.scale, b.scale)
    zscale = div(a.scale*b.scale, s)
    ainf = div(a.scale, s)
    binf = div(b.scale, s)
    z = parent(a)(divexact(inflate(a.data, binf), inflate(b.data, ainf)), zscale)
    z = rescale!(z)
    return z
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

@doc Markdown.doc"""
    divexact(x::Generic.PuiseuxSeriesElem, y::Union{Integer, Rational, AbstractFloat})
> Return $a/b$ where the quotient is expected to be exact.
"""
function divexact(x::PuiseuxSeriesElem, y::Union{Integer, Rational, AbstractFloat})
   return parent(x)(divexact(x.data, y), x.scale)
end

@doc Markdown.doc"""
    divexact(x::Generic.PuiseuxSeriesElem{T}, y::T) where {T <: RingElem}
> Return $a/b$ where the quotient is expected to be exact.
"""
function divexact(x::PuiseuxSeriesElem{T}, y::T) where {T <: RingElem}
   return parent(x)(divexact(x.data, y), x.scale)
end

###############################################################################
#
#   Inversion
#
###############################################################################

function inv(a::PuiseuxSeriesElem{T}) where T <: RingElement
   z = parent(a)(inv(a.data), a.scale)
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
   if iszero(a.data)
      return parent(a)(a.data^b, a.scale)
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
   z = parent(a)(a.data^numerator(b), a.scale*denominator(b))
   z = rescale!(z)
   return z
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::PuiseuxSeriesElem{T}, b::PuiseuxSeriesElem{T}) where T <: RingElement
   fl = check_parent(a, b, false)
   !fl && return false
   s = gcd(a.scale, b.scale)
   zscale = div(a.scale*b.scale, s)
   ainf = div(a.scale, s)
   binf = div(b.scale, s)
   return inflate(a.data, binf) == inflate(b.data, ainf)
end

function isequal(a::PuiseuxSeriesElem{T}, b::PuiseuxSeriesElem{T}) where T <: RingElement
   return a.scale == b.scale && isequal(a.data, b.data)
end

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

==(x::PuiseuxSeriesElem{T}, y::T) where T <: RingElem = x.data == y

==(x::T, y::PuiseuxSeriesElem{T}) where T <: RingElem = y == x

==(x::PuiseuxSeriesElem, y::Union{Integer, Rational, AbstractFloat}) = x.data == y

==(x::Union{Integer, Rational, AbstractFloat}, y::PuiseuxSeriesElem) = y == x

###############################################################################
#
#   Square root
#
###############################################################################

@doc Markdown.doc"""
    sqrt(a::Generic.PuiseuxSeriesElem{T}) where T <: RingElement
> Return the square root of the given Puiseux series $a$.
"""
function Base.sqrt(a::PuiseuxSeriesElem{T}) where T <: RingElement
   val = valuation(a.data)
   S = parent(a)
   if mod(val, 2) != 0
      return S(sqrt(inflate(a.data, 2)), a.scale*2)
   else
      return S(sqrt(a.data), a.scale)
   end
end

###############################################################################
#
#   Exponential
#
###############################################################################

@doc Markdown.doc"""
    exp(a::Generic.PuiseuxSeriesElem{T}) where T <: RingElement
> Return the exponential of the given Puiseux series $a$.
"""
function Base.exp(a::PuiseuxSeriesElem{T}) where T <: RingElement
   z = parent(a)(exp(a.data), a.scale)
   z = rescale!(z)
   return z
end

###############################################################################
#
#   Random elements
#
###############################################################################

function rand(rng::AbstractRNG, S::PuiseuxSeriesRing, val_range::UnitRange{Int}, scale_range::UnitRange{Int}, v...)
   (first(scale_range) <= 0 || last(scale_range) <= 0) && error("Scale must be positive")
   return S(rand(rng, laurent_ring(S), val_range, v...), rand(rng, scale_range))
end

function rand(rng::AbstractRNG, S::PuiseuxSeriesField, val_range::UnitRange{Int}, scale_range::UnitRange{Int}, v...)
   (first(scale_range) <= 0 || last(scale_range) <= 0) && error("Scale must be positive")
   return S(rand(rng, laurent_ring(S), val_range, v...), rand(rng, scale_range))
end

function rand(S::Union{PuiseuxSeriesRing,PuiseuxSeriesField}, val_range, scale_range, v...)
   rand(Random.GLOBAL_RNG, S, val_range, scale_range, v...)
end

###############################################################################
#
#   Unsafe operations
#
###############################################################################

function zero!(a::PuiseuxSeriesElem{T}) where T <: RingElement
   zero!(a.data)
   set_scale(a, 1)
   return a
end

function mul!(c::PuiseuxSeriesElem{T}, a::PuiseuxSeriesElem{T}, b::PuiseuxSeriesElem{T}) where T <: RingElement
    s = gcd(a.scale, b.scale)
    zscale = div(a.scale*b.scale, s)
    ainf = div(a.scale, s)
    binf = div(b.scale, s)
    c.data = mul!(c.data, inflate(a.data, binf), inflate(b.data, ainf))
    c.scale = zscale
    c = rescale!(c)
    return c
end

function add!(c::PuiseuxSeriesElem{T}, a::PuiseuxSeriesElem{T}, b::PuiseuxSeriesElem{T}) where T <: RingElement
    s = gcd(a.scale, b.scale)
    zscale = div(a.scale*b.scale, s)
    ainf = div(a.scale, s)
    binf = div(b.scale, s)
    c.data = add!(c.data, inflate(a.data, binf), inflate(b.data, ainf))
    c.scale = zscale
    c = rescale!(c)
    return c
end

function addeq!(c::PuiseuxSeriesElem{T}, a::PuiseuxSeriesElem{T}) where T <: RingElement
    s = gcd(c.scale, a.scale)
    zscale = div(c.scale*a.scale, s)
    ainf = div(a.scale, s)
    cinf = div(c.scale, s)
    cnew = inflate(c.data, ainf)
    c.data = addeq!(cnew, inflate(a.data, cinf))
    c.scale = zscale
    c = rescale!(c)
    return c
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{PuiseuxSeriesRingElem{T}}, ::Type{PuiseuxSeriesRingElem{T}}) where T <: RingElement = PuiseuxSeriesRingElem{T}

promote_rule(::Type{PuiseuxSeriesFieldElem{T}}, ::Type{PuiseuxSeriesFieldElem{T}}) where T <: FieldElement = PuiseuxSeriesFieldElem{T}

function promote_rule(::Type{PuiseuxSeriesRingElem{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? PuiseuxSeriesRingElem{T} : Union{}
end

function promote_rule(::Type{PuiseuxSeriesFieldElem{T}}, ::Type{U}) where {T <: FieldElement, U <: RingElement}
   promote_rule(T, U) == T ? PuiseuxSeriesFieldElem{T} : Union{}
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
   z = PuiseuxSeriesRingElem{T}(laurent_ring(R)(), 1)
   z.parent = R
   return z
end

function (R::PuiseuxSeriesField{T})() where {T <: FieldElement}
   z = PuiseuxSeriesFieldElem{T}(laurent_ring(R)(), 1)
   z.parent = R
   return z
end

function (R::PuiseuxSeriesRing{T})(b::LaurentSeriesRingElem{T}, scale::Int) where T <: RingElement
   z = PuiseuxSeriesRingElem{T}(b, scale)
   z.parent = R
   z = rescale!(z)
   return z
end

function (R::PuiseuxSeriesField{T})(b::LaurentSeriesFieldElem{T}, scale::Int) where T <: FieldElement
   z = PuiseuxSeriesFieldElem{T}(b, scale)
   z.parent = R
   z = rescale!(z)
   return z
end

function (R::PuiseuxSeriesRing{T})(b::Union{Integer, Rational, AbstractFloat}) where T <: RingElement
   z = PuiseuxSeriesRingElem{T}(laurent_ring(R)(b), 1)
   z.parent = R
   return z
end

function (R::PuiseuxSeriesField{T})(b::Union{Rational, AbstractFloat}) where T <: FieldElement
   z = PuiseuxSeriesFieldElem{T}(laurent_ring(R)(b), 1)
   z.parent = R
   return z
end

function (R::PuiseuxSeriesRing{T})(b::T) where T <: RingElem
   parent(b) != base_ring(R) && error("Unable to coerce to Puiseux series")
   z = PuiseuxSeriesRingElem{T}(laurent_ring(R)(b), 1)
   z.parent = R
   return z
end

function (R::PuiseuxSeriesField{T})(b::T) where T <: FieldElem
   parent(b) != base_ring(R) && error("Unable to coerce to Puiseux series")
   z = PuiseuxSeriesFieldElem{T}(laurent_ring(R)(b), 1)
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

@doc Markdown.doc"""
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
