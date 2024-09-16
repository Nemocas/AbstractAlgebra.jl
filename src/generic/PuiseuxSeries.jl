###############################################################################
#
#   PuiseuxSeries.jl : Generic Puiseux series over rings and fields
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

@doc raw"""
    laurent_ring(R::PuiseuxSeriesRing{T}) where T <: RingElement

Return the `LaurentSeriesRing` underlying the given `PuiseuxSeriesRing`.
"""
laurent_ring(R::PuiseuxSeriesRing{T}) where T <: RingElement = R.laurent_ring::LaurentSeriesRing{T}

@doc raw"""
    laurent_ring(R::PuiseuxSeriesField{T}) where T <: FieldElement

Return the `LaurentSeriesField` underlying the given `PuiseuxSeriesField`.
"""
laurent_ring(R::PuiseuxSeriesField{T}) where T <: FieldElement = R.laurent_ring::LaurentSeriesField{T}

@doc raw"""
    O(a::Generic.PuiseuxSeriesElem{T}) where T <: RingElement

Return $0 + O(x^\mathrm{val}(a))$. Usually this function is called with $x^n$
as parameter for some rational $n$. Then the function returns the Puiseux series
$0 + O(x^n)$, which can be used to set the precision of a Puiseux series when
constructing it.
"""
function O(a::PuiseuxSeriesElem{T}) where T <: RingElement
   val = valuation(a)
   par = parent(a)
   laur = laurent_ring(par)(Vector{T}(undef, 0), 0, numerator(val), numerator(val), 1)
   return parent(a)(laur, denominator(val))
end

parent_type(::Type{PuiseuxSeriesRingElem{T}}) where T <: RingElement = PuiseuxSeriesRing{T}

parent_type(::Type{PuiseuxSeriesFieldElem{T}}) where T <: FieldElement = PuiseuxSeriesField{T}

parent(a::PuiseuxSeriesElem) = a.parent

elem_type(::Type{PuiseuxSeriesRing{T}}) where T <: RingElement = PuiseuxSeriesRingElem{T}

elem_type(::Type{PuiseuxSeriesField{T}}) where T <: FieldElement = PuiseuxSeriesFieldElem{T}

base_ring_type(::Type{PuiseuxSeriesRing{T}}) where T <: RingElement = parent_type(T)

base_ring_type(::Type{PuiseuxSeriesField{T}}) where T <: RingElement = parent_type(T)

base_ring(R::PuiseuxSeriesRing{T}) where T <: RingElement = base_ring(laurent_ring(R))

base_ring(R::PuiseuxSeriesField{T}) where T <: FieldElement = base_ring(laurent_ring(R))

@doc raw"""
    max_precision(R::PuiseuxSeriesRing{T}) where T <: RingElement

Return the maximum precision of the underlying Laurent series ring.
"""
max_precision(R::PuiseuxSeriesRing{T}) where T <: RingElement = max_precision(laurent_ring(R))

@doc raw"""
    max_precision(R::PuiseuxSeriesField{T}) where T <: FieldElement

Return the maximum precision of the underlying Laurent series field.
"""
max_precision(R::PuiseuxSeriesField{T}) where T <: FieldElement = max_precision(laurent_ring(R))

@doc raw"""
    var(R::PuiseuxSeriesRing{T}) where T <: RingElement

Return a symbol representing the variable of the given Puiseux series ring.
"""
var(R::PuiseuxSeriesRing{T}) where T <: RingElement = var(laurent_ring(R))

@doc raw"""
    var(R::PuiseuxSeriesField{T}) where T <: FieldElement

Return a symbol representing the variable of the given Puiseux series field.
"""
var(R::PuiseuxSeriesField{T}) where T <: FieldElement = var(laurent_ring(R))

function is_domain_type(::Type{T}) where {S <: RingElement, T <: PuiseuxSeriesElem{S}}
   return is_domain_type(S)
end

is_exact_type(a::Type{T}) where T <: PuiseuxSeriesElem = false

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

@doc raw"""
    precision(a::Generic.PuiseuxSeriesElem)

Return the precision of the given Puiseux series in absolute terms.
"""
precision(a::PuiseuxSeriesElem) = precision(a.data)//a.scale

@doc raw"""
    valuation(a::Generic.PuiseuxSeriesElem)

Return the valuation of the given Puiseux series, i.e. the exponent of the first
nonzero term (or the precision if it is arithmetically zero).
"""
valuation(a::PuiseuxSeriesElem) = valuation(a.data)//a.scale

scale(a::PuiseuxSeriesElem) = a.scale

function set_precision!(a::PuiseuxSeriesElem, prec::Rational{Int})
   s = scale(a)
   n = numerator(prec)
   d = denominator(prec)
   sa = lcm(s, d)
   a.data = inflate(a.data, div(sa, s))
   a.data = set_precision!(a.data, n*div(sa, d))
   a.scale = sa
   return rescale!(a)
end

set_precision!(a::PuiseuxSeriesElem, prec::Int) = set_precision!(a, prec//1)

function set_valuation!(a::PuiseuxSeriesElem, val::Rational{Int})
   s = scale(a)
   n = numerator(val)
   d = denominator(val)
   sa = lcm(s, d)
   a.data = inflate(a.data, div(sa, s))
   a.data = set_valuation!(a.data, n*div(sa, d))
   a.scale = sa
   return rescale!(a)
end

set_valuation!(a::PuiseuxSeriesElem, val::Int) = set_valuation!(a, val//1)

@doc raw"""
    coeff(a::Generic.PuiseuxSeriesElem, n::Int)

Return the coefficient of the term of exponent $n$ of the given Puiseux series.
"""
function coeff(a::PuiseuxSeriesElem, n::Int)
   s = scale(a)
   return coeff(a.data, n*s)
end

@doc raw"""
    coeff(a::Generic.PuiseuxSeriesElem, r::Rational{Int})

Return the coefficient of the term of exponent $r$ of the given Puiseux series.
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

zero(R::PuiseuxSeriesRing) = R(0)

zero(R::PuiseuxSeriesField) = R(0)

one(R::PuiseuxSeriesField) = R(1)

one(R::PuiseuxSeriesRing) = R(1)

@doc raw"""
    gen(R::PuiseuxSeriesRing)

Return the generator of the Puiseux series ring, i.e. $x + O(x^{n + 1})$ where
$n$ is the maximum precision of the Puiseux series ring $R$.
"""
function gen(R::PuiseuxSeriesRing)
   S = laurent_ring(R)
   return R(gen(S), 1)
end

@doc raw"""
    gen(R::PuiseuxSeriesField)

Return the generator of the Puiseux series ring, i.e. $x + O(x^{n + 1})$ where
$n$ is the maximum precision of the Puiseux series ring $R$.
"""
function gen(R::PuiseuxSeriesField)
   S = laurent_ring(R)
   return R(gen(S), 1)
end

iszero(a::PuiseuxSeriesElem) = iszero(a.data)

function isone(a::PuiseuxSeriesElem)
   return isone(a.data)
end

@doc raw"""
    is_gen(a::Generic.PuiseuxSeriesElem)

Return `true` if the given Puiseux series is arithmetically equal to the
generator of its Puiseux series ring to its current precision, otherwise return
`false`.
"""
function is_gen(a::PuiseuxSeriesElem)
   return valuation(a) == 1 && pol_length(a.data) == 1 && isone(polcoeff(a.data, 0))
end

is_unit(a::PuiseuxSeriesElem) = valuation(a) == 0 && is_unit(polcoeff(a.data, 0))

@doc raw"""
    modulus(a::Generic.PuiseuxSeriesElem{T}) where {T <: ResElem}

Return the modulus of the coefficients of the given Puiseux series.
"""
modulus(a::PuiseuxSeriesElem{T}) where {T <: ResElem} = modulus(base_ring(a))

@doc raw"""
    rescale!(a::Generic.PuiseuxSeriesElem)

Rescale so that the scale of the given Puiseux series and the scale of the underlying
Laurent series are coprime. This function is used internally, as all user facing
functions are assumed to rescale their output.
"""
function rescale!(a::PuiseuxSeriesElem)
   if !iszero(a)
      d = gcd(a.scale, gcd(scale(a.data), gcd(valuation(a.data), precision(a.data))))
      if d != 1
         a.data = set_scale!(a.data, div(scale(a.data), d))
         a.data = set_precision!(a.data, AbstractAlgebra.div(precision(a.data), d))
         a.data = set_valuation!(a.data, AbstractAlgebra.div(valuation(a.data), d))
         a.scale = div(a.scale, d)
      end
   else
      d = gcd(precision(a.data), a.scale)
      if d != 1
         a.data = set_precision!(a.data, AbstractAlgebra.div(precision(a.data), d))
         a.data = set_valuation!(a.data, AbstractAlgebra.div(valuation(a.data), d))
         a.scale = div(a.scale, d)
      end
   end
   return a
end

function deepcopy_internal(a::PuiseuxSeriesElem{T}, dict::IdDict) where {T <: RingElement}
    return parent(a)(deepcopy_internal(a.data, dict), a.scale)
end

function characteristic(a::PuiseuxSeriesRing{T}) where T <: RingElement
   return characteristic(base_ring(a))
end

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function AbstractAlgebra.expressify(a::PuiseuxSeriesElem,
                                    x = var(parent(a)); context = nothing)
    b = a.data
    v = valuation(b)
    len = pol_length(b)
    den = a.scale
    sc = scale(b)

    sum = Expr(:call, :+)

    for i in 0:len - 1
        c = polcoeff(b, i)
        expo = (i * sc + v)//den
        if !iszero(c)
            if iszero(expo)
                xk = 1
            elseif isone(expo)
                xk = x
            else
                xk = Expr(:call, :^, x, expressify(expo, context = context))
            end
            if isone(c)
                push!(sum.args, Expr(:call, :*, xk))
            else
                push!(sum.args, Expr(:call, :*, expressify(c, context = context), xk))
            end
        end
    end
    expo = precision(b)//den
    push!(sum.args, Expr(:call, :O, Expr(:call, :^, x, expressify(expo; context))))
    return sum
end

function Base.show(io::IO, ::MIME"text/plain", a::PuiseuxSeriesElem)
  print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

function Base.show(io::IO, a::PuiseuxSeriesElem)
  print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

function show(io::IO, p::PuiseuxSeriesRing)
   @show_name(io, p)
   @show_special(io, p)
   if is_terse(io)
      print(io, "Puiseux series ring")
   else
      io = pretty(io)
      print(io, "Puiseux series ring in ", var(laurent_ring(p)), " over ")
      print(terse(io), Lowercase(), base_ring(p))
   end
end

function show(io::IO, p::PuiseuxSeriesField)
   @show_name(io, p)
   @show_special(io, p)
   if is_terse(io)
      print(io, "Puiseux series field")
   else
      io = pretty(io)
      print(io, "Puiseux series field in ", var(laurent_ring(p)), " over ")
      print(terse(io), Lowercase(), base_ring(p))
   end
end

###############################################################################
#
#   Map coefficients
#
###############################################################################

function _make_parent(g::T, p::PuiseuxSeriesElem, cached::Bool) where T
   R = parent(g(zero(base_ring(p))))
   S = parent(p)
   sym = var(S)
   max_prec = max_precision(S)
   return AbstractAlgebra.puiseux_series_ring(R, max_prec, sym; cached=cached)[1]
end

function map_coefficients(g::T, p::PuiseuxSeriesElem{<:RingElement};
                    cached::Bool = true,
                    parent::Ring = _make_parent(g, p, cached)) where T
   return _map(g, p, parent)
end

function _map(g::T, p::PuiseuxSeriesElem, Rx) where T
   R = base_ring(Rx)
   res = Rx(map_coefficients(g, p.data), scale(p))
   res = rescale!(res)
   return res
end

################################################################################
#
#  Change base ring
#
################################################################################

function _change_puiseux_series_ring(R, Rx, cached)
   P, _ = AbstractAlgebra.puiseux_series_ring(R, max_precision(Rx),
                                               var(Rx), cached = cached)
   return P
end

function change_base_ring(R::Ring, p::PuiseuxSeriesElem{T};
                    cached::Bool = true, parent::Ring =
          _change_puiseux_series_ring(R, parent(p), cached)) where T <: RingElement
   return _map(R, p, parent)
end

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

function divexact(a::PuiseuxSeriesElem{T}, b::PuiseuxSeriesElem{T}; check::Bool=true) where T <: RingElement
    s = gcd(a.scale, b.scale)
    zscale = div(a.scale*b.scale, s)
    ainf = div(a.scale, s)
    binf = div(b.scale, s)
    z = parent(a)(divexact(inflate(a.data, binf), inflate(b.data, ainf); check=check), zscale)
    z = rescale!(z)
    return z
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(x::PuiseuxSeriesElem, y::Union{Integer, Rational, AbstractFloat}; check::Bool=true)
   return parent(x)(divexact(x.data, y; check=check), x.scale)
end

function divexact(x::PuiseuxSeriesElem{T}, y::T; check::Bool=true) where {T <: RingElem}
   return parent(x)(divexact(x.data, y; check=check), x.scale)
end

###############################################################################
#
#   Inversion
#
###############################################################################

@doc raw"""
    Base.inv(a::PuiseuxSeriesElem{T}) where T <: RingElement

Return the inverse of the power series $a$, i.e. $1/a$, if it exists.
Otherwise an exception is raised.
"""
function Base.inv(a::PuiseuxSeriesElem{T}) where T <: RingElement
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
   if b == 0
      # in fact, the result would be exact 1 if we had exact series
      return one(parent(a))
   elseif iszero(a.data)
      return parent(a)(a.data^b, a.scale)
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
   (pol_length(a.data) != 1 || !isone(polcoeff(a.data, 0))) && error("Rational power not implemented")
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

@doc raw"""
    sqrt(a::Generic.PuiseuxSeriesElem{T}; check::Bool=true) where T <: RingElement

Return the square root of the given Puiseux series $a$. By default the function
will throw an exception if the input is not square. If `check=false` this test
is omitted.
"""
function sqrt_classical(a::PuiseuxSeriesElem{T}; check::Bool=true) where T <: RingElement
   val = valuation(a.data)
   S = parent(a)
   R = base_ring(S)
   if mod(val, 2) != 0 || (characteristic(R) == 2 && isodd(scale(a.data)))
      d = inflate(a.data, 2)
      sscale = a.scale*2
   else
      d = a.data
      sscale = a.scale
   end
   flag, s = sqrt_classical(d; check=check)
   if check && !flag
      return false, zero(S)
   end
   return true, S(s, sscale)
end

@doc raw"""
    sqrt(a::Generic.PuiseuxSeriesElem{T}; check::Bool=true) where T <: RingElement

Return the square root of the given Puiseux series $a$. By default the function
will throw an exception if the input is not square. If `check=false` this test
is omitted.
"""
function Base.sqrt(a::PuiseuxSeriesElem{T}; check::Bool=true) where T <: RingElement
   flag, s = sqrt_classical(a; check=check)
   check && !flag && error("Not a square in sqrt")
   return s
end

function is_square(a::PuiseuxSeriesElem{T}) where T <: RingElement
   flag, s = sqrt_classical(a; check=true)
   return flag
end

function is_square_with_sqrt(a::PuiseuxSeriesElem{T}) where T <: RingElement
   return sqrt_classical(a; check=true)
end

########################################################
#
#   Derivative and integral
#
###############################################################################

@doc raw"""
    derivative(a::Generic.PuiseuxSeriesElem{T}) where T <: RingElement

Return the derivative of the given Puiseux series $a$.
"""
function derivative(a::PuiseuxSeriesElem{T}) where T <: RingElement
   S = parent(a)
   s = scale(a)
   z = derivative(a.data)
   z = set_valuation!(z, valuation(z) - s + 1)
   z = set_precision!(z, precision(z) - s + 1)
   r = divexact(S(z, s), s)
   return rescale!(r)
end

@doc raw"""
    integral(a::Generic.PuiseuxSeriesElem{T}) where T <: RingElement

Return the integral of the given Puiseux series $a$.
"""
function integral(a::PuiseuxSeriesElem{T}) where T <: RingElement
   S = parent(a)
   s = scale(a)
   z = s*a.data
   z = set_valuation!(z, valuation(z) + s - 1)
   z = set_precision!(z, precision(z) + s - 1)
   r = S(integral(z), s)
   return rescale!(r)
end

###############################################################################
#
#   Exponential
#
###############################################################################

@doc raw"""
    exp(a::Generic.PuiseuxSeriesElem{T}) where T <: RingElement

Return the exponential of the given Puiseux series $a$.
"""
function Base.exp(a::PuiseuxSeriesElem{T}) where T <: RingElement
   z = parent(a)(exp(a.data), a.scale)
   z = rescale!(z)
   return z
end

@doc raw"""
    log(a::Generic.PuiseuxSeriesElem{T}) where T <: RingElement

Return the logarithm of the given Puiseux series $a$.
"""
function Base.log(a::PuiseuxSeriesElem{T}) where T <: RingElement
   z = parent(a)(log(a.data), a.scale)
   z = rescale!(z)
   return z
end

###############################################################################
#
#   Random elements
#
###############################################################################

const PuiseuxSeriesRingOrField = Union{PuiseuxSeriesRing,PuiseuxSeriesField}

RandomExtensions.maketype(S::PuiseuxSeriesRingOrField, ::AbstractUnitRange{Int}, _) = elem_type(S)

RandomExtensions.make(S::PuiseuxSeriesRingOrField, val_range::AbstractUnitRange{Int},
                      scale_range::AbstractUnitRange{Int}, vs...) =
     make(S, scale_range, make(laurent_ring(S), val_range, vs...))

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make3{<:RingElement,
                                         <:PuiseuxSeriesRingOrField,
                                         <:AbstractUnitRange{Int}}})
   S, scale_range, v = sp[][1:end]
   (first(scale_range) <= 0 || last(scale_range) <= 0) && error("Scale must be positive")
   return S(rand(rng, v), rand(rng, scale_range))
end

rand(rng::AbstractRNG, S::PuiseuxSeriesRingOrField, val_range::AbstractUnitRange{Int},
     scale_range::AbstractUnitRange{Int}, v...) =
        rand(rng, make(S, val_range, scale_range, v...))

rand(S::PuiseuxSeriesRingOrField, val_range, scale_range, v...) =
   rand(Random.GLOBAL_RNG, S, val_range, scale_range, v...)


###############################################################################
#
#   Unsafe operations
#
###############################################################################

function zero!(a::PuiseuxSeriesElem{T}) where T <: RingElement
   zero!(a.data)
   a.scale = 1
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


function PuiseuxSeriesRing(R::AbstractAlgebra.Ring, prec::Int, s::Symbol; cached::Bool=true)
   S, x = AbstractAlgebra.laurent_series_ring(R, prec, s; cached=cached)
   T = elem_type(R)

   parent_obj = PuiseuxSeriesRing{T}(S, cached)

   return parent_obj, gen(parent_obj)
end

function PuiseuxSeriesRing(R::AbstractAlgebra.Field, prec::Int, s::Symbol; cached::Bool=true)
   S, x = AbstractAlgebra.laurent_series_field(R, prec, s; cached=cached)
   T = elem_type(R)

   parent_obj = PuiseuxSeriesField{T}(S, cached)

   return parent_obj, gen(parent_obj)
end

function PuiseuxSeriesField(R::AbstractAlgebra.Field, prec::Int, s::Symbol; cached::Bool=true)
   S, x = AbstractAlgebra.laurent_series_field(R, prec, s; cached=cached)
   T = elem_type(R)

   parent_obj = PuiseuxSeriesField{T}(S, cached)

   return parent_obj, gen(parent_obj)
end
