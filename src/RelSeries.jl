###############################################################################
#
#   RelSeries.jl : Power series over rings, capped relative precision
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

@doc raw"""
    O(a::RelPowerSeriesRingElem{T}) where T <: RingElement

Return $0 + O(x^\mathrm{deg}(a))$. Usually this function is called with $x^n$
as parameter. Then the function returns the power series $0 + O(x^n)$, which
can be used to set the precision of a power series when constructing it.
"""
function O(a::RelPowerSeriesRingElem{T}) where T <: RingElement
   val = pol_length(a) + valuation(a) - 1
   val < 0 && throw(DomainError(a, "pol_length(a) + valuation(a) must be >= 1"))
   return parent(a)(Vector{T}(undef, 0), 0, val, val)
end

parent(a::SeriesElem) = a.parent

base_ring_type(::Type{SeriesRing{T}}) where T <: RingElement = parent_type(T)

base_ring(R::SeriesRing{T}) where T <: RingElement = R.base_ring::parent_type(T)

function is_domain_type(::Type{T}) where {S <: RingElement, T <: SeriesElem{S}}
   return is_domain_type(S)
end

is_exact_type(a::Type{T}) where T <: SeriesElem = false

@doc raw"""
    var(a::SeriesRing)

Return the internal name of the generator of the power series ring. Note that
this is returned as a `Symbol` not a `String`.
"""
var(a::SeriesRing) = a.S

function check_parent(a::SeriesElem, b::SeriesElem, throw::Bool = true)
   b = parent(a) != parent(b)
   b && throw && error("Incompatible power series rings in power series operation")
   return !b
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::RelPowerSeriesRingElem, h::UInt)
   b = 0xb44d6896204881f3%UInt
   for i in 0:pol_length(a) - 1
      b = xor(b, hash(polcoeff(a, i), h), h)
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

@doc raw"""
    pol_length(a::RelPowerSeriesRingElem)

Return the length of the polynomial underlying the given power series. This
will be zero if the power series has no nonzero terms.
"""
pol_length(a::RelPowerSeriesRingElem) = a.length

@doc raw"""
    precision(a::RelPowerSeriesRingElem)

Return the precision of the given power series in absolute terms. This will
be the sum of the valuation and the length of the underlying polynomial.
"""
precision(a::RelPowerSeriesRingElem) = a.prec

@doc raw"""
    valuation(a::RelPowerSeriesRingElem)

Return the valuation of the given power series, i.e. the degree of the first
nonzero term (or the precision if it is arithmetically zero).
"""
valuation(a::RelPowerSeriesRingElem) = a.val

@doc raw"""
    max_precision(R::SeriesRing)

Return the maximum relative precision of power series in the given power
series ring.
"""
max_precision(R::SeriesRing) = R.prec_max

function set_length!(a::SeriesElem, len::Int)
   a.length = len
   return a
end

function set_precision!(a::SeriesElem, prec::Int)
   a.prec = prec
   return a
end

function set_valuation!(a::RelPowerSeriesRingElem, val::Int)
   a.val = val
   return a
end

function coeff(a::RelPowerSeriesRingElem, n::Int)
   if n < valuation(a)
      return base_ring(a)()
   else
      return polcoeff(a, n - valuation(a))
   end
end

zero(R::SeriesRing) = R(0)

one(R::SeriesRing) = R(1)

iszero(a::RelPowerSeriesRingElem) = pol_length(a) == 0

function isone(a::RelPowerSeriesRingElem)
   return valuation(a) == 0 && pol_length(a) == 1 && isone(polcoeff(a, 0))
end

@doc raw"""
    is_gen(a::RelPowerSeriesRingElem)

Return `true` if the given power series is arithmetically equal to the
generator of its power series ring to its current precision, otherwise return
`false`.
"""
function is_gen(a::RelPowerSeriesRingElem)
   return valuation(a) == 1 && pol_length(a) == 1 && isone(polcoeff(a, 0))
end

is_unit(a::RelPowerSeriesRingElem) = valuation(a) == 0 && is_unit(polcoeff(a, 0))

@doc raw"""
    modulus(a::SeriesElem{T}) where {T <: ResElem}

Return the modulus of the coefficients of the given power series.
"""
modulus(a::SeriesElem{T}) where {T <: Union{ResElem, FinFieldElem}} = modulus(base_ring(a))

function renormalize!(z::RelPowerSeriesRingElem)
   i = 0
   zlen = pol_length(z)
   zval = valuation(z)
   zprec = precision(z)
   while i < zlen && iszero(polcoeff(z, i))
      i += 1
   end
   z = set_precision!(z, zprec)
   if i == zlen
      z = set_length!(z, 0)
      z = set_valuation!(z, zprec)
   else
      z = set_valuation!(z, zval + i)
      for j = 1:zlen - i
         z = setcoeff!(z, j - 1, polcoeff(z, j + i - 1))
      end
      z = set_length!(z, zlen - i)
   end
   return nothing
end

###############################################################################
#
#   Similar and zero
#
###############################################################################

function similar(x::RelPowerSeriesRingElem, R::Ring, max_prec::Int,
                                   s::VarName=var(parent(x)); cached::Bool=true)
   TT = elem_type(R)
   V = Vector{TT}(undef, 0)
   p = Generic.RelSeries{TT}(V, 0, max_prec, max_prec)
   # Default similar is supposed to return a Generic series
   if base_ring(x) === R && Symbol(s) == var(parent(x)) &&
            x isa Generic.RelSeries{TT} &&
            max_precision(parent(x)) == max_prec
       # steal parent in case it is not cached
       p.parent = parent(x)
   else
       p.parent = Generic.RelPowerSeriesRing{TT}(R, max_prec, Symbol(s), cached)
   end
   return p
end

similar(x::RelPowerSeriesRingElem, R::Ring,
                                   var::VarName=var(parent(x)); cached::Bool=true) =
   similar(x, R, max_precision(parent(x)), Symbol(var); cached)


similar(x::RelPowerSeriesRingElem, max_prec::Int,
                                   var::VarName=var(parent(x)); cached::Bool=true) =
   similar(x, base_ring(x), max_prec, Symbol(var); cached)


similar(x::RelPowerSeriesRingElem, var::VarName=var(parent(x)); cached::Bool=true) =
   similar(x, base_ring(x), max_precision(parent(x)), Symbol(var); cached)


zero(a::RelPowerSeriesRingElem, R::Ring, max_prec::Int,
                                   var::VarName=var(parent(a)); cached::Bool=true) =
   similar(a, R, max_prec, Symbol(var); cached=cached)


zero(a::RelPowerSeriesRingElem, R::Ring,
                                   var::VarName=var(parent(a)); cached::Bool=true) =
   similar(a, R, Symbol(var); cached)


zero(a::RelPowerSeriesRingElem, max_prec::Int,
                                   var::VarName=var(parent(a)); cached::Bool=true) =
   similar(a, max_prec, Symbol(var); cached)


zero(a::RelPowerSeriesRingElem, var::VarName=var(parent(a)); cached::Bool=true) =
   similar(a, Symbol(var); cached)

###############################################################################
#
#   rel_series constructor
#
###############################################################################

function rel_series(R::Ring, arr::Vector{T}, len::Int, prec::Int, val::Int, var::VarName=:x; max_precision::Int=prec, cached::Bool=true) where T
   prec < len + val && error("Precision too small for given data")
   TT = elem_type(R)
   coeffs = T == Any && length(arr) == 0 ? elem_type(R)[] : map(R, arr)
   p = Generic.RelSeries{TT}(coeffs, len, prec, val)
   # Default is supposed to return a Generic polynomial
   p.parent = Generic.RelPowerSeriesRing{TT}(R, max_precision, Symbol(var), cached)
   return p
end

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function expressify(a::RelPowerSeriesRingElem, x = var(parent(a)); context = nothing)
    sum = Expr(:call, :+)
    v = valuation(a)
    for i in 0:pol_length(a) - 1
        k = i + v
        c = polcoeff(a, i)
        if !iszero(c)
            if k == 0
                xk = 1
            elseif k == 1
                xk = x
            else
                xk = Expr(:call, :^, x, k)
            end
            if isone(c)
                push!(sum.args, Expr(:call, :*, xk))
            else
                push!(sum.args, Expr(:call, :*, expressify(c, context = context), xk))
            end
        end
    end
    push!(sum.args, Expr(:call, :O, Expr(:call, :^, x, precision(a))))
    return sum
end

@enable_all_show_via_expressify SeriesElem

function show(io::IO, ::MIME"text/plain", a::SeriesRing)
  print(io, "Univariate power series ring in ", var(a), " with precision ", a.prec_max)
  println(io)
  io = pretty(io)
  print(io, Indent(), "over ", Lowercase(), base_ring(a))
  print(io, Dedent())
end

function show(io::IO, a::SeriesRing)
  if get(io, :supercompact, false)
    print(io, "Univariate power series ring")
  else
    io = pretty(io)
    print(io, "Univariate power series ring over " )
    print(IOContext(io, :supercompact => true), Lowercase(), base_ring(a))
  end
end
###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::RelPowerSeriesRingElem)
   len = pol_length(a)
   z = parent(a)()
   z = set_precision!(z, precision(a))
   z = set_valuation!(z, valuation(a))
   fit!(z, len)
   for i = 1:len
      z = setcoeff!(z, i - 1, -polcoeff(a, i - 1))
   end
   return z
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::RelPowerSeriesRingElem{T}, b::RelPowerSeriesRingElem{T}) where T <: RingElement
   check_parent(a, b)
   lena = pol_length(a)
   lenb = pol_length(b)
   vala = valuation(a)
   valb = valuation(b)
   valz = min(vala, valb)
   prec = min(precision(a), precision(b))
   mina = min(vala + lena, prec)
   minb = min(valb + lenb, prec)
   lenz = max(mina, minb) - valz
   R = base_ring(a)
   z = parent(a)()
   fit!(z, lenz)
   z = set_precision!(z, prec)
   z = set_valuation!(z, valz)
   if vala >= valb
      for i = 1:min(lenb, vala - valb)
         z = setcoeff!(z, i - 1, deepcopy(polcoeff(b, i - 1)))
      end
      for i = lenb + 1:min(vala - valb, lenz)
         z = setcoeff!(z, i - 1, R())
      end
      for i = vala - valb + 1:lenb
         z = setcoeff!(z, i - 1, polcoeff(a, i - vala + valb - 1) + polcoeff(b, i - 1))
      end
      for i = max(lenb, vala - valb) + 1:lena + vala - valb
         z = setcoeff!(z, i - 1, deepcopy(polcoeff(a, i - vala + valb - 1)))
      end
      for i = lena + vala - valb + 1:lenb
         z = setcoeff!(z, i - 1, deepcopy(polcoeff(b, i - 1)))
      end
   else
      for i = 1:min(lena, valb - vala)
         z = setcoeff!(z, i - 1, deepcopy(polcoeff(a, i - 1)))
      end
      for i = lena + 1:min(valb - vala, lenz)
         z = setcoeff!(z, i - 1, R())
      end
      for i = valb - vala + 1:lena
         z = setcoeff!(z, i - 1, polcoeff(a, i - 1) + polcoeff(b, i - valb + vala - 1))
      end
      for i = max(lena, valb - vala) + 1:lenb + valb - vala
         z = setcoeff!(z, i - 1, deepcopy(polcoeff(b, i - valb + vala - 1)))
      end
      for i = lenb + valb - vala + 1:lena
         z = setcoeff!(z, i - 1, deepcopy(polcoeff(a, i - 1)))
      end
   end
   z = set_length!(z, normalise(z, lenz))
   renormalize!(z)
   return z
end

function -(a::RelPowerSeriesRingElem{T}, b::RelPowerSeriesRingElem{T}) where T <: RingElement
   check_parent(a, b)
   lena = pol_length(a)
   lenb = pol_length(b)
   vala = valuation(a)
   valb = valuation(b)
   valz = min(vala, valb)
   prec = min(precision(a), precision(b))
   mina = min(vala + lena, prec)
   minb = min(valb + lenb, prec)
   lenz = max(mina, minb) - valz
   R = base_ring(a)
   z = parent(a)()
   fit!(z, lenz)
   z = set_precision!(z, prec)
   z = set_valuation!(z, valz)
   if vala >= valb
      for i = 1:min(lenb, vala - valb)
         z = setcoeff!(z, i - 1, -polcoeff(b, i - 1))
      end
      for i = lenb + 1:min(vala - valb, lenz)
         z = setcoeff!(z, i - 1, R())
      end
      for i = vala - valb + 1:lenb
         z = setcoeff!(z, i - 1, polcoeff(a, i - vala + valb - 1) - polcoeff(b, i - 1))
      end
      for i = max(lenb, vala - valb) + 1:lena + vala - valb
         z = setcoeff!(z, i - 1, deepcopy(polcoeff(a, i - vala + valb - 1)))
      end
      for i = lena + vala - valb + 1:lenb
         z = setcoeff!(z, i - 1, -polcoeff(b, i - 1))
      end
   else
      for i = 1:min(lena, valb - vala)
         z = setcoeff!(z, i - 1, deepcopy(polcoeff(a, i - 1)))
      end
      for i = lena + 1:min(valb - vala, lenz)
         z = setcoeff!(z, i - 1, R())
      end
      for i = valb - vala + 1:lena
         z = setcoeff!(z, i - 1, polcoeff(a, i - 1) - polcoeff(b, i - valb + vala - 1))
      end
      for i = max(lena, valb - vala) + 1:lenb + valb - vala
         z = setcoeff!(z, i - 1, -polcoeff(b, i - valb + vala - 1))
      end
      for i = lenb + valb - vala + 1:lena
         z = setcoeff!(z, i - 1, deepcopy(polcoeff(a, i - 1)))
      end
   end
   z = set_length!(z, normalise(z, lenz))
   renormalize!(z)
   return z
end

function *(a::RelPowerSeriesRingElem{T}, b::RelPowerSeriesRingElem{T}) where T <: RingElement
   check_parent(a, b)
   lena = pol_length(a)
   lenb = pol_length(b)
   aval = valuation(a)
   bval = valuation(b)
   zval = aval + bval
   prec = min(precision(a) - aval, precision(b) - bval)
   lena = min(lena, prec)
   lenb = min(lenb, prec)
   if lena == 0 || lenb == 0
      return parent(a)(Vector{T}(undef, 0), 0, prec + zval, zval)
   end
   t = base_ring(a)()
   lenz = min(lena + lenb - 1, prec)
   d = Vector{T}(undef, lenz)
   for i = 1:min(lena, lenz)
      d[i] = polcoeff(a, i - 1)*polcoeff(b, 0)
   end
   if lenz > lena
      for j = 2:min(lenb, lenz - lena + 1)
          d[lena + j - 1] = polcoeff(a, lena - 1)*polcoeff(b, j - 1)
      end
   end
   for i = 1:lena - 1
      if lenz > i
         for j = 2:min(lenb, lenz - i + 1)
            t = mul!(t, polcoeff(a, i - 1), polcoeff(b, j - 1))
            d[i + j - 1] = addeq!(d[i + j - 1], t)
         end
      end
   end
   z = parent(a)(d, lenz, prec + zval, zval)
   z = set_length!(z, normalise(z, lenz))
   renormalize!(z)
   return z
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(a::T, b::RelPowerSeriesRingElem{T}) where {T <: RingElem}
   len = pol_length(b)
   z = parent(b)()
   fit!(z, len)
   z = set_precision!(z, precision(b))
   z = set_valuation!(z, valuation(b))
   for i = 1:len
      z = setcoeff!(z, i - 1, a*polcoeff(b, i - 1))
   end
   z = set_length!(z, normalise(z, len))
   renormalize!(z)
   return z
end

function *(a::Union{Integer, Rational, AbstractFloat}, b::RelPowerSeriesRingElem)
   len = pol_length(b)
   z = parent(b)()
   fit!(z, len)
   z = set_precision!(z, precision(b))
   z = set_valuation!(z, valuation(b))
   for i = 1:len
      z = setcoeff!(z, i - 1, a*polcoeff(b, i - 1))
   end
   z = set_length!(z, normalise(z, len))
   renormalize!(z)
   return z
end

*(a::RelPowerSeriesRingElem{T}, b::T) where {T <: RingElem} = b*a

*(a::RelPowerSeriesRingElem, b::Union{Integer, Rational, AbstractFloat}) = b*a

###############################################################################
#
#   Shifting
#
###############################################################################

@doc raw"""
    shift_left(x::RelPowerSeriesRingElem{T}, n::Int) where T <: RingElement

Return the power series $x$ shifted left by $n$ terms, i.e. multiplied by
$x^n$.
"""
function shift_left(x::RelPowerSeriesRingElem{T}, n::Int) where T <: RingElement
   n < 0 && throw(DomainError(n, "n must be >= 0"))
   xlen = pol_length(x)
   if xlen == 0
      z = zero(parent(x))
      z = set_precision!(z, precision(x) + n)
      z = set_valuation!(z, valuation(x) + n)
      return z
   end
   z = parent(x)()
   fit!(z, xlen)
   z = set_precision!(z, precision(x) + n)
   z = set_valuation!(z, valuation(x) + n)
   for i = 1:xlen
      z = setcoeff!(z, i - 1, polcoeff(x, i - 1))
   end
   return z
end

@doc raw"""
    shift_right(x::RelPowerSeriesRingElem{T}, n::Int) where T <: RingElement

Return the power series $x$ shifted right by $n$ terms, i.e. divided by
$x^n$.
"""
function shift_right(x::RelPowerSeriesRingElem{T}, n::Int) where T <: RingElement
   n < 0 && throw(DomainError(n, "n must be >= 0"))
   xlen = pol_length(x)
   xval = valuation(x)
   xprec = precision(x)
   z = parent(x)()
   if n >= xlen + xval
      z = set_precision!(z, max(0, xprec - n))
      z = set_valuation!(z, max(0, xprec - n))
   else
      zlen = min(xlen + xval - n, xlen)
      fit!(z, zlen)
      z = set_precision!(z, max(0, xprec - n))
      z = set_valuation!(z, max(0, xval - n))
      for i = 1:zlen
         z = setcoeff!(z, i - 1, polcoeff(x, i + xlen  - zlen - 1))
      end
      renormalize!(z)
   end
   return z
end

###############################################################################
#
#   Truncation
#
###############################################################################

@doc raw"""
    truncate(a::RelPowerSeriesRingElem{T}, n::Int) where T <: RingElement

Return $a$ truncated to (absolute) precision $n$.
"""
function truncate(a::RelPowerSeriesRingElem{T}, n::Int) where T <: RingElement
   n < 0 && throw(DomainError(n, "n must be >= 0"))
   alen = pol_length(a)
   aprec = precision(a)
   aval = valuation(a)
   if aprec <= n
      return a
   end
   z = parent(a)()
   z = set_precision!(z, n)
   if n <= aval
      z = set_length!(z, 0)
      z = set_valuation!(z, n)
   else
      fit!(z, n - aval)
      z = set_valuation!(z, aval)
      for i = 1:min(n - aval, alen)
         z = setcoeff!(z, i - 1, polcoeff(a, i - 1))
      end
      z = set_length!(z, normalise(z, n - aval))
   end
   return z
end

# Intended only for internal use, does not renormalize, assumes n >= 0
# Only efficient if valuation(a) == valuation(b) == 0
function mullow(a::RelPowerSeriesRingElem{T}, b::RelPowerSeriesRingElem{T}, n::Int) where T <: RingElement
   lena = pol_length(a)
   lenb = pol_length(b)
   if lena == 0 || lenb == 0
      return zero(parent(a))
   end
   prec = min(precision(a), precision(b))
   t = base_ring(a)()
   lenz = min(lena + lenb - 1, n)
   d = Vector{T}(undef, lenz)
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
   z = parent(a)(d, lenz, prec, 0)
   z = set_length!(z, normalise(z, lenz))
   return z
end

###############################################################################
#
#   Powering
#
###############################################################################

@doc raw"""
    ^(a::RelPowerSeriesRingElem{T}, b::Int) where T <: RingElement

Return $a^b$. We require $b \geq 0$.
"""
function ^(a::RelPowerSeriesRingElem{T}, b::Int) where T <: RingElement
   b < 0 && throw(DomainError(b, "exponent must be >= 0"))
   # special case powers of x for constructing power series efficiently
   if b == 0
      # in fact, the result would be exact 1 if we had exact series
      z = one(parent(a))
      return z
   elseif pol_length(a) == 0
      z = parent(a)()
      z = set_precision!(z, b*valuation(a))
      z = set_valuation!(z, b*valuation(a))
      return z
   elseif is_gen(a)
      z = parent(a)()
      fit!(z, 1)
      z = set_precision!(z, b + precision(a) - 1)
      z = set_valuation!(z, b)
      z = setcoeff!(z, 0, deepcopy(polcoeff(a, 0)))
      z = set_length!(z, 1)
      return z
   elseif pol_length(a) == 1
      c = polcoeff(a, 0)^b
      z = parent(a)(c)
      z = set_precision!(z, (b - 1)*valuation(a) + precision(a))
      z = set_valuation!(z, iszero(c) ? precision(z) : b*valuation(a))
      return z
   elseif b == 1
      return deepcopy(a)
   else
      bit = ~((~UInt(0)) >> 1)
      while (UInt(bit) & b) == 0
         bit >>= 1
      end
      val = valuation(a)
      a = shift_right(a, val)
      prec = precision(a)
      z = a
      bit >>= 1
      while bit !=0
         z = mullow(z, z, prec)
         if (UInt(bit) & b) != 0
            z = mullow(z, a, prec)
         end
         bit >>= 1
      end
      z = set_valuation!(z, b*val)
      z = set_precision!(z, b*val + prec)
      renormalize!(z)
      return z
   end
end

###############################################################################
#
#   Comparison
#
###############################################################################

@doc raw"""
    ==(x::RelPowerSeriesRingElem{T}, y::RelPowerSeriesRingElem{T}) where T <: RingElement

Return `true` if $x == y$ arithmetically, otherwise return `false`. Recall
that power series to different precisions may still be arithmetically
equal to the minimum of the two precisions.
"""
function ==(x::RelPowerSeriesRingElem{T}, y::RelPowerSeriesRingElem{T}) where T <: RingElement
   b = check_parent(x, y, false)
   !b && return false

   xval = valuation(x)
   xprec = precision(x)
   yval = valuation(y)
   yprec = precision(y)
   prec = min(xprec, yprec)
   if prec <= xval && prec <= yval
      return true
   end
   if xval != yval
      return false
   end
   xlen = normalise(x, min(pol_length(x), prec - xval))
   ylen = normalise(y, min(pol_length(y), prec - yval))
   if xlen != ylen
      return false
   end
   for i = 1:xlen
      if polcoeff(x, i - 1) != polcoeff(y, i - 1)
         return false
      end
   end
   return true
end

@doc raw"""
    isequal(x::RelPowerSeriesRingElem{T}, y::RelPowerSeriesRingElem{T}) where T <: RingElement

Return `true` if $x == y$ exactly, otherwise return `false`. Only if the
power series are precisely the same, to the same precision, are they declared
equal by this function.
"""
function isequal(x::RelPowerSeriesRingElem{T}, y::RelPowerSeriesRingElem{T}) where T <: RingElement
   if parent(x) != parent(y)
      return false
   end
   if precision(x) != precision(y) || pol_length(x) != pol_length(y) ||
      valuation(x) != valuation(y)
      return false
   end
   for i = 1:pol_length(x)
      if !isequal(polcoeff(x, i - 1), polcoeff(y, i - 1))
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

@doc raw"""
    ==(x::RelPowerSeriesRingElem{T}, y::T) where {T <: RingElem}

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::RelPowerSeriesRingElem{T}, y::T) where {T <: RingElem} = precision(x) == 0 ||
           ((pol_length(x) == 0 && iszero(y)) || (pol_length(x) == 1 &&
             valuation(x) == 0 && polcoeff(x, 0) == y))

@doc raw"""
    ==(x::T, y::RelPowerSeriesRingElem{T}) where {T <: RingElem}

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::T, y::RelPowerSeriesRingElem{T}) where {T <: RingElem} = y == x

@doc raw"""
    ==(x::RelPowerSeriesRingElem, y::Union{Integer, Rational, AbstractFloat})

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::RelPowerSeriesRingElem, y::Union{Integer, Rational, AbstractFloat}) = precision(x) == 0 ||
                  ((pol_length(x) == 0 && iszero(y)) || (pol_length(x) == 1 &&
                    valuation(x) == 0 && polcoeff(x, 0) == y))

@doc raw"""
    ==(x::Union{Integer, Rational, AbstractFloat}, y::RelPowerSeriesRingElem)

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::Union{Integer, Rational, AbstractFloat}, y::RelPowerSeriesRingElem) = y == x

###############################################################################
#
#   Approximation
#
###############################################################################

function Base.isapprox(f::RelPowerSeriesRingElem, g::RelPowerSeriesRingElem; atol::Real=sqrt(eps()))
   check_parent(f, g)
   nmin = min(precision(f), precision(g))
   i = 1
   while i <= nmin
      if !isapprox(coeff(f, i - 1), coeff(g, i - 1); atol=atol)
         return false
      end
      i += 1
   end
   return true
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(x::RelPowerSeriesRingElem{T}, y::RelPowerSeriesRingElem{T}; check::Bool=true) where T <: RingElement
   check_parent(x, y)
   iszero(y) && throw(DivideError())
   v2 = valuation(y)
   if v2 != 0
      v1 = valuation(x)
      if check && v1 < v2
         error("Not an exact division")
      end
      x = shift_right(x, v2)
      y = shift_right(y, v2)
   else
      x = deepcopy(x)
   end
   y = truncate(y, precision(x))
   res = parent(x)()
   res = set_precision!(res, min(precision(x), valuation(x) + precision(y)))
   res = set_valuation!(res, valuation(x))
   lc = coeff(y, 0)
   check && lc == 0 && error("Not an exact division")
   lenr = precision(x) - valuation(x)
   for i = 0:lenr - 1
      q = divexact(polcoeff(x, i), lc; check=check)
      res = setcoeff!(res, i, q)
      for j = 0:min(precision(y) - 1, lenr - i - 1)
         x = setcoeff!(x, i + j, polcoeff(x, i + j) - polcoeff(y, j)*q)
      end
   end
   res = set_length!(res, normalise(res, pol_length(res)))
   return res
end

function divexact(x::RelPowerSeriesRingElem{T}, y::RelPowerSeriesRingElem{T}; check::Bool=true) where T <: FieldElement
   check_parent(x, y)
   iszero(y) && throw(DivideError())
   v2 = valuation(y)
   if v2 != 0
      v1 = valuation(x)
      if check && v1 < v2
         error("Not an exact division")
      end
      x = shift_right(x, v2)
      y = shift_right(y, v2)
   else
      x = deepcopy(x)
   end
   y = truncate(y, precision(x))
   return x*inv(y)
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(x::RelPowerSeriesRingElem, y::Union{Integer, Rational, AbstractFloat}; check::Bool=true)
   y == 0 && throw(DivideError())
   lenx = pol_length(x)
   z = parent(x)()
   fit!(z, lenx)
   z = set_precision!(z, precision(x))
   z = set_valuation!(z, valuation(x))
   for i = 1:lenx
      z = setcoeff!(z, i - 1, divexact(polcoeff(x, i - 1), y; check=check))
   end
   return z
end

function divexact(x::RelPowerSeriesRingElem{T}, y::T; check::Bool=true) where {T <: RingElem}
   iszero(y) && throw(DivideError())
   lenx = pol_length(x)
   z = parent(x)()
   fit!(z, lenx)
   z = set_precision!(z, precision(x))
   z = set_valuation!(z, valuation(x))
   for i = 1:lenx
      z = setcoeff!(z, i - 1, divexact(polcoeff(x, i - 1), y; check=check))
   end
   return z
end

###############################################################################
#
#   Inversion
#
###############################################################################

@doc raw"""
    Base.inv(a::RelPowerSeriesRingElem)

Return the inverse of the power series $a$, i.e. $1/a$.
"""
function Base.inv(a::RelPowerSeriesRingElem)
   iszero(a) && throw(DivideError())
   !is_unit(a) && error("Unable to invert power series")
   R = base_ring(a)
   a1 = polcoeff(a, 0)
   ainv = parent(a)()
   fit!(ainv, precision(a))
   ainv = set_precision!(ainv, precision(a))
   ainv = set_valuation!(ainv, 0)
   if precision(a) != 0
      ainv = setcoeff!(ainv, 0, divexact(one(R), a1))
   end
   a1 = -a1
   s = R()
   t = R()
   for n = 2:precision(a)
      s = mul_red!(s, polcoeff(a, 1), polcoeff(ainv, n - 2), false)
      for i = 2:min(n, pol_length(a)) - 1
         s = addmul_delayed_reduction!(s, polcoeff(a, i), polcoeff(ainv, n - i - 1), t)
      end
      s = reduce!(s)
      ainv = setcoeff!(ainv, n - 1, divexact(s, a1))
   end
   ainv = set_length!(ainv, normalise(ainv, precision(a)))
   return ainv
end

function Base.inv(a::RelPowerSeriesRingElem{T}) where T <: FieldElement
    prec = precision(a)
    @assert valuation(a) == 0
    @assert prec != 0
    R = parent(a)
    x = R(inv(coeff(a, 0)))
    x = set_precision!(x, 1)
    la = [prec]
    while la[end] > 1
        push!(la, div(la[end] + 1, 2))
    end 
    two = R(2)
    two = set_precision!(two, prec)
    n = length(la) - 1
    y = R()
    minus_a = -a
    while n > 0
        # x -> x*(2 - xa) is the lifting recursion
        x = set_precision!(x, la[n])
        y = set_precision!(y, la[n])
        y = mul!(y, minus_a, x)
        y = addeq!(y, two)
        x = mul!(x, x, y)
        n -= 1 
    end
    return x
end

###############################################################################
#
#   Division with remainder
#
###############################################################################

function Base.divrem(a::RelPowerSeriesRingElem{T}, b::RelPowerSeriesRingElem{T}) where {T <: FieldElement}
   check_parent(a, b)
   if pol_length(b) == 0
      throw(DivideError())
   end
   if valuation(a) < valuation(b)
      return zero(parent(a)), a
   end
   # valuation(a) >= valuation(b), so the exact division works
   q = divexact(a, b)
   return q, a - q*b
end

###############################################################################
#
#   Composition
#
###############################################################################

@doc raw"""
    compose(a::RelPowerSeriesRingElem, b::RelPowerSeriesRingElem)

Compose the series $a$ with the series $b$ and return the result,
i.e. return $a\circ b$. The two series do not need to be in the same ring,
however the series $b$ must have positive valuation or an exception is raised.
"""
function compose(a::RelPowerSeriesRingElem, b::RelPowerSeriesRingElem)
   valuation(b) == 0 && error("Series being substituted must have positive valuation")
   i = pol_length(a)
   R = base_ring(a)
   S = parent(b)
   if i == 0
      return zero(R) + zero(S)
   end
   z = polcoeff(a, i - 1) * one(S)
   while i > 1
      i -= 1
      c = S(polcoeff(a, i - 1))
      z = z*b
      if !iszero(c)
         c = set_precision!(c, precision(z))
         z += c
      end
   end
   z *= b^valuation(a)
   zprec = min(precision(z), valuation(b)*precision(a))
   z = set_precision!(z, zprec)
   z = set_valuation!(z, min(valuation(z), zprec))
   zlen = max(0, precision(z) - valuation(z))
   z = set_length!(z, min(zlen, pol_length(z)))
   return z
end

# General substitution is not well-defined
function subst(a::SeriesElem, b::SeriesElem)
   return compose(a, b)
end

###############################################################################
#
#   Square root
#
###############################################################################

function sqrt_classical_char2(a::RelPowerSeriesRingElem; check::Bool=true)
   S = parent(a)
   R = base_ring(a)
   prec = div(precision(a) + 1, 2)
   if iszero(a)
      asqrt = parent(a)()
      asqrt = set_precision!(asqrt, prec)
      asqrt = set_valuation!(asqrt, prec)
      return true, asqrt
   end
   aval = valuation(a)
   if check && !iseven(aval)
      return false, S()
   end
   aval2 = div(aval, 2)
   asqrt = parent(a)()
   fit!(asqrt, prec)
   asqrt = set_precision!(asqrt, prec)
   asqrt = set_valuation!(asqrt, aval2)
   if check
      for i = 1:2:precision(a) - aval - 1 # series must have even exponents
         if !iszero(polcoeff(a, i))
            return false, S()
         end
      end
   end
   for i = 0:prec - aval2 - 1
      c = polcoeff(a, 2*i)
      if check && !is_square(c)
         return false, S()
      end
      asqrt = setcoeff!(asqrt, i, sqrt(c; check=false))
   end
   asqrt = set_length!(asqrt, normalise(asqrt, prec))
   return true, asqrt
end

function sqrt_classical(a::RelPowerSeriesRingElem; check::Bool=true)
   R = base_ring(a)
   S = parent(a)
   aval = valuation(a)
   if check && !iseven(aval)
      return false, S()
   end
   !is_domain_type(elem_type(R)) && error("Sqrt not implemented over non-integral domains")
   if characteristic(R) == 2
      return sqrt_classical_char2(a, check=check)
   end
   aval2 = div(aval, 2)
   prec = precision(a) - aval
   if prec == 0
      asqrt = parent(a)()
      asqrt = set_precision!(asqrt, aval2)
      asqrt = set_valuation!(asqrt, aval2)
      return true, asqrt
   end
   asqrt = parent(a)()
   fit!(asqrt, prec)
   asqrt = set_precision!(asqrt, prec + aval2)
   asqrt = set_valuation!(asqrt, aval2)
   if prec > 0
      c = polcoeff(a, 0)
      if check && !is_square(c)
         return false, zero(S)
      end
      g = sqrt(c; check=check)
      asqrt = setcoeff!(asqrt, 0, g)
      g2 = g + g
   end
   p = R()
   for n = 1:prec - 1
      c = R()
      for i = 1:div(n - 1, 2)
         j = n - i
         p = mul!(p, polcoeff(asqrt, i), polcoeff(asqrt, j))
         c = addeq!(c, p)
      end
      c *= 2
      if (n % 2) == 0
         i = div(n, 2)
         p = mul!(p, polcoeff(asqrt, i), polcoeff(asqrt, i))
         c = addeq!(c, p)
      end
      c = polcoeff(a, n) - c
      if check
         flag, c = divides(c, g2)
         if !flag
            return false, zero(S)
         end
      else
         c = divexact(c, g2; check=check)
      end
      asqrt = setcoeff!(asqrt, n, c)
   end
   asqrt = set_length!(asqrt, normalise(asqrt, prec))
   return true, asqrt
end

@doc raw"""
    sqrt(a::RelPowerSeriesRingElem)

Return the square root of the power series $a$. By default the function raises
an exception if the input is not a square. If `check=false` this check is
omitted.
"""
function Base.sqrt(a::RelPowerSeriesRingElem; check::Bool=true)
   flag, q = sqrt_classical(a; check=check)
   if check && !flag
      error("Not a square in sqrt")
   end
   return q
end

function is_square(a::RelPowerSeriesRingElem)
   flag, q = sqrt_classical(a; check=true)
   return flag
end

function is_square_with_sqrt(a::RelPowerSeriesRingElem)
   return sqrt_classical(a; check=true)
end

###############################################################################
#
#  Derivative and Integral
#
###############################################################################

@doc raw"""
    derivative(f::RelPowerSeriesRingElem{T})

Return the derivative of the power series $f$.

```
julia> R, x = power_series_ring(QQ, 10, "x")
(Univariate power series ring in x over Rationals, x + O(x^11))

julia> f = 2 + x + 3x^3
2 + x + 3*x^3 + O(x^10)

julia> derivative(f)
1 + 9*x^2 + O(x^9)
```
"""
function derivative(f::RelPowerSeriesRingElem{T}) where T <: RingElement
   g = parent(f)()
   g = set_precision!(g, precision(f) - 1)
   fit!(g, pol_length(f))
   v = valuation(f)
   g = set_valuation!(g, 0)
   if v == 0
      for i = 1:pol_length(f) - 1
         g = setcoeff!(g, i - 1, i*polcoeff(f, i))
      end
   else
      for i = 0:pol_length(f) - 1
         g = setcoeff!(g, i, (i + v)*polcoeff(f, i))
      end
      g = set_valuation!(g, v - 1)
   end  
   g = set_length!(g, normalise(g, pol_length(f)))
   renormalize!(g)
   return g
end

@doc raw"""
    integral(f::RelPowerSeriesRingElem{T})

Return the integral of the power series $f$.

```
julia> R, x = power_series_ring(QQ, 10, "x")
(Univariate power series ring in x over Rationals, x + O(x^11))

julia> f = 2 + x + 3x^3
2 + x + 3*x^3 + O(x^10)

julia> integral(f)
2*x + 1//2*x^2 + 3//4*x^4 + O(x^11)
```
"""
function integral(f::RelPowerSeriesRingElem{T}) where T <: RingElement
   g = parent(f)()
   fit!(g, pol_length(f))
   g = set_precision!(g, precision(f) + 1)
   v = valuation(f)
   g = set_valuation!(g, v + 1)
   for i = 1:pol_length(f)
      c = polcoeff(f, i - 1)
      if !iszero(c)
         g = setcoeff!(g, i - 1, divexact(c, i + v))
      end
   end
   g = set_length!(g, normalise(g, pol_length(f)))
   renormalize!(g)
   return g
end

###############################################################################
#
#   Special functions
#
###############################################################################

@doc raw"""
    log(a::SeriesElem{T}) where T <: FieldElement

Return the logarithm of the power series $a$.
"""
function Base.log(a::SeriesElem{T}) where T <: FieldElement
   @assert valuation(a) == 0 
   if isone(coeff(a, 0))
      return integral(derivative(a)*inv(a))
   else
      # Definition only works if series is monic, so divide through by constant
      c = coeff(a, 0)
      clog = log(c)
      adivc = divexact(a, c)
      return integral(derivative(adivc)*inv(adivc)) + clog
   end
end

@doc raw"""
    exp(a::RelPowerSeriesRingElem)

Return the exponential of the power series $a$.
"""
function Base.exp(a::RelPowerSeriesRingElem{T}) where T <: RingElement
   if iszero(a)
      z = one(parent(a))
      z = set_precision!(z, precision(a))
      return z
   end
   vala = valuation(a)
   preca = precision(a)
   z = parent(a)()
   R = base_ring(a)
   fit!(z, preca)
   z = set_precision!(z, preca)
   z = set_valuation!(z, 0)
   c = vala == 0 ? polcoeff(a, 0) : R()
   z = setcoeff!(z, 0, exp(c))
   len = pol_length(a) + vala
   C = R()
   d = derivative(a)
   vald = valuation(d)
   for k = 1 : preca - 1
      s = R()
      for j = 1 : min(k + 1, len) - 1
         c = j > vald ? polcoeff(d, j - vald - 1) : R()
         s = addmul_delayed_reduction!(s, c, polcoeff(z, k - j), C)
      end
      s = reduce!(s)
      !is_unit(R(k)) && error("Unable to divide in exp")
      z = setcoeff!(z, k, divexact(s, k))
   end
   z = set_length!(z, normalise(z, preca))
   return z
end

function Base.exp(a::RelPowerSeriesRingElem{T}) where T <: FieldElement
   if iszero(a)
      b = one(parent(a))
      b = set_precision!(b, precision(a))
      return b
   end
   R = base_ring(a)
   c = one(R)
   if valuation(a) == 0
      a = deepcopy(a)
      c = exp(coeff(a, 0))
      a = setcoeff!(a, 0, R())
   end
   x = parent(a)([one(R)], 1, min(2, precision(a)), 0)
   prec = precision(a)
   la = [prec]
   while la[end] > 1
      push!(la, div(la[end] + 1, 2))
   end
   one1 = parent(a)([one(R)], 1, 2, 0)
   n = length(la) - 1
   # x -> x*(1 - log(a) + a) is the recursion
   while n > 0
      x = set_precision!(x, la[n])
      one1 = set_precision!(one1, la[n])
      t = -log(x)
      t = addeq!(t, one1)
      t = addeq!(t, a)
      x = mul!(x, x, t)
      n -= 1 
   end
   if !isone(c)
      x *= c
   end
   return x
end

################################################################################
#
#  Map
#
################################################################################

function _make_parent(g, p::RelPowerSeriesRingElem, cached::Bool)
   R = parent(g(zero(base_ring(p))))
   S = parent(p)
   sym = var(S)
   max_prec = max_precision(S)
   return power_series_ring(R, max_prec, sym; cached=cached)[1]
end

function map_coefficients(g, p::RelPowerSeriesRingElem{<:RingElement};
                    cached::Bool = true,
                    parent::Ring = _make_parent(g, p, cached))
   return _map(g, p, parent)
end

function _map(g, p::RelPowerSeriesRingElem, Rx)
   R = base_ring(Rx)
   new_coefficients = elem_type(R)[let c = polcoeff(p, i)
                                     iszero(c) ? zero(R) : R(g(c))
                                   end for i in 0:pol_length(p) - 1]
   res = Rx(new_coefficients, pol_length(p), precision(p), valuation(p))
   res = set_length!(res, normalise(res, pol_length(res)))
   renormalize!(res)
   return res
end

################################################################################
#
#  Change base ring
#
################################################################################

function _change_rel_series_ring(R, Rx, cached)
   P, _ = power_series_ring(R, max_precision(Rx), var(Rx), cached = cached)
   return P
end

function change_base_ring(R::Ring, p::RelPowerSeriesRingElem{T};
                    cached::Bool = true, parent::Ring =
          _change_rel_series_ring(R, parent(p), cached)) where T <: RingElement
   return _map(R, p, parent)
end

###############################################################################
#
#   Random elements
#
###############################################################################

RandomExtensions.maketype(S::SeriesRing, ::AbstractUnitRange{Int}, _) = elem_type(S)

function RandomExtensions.make(S::SeriesRing, val_range::AbstractUnitRange{Int}, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, val_range, vs[1]) # forward to default Make constructor
   else
      Make(S, val_range, make(R, vs...))
   end
end

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make3{<:RingElement, <:SeriesRing, <:AbstractUnitRange{Int}}})
   S, val_range, v = sp[][1:end]
   R = base_ring(S)
   f = S()
   x = gen(S)
   for i = 0:S.prec_max - 1
      f += rand(rng, v)*x^i
   end
   return shift_left(f, rand(rng, val_range))
end

rand(rng::AbstractRNG, S::SeriesRing, val_range::AbstractUnitRange{Int}, v...) =
   rand(rng, make(S, val_range, v...))

rand(S::SeriesRing, val_range, v...) = rand(Random.GLOBAL_RNG, S, val_range, v...)

###############################################################################
#
#   power_series_ring constructor
#
###############################################################################

@doc raw"""
    power_series_ring(R::Ring, prec::Int, s::VarName; cached::Bool=true, model=:capped_relative)

Return a tuple $(S, x)$ consisting of the parent object `S` of a power series
ring over the given base ring and a generator `x` for the power series ring.
The maximum precision of power series in the ring is set to `prec`. If the
model is set to `:capped_relative` this is taken as a maximum relative
precision, and if it is set to `:capped_absolute` this is take to be a
maximum absolute precision. The supplied string `s` specifies the way the
generator of the power series ring will be printed. By default, the parent
object `S` will be cached so that supplying the same base ring, string and
precision in future will return the same parent object and generator. If
caching of the parent object is not required, `cached` can be set to `false`.
"""
power_series_ring(R::Ring, prec::Int, s::VarName; cached::Bool=true, model=:capped_relative) =
   Generic.power_series_ring(R, prec, Symbol(s); cached, model)

function AbsPowerSeriesRing(R::Ring, prec::Int)
   T = elem_type(R)
   return Generic.AbsPowerSeriesRing{T}(R, prec, :x, false)
end

function RelPowerSeriesRing(R::Ring, prec::Int)
   T = elem_type(R)
   return Generic.RelPowerSeriesRing{T}(R, prec, :x, false)
end

