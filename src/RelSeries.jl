###############################################################################
#
#   RelSeries.jl : Power series over rings, capped relative precision
#
###############################################################################

export PowerSeriesRing, coeff, polcoeff, rel_series, rel_series_type,
       renormalize!, set_length!, set_precision!, set_valuation!

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

@doc Markdown.doc"""
    O(a::RelSeriesElem{T}) where T <: RingElement

Return $0 + O(x^\mathrm{deg}(a))$. Usually this function is called with $x^n$
as parameter. Then the function returns the power series $0 + O(x^n)$, which
can be used to set the precision of a power series when constructing it.
"""
function O(a::RelSeriesElem{T}) where T <: RingElement
   val = pol_length(a) + valuation(a) - 1
   val < 0 && throw(DomainError(a, "pol_length(a) + valuation(a) must be >= 1"))
   return parent(a)(Array{T}(undef, 0), 0, val, val)
end

parent(a::SeriesElem) = a.parent

base_ring(R::SeriesRing{T}) where T <: RingElement = R.base_ring::parent_type(T)

base_ring(a::SeriesElem) = base_ring(parent(a))

function isdomain_type(::Type{T}) where {S <: RingElement, T <: SeriesElem{S}}
   return isdomain_type(S)
end

isexact_type(a::Type{T}) where T <: SeriesElem = false

@doc Markdown.doc"""
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

function Base.hash(a::RelSeriesElem, h::UInt)
   b = 0xb44d6896204881f3%UInt
   for i in 0:pol_length(a) - 1
      b = xor(b, hash(polcoeff(a, i), h), h)
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   return b
end

@doc Markdown.doc"""
    pol_length(a::RelSeriesElem)

Return the length of the polynomial underlying the given power series. This
will be zero if the power series has no nonzero terms.
"""
pol_length(a::RelSeriesElem) = a.length

@doc Markdown.doc"""
    precision(a::RelSeriesElem)

Return the precision of the given power series in absolute terms. This will
be the sum of the valuation and the length of the underlying polynomial.
"""
precision(a::RelSeriesElem) = a.prec

@doc Markdown.doc"""
    valuation(a::RelSeriesElem)

Return the valuation of the given power series, i.e. the degree of the first
nonzero term (or the precision if it is arithmetically zero).
"""
valuation(a::RelSeriesElem) = a.val

@doc Markdown.doc"""
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

function set_valuation!(a::RelSeriesElem, val::Int)
   a.val = val
   return a
end

function coeff(a::RelSeriesElem, n::Int)
   if n < valuation(a)
      return base_ring(a)()
   else
      return polcoeff(a, n - valuation(a))
   end
end

zero(R::SeriesRing) = R(0)

one(R::SeriesRing) = R(1)

iszero(a::RelSeriesElem) = pol_length(a) == 0

function isone(a::RelSeriesElem)
   return valuation(a) == 0 && pol_length(a) == 1 && isone(polcoeff(a, 0))
end

@doc Markdown.doc"""
    isgen(a::RelSeriesElem)

Return `true` if the given power series is arithmetically equal to the
generator of its power series ring to its current precision, otherwise return
`false`.
"""
function isgen(a::RelSeriesElem)
   return valuation(a) == 1 && pol_length(a) == 1 && isone(polcoeff(a, 0))
end

isunit(a::RelSeriesElem) = valuation(a) == 0 && isunit(polcoeff(a, 0))

@doc Markdown.doc"""
    modulus(a::SeriesElem{T}) where {T <: ResElem}

Return the modulus of the coefficients of the given power series.
"""
modulus(a::SeriesElem{T}) where {T <: Union{ResElem, FinFieldElem}} = modulus(base_ring(a))

function renormalize!(z::RelSeriesElem)
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

function similar(x::RelSeriesElem, R::Ring, max_prec::Int,
                                 var::Symbol=var(parent(x)); cached::Bool=true)
   TT = elem_type(R)
   V = Vector{TT}(undef, 0)
   p = Generic.RelSeries{TT}(V, 0, max_prec, max_prec)
   # Default similar is supposed to return a Generic series
   p.parent = Generic.RelSeriesRing{TT}(R, max_prec, var, cached)
   return p
end

function similar(x::RelSeriesElem, R::Ring,
                                 var::Symbol=var(parent(x)); cached::Bool=true)
   return similar(x, R, max_precision(parent(x)), var; cached = cached)
end

function similar(x::RelSeriesElem, max_prec::Int,
                                 var::Symbol=var(parent(x)); cached::Bool=true)
   return similar(x, base_ring(x), max_prec, var; cached=cached)
end

function similar(x::RelSeriesElem,
                                 var::Symbol=var(parent(x)); cached::Bool=true)
   return similar(x, base_ring(x),
                  max_precision(parent(x)), var; cached=cached)
end

function similar(x::RelSeriesElem, R::Ring, max_prec::Int,
                                                var::String; cached::Bool=true)
   return similar(x, R, max_prec, Symbol(var); cached=cached)
end

function similar(x::RelSeriesElem, R::Ring, var::String; cached::Bool=true)
   return similar(x, R, max_precision(parent(x)), Symbol(var); cached=cached)
end

function similar(x::RelSeriesElem, max_prec::Int,
                                                var::String; cached::Bool=true)
   return similar(x, base_ring(x), max_prec, Symbol(var); cached=cached)
end

function similar(x::RelSeriesElem, var::String; cached::Bool=true)
   return similar(x, base_ring(x), max_precision(parent(x)),
                  Symbol(var); cached=cached)
end

function zero(a::RelSeriesElem, R::Ring, max_prec::Int,
                                 var::Symbol=var(parent(a)); cached::Bool=true)
   return similar(a, R, max_prec, var; cached=cached)
end

function zero(a::RelSeriesElem, R::Ring,
                                 var::Symbol=var(parent(a)); cached::Bool=true)
   return similar(a, R, max_precision(parent(a)), var; cached=cached)
end

function zero(a::RelSeriesElem, max_prec::Int,
                                 var::Symbol=var(parent(a)); cached::Bool=true)
   return similar(a, base_ring(a), max_prec, var; cached=cached)
end

zero(a::RelSeriesElem, var::Symbol=var(parent(a)); cached::Bool=true) =
   similar(a, base_ring(a), max_precision(parent(a)), var; cached=cached)

function zero(a::RelSeriesElem, R::Ring, max_prec::Int,
                                                var::String; cached::Bool=true)
   return zero(a, R, max_prec, Symbol(var); cached=cached)
end

zero(a::RelSeriesElem, R::Ring, var::String; cached::Bool=true) =
   zero(a, R, max_precision(parent(a)), Symbol(var); cached=cached)

zero(a::RelSeriesElem, max_prec::Int, var::String; cached::Bool=true) =
   zero(a, base_ring(a), max_prec, Symbol(var); cached=cached)

zero(a::RelSeriesElem, var::String; cached::Bool=true) =
   zero(a, base_ring(a), max_precision(parent(a)), Symbol(var); cached=cached)

###############################################################################
#
#   rel_series constructor
#
###############################################################################

function rel_series(R::Ring, arr::Vector{T}, len::Int, prec::Int, val::Int, var::AbstractString="x"; max_precision::Int=prec, cached::Bool=true) where T
   prec < len + val && error("Precision too small for given data")
   TT = elem_type(R)
   coeffs = T == Any && length(arr) == 0 ? elem_type(R)[] : map(R, arr)
   p = Generic.RelSeries{TT}(coeffs, len, prec, val)
   # Default is supposed to return a Generic polynomial
   p.parent = Generic.RelSeriesRing{TT}(R, max_precision, Symbol(var), cached)
   return p
end

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function expressify(a::RelSeriesElem,
                                    x = var(parent(a)); context = nothing)
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

function Base.show(io::IO, a::SeriesElem)
  print(io, obj_to_string(a, context = io))
end

function Base.show(io::IO, ::MIME"text/plain", a::SeriesElem)
  print(io, obj_to_string(a, context = io))
end

function show(io::IO, a::SeriesRing)
   print(io, "Univariate power series ring in ", var(a), " over ")
   print(IOContext(io, :compact => true), base_ring(a))
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::RelSeriesElem)
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

function +(a::RelSeriesElem{T}, b::RelSeriesElem{T}) where T <: RingElement
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

function -(a::RelSeriesElem{T}, b::RelSeriesElem{T}) where T <: RingElement
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

function *(a::RelSeriesElem{T}, b::RelSeriesElem{T}) where T <: RingElement
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
      return parent(a)(Array{T}(undef, 0), 0, prec + zval, zval)
   end
   t = base_ring(a)()
   lenz = min(lena + lenb - 1, prec)
   d = Array{T}(undef, lenz)
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

function *(a::T, b::RelSeriesElem{T}) where {T <: RingElem}
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

function *(a::Union{Integer, Rational, AbstractFloat}, b::RelSeriesElem)
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

*(a::RelSeriesElem{T}, b::T) where {T <: RingElem} = b*a

*(a::RelSeriesElem, b::Union{Integer, Rational, AbstractFloat}) = b*a

###############################################################################
#
#   Shifting
#
###############################################################################

@doc Markdown.doc"""
    shift_left(x::RelSeriesElem{T}, n::Int) where T <: RingElement

Return the power series $x$ shifted left by $n$ terms, i.e. multiplied by
$x^n$.
"""
function shift_left(x::RelSeriesElem{T}, n::Int) where T <: RingElement
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

@doc Markdown.doc"""
    shift_right(x::RelSeriesElem{T}, n::Int) where T <: RingElement

Return the power series $x$ shifted right by $n$ terms, i.e. divided by
$x^n$.
"""
function shift_right(x::RelSeriesElem{T}, n::Int) where T <: RingElement
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

@doc Markdown.doc"""
    truncate(a::RelSeriesElem{T}, n::Int) where T <: RingElement

Return $a$ truncated to (absolute) precision $n$.
"""
function truncate(a::RelSeriesElem{T}, n::Int) where T <: RingElement
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
function mullow(a::RelSeriesElem{T}, b::RelSeriesElem{T}, n::Int) where T <: RingElement
   lena = pol_length(a)
   lenb = pol_length(b)
   if lena == 0 || lenb == 0
      return zero(parent(a))
   end
   prec = min(precision(a), precision(b))
   t = base_ring(a)()
   lenz = min(lena + lenb - 1, n)
   d = Array{T}(undef, lenz)
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

@doc Markdown.doc"""
    ^(a::RelSeriesElem{T}, b::Int) where T <: RingElement

Return $a^b$. We require $b \geq 0$.
"""
function ^(a::RelSeriesElem{T}, b::Int) where T <: RingElement
   b < 0 && throw(DomainError(b, "exponent must be >= 0"))
   # special case powers of x for constructing power series efficiently
   if pol_length(a) == 0
      z = parent(a)()
      z = set_precision!(z, b*valuation(a))
      z = set_valuation!(z, b*valuation(a))
      return z
   elseif b == 0
      # in fact, the result would be exact 1 if we had exact series
      z = one(parent(a))
      return z
   elseif isgen(a)
      z = parent(a)()
      fit!(z, 1)
      z = set_precision!(z, b + precision(a) - 1)
      z = set_valuation!(z, b)
      z = setcoeff!(z, 0, deepcopy(polcoeff(a, 0)))
      z = set_length!(z, 1)
      return z
   elseif pol_length(a) == 1
      z = parent(a)(polcoeff(a, 0)^b)
      z = set_precision!(z, (b - 1)*valuation(a) + precision(a))
      z = set_valuation!(z, b*valuation(a))
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

@doc Markdown.doc"""
    ==(x::RelSeriesElem{T}, y::RelSeriesElem{T}) where T <: RingElement

Return `true` if $x == y$ arithmetically, otherwise return `false`. Recall
that power series to different precisions may still be arithmetically
equal to the minimum of the two precisions.
"""
function ==(x::RelSeriesElem{T}, y::RelSeriesElem{T}) where T <: RingElement
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

@doc Markdown.doc"""
    isequal(x::RelSeriesElem{T}, y::RelSeriesElem{T}) where T <: RingElement

Return `true` if $x == y$ exactly, otherwise return `false`. Only if the
power series are precisely the same, to the same precision, are they declared
equal by this function.
"""
function isequal(x::RelSeriesElem{T}, y::RelSeriesElem{T}) where T <: RingElement
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

@doc Markdown.doc"""
    ==(x::RelSeriesElem{T}, y::T) where {T <: RingElem}

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::RelSeriesElem{T}, y::T) where {T <: RingElem} = precision(x) == 0 ||
           ((pol_length(x) == 0 && iszero(y)) || (pol_length(x) == 1 &&
             valuation(x) == 0 && polcoeff(x, 0) == y))

@doc Markdown.doc"""
    ==(x::T, y::RelSeriesElem{T}) where {T <: RingElem}

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::T, y::RelSeriesElem{T}) where {T <: RingElem} = y == x

@doc Markdown.doc"""
    ==(x::RelSeriesElem, y::Union{Integer, Rational, AbstractFloat})

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::RelSeriesElem, y::Union{Integer, Rational, AbstractFloat}) = precision(x) == 0 ||
                  ((pol_length(x) == 0 && iszero(y)) || (pol_length(x) == 1 &&
                    valuation(x) == 0 && polcoeff(x, 0) == y))

@doc Markdown.doc"""
    ==(x::Union{Integer, Rational, AbstractFloat}, y::RelSeriesElem)

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::Union{Integer, Rational, AbstractFloat}, y::RelSeriesElem) = y == x

###############################################################################
#
#   Approximation
#
###############################################################################

function Base.isapprox(f::RelSeriesElem, g::RelSeriesElem; atol::Real=sqrt(eps()))
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

function divexact(x::RelSeriesElem{T}, y::RelSeriesElem{T}) where T <: RingElement
   check_parent(x, y)
   iszero(y) && throw(DivideError())
   v2 = valuation(y)
   if v2 != 0
      v1 = valuation(x)
      if v1 >= v2
         x = shift_right(x, v2)
         y = shift_right(y, v2)
      else
         error("Not an exact division")
      end
   else
      x = deepcopy(x)
   end
   y = truncate(y, precision(x))
   res = parent(x)()
   res = set_precision!(res, min(precision(x), valuation(x) + precision(y)))
   res = set_valuation!(res, valuation(x))
   lc = coeff(y, 0)
   lc == 0 && error("Not an exact division")
   lenr = precision(x) - valuation(x)
   for i = 0:lenr - 1
      flag, q = divides(polcoeff(x, i), lc)
      !flag && error("Not an exact division")
      res = setcoeff!(res, i, q)
      for j = 0:min(precision(y) - 1, lenr - i - 1)
         x = setcoeff!(x, i + j, polcoeff(x, i + j) - polcoeff(y, j)*q)
      end
   end
   res = set_length!(res, normalise(res, pol_length(res)))
   return res
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(x::RelSeriesElem, y::Union{Integer, Rational, AbstractFloat})
   y == 0 && throw(DivideError())
   lenx = pol_length(x)
   z = parent(x)()
   fit!(z, lenx)
   z = set_precision!(z, precision(x))
   z = set_valuation!(z, valuation(x))
   for i = 1:lenx
      z = setcoeff!(z, i - 1, divexact(polcoeff(x, i - 1), y))
   end
   return z
end

function divexact(x::RelSeriesElem{T}, y::T) where {T <: RingElem}
   iszero(y) && throw(DivideError())
   lenx = pol_length(x)
   z = parent(x)()
   fit!(z, lenx)
   z = set_precision!(z, precision(x))
   z = set_valuation!(z, valuation(x))
   for i = 1:lenx
      z = setcoeff!(z, i - 1, divexact(polcoeff(x, i - 1), y))
   end
   return z
end

###############################################################################
#
#   Inversion
#
###############################################################################

@doc Markdown.doc"""
    Base.inv(a::RelSeriesElem)

Return the inverse of the power series $a$, i.e. $1/a$.
"""
function Base.inv(a::RelSeriesElem)
   iszero(a) && throw(DivideError())
   !isunit(a) && error("Unable to invert power series")
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

function Base.inv(a::RelSeriesElem{T}) where T <: FieldElement
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
#   Square root
#
###############################################################################

@doc Markdown.doc"""
    sqrt(a::RelSeriesElem)

Return the square root of the power series $a$.
"""
function Base.sqrt(a::RelSeriesElem)
   aval = valuation(a)
   !iseven(aval) && error("Not a square in sqrt")
   R = base_ring(a)
   !isdomain_type(elem_type(R)) && error("Sqrt not implemented over non-integral domains")
   aval2 = div(aval, 2)
   prec = precision(a) - aval
   if prec == 0
      asqrt = parent(a)()
      asqrt = set_precision!(asqrt, aval2)
      asqrt = set_valuation!(asqrt, aval2)
      return asqrt
   end
   asqrt = parent(a)()
   fit!(asqrt, prec)
   asqrt = set_precision!(asqrt, prec + aval2)
   asqrt = set_valuation!(asqrt, aval2)
   if prec > 0
      g = sqrt(polcoeff(a, 0))
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
      c = divexact(c, g2)
      asqrt = setcoeff!(asqrt, n, c)
   end
   asqrt = set_length!(asqrt, normalise(asqrt, prec))
   return asqrt
end

###############################################################################
#
#  Derivative and Integral
#
###############################################################################

@doc Markdown.doc"""
    derivative(f::RelSeriesElem{T})

Return the derivative of the power series $f$.

```
julia> R, x = PowerSeriesRing(QQ, 10, "x")
(Univariate power series ring in x over Rationals, x + O(x^11))

julia> f = 2 + x + 3x^3
2 + x + 3*x^3 + O(x^10)

julia> derivative(f)
1 + 9*x^2 + O(x^9)
```
"""
function derivative(f::RelSeriesElem{T}) where T <: RingElement
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

@doc Markdown.doc"""
    integral(f::RelSeriesElem{T})

Return the integral of the power series $f$.

```
julia> R, x = PowerSeriesRing(QQ, 10, "x")
(Univariate power series ring in x over Rationals, x + O(x^11))

julia> f = 2 + x + 3x^3
2 + x + 3*x^3 + O(x^10)

julia> integral(f)
2*x + 1//2*x^2 + 3//4*x^4 + O(x^11)
```
"""
function integral(f::RelSeriesElem{T}) where T <: RingElement
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

@doc Markdown.doc"""
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

@doc Markdown.doc"""
    exp(a::RelSeriesElem)

Return the exponential of the power series $a$.
"""
function Base.exp(a::RelSeriesElem{T}) where T <: RingElement
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
      !isunit(R(k)) && error("Unable to divide in exp")
      z = setcoeff!(z, k, divexact(s, k))
   end
   z = set_length!(z, normalise(z, preca))
   return z
end

function Base.exp(a::RelSeriesElem{T}) where T <: FieldElement
   if iszero(a)
      b = parent(a)(1)
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
   x = parent(a)([R(1)], 1, min(2, precision(a)), 0)
   prec = precision(a)
   la = [prec]
   while la[end] > 1
      push!(la, div(la[end] + 1, 2))
   end
   one1 = parent(a)([R(1)], 1, 2, 0)
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

function _make_parent(g, p::RelSeriesElem, cached::Bool)
   R = parent(g(zero(base_ring(p))))
   S = parent(p)
   sym = String(var(S))
   max_prec = max_precision(S)
   return PowerSeriesRing(R, max_prec, sym; cached=cached)[1]
end

function map_coefficients(g, p::RelSeriesElem{<:RingElement};
                    cached::Bool = true,
                    parent::Ring = _make_parent(g, p, cached))
   return _map(g, p, parent)
end

function _map(g, p::RelSeriesElem, Rx)
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
   P, _ = PowerSeriesRing(R, max_precision(Rx),
                                               string(var(Rx)), cached = cached)
   return P
end

function change_base_ring(R::Ring, p::RelSeriesElem{T};
                    cached::Bool = true, parent::Ring =
          _change_rel_series_ring(R, parent(p), cached)) where T <: RingElement
   return _map(R, p, parent)
end

###############################################################################
#
#   Random elements
#
###############################################################################

RandomExtensions.maketype(S::SeriesRing, ::UnitRange{Int}, _) = elem_type(S)

function RandomExtensions.make(S::SeriesRing, val_range::UnitRange{Int}, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, val_range, vs[1]) # forward to default Make constructor
   else
      make(S, val_range, make(R, vs...))
   end
end

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make3{<:RingElement, <:SeriesRing, UnitRange{Int}}})
   S, val_range, v = sp[][1:end]
   R = base_ring(S)
   f = S()
   x = gen(S)
   for i = 0:S.prec_max - 1
      f += rand(rng, v)*x^i
   end
   return shift_left(f, rand(rng, val_range))
end

rand(rng::AbstractRNG, S::SeriesRing, val_range::UnitRange{Int}, v...) =
   rand(rng, make(S, val_range, v...))

rand(S::SeriesRing, val_range, v...) = rand(Random.GLOBAL_RNG, S, val_range, v...)

###############################################################################
#
#   PowerSeriesRing constructor
#
###############################################################################

@doc Markdown.doc"""
    PowerSeriesRing(R::Ring, prec::Int, s::Union{AbstractString, char, Symbol}; cached=true, model=:capped_relative)

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
PowerSeriesRing(R::Ring, prec::Int, s::Union{AbstractString, Char, Symbol}; cached=true, model=:capped_relative)

function PowerSeriesRing(R::Ring, prec::Int, s::Symbol; cached=true, model=:capped_relative)
   return Generic.PowerSeriesRing(R, prec, s; cached=cached, model=model)
end

function PowerSeriesRing(R::Ring, prec::Int, s::Char; cached=true, model=:capped_relative)
   return Generic.PowerSeriesRing(R, prec, Symbol(s); cached=cached, model=model)
end

function PowerSeriesRing(R::Ring, prec::Int, s::AbstractString; cached=true, model=:capped_relative)
   return Generic.PowerSeriesRing(R, prec, Symbol(s); cached=cached, model=model)
end
