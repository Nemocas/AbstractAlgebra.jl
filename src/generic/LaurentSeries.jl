###############################################################################
#
#   LaurentSeries.jl : Generic Laurent series over rings and fields,
#                      capped relative precision
#
###############################################################################

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

@doc raw"""
    O(a::Generic.LaurentSeriesElem{T}) where T <: RingElement

Return $0 + O(x^\mathrm{val}(a))$. Usually this function is called with $x^n$
as parameter. Then the function returns the power series $0 + O(x^n)$, which
can be used to set the precision of a power series when constructing it.
"""
function O(a::LaurentSeriesElem{T}) where T <: RingElement
   val = valuation(a)
   return parent(a)(Vector{T}(undef, 0), 0, val, val, 1)
end

parent_type(::Type{LaurentSeriesRingElem{T}}) where T <: RingElement = LaurentSeriesRing{T}

parent_type(::Type{LaurentSeriesFieldElem{T}}) where T <: FieldElement = LaurentSeriesField{T}

parent(a::LaurentSeriesElem) = a.parent

elem_type(::Type{LaurentSeriesRing{T}}) where T <: RingElement = LaurentSeriesRingElem{T}

elem_type(::Type{LaurentSeriesField{T}}) where T <: FieldElement = LaurentSeriesFieldElem{T}

base_ring_type(::Type{LaurentSeriesRing{T}}) where T <: RingElement = parent_type(T)

base_ring_type(::Type{LaurentSeriesField{T}}) where T <: FieldElement = parent_type(T)

base_ring(R::LaurentSeriesRing{T}) where T <: RingElement = R.base_ring::parent_type(T)

base_ring(R::LaurentSeriesField{T}) where T <: FieldElement = R.base_ring::parent_type(T)

function is_domain_type(::Type{T}) where {S <: RingElement, T <: LaurentSeriesElem{S}}
   return is_domain_type(S)
end

is_exact_type(a::Type{T}) where T <: LaurentSeriesElem = false

@doc raw"""
    var(a::LaurentSeriesRing)

Return the internal name of the generator of the power series ring. Note that
this is returned as a `Symbol` not a `String`.
"""
var(a::LaurentSeriesRing) = a.S

@doc raw"""
    var(a::LaurentSeriesField)

Return the internal name of the generator of the power series ring. Note that
this is returned as a `Symbol` not a `String`.
"""
var(a::LaurentSeriesField) = a.S

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function Base.hash(a::LaurentSeriesElem, h::UInt)
   b = 0xb163af5694734274%UInt
   for i in 0:pol_length(a) - 1
      b = xor(b, hash(polcoeff(a, i), h))
      b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
   end
   b = xor(b, hash(scale(a), h))
   return b
end

@doc raw"""
    pol_length(a::Generic.LaurentSeriesElem)

Return the length of the polynomial underlying the given power series. This
will be zero if the power series has no nonzero terms.
"""
pol_length(a::LaurentSeriesElem) = a.length

@doc raw"""
    precision(a::Generic.LaurentSeriesElem)

Return the precision of the given power series in absolute terms. This will
be the sum of the valuation and the length of the underlying polynomial.
"""
precision(a::LaurentSeriesElem) = a.prec

@doc raw"""
    valuation(a::Generic.LaurentSeriesElem)

Return the valuation of the given power series, i.e. the degree of the first
nonzero term (or the precision if it is arithmetically zero).
"""
valuation(a::LaurentSeriesElem) = a.val

@doc raw"""
    scale(a::Generic.LaurentSeriesElem)

Return the scale factor of the polynomial underlying the given power series.
"""
scale(a::LaurentSeriesElem) = a.scale

@doc raw"""
    max_precision(R::LaurentSeriesRing)

Return the maximum relative precision of power series in the given power
series ring.
"""
max_precision(R::LaurentSeriesRing) = R.prec_max

@doc raw"""
    max_precision(R::LaurentSeriesField)

Return the maximum relative precision of power series in the given power
series ring.
"""
max_precision(R::LaurentSeriesField) = R.prec_max

@doc raw"""
    exp_gcd(a::Generic.LaurentSeriesElem)

Return the GCD of the exponents of the polynomial underlying the given Laurent series.
"""
function exp_gcd(a::LaurentSeriesElem)
   n = 0
   s = scale(a)
   for i = 1:pol_length(a) - 1
      if n == 1
         return n
      end
      if polcoeff(a, i) != 0
         n = gcd(n, i)
      end
   end
   return n
end

function normalise(a::LaurentSeriesElem, len::Int)
   while len > 0 && iszero(a.coeffs[len])
      len -= 1
   end
   return len
end

function set_length!(a::LaurentSeriesElem, len::Int)
   a.length = len
   return a
end

function set_precision!(a::LaurentSeriesElem, prec::Int)
   a.prec = prec
   return a
end

function set_valuation!(a::LaurentSeriesElem, val::Int)
   a.val = val
   return a
end

@doc raw"""
    set_scale!(a::Generic.LaurentSeriesElem, scale::Int)

Set the scale factor of the polynomial underlying the given series to the given value.
"""
function set_scale!(a::LaurentSeriesElem, scale::Int)
   a.scale = scale
   return a
end

function polcoeff(a::LaurentSeriesElem, n::Int)
   n < 0  && throw(DomainError(n, "n must be >= 0"))
   return n >= pol_length(a) ? zero(base_ring(a)) : a.coeffs[n + 1]
end

function coeff(a::LaurentSeriesElem, n::Int)
   if n < valuation(a)
      return base_ring(a)()
   else
      i = n - valuation(a)
      if mod(i, scale(a)) != 0
         return base_ring(a)()
      else
         return polcoeff(a, div(i, scale(a)))
      end
   end
end

@doc raw"""
    rescale!(a::Generic.LaurentSeriesElem)

Rescale the polynomial underlying the series so that the GCD of its exponents is 1.
This is only used internally, since the result of every user facing function is a
rescaled series.
"""
function rescale!(a::LaurentSeriesElem)
   s = exp_gcd(a)
   if s > 1
      zlen = div(pol_length(a) - 1, s) + 1
      for i = 1:zlen - 1
         t = polcoeff(a, i)
         a = setcoeff!(a, i, polcoeff(a, i*s))
         a = setcoeff!(a, i*s, t)
      end
      a = set_scale!(a, s*scale(a))
      a = set_length!(a, zlen)
   elseif pol_length(a) <= 1
      a = set_scale!(a, 1)
   end
   return a
end

@doc raw"""
    downscale(a::Generic.LaurentSeriesElem{T}, n::Int) where T <: RingElement

Inflate the underlying polynomial by a factor of $n$. This inserts zero coefficients
for padding. It is assumed that the scale factor of $a$ is divisible by $n$.
"""
function downscale(a::LaurentSeriesElem{T}, n::Int) where T <: RingElement
   n <= 0 && throw(DomainError(n, "n must be > 0"))
   lena = pol_length(a)
   if n == 1 || lena == 0
      return a
   end
   R = base_ring(a)
   lenz = (lena - 1)*n + 1
   d = Vector{T}(undef, lenz)
   j = 0
   pn = 0
   for i = 0:lenz - 1
      if i == pn
         d[i + 1] = polcoeff(a, j)
         j += 1
         pn += n
      else
         d[i + 1] = R()
      end
   end
   S = typeof(a)
   z = S(d, lenz, precision(a), valuation(a), div(scale(a), n))
   z.parent = parent(a)
   return z
end

@doc raw"""
    upscale(a::Generic.LaurentSeriesElem{T}, n::Int) where T <: RingElement

Deflate the underlying polynomial by a factor of $n$. This removes zero coefficients
that existed for padding. It is assumed that the spacing of nonzero coefficients of
$a$ is divisible by $n$.
"""
function upscale(a::LaurentSeriesElem{T}, n::Int) where T <: RingElement
   n <= 0 && throw(DomainError(n, "n must be > 0"))
   lena = pol_length(a)
   if n == 1 || lena == 0
      return a
   end
   R = base_ring(a)
   lenz = div(lena - 1, n) + 1
   d = Vector{T}(undef, lenz)
   j = 0
   for i = 1:lenz
      d[i] = polcoeff(a, j)
      j += n
   end
   S = typeof(a)
   z = S(d, lenz, precision(a), valuation(a), scale(a)*n)
   z.parent = parent(a)
   return z
end

zero(R::LaurentSeriesRing) = R(0)

zero(R::LaurentSeriesField) = R(0)

one(R::LaurentSeriesField) = R(1)

one(R::LaurentSeriesRing) = R(1)

@doc raw"""
    gen(R::LaurentSeriesRing)

Return the generator of the power series ring, i.e. $x + O(x^{n + 1})$ where
$n$ is the maximum precision of the power series ring $R$.
"""
function gen(R::LaurentSeriesRing)
   S = base_ring(R)
   return R([one(S)], 1, max_precision(R) + 1, 1, 1)
end

@doc raw"""
    gen(R::LaurentSeriesField)

Return the generator of the power series ring, i.e. $x + O(x^{n + 1})$ where
$n$ is the maximum precision of the power series ring $R$.
"""
function gen(R::LaurentSeriesField)
   S = base_ring(R)
   return R([one(S)], 1, max_precision(R) + 1, 1, 1)
end

iszero(a::LaurentSeriesElem) = pol_length(a) == 0

function isone(a::LaurentSeriesElem)
   return valuation(a) == 0 && pol_length(a) == 1 && isone(polcoeff(a, 0))
end

@doc raw"""
    is_gen(a::Generic.LaurentSeriesElem)

Return `true` if the given power series is arithmetically equal to the
generator of its power series ring to its current precision, otherwise return
`false`.
"""
function is_gen(a::LaurentSeriesElem)
   return valuation(a) == 1 && pol_length(a) == 1 && isone(polcoeff(a, 0))
end

is_unit(a::LaurentSeriesElem) = valuation(a) == 0 && is_unit(polcoeff(a, 0))

@doc raw"""
    modulus(a::Generic.LaurentSeriesElem{T}) where {T <: ResElem}

Return the modulus of the coefficients of the given power series.
"""
modulus(a::LaurentSeriesElem{T}) where {T <: ResElem} = modulus(base_ring(a))

function deepcopy_internal(a::LaurentSeriesElem{T}, dict::IdDict) where {T <: RingElement}
   coeffs = Vector{T}(undef, pol_length(a))
   for i = 1:pol_length(a)
      coeffs[i] = deepcopy_internal(polcoeff(a, i - 1), dict)
   end
   return parent(a)(coeffs, pol_length(a), precision(a), valuation(a), scale(a))
end

function renormalize!(z::LaurentSeriesElem)
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
      z = set_scale!(z, 1)
   elseif i != 0
      z = set_valuation!(z, zval + i*scale(z))
      for j = 1:zlen - i
         z = setcoeff!(z, j - 1, polcoeff(z, j + i - 1))
      end
      z = set_length!(z, zlen - i)
   end
   return nothing
end

function characteristic(a::LaurentSeriesRing{T}) where T <: RingElement
   return characteristic(base_ring(a))
end

###############################################################################
#
#   Similar and zero
#
###############################################################################

function similar(x::LaurentSeriesElem, R::Ring, max_prec::Int,
      s::Symbol; cached::Bool=true)
   TT = elem_type(R)
   V = Vector{TT}(undef, 0)
   p = Generic.LaurentSeriesRingElem{TT}(V, 0, max_prec, max_prec, 1)
   # Default similar is supposed to return a Generic series
   if base_ring(x) === R && s == var(parent(x)) &&
         x isa Generic.LaurentSeriesRingElem{TT} &&
      max_precision(parent(x)) == max_prec
      # steal parent in case it is not cached
      p.parent = parent(x)
   else
      p.parent = Generic.LaurentSeriesRing{TT}(R, max_prec, s, cached)
   end
   return p
end

function similar(x::LaurentSeriesElem, R::Field, max_prec::Int,
      s::Symbol; cached::Bool=true)
   TT = elem_type(R)
   V = Vector{TT}(undef, 0)
   p = Generic.LaurentSeriesFieldElem{TT}(V, 0, max_prec, max_prec, 1)
   # Default similar is supposed to return a Generic series
   if base_ring(x) === R && s == var(parent(x)) &&
         x isa Generic.LaurentSeriesFieldElem{TT} &&
      max_precision(parent(x)) == max_prec
      # steal parent in case it is not cached
      p.parent = parent(x)
   else
      p.parent = Generic.LaurentSeriesField{TT}(R, max_prec, s, cached)
   end
   return p
end

similar(x::LaurentSeriesElem, R::Ring, max_prec::Int,
                              var::VarName=var(parent(x)); cached::Bool=true) =
   similar(x, R, max_prec, Symbol(var); cached)

similar(x::LaurentSeriesElem, R::Ring,
                              var::VarName=var(parent(x)); cached::Bool=true) =
   similar(x, R, max_precision(parent(x)), Symbol(var); cached)

similar(x::LaurentSeriesElem, max_prec::Int,
                              var::VarName=var(parent(x)); cached::Bool=true) =
   similar(x, base_ring(x), max_prec, Symbol(var); cached)

similar(x::LaurentSeriesElem, var::VarName=var(parent(x)); cached::Bool=true) =
   similar(x, base_ring(x), max_precision(parent(x)), Symbol(var); cached)

zero(a::LaurentSeriesElem, R::Ring, max_prec::Int,
                              var::VarName=var(parent(a)); cached::Bool=true) =
   similar(a, R, max_prec, Symbol(var); cached)

zero(a::LaurentSeriesElem, R::Ring,
                              var::VarName=var(parent(a)); cached::Bool=true) =
   similar(a, R, Symbol(var); cached)

zero(a::LaurentSeriesElem, max_prec::Int,
                              var::VarName=var(parent(a)); cached::Bool=true) =
   similar(a, max_prec, Symbol(var); cached)

zero(a::LaurentSeriesElem,    var::VarName=var(parent(a)); cached::Bool=true) =
   similar(a, Symbol(var); cached)

###############################################################################
#
#   laurent_series constructor
#
###############################################################################

function laurent_series(R::Ring, arr::Vector{T}, len::Int, prec::Int, val::Int, scale::Int, var::VarName=:x; max_precision::Int=prec, cached::Bool=true) where T
   scale <= 0 && error("Scale must be positive")
   prec < (len - 1)*scale + val + 1 && error("Precision too small for given data")
   TT = elem_type(R)
   coeffs = T === Any && length(arr) == 0 ? elem_type(R)[] : map(R, arr)
   p = Generic.LaurentSeriesRingElem{TT}(coeffs, len, prec, val, scale)
   # Default is supposed to return a Generic Laurent series
   p.parent = Generic.LaurentSeriesRing{TT}(R, max_precision, Symbol(var), cached)
   return p
end

function laurent_series(R::Field, arr::Vector{T}, len::Int, prec::Int, val::Int, scale::Int, var::VarName=:x; max_precision::Int=prec, cached::Bool=true) where T
   scale <= 0 && error("Scale must be positive")
   prec < (len - 1)*scale + val + 1 && error("Precision too small for given data")
   TT = elem_type(R)
   coeffs = T === Any && length(arr) == 0 ? elem_type(R)[] : map(R, arr)
   p = Generic.LaurentSeriesFieldElem{TT}(coeffs, len, prec, val, scale)
   # Default is supposed to return a Generic Laurent series
   p.parent = Generic.LaurentSeriesField{TT}(R, max_precision, Symbol(var), cached)
   return p
end

###############################################################################
#
#   AbstractString I/O
#
###############################################################################

function AbstractAlgebra.expressify(a::LaurentSeriesElem,
                                    x = var(parent(a)); context = nothing)
    sum = Expr(:call, :+)
    v = valuation(a)
    sc = scale(a)
    len = pol_length(a)

    for i in 0:len - 1
        c = polcoeff(a, i)
        expo = i * sc + v
        if !iszero(c)
            if expo == 0
                xk = 1
            elseif expo == 1
                xk = x
            else
                xk = Expr(:call, :^, x, expo)
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

function Base.show(io::IO, ::MIME"text/plain", a::LaurentSeriesElem)
  print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

function Base.show(io::IO, a::LaurentSeriesElem)
  print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

function show(io::IO, p::LaurentSeriesRing)
   @show_name(io, p)
   @show_special(io, p)
   if is_terse(io)
      print(io, "Laurent series ring")
   else
      io = pretty(io)
      print(io, "Laurent series ring in ", var(p), " over ")
      print(terse(io), Lowercase(), base_ring(p))
   end
end

function show(io::IO, p::LaurentSeriesField)
   @show_name(io, p)
   @show_special(io, p)
   if is_terse(io)
      print(io, "Laurent series field")
   else
      io = pretty(io)
      print(io, "Laurent series field in ", var(p), " over ")
      print(terse(io), Lowercase(), base_ring(p))
   end
end

###############################################################################
#
#   Map coefficients
#
###############################################################################

function _make_parent(g::T, p::LaurentSeriesElem, cached::Bool) where T
   R = parent(g(zero(base_ring(p))))
   S = parent(p)
   sym = var(S)
   max_prec = max_precision(S)
   return AbstractAlgebra.laurent_series_ring(R, max_prec, sym; cached=cached)[1]
end

function map_coefficients(g::T, p::LaurentSeriesElem{<:RingElement};
                    cached::Bool = true,
		    parent::Ring = _make_parent(g, p, cached)) where {T}
   return _map(g, p, parent)
end

function _map(g::T, p::LaurentSeriesElem, Rx) where T
   R = base_ring(Rx)
   new_coefficients = elem_type(R)[let c = polcoeff(p, i)
                                     iszero(c) ? zero(R) : R(g(c))
                                   end for i in 0:pol_length(p) - 1]
   res = Rx(new_coefficients, pol_length(p), precision(p), valuation(p), scale(p), false)
   res = set_length!(res, normalise(res, pol_length(res)))
   renormalize!(res)
   res = rescale!(res)
   return res
end

################################################################################
#
#  Change base ring
#
################################################################################

function _change_laurent_series_ring(R, Rx, cached)
   P, _ = AbstractAlgebra.laurent_series_ring(R, max_precision(Rx),
                                               var(Rx), cached = cached)
   return P
end

function change_base_ring(R::Ring, p::LaurentSeriesElem{T};
                    cached::Bool = true, parent::Ring =
          _change_laurent_series_ring(R, parent(p), cached)) where T <: RingElement
   return _map(R, p, parent)
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::LaurentSeriesElem)
   len = pol_length(a)
   z = parent(a)()
   z = set_precision!(z, precision(a))
   z = set_valuation!(z, valuation(a))
   z = set_scale!(z, scale(a))
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

function +(a::LaurentSeriesElem{T}, b::LaurentSeriesElem{T}) where {T <: RingElement}
   check_parent(a, b)
   lena = pol_length(a)
   lenb = pol_length(b)
   vala = valuation(a)
   valb = valuation(b)
   valz = min(vala, valb)
   prec = min(precision(a), precision(b))
   sa = scale(a)
   sb = scale(b)
   if lena == 1
      sa = sb
   elseif lenb == 1
      sb = sa
   end
   sz = gcd(gcd(sa, sb), abs(vala - valb))
   mina = min(vala + lena*sa, prec)
   minb = min(valb + lenb*sb, prec)
   lenz = max(mina, minb) - valz
   lenz = div(lenz + sz - 1, sz)
   R = base_ring(a)
   z = parent(a)()
   fit!(z, lenz)
   z = set_precision!(z, prec)
   z = set_valuation!(z, valz)
   z = set_scale!(z, sz)
   pa = vala
   pb = valb
   j = 0
   k = 0
   for i = 0: lenz - 1
      pi = valz + sz*i
      if pi == pa && pi < mina
         if pi == pb && pi < minb
            z = setcoeff!(z, i, polcoeff(a, j) + polcoeff(b, k))
            pb += sb
            k += 1
         else
            z = setcoeff!(z, i, deepcopy(polcoeff(a, j)))
         end
         j += 1
         pa += sa
      elseif pi == pb && pi < minb
         z = setcoeff!(z, i, deepcopy(polcoeff(b, k)))
         k += 1
         pb += sb
      else
         z = setcoeff!(z, i, R())
      end
   end
   z = set_length!(z, normalise(z, lenz))
   renormalize!(z)
   z = rescale!(z)
   return z
end

function -(a::LaurentSeriesElem{T}, b::LaurentSeriesElem{T}) where {T <: RingElement}
   check_parent(a, b)
   lena = pol_length(a)
   lenb = pol_length(b)
   vala = valuation(a)
   valb = valuation(b)
   valz = min(vala, valb)
   prec = min(precision(a), precision(b))
   sa = scale(a)
   sb = scale(b)
   if lena == 1
      sa = sb
   elseif lenb == 1
      sb = sa
   end
   sz = gcd(gcd(sa, sb), abs(vala - valb))
   mina = min(vala + lena*sa, prec)
   minb = min(valb + lenb*sb, prec)
   lenz = max(mina, minb) - valz
   lenz = div(lenz + sz - 1, sz)
   R = base_ring(a)
   z = parent(a)()
   fit!(z, lenz)
   z = set_precision!(z, prec)
   z = set_valuation!(z, valz)
   z = set_scale!(z, sz)
   pa = vala
   pb = valb
   j = 0
   k = 0
   for i = 0: lenz - 1
      pi = valz + sz*i
      if pi == pa && pi < mina
         if pi == pb && pi < minb
            z = setcoeff!(z, i, polcoeff(a, j) - polcoeff(b, k))
            pb += sb
            k += 1
         else
            z = setcoeff!(z, i, deepcopy(polcoeff(a, j)))
         end
         j += 1
         pa += sa
      elseif pi == pb && pi < minb
         z = setcoeff!(z, i, -polcoeff(b, k))
         k += 1
         pb += sb
      else
         z = setcoeff!(z, i, R())
      end
   end
   z = set_length!(z, normalise(z, lenz))
   renormalize!(z)
   z = rescale!(z)
   return z
end

function *(a::LaurentSeriesElem{T}, b::LaurentSeriesElem{T}) where {T <: RingElement}
   check_parent(a, b)
   lena = pol_length(a)
   lenb = pol_length(b)
   if lena > lenb
      return b*a
   end
   aval = valuation(a)
   bval = valuation(b)
   zval = aval + bval
   prec = min(precision(a) - aval, precision(b) - bval)
   sa = scale(a)
   sb = scale(b)
   if lena == 1
      sa = sb
   elseif lenb == 1
      sb = sa
   end
   sz = gcd(sa, sb)
   lena = min(lena*sa, prec)
   lenb = min(lenb*sb, prec)
   if lena == 0 || lenb == 0
      return parent(a)(Vector{T}(undef, 0), 0, prec + zval, zval, 1)
   end
   t = base_ring(a)()
   da = div(sa, sz)
   db = div(sb, sz)
   a = downscale(a, da)
   b = downscale(b, db)
   lena = pol_length(a)
   lenb = pol_length(b)
   lenz = min(lena + lenb - 1, div(prec + sz - 1, sz))
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
         ai = polcoeff(a, i - 1)
         if ai != 0
            for j = 2:min(lenb, lenz - i + 1)
               t = mul!(t, ai, polcoeff(b, j - 1))
               d[i + j - 1] = add!(d[i + j - 1], t)
            end
         end
      end
   end
   z = parent(a)(d, lenz, prec + zval, zval, sz)
   z = set_length!(z, normalise(z, lenz))
   renormalize!(z)
   z = rescale!(z)
   return z
end

###############################################################################
#
#   Ad hoc binary operators
#
###############################################################################

function *(a::T, b::LaurentSeriesElem{T}) where {T <: RingElem}
   len = pol_length(b)
   z = parent(b)()
   fit!(z, len)
   z = set_precision!(z, precision(b))
   z = set_valuation!(z, valuation(b))
   z = set_scale!(z, scale(b))
   for i = 1:len
      z = setcoeff!(z, i - 1, a*polcoeff(b, i - 1))
   end
   z = set_length!(z, normalise(z, len))
   renormalize!(z)
   z = rescale!(z)
   return z
end

function *(a::Union{Integer, Rational, AbstractFloat}, b::LaurentSeriesElem)
   len = pol_length(b)
   z = parent(b)()
   fit!(z, len)
   z = set_precision!(z, precision(b))
   z = set_valuation!(z, valuation(b))
   z = set_scale!(z, scale(b))
   for i = 1:len
      z = setcoeff!(z, i - 1, a*polcoeff(b, i - 1))
   end
   z = set_length!(z, normalise(z, len))
   renormalize!(z)
   z = rescale!(z)
   return z
end

*(a::LaurentSeriesElem{T}, b::T) where {T <: RingElem} = b*a

*(a::LaurentSeriesElem, b::Union{Integer, Rational, AbstractFloat}) = b*a

###############################################################################
#
#   Shifting
#
###############################################################################

@doc raw"""
    shift_left(x::Generic.LaurentSeriesElem{T}, n::Int) where {T <: RingElement}

Return the power series $x$ shifted left by $n$ terms, i.e. multiplied by
$x^n$.
"""
function shift_left(x::LaurentSeriesElem{T}, n::Int) where {T <: RingElement}
   z = deepcopy(x)
   z = set_precision!(z, precision(x) + n)
   z = set_valuation!(z, valuation(x) + n)
   return z
end

@doc raw"""
    shift_right(x::Generic.LaurentSeriesElem{T}, n::Int) where {T <: RingElement}

Return the power series $x$ shifted right by $n$ terms, i.e. divided by
$x^n$.
"""
function shift_right(x::LaurentSeriesElem{T}, n::Int) where {T <: RingElement}
   z = deepcopy(x)
   z = set_precision!(z, precision(x) - n)
   z = set_valuation!(z, valuation(x) - n)
   return z
end

###############################################################################
#
#   Truncation
#
###############################################################################

@doc raw"""
    truncate(a::Generic.LaurentSeriesElem{T}, n::Int) where {T <: RingElement}

Return $a$ truncated to (absolute) precision $n$.
"""
function truncate(a::LaurentSeriesElem{T}, n::Int) where {T <: RingElement}
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
      z = set_scale!(z, 1)
   else
      sa = scale(a)
      zlen = div(n - aval + sa - 1, sa)
      zlen = min(zlen, alen)
      fit!(z, zlen)
      z = set_valuation!(z, aval)
      for i = 0:zlen - 1
         z = setcoeff!(z, i, polcoeff(a, i))
      end
      z = set_length!(z, normalise(z, zlen))
      z = set_scale!(z, sa)
      z = rescale!(z)
   end
   return z
end

# Intended only for internal use, does not renormalize, assumes n >= 0
# Requires valuation(a) == valuation(b) == 0 and scale(a) == scale(b)
function mullow(a::LaurentSeriesElem{T}, b::LaurentSeriesElem{T}, n::Int) where {T <: RingElement}
   lena = pol_length(a)
   lenb = pol_length(b)
   if lena == 0 || lenb == 0
      z = zero(parent(a))
      zprec = valuation(a) + valuation(b)
      z = set_valuation!(z, zprec)
      z = set_precision!(z, zprec)
      z = set_scale!(z, scale(a))
      return z
   end
   s = scale(a)
   prec = min(precision(a), precision(b))
   t = base_ring(a)()
   lenz = min(lena + lenb - 1, div(n + s - 1, s))
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
            d[i + j - 1] = add!(d[i + j - 1], t)
         end
      end
   end
   z = parent(a)(d, lenz, prec, 0, s, false)
   z = set_length!(z, normalise(z, lenz))
   return z
end

###############################################################################
#
#   Inflation/deflation
#
###############################################################################

function inflate(a::LaurentSeriesElem{T}, b::Int) where {T <: RingElement}
    return parent(a)(deepcopy(a.coeffs), pol_length(a), b*a.prec, b*a.val, b*a.scale)
end

function deflate(a::LaurentSeriesElem{T}, b::Int) where {T <: RingElement}
    return parent(a)(deepcopy(a.coeffs), pol_length(a), div(a.prec, b), div(a.val, b), div(a.scale, b))
end

###############################################################################
#
#   Powering
#
###############################################################################

@doc raw"""
    ^(a::Generic.LaurentSeriesElem{T}, b::Int) where {T <: RingElement}

Return $a^b$. We require $b \geq 0$.
"""
function ^(a::LaurentSeriesElem{T}, b::Int) where {T <: RingElement}
   # special case powers of x for constructing power series efficiently
   if b == 0
      # in fact, the result would be exact 1 if we had exact series
      z = one(parent(a))
      return z
   elseif pol_length(a) == 0
      z = parent(a)()
      z = set_precision!(z, b*valuation(a))
      z = set_valuation!(z, b*valuation(a))
      z = set_scale!(z, 1)
      return z
   elseif is_gen(a)
      z = parent(a)()
      fit!(z, 1)
      z = set_precision!(z, b + precision(a) - 1)
      z = set_valuation!(z, b)
      z = setcoeff!(z, 0, deepcopy(polcoeff(a, 0)))
      z = set_scale!(z, 1)
      z = set_length!(z, 1)
      return z
   elseif pol_length(a) == 1
      c = polcoeff(a, 0)^b
      z = parent(a)(c)
      z = set_precision!(z, (b - 1)*valuation(a) + precision(a))
      z = set_valuation!(z, iszero(c) ? precision(z) : b*valuation(a))
      z = set_scale!(z, 1)
      return z
   elseif b == 1
      return deepcopy(a)
   elseif b == -1
      return inv(a)
   end

   if b < 0
      a = inv(a)
      b = -b
   end

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
   if pol_length(z) <= 1
      z = set_scale!(z, 1)
   end
   renormalize!(z)
   z = rescale!(z)
   return z
end

###############################################################################
#
#   Comparison
#
###############################################################################

@doc raw"""
    ==(x::Generic.LaurentSeriesElem{T}, y::Generic.LaurentSeriesElem{T}) where {T <: RingElement}

Return `true` if $x == y$ arithmetically, otherwise return `false`. Recall
that power series to different precisions may still be arithmetically
equal to the minimum of the two precisions.
"""
function ==(x::LaurentSeriesElem{T}, y::LaurentSeriesElem{T}) where {T <: RingElement}
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
   sx = scale(x)
   sy = scale(y)
   xlen = min(pol_length(x), div(prec - xval + sx - 1, sx))
   ylen = min(pol_length(y), div(prec - yval + sy - 1, sy))
   i = 0
   j = 0
   while i < xlen && j < ylen
      while polcoeff(x, i) == 0 && i < xlen
         i += 1
      end
      while polcoeff(y, j) == 0 && j < ylen
         j += 1
      end
      if i < xlen && j < ylen
         if i*sx != j*sy || polcoeff(x, i) != polcoeff(y, j)
            return false
         end
         i += 1
         j += 1
      end
   end
   while i < xlen
      if polcoeff(x, i) != 0
         return false
      end
      i += 1
   end
   while j < ylen
      if polcoeff(y, j) != 0
         return false
      end
      j += 1
   end
   return true
end

@doc raw"""
    isequal(x::Generic.LaurentSeriesElem{T}, y::Generic.LaurentSeriesElem{T}) where {T <: RingElement}

Return `true` if $x == y$ exactly, otherwise return `false`. Only if the
power series are precisely the same, to the same precision, are they declared
equal by this function.
"""
function isequal(x::LaurentSeriesElem{T}, y::LaurentSeriesElem{T}) where {T <: RingElement}
   if parent(x) != parent(y)
      return false
   end
   if precision(x) != precision(y) || pol_length(x) != pol_length(y) ||
      valuation(x) != valuation(y) || scale(x) != scale(y)
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
    ==(x::Generic.LaurentSeriesElem{T}, y::T) where {T <: RingElem}

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::LaurentSeriesElem{T}, y::T) where {T <: RingElem} = precision(x) == 0 ||
           ((pol_length(x) == 0 && iszero(y)) || (pol_length(x) == 1 &&
             valuation(x) == 0 && polcoeff(x, 0) == y))

@doc raw"""
    ==(x::T, y::Generic.LaurentSeriesElem{T}) where {T <: RingElem}

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::T, y::LaurentSeriesElem{T}) where {T <: RingElem} = y == x

@doc raw"""
    ==(x::Generic.LaurentSeriesElem, y::Union{Integer, Rational, AbstractFloat})

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::LaurentSeriesElem, y::Union{Integer, Rational, AbstractFloat}) = precision(x) == 0 ||
                  ((pol_length(x) == 0 && iszero(y)) || (pol_length(x) == 1 &&
                    valuation(x) == 0 && polcoeff(x, 0) == y))

@doc raw"""
    ==(x::Union{Integer, Rational, AbstractFloat}, y::Generic.LaurentSeriesElem)

Return `true` if $x == y$ arithmetically, otherwise return `false`.
"""
==(x::Union{Integer, Rational, AbstractFloat}, y::LaurentSeriesElem) = y == x

###############################################################################
#
#   Approximation
#
###############################################################################

function Base.isapprox(f::LaurentSeriesElem, g::LaurentSeriesElem; atol::Real=sqrt(eps()))
   check_parent(f, g)
   pmin = min(precision(f), precision(g))
   vmin = min(valuation(f), valuation(g))
   for i = vmin:pmin - 1
      if !isapprox(coeff(f, i), coeff(g, i); atol=atol)
         return false
      end
   end
   return true
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(x::LaurentSeriesElem{T},
           y::LaurentSeriesElem{T}; check::Bool=true) where {T <: RingElement}
   check_parent(x, y)
   iszero(y) && throw(DivideError())
   v2 = valuation(y)
   if v2 != 0
      x = shift_right(x, v2)
      y = shift_right(y, v2)
   else
      x = deepcopy(x)
   end
   res = parent(x)()
   y = truncate(y, precision(x) - valuation(x))
   res = set_precision!(res, min(precision(x), valuation(x) + precision(y)))
   res = set_valuation!(res, valuation(x))
   sx = scale(x)
   sy = scale(y)
   sr = gcd(sx, sy)
   dx = div(sx, sr)
   dy = div(sy, sr)
   x = downscale(x, dx)
   y = downscale(y, dy)
   res = set_scale!(res, sr)
   lc = coeff(y, 0)
   lenr = div(precision(res) - valuation(res) + sr - 1, sr)
   leny = div(precision(y) + sr - 1, sr)
   for i = 0:lenr - 1
      q = divexact(polcoeff(x, i), lc; check=check)
      res = setcoeff!(res, i, q)
      for j = 0:min(leny - 1, lenr - i - 1)
         x = setcoeff!(x, i + j, polcoeff(x, i + j) - polcoeff(y, j)*q)
      end
   end
   res = set_length!(res, normalise(res, pol_length(res)))
   res = rescale!(res)
   return res
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(x::LaurentSeriesElem, y::Union{Integer, Rational, AbstractFloat}; check::Bool=true)
   y == 0 && throw(DivideError())
   lenx = pol_length(x)
   z = parent(x)()
   fit!(z, lenx)
   z = set_precision!(z, precision(x))
   z = set_valuation!(z, valuation(x))
   z = set_scale!(z, scale(x))
   for i = 1:lenx
      z = setcoeff!(z, i - 1, divexact(polcoeff(x, i - 1), y; check=check))
   end
   return z
end

function divexact(x::LaurentSeriesElem{T}, y::T; check::Bool=true) where {T <: RingElem}
   iszero(y) && throw(DivideError())
   lenx = pol_length(x)
   z = parent(x)()
   fit!(z, lenx)
   z = set_precision!(z, precision(x))
   z = set_valuation!(z, valuation(x))
   z = set_scale!(z, scale(x))
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
    Base.inv(a::Generic.LaurentSeriesElem)

Return the inverse of the power series $a$, i.e. $1/a$.
"""
function Base.inv(a::LaurentSeriesElem)
   iszero(a) && throw(DivideError())
   a1 = polcoeff(a, 0)
   ainv = parent(a)()
   sa = scale(a)
   lenz = div(precision(a) - valuation(a) + sa - 1, sa)
   fit!(ainv, lenz)
   ainv = set_precision!(ainv, precision(a) - 2*valuation(a))
   ainv = set_valuation!(ainv, -valuation(a))
   ainv = set_scale!(ainv, sa)
   !is_unit(a1) && error("Unable to invert power series")
   if lenz != 0
      ainv = setcoeff!(ainv, 0, divexact(one(base_ring(a)), a1))
   end
   a1 = -a1
   for n = 2:lenz
      s = polcoeff(a, 1)*polcoeff(ainv, n - 2)
      for i = 2:min(n, pol_length(a)) - 1
         s += polcoeff(a, i)*polcoeff(ainv, n - i - 1)
      end
      ainv = setcoeff!(ainv, n - 1, divexact(s, a1))
   end
   ainv = set_length!(ainv, normalise(ainv, lenz))
   ainv = rescale!(ainv)
   return ainv
end

###############################################################################
#
#   Square root
#
###############################################################################

function sqrt_classical_char2(a::LaurentSeriesElem; check::Bool=true)
   S = parent(a)
   R = base_ring(a)
   prec = div(precision(a) + 1, 2)
   if iszero(a)
      asqrt = parent(a)()
      asqrt = set_precision!(asqrt, prec)
      asqrt = set_valuation!(asqrt, prec)
      asqrt = set_scale!(asqrt, 1)
      return true, asqrt
   end
   aval = valuation(a)
   if check && !iseven(aval)
      return false, S()
   end
   s = scale(a)
   zlen = pol_length(a)
   if check
      if isodd(s) && zlen > 1
         # valuation is even, scale of square must be even
         # unless polynomial length is one
         return false, S()
      end
   end
   asqrt = parent(a)()
   zlen = pol_length(a)
   fit!(asqrt, zlen)
   aval2 = div(aval, 2)
   asqrt = set_precision!(asqrt, prec)
   asqrt = set_valuation!(asqrt, aval2)
   asqrt = set_scale!(asqrt, div(s + 1, 2)) # deals with s = 1
   for i = 0:zlen - 1
      c = polcoeff(a, i)
      if check && !is_square(c)
         return false, S()
      end
      asqrt = setcoeff!(asqrt, i, sqrt(c; check=false))
   end
   asqrt = set_length!(asqrt, zlen)
   return true, asqrt
end

function sqrt_classical(a::LaurentSeriesElem; check::Bool=true)
   S = parent(a)
   R = base_ring(a)
   !is_domain_type(elem_type(R)) && error("Sqrt not implemented over non-integral domains")
   if characteristic(R) == 2
      return sqrt_classical_char2(a, check=check)
   end
   aval = valuation(a)
   if check && !iseven(aval)
      return false, S()
   end
   aval2 = div(aval, 2)
   prec = precision(a) - aval
   if prec == 0
      asqrt = parent(a)()
      asqrt = set_precision!(asqrt, aval2)
      asqrt = set_valuation!(asqrt, aval2)
      asqrt = set_scale!(asqrt, 1)
      return true, asqrt
   end
   asqrt = parent(a)()
   s = scale(a)
   zlen = div(prec + s - 1, s)
   fit!(asqrt, zlen)
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
   for n = 1:zlen - 1
      c = R()
      for i = 1:div(n - 1, 2)
         j = n - i
         p = mul!(p, polcoeff(asqrt, i), polcoeff(asqrt, j))
         c = add!(c, p)
      end
      c *= 2
      if (n % 2) == 0
         i = div(n, 2)
         p = mul!(p, polcoeff(asqrt, i), polcoeff(asqrt, i))
         c = add!(c, p)
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
    asqrt = set_scale!(asqrt, s)
    asqrt = set_length!(asqrt, normalise(asqrt, zlen))
    asqrt = rescale!(asqrt)
    return true, asqrt
end

@doc raw"""
    sqrt(a::Generic.LaurentSeriesElem; check::Bool=true)

Return the square root of the power series $a$. By default the function will
throw an exception if the input is not square. If `check=false` this test is
omitted.
"""
function Base.sqrt(a::LaurentSeriesElem; check::Bool=true)
   flag, s = sqrt_classical(a, check=check)
   check && !flag && error("Not a square in sqrt")
   return s
end

function is_square(a::LaurentSeriesElem)
   flag, q = sqrt_classical(a; check=true)
   return flag
end

function is_square_with_sqrt(a::LaurentSeriesElem)
   return sqrt_classical(a; check=true)
end

###############################################################################
#
#   Derivative and integral
#
###############################################################################

function derivative(f::LaurentSeriesElem{T}) where T <: RingElement
   g = parent(f)()
   g = set_precision!(g, precision(f) - 1)
   sf = scale(f)
   g = set_scale!(g, sf)
   fit!(g, pol_length(f))
   v = valuation(f)
   if v == 0
      g = set_valuation!(g, sf - 1)
      for i = 1:pol_length(f) - 1
         g = setcoeff!(g, i - 1, sf*i*polcoeff(f, i))
      end
   else
      g = set_valuation!(g, v - 1)
      for i = 0:pol_length(f) - 1
         g = setcoeff!(g, i, (sf*i + v)*polcoeff(f, i))
      end
   end
   g = set_length!(g, normalise(g, pol_length(f)))
   renormalize!(g)
   g = rescale!(g)
   return g
end

function integral(f::LaurentSeriesElem{T}) where T <: RingElement
   g = parent(f)()
   fit!(g, pol_length(f))
   g = set_precision!(g, precision(f) + 1)
   v = valuation(f)
   g = set_valuation!(g, v + 1)
   sf = scale(f)
   g = set_scale!(g, sf)
   for i = 1:pol_length(f)
      c = polcoeff(f, i - 1)
      if !iszero(c)
         g = setcoeff!(g, i - 1, divexact(c, (i - 1)*sf + 1 + v))
      end
   end
   g = set_length!(g, normalise(g, pol_length(f)))
   renormalize!(g)
   g = rescale!(g)
   return g
end

###############################################################################
#
#   Special functions
#
###############################################################################

function Base.log(a::LaurentSeriesElem{T}) where T <: FieldElement
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
    exp(a::Generic.LaurentSeriesElem)

Return the exponential of the power series $a$.
"""
function Base.exp(a::LaurentSeriesElem{T}) where T <: RingElement
   if iszero(a)
      z = one(parent(a))
      z = set_precision!(z, precision(a))
      return z
   end
   vala = valuation(a)
   preca = precision(a)
   vala < 0 && error("Valuation must be non-negative in exp")
   sc = scale(a)
   gs = gcd(sc, gcd(vala, preca))
   if sc != gs
      a = downscale(a, div(sc, gs))
      sc = gs
   end
   if sc != 1
      vala = div(vala, sc)
      preca = div(preca, sc)
      a = parent(a)(a.coeffs, pol_length(a), preca, vala, 1, false)
   end
   z = parent(a)()
   R = base_ring(a)
   fit!(z, preca)
   z = set_precision!(z, preca)
   z = set_valuation!(z, 0)
   c = vala == 0 ? polcoeff(a, 0) : R()
   z = setcoeff!(z, 0, exp(c))
   len = pol_length(a) + vala
   for k = 1 : preca - 1
      s = R()
      for j = 1 : min(k + 1, len) - 1
         c = j >= vala ? polcoeff(a, j - vala) : R()
         s += j * c * polcoeff(z, k - j)
      end
      flag, q = divides(s, parent(s)(k))
      !flag && error("Unable to divide in exp")
      z = setcoeff!(z, k, q)
   end
   z = set_length!(z, normalise(z, preca))
   z = inflate(z, sc)
   z = rescale!(z)
   return z
end

function Base.exp(a::LaurentSeriesElem{T}) where T <: FieldElement
   if iszero(a)
      z = one(parent(a))
      z = set_precision!(z, precision(a))
      return z
   end
   vala = valuation(a)
   preca = precision(a)
   vala < 0 && error("Valuation must be non-negative in exp")
   R = base_ring(a)
   c = one(R)
   if valuation(a) == 0
      a = deepcopy(a)
      c = exp(coeff(a, 0))
      a = setcoeff!(a, 0, R())
   end
   z = parent(a)([one(R)], 1, min(2, preca), 0, 1, false)
   la = [preca]
   while la[end] > 1
      push!(la, div(la[end] + 1, 2))
   end
   one1 = parent(a)([one(R)], 1, 2, 0, 1, false)
   n = length(la) - 1
   # z -> z*(1 - log(a) + a) is the recursion
   while n > 0
      z = set_precision!(z, la[n])
      one1 = set_precision!(one1, la[n])
      t = -log(z)
      t = add!(t, one1)
      t = add!(t, a)
      z = mul!(z, z, t)
      n -= 1
   end
   if !isone(c)
      z *= c
   end
   renormalize!(z)
   z = rescale!(z)
   return z
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function zero!(a::LaurentSeriesElem)
   a.length = 0
   a.prec = parent(a).prec_max
   a.val = a.prec
   a.scale = 1
   return a
end

function fit!(c::LaurentSeriesElem{T}, n::Int) where {T <: RingElement}
   if length(c.coeffs) < n
      resize!(c.coeffs, n)
      for i = pol_length(c) + 1:n
         c.coeffs[i] = zero(base_ring(c))
      end
   end
   return nothing
end

function setcoeff!(c::LaurentSeriesElem{T}, n::Int, a::T) where {T <: RingElement}
   s = c.scale
   if (a != 0 && div(precision(c) - valuation(c) + s - 1, s) > n) || n + 1 <= c.length
      fit!(c, n + 1)
      c.coeffs[n + 1] = a
      c.length = max(pol_length(c), n + 1)
      # don't normalise
   end
   return c
end

function mul!(c::LaurentSeriesElem{T}, a::LaurentSeriesElem{T}, b::LaurentSeriesElem{T}) where {T <: RingElement}
   lena = pol_length(a)
   lenb = pol_length(b)
   if lena > lenb
      return mul!(c, b, a)
   end
   aval = valuation(a)
   bval = valuation(b)
   prec = min(precision(a) - aval, precision(b) - bval)
   sa = scale(a)
   sb = scale(b)
   if lena == 1
      sa = sb
   elseif lenb == 1
      sb = sa
   end
   sz = gcd(sa, sb)
   lena = min(lena*sa, prec)
   lenb = min(lenb*sb, prec)
   if lena == 0 || lenb == 0
      c.length = 0
      c.prec = prec + aval + bval
      c.val = aval + bval
      c.scale = 1
      return c
   end
   t = base_ring(a)()
   da = div(sa, sz)
   db = div(sb, sz)
   a2 = downscale(a, da)
   b2 = downscale(b, db)
   a = a === c && a2 === a ? deepcopy(a2) : a2
   b = b === c && b2 === b ? deepcopy(b2) : b2
   lena = pol_length(a)
   lenb = pol_length(b)
   lenc = min(lena + lenb - 1, div(prec + sz - 1, sz))
   fit!(c, lenc)
   for i = 1:min(lena, lenc)
      c.coeffs[i] = polcoeff(a, i - 1) * polcoeff(b, 0)
   end
   if lenc > lena
      for i = 2:min(lenb, lenc - lena + 1)
         c.coeffs[lena + i - 1] = polcoeff(a, lena - 1) * polcoeff(b, i - 1)
      end
   end
   for i = 1:lena - 1
      if lenc > i
         ai = polcoeff(a, i - 1)
         if ai != 0
            for j = 2:min(lenb, lenc - i + 1)
               t = mul!(t, polcoeff(a, i - 1), polcoeff(b, j - 1))
               c.coeffs[i + j - 1] = add!(c.coeffs[i + j - 1], t)
            end
         end
      end
   end
   c.length = normalise(c, lenc)
   c.val = aval + bval
   c.prec = prec + c.val
   c.scale = sz
   renormalize!(c)
   c = rescale!(c)
   return c
end

function add!(c::LaurentSeriesElem{T}, a::LaurentSeriesElem{T}) where {T <: RingElement}
   # TODO: write a version which doesn't make a copy
   b = deepcopy(c)
   return add!(c, b, a)
end

function add!(c::LaurentSeriesElem{T}, a::LaurentSeriesElem{T}, b::LaurentSeriesElem{T}) where {T <: RingElement}
   if c === a
      return add!(c, b)
   elseif c === b
      return add!(c, a)
   end
   lena = pol_length(a)
   lenb = pol_length(b)
   valb = valuation(b)
   vala = valuation(a)
   valr = min(vala, valb)
   precb = precision(b)
   preca = precision(a)
   prec = min(precb, preca)
   sa = scale(a)
   sb = scale(b)
   if lena == 1
      sa = sb
   elseif lenb == 1
      sb = sa
   end
   sc = gcd(gcd(sa, sb), abs(vala - valb))
   mina = min(vala + lena*sa, prec)
   minb = min(valb + lenb*sb, prec)
   lenr = max(mina, minb) - valr
   lenr = div(lenr + sc - 1, sc)
   R = base_ring(c)
   fit!(c, lenr)
   c.prec = prec
   c.val = valr
   c.scale = sc
   pa = vala
   pb = valb
   j = 0
   k = 0
   for i = 0: lenr - 1
      pi = valr + sc*i
      if pi == pa && pi < mina
         if pi == pb && pi < minb
            c.coeffs[i + 1] = add!(c.coeffs[i + 1], polcoeff(a, j), polcoeff(b, k))
            pb += sb
            k += 1
         else
            c.coeffs[i + 1] = deepcopy(polcoeff(a, j))
         end
         j += 1
         pa += sa
      elseif pi == pb && pi < minb
         c.coeffs[i + 1] = deepcopy(polcoeff(b, k))
         k += 1
         pb += sb
      else
         c.coeffs[i + 1] = R()
      end
   end
   c = set_length!(c, normalise(c, lenr))
   renormalize!(c)
   c = rescale!(c)
   return c
end

###############################################################################
#
#   Random elements
#
###############################################################################

const LaurentSeriesRingOrField = Union{LaurentSeriesRing,LaurentSeriesField}

RandomExtensions.maketype(S::LaurentSeriesRingOrField, ::AbstractUnitRange{Int}, _) = elem_type(S)

function RandomExtensions.make(S::LaurentSeriesRingOrField, val_range::AbstractUnitRange{Int}, vs...)
   R = base_ring(S)
   if length(vs) == 1 && elem_type(R) == Random.gentype(vs[1])
      Make(S, val_range, vs[1]) # forward to default Make constructor
   else
      Make(S, val_range, make(R, vs...))
   end
end

function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make3{<:RingElement,
                                         <:LaurentSeriesRingOrField,
                                         <:AbstractUnitRange{Int}}})
   S, val_range, v = sp[][1:end]
   R = base_ring(S)
   f = S()
   x = gen(S)
   for i = 0:S.prec_max - 1
      f += rand(rng, v)*x^i
   end
   return shift_left(f, rand(rng, val_range))
end

rand(rng::AbstractRNG, S::LaurentSeriesRingOrField, val_range::AbstractUnitRange{Int}, v...) =
   rand(rng, make(S, val_range, v...))

rand(S::LaurentSeriesRingOrField, val_range, v...) =
   rand(Random.default_rng(), S, val_range, v...)

###############################################################################
#
#   Conformance test element generation
#
###############################################################################

function ConformanceTests.generate_element(R::LaurentSeriesRing{BigInt})
  rand(R, 0:12, -10:10)
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

promote_rule(::Type{LaurentSeriesRingElem{T}}, ::Type{LaurentSeriesRingElem{T}}) where T <: RingElement = LaurentSeriesRingElem{T}

promote_rule(::Type{LaurentSeriesFieldElem{T}}, ::Type{LaurentSeriesFieldElem{T}}) where T <: FieldElement = LaurentSeriesFieldElem{T}

function promote_rule(::Type{LaurentSeriesRingElem{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
   promote_rule(T, U) == T ? LaurentSeriesRingElem{T} : Union{}
end

function promote_rule(::Type{LaurentSeriesFieldElem{T}}, ::Type{U}) where {T <: FieldElement, U <: RingElement}
   promote_rule(T, U) == T ? LaurentSeriesFieldElem{T} : Union{}
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (R::LaurentSeriesRing{T})(b::RingElement) where {T <: RingElement}
   return R(base_ring(R)(b))
end

function (R::LaurentSeriesField{T})(b::RingElement) where {T <: FieldElement}
   return R(base_ring(R)(b))
end

function (R::LaurentSeriesRing{T})() where {T <: RingElement}
   z = LaurentSeriesRingElem{T}(Vector{T}(undef, 0), 0, R.prec_max, R.prec_max, 1)
   z.parent = R
   return z
end

function (R::LaurentSeriesField{T})() where {T <: FieldElement}
   z = LaurentSeriesFieldElem{T}(Vector{T}(undef, 0), 0, R.prec_max, R.prec_max, 1)
   z.parent = R
   return z
end

function (R::LaurentSeriesRing{T})(b::Union{Integer, Rational, AbstractFloat}) where {T <: RingElement}
   if b == 0
      z = LaurentSeriesRingElem{T}(Vector{T}(undef, 0), 0, R.prec_max, R.prec_max, 1)
   else
      z = LaurentSeriesRingElem{T}([base_ring(R)(b)], 1, R.prec_max, 0, 1)
   end
   z.parent = R
   return z
end

function (R::LaurentSeriesField{T})(b::Union{Rational, AbstractFloat}) where {T <: FieldElement}
   if b == 0
      z = LaurentSeriesFieldElem{T}(Vector{T}(undef, 0), 0, R.prec_max, R.prec_max, 1)
   else
      z = LaurentSeriesFieldElem{T}([base_ring(R)(b)], 1, R.prec_max, 0, 1)
   end
   z.parent = R
   return z
end

function (R::LaurentSeriesRing{T})(b::T) where {T <: RingElem}
   parent(b) != base_ring(R) && error("Unable to coerce to power series")
   if iszero(b)
      z = LaurentSeriesRingElem{T}(Vector{T}(undef, 0), 0, R.prec_max, R.prec_max, 1)
   else
      z = LaurentSeriesRingElem{T}([b], 1, R.prec_max, 0, 1)
   end
   z.parent = R
   return z
end

function (R::LaurentSeriesField{T})(b::T) where {T <: FieldElem}
   parent(b) != base_ring(R) && error("Unable to coerce to power series")
   if iszero(b)
      z = LaurentSeriesFieldElem{T}(Vector{T}(undef, 0), 0, R.prec_max, R.prec_max, 1)
   else
      z = LaurentSeriesFieldElem{T}([b], 1, R.prec_max, 0, 1)
   end
   z.parent = R
   return z
end

function (R::LaurentSeriesRing{T})(b::LaurentSeriesElem{T}) where {T <: RingElement}
   parent(b) != R && error("Unable to coerce power series")
   return b
end

function (R::LaurentSeriesField{T})(b::LaurentSeriesElem{T}) where {T <: FieldElement}
   parent(b) != R && error("Unable to coerce power series")
   return b
end

function (R::LaurentSeriesRing{T})(b::Vector{T}, len::Int, prec::Int, val::Int, scale::Int, rescale::Bool=true) where {T <: RingElement}
   if length(b) > 0
      parent(b[1]) != base_ring(R) && error("Unable to coerce to power series")
   end
   z = LaurentSeriesRingElem{T}(b, len, prec, val, scale)
   z.parent = R
   if rescale
      z = rescale!(z)
   end
   return z
end

function (R::LaurentSeriesField{T})(b::Vector{T}, len::Int, prec::Int, val::Int, scale::Int, rescale::Bool=true) where {T <: RingElement}
   if length(b) > 0
      parent(b[1]) != base_ring(R) && error("Unable to coerce to power series")
   end
   z = LaurentSeriesFieldElem{T}(b, len, prec, val, scale)
   z.parent = R
   if rescale
      z = rescale!(z)
   end
   return z
end

###############################################################################
#
#   power_series_ring constructor
#
###############################################################################

function laurent_series_ring(R::AbstractAlgebra.Ring, prec::Int, s::VarName; cached::Bool=true)
   T = elem_type(R)

   parent_obj = LaurentSeriesRing{T}(R, prec, Symbol(s), cached)

   return parent_obj, gen(parent_obj)
end

function laurent_series_ring(R::AbstractAlgebra.Field, prec::Int, s::VarName; cached::Bool=true)
   T = elem_type(R)

   parent_obj = LaurentSeriesField{T}(R, prec, Symbol(s), cached)

   return parent_obj, gen(parent_obj)
end

function laurent_series_field(R::AbstractAlgebra.Field, prec::Int, s::VarName; cached::Bool=true)
   T = elem_type(R)

   parent_obj = LaurentSeriesField{T}(R, prec, Symbol(s), cached)

   return parent_obj, gen(parent_obj)
end
